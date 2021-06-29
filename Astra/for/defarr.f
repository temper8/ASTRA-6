C Modules:	>>>>>>> DEFARR, IFSTEP, OLDNEW <<<<<<<
C			TRANSF, CHEBFT, CHEBPC, PFITN
C			RUNN,	RUNT,	RUNF,
C			ISTORY,	APPST,	APPSR
C----------------------------------------------------------------------|
	subroutine	DEFARR
C----------------------------------------------------------------------|
C 1) Check positiveness of Z_eff, n_e, n_i, T_e, T_i.  Stop if negative.
C 2) Extend definition of all standard arrays beyond NA1.
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer j,js,jj,length
	character str*2
	double precision	YV,YF,YMU,YN,YNE,YNI,YTE,YTI,YZF
	call	markloc("DEFARR"//char(0))
	YZF = 0.
	YNE = 0.
	YNI = 0.
	YTE = 0.
	YTI = 0.
	YV = GP2*BTOR*HRO
	do	j=1,NA1
	   SQEPS(j) = SQRT(AMETR(j)/(RTOR+SHIF(j)))
! VP = c<E_||/h>/B_p	in Hinton & Hazeltine Rev.Mod.Phys.48(1976) p.297
	   if (MU(j) .le. 0.d0)	then
	      write(*,'(/A,F10.6,A)')'>>> ERROR >>>  Time =',TIME,' sec'
	      write(*,'(2A,F10.6)')
     >	'               The rotational transform is less or equal zero '
     >,			'at rho_N =',RHO(j)/ROC
	      if (TIME .le. TSTART+TAU/2.)	write(*,'(A/2A,1H"/)')
     >	'               Most probably it has not been properly defined'
     >,	'               by the data file "',RDNAME(1:length(RDNAME))
	      call	a_stop
	   endif
	   VP(j) = ULON(j)/(YV*j*MU(j))
	   YZF = max(YZF,ZEF(j))
	   YNE = max(YNE,NE(j))
	   YNI = max(YNI,NI(j))
	   YTE = max(YTE,TE(j))
	   YTI = max(YTI,TI(j))
	enddo
	str = '3M'
	if (min(YZF,YNE,YNI,YTE,YTI) .gt. 0.)	then
	   goto	11
	elseif (YZF .le. 0.)	then
	   str = "ZF"
	elseif (YNE .le. 0.)	then
	   str = "NE"
	elseif (YNI .le. 0.)	then
	   str = "NI"
	elseif (YTE .le. 0.)	then
	   str = "TE"
	elseif (YTI .le. 0.)	then
	   str = "TI"
	endif
	write(*,'(/A,F6.3)')'>>> ERROR >>>  Time =',TIME
	if (str .ne. '3M' .and. str .ne. 'ZF') write(*,'(3A/)')
     >	'               The variable "',str,'" is less or equal zero'
	if (str .eq. '3M') write(*,'(A/)')
     >	'               The 3M equilibrium solver does not converge'
	if (str .eq. 'ZF') write(*,'(A/)')
     >	'               The variable "ZEF" is less or equal zero'
	call	a_stop

 11	continue
	if (NA1 .eq. NB1)	goto	12

C Determine arrays beyond ABC
	YF = GP2*BTOR
	YMU = (MU(NA1)-MV(NA1))*ROC**2
	js = 1
	do	j=1,NSBR
	   if (DTNAME(24+4*j) .eq. 'NEUTAB')	js = 3
	enddo
	do	j=NA1+1,NB1
	   do	jj=10,80		! Extended to 80
		STAARR(j,jj)=STAARR(NA1,jj)
	   enddo
C   Common /A_IONS/: NN -> VIMP3	Exclude NN,TN if NEUTAB is active
	   do	jj=js,31
		EXTARR(j,jj)=EXTARR(NA1,jj)
	   enddo
C   Common /A_CARS/: CAR1 -> CAR32
	   do	jj=1,32
		CAR(j,jj)=CAR(NA1,jj)
	   enddo

	   if (js.eq.3 .and. j.gt.NAB)	then
	      NN(j) = NN(NAB)
	      TN(j) = TN(NAB)
	   endif
	   YN = exp((ABC-AMETR(j))/WNE)
	   NE(j) = NE(NA1)*YN
	   NI(j) = NI(NA1)*YN
	   NALF(j) = NALF(NA1)*YN
	   NHE3(j) = NHE3(NA1)*YN
	   NHYDR(j) = NHYDR(NA1)*YN
	   NDEUT(j) = NDEUT(NA1)*YN
	   NTRIT(j) = NTRIT(NA1)*YN
	   TE(j) = TE(NA1)*exp((ABC-AMETR(j))/WTE)
	   TI(j) = TI(NA1)*exp((ABC-AMETR(j))/WTI)
C Should be changed when a SOL equilibrium is implemented
	   MV(j) = MV(NA1)
	   FV(j) = FV(j-1)+YF*RHO(j)*MV(j)*(RHO(j)-RHO(j-1))
	   MU(j) = MV(j)+YMU/RHO(j)**2
	   FP(j) = FP(j-1)+YF*RHO(j)*MU(j)*(RHO(j)-RHO(j-1))
	   CU(j) = 0.
	   CUBS(j) = 0.
	   CV(j) = 0.
	   CD(j) = 0.
	   CUTOR(j) = 0.
	enddo
 12	continue
	end
C======================================================================|
	double precision function ARRNA1(ARR,H)
C----------------------------------------------------------------------|
C Compute edge value of the array ARR(NA1) assuming that d^3(ARR)/dr^3=0
C----------------------------------------------------------------------|
	double precision ARR(*),H
	print *, ARR
	print *, ARR(-2),  ARR(-1), ARR(0), ARR(1)
	ARRNA1 = ARR(1)+H*(2.*ARR(1)+ARR(-2)-3.*ARR(-1))
	end
C======================================================================|
	integer	function	IFSTEP(IFCONV)
C	IFCONV dummy parameter
C Input
C       LEQ(1)  LEQ(2)  LEQ(3)  LEQ(4)  LEQ(5)  LEQ(6-9)
C	 NE	 TE	 TI	 CU	equil    free
C       LEQ(10) LEQ(11) LEQ(12), etc
C        F0	 F1	 F2
C       LEQ(j) - usage of j-th equation in the model
C	       =-1 no equation (default)
C	       = 0 type: AS
C	       = 1 type: EQ
C	       = 2 type: FU
C	       = 3 type: heat conductivity flux + heat convection
C Returned value
C	 1 - time step made, normal exit
C	 0 - time step will be repeated
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	double precision	CTAU,TAUO,YY,YDTIME
	integer IFCONV,j,jj,ICALL
	save	YDTIME,ICALL
	data	ICALL/0/
	call	markloc("IFSTEP"//char(0))
	if (ICALL .eq. 0)	then
	   ICALL = 1
	   YDTIME = TIME
	endif
	CTAU = 1./TAUINC
	do	j = 1,NA
	   if(LEQ(1).gt.0)  CTAU = MAX(CTAU,ABS(NEO(j)/NE(j)-1.)/DELVAR)
	   if(LEQ(2).gt.0)  CTAU = MAX(CTAU,ABS(TEO(j)/TE(j)-1.)/DELVAR)
	   if(LEQ(3).gt.0)  CTAU = MAX(CTAU,ABS(TIO(j)/TI(j)-1.)/DELVAR)
C	   if(LEQ(4).gt.0)  CTAU = MAX(CTAU,ABS(FPO(j)/FP(j)-1.)/DELVAR)
	   do	jj=0,9
	      if (LEQ(jj+10) .le. 0)	goto	1
	      YY = 0.5*(abs(FJO(j,jj))+abs(FJ(j,jj)))
	      if (YY .lt. 1.E-6)	then		! Allow zero FJ
		 YY = abs(FJO(j,jj)-FJ(j,jj))
	      else
		 YY = abs(FJO(j,jj)-FJ(j,jj))/YY
	      endif
	      CTAU = MAX(CTAU,YY/DELVAR)
 1	      continue
	   enddo
	enddo
	TAUO = TAU
	TAU = MIN(TAUMAX,TAU/CTAU,DTOUT,DPOUT)
	TAU = MAX(TAUMIN,TAU)
C The condition CTAU >= TAUINC**2 can be fulfilled if at least 2 of unknowns
C have relative variation above DELVAR
	if (CTAU.le.TAUINC**2 .or. TAUO.le.TAUMIN*TAUINC)	then
	   YDTIME = YDTIME+TAU
	   TIME = YDTIME
	   IFSTEP = 1
	   NSTEPS = NSTEPS+1
	   return		! -> proceed to the next time step
	endif

C repeat internal loop in "equftn.tmp":
	do	j=1,NB1
	   NE(j) = NEO(j)
	   NI(j) = NIO(j)
	   TE(j) = TEO(j)
	   TI(j) = TIO(j)
	   FP(j) = FPO(j)
	   VR(j) = VRO(j)
	   do	jj=0,9
	      FJ(j,jj) = FJO(j,jj)
	   enddo
	enddo
	IFSTEP = 0
	end
C======================================================================|
	subroutine	OLDNEW
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	integer j,jj
C----------------------------------------------------------------------|
	call	markloc("OLDNEW"//char(0))
	do	j=1,NB1
	   NEO(j) = NE(j)
	   NIO(j) = NI(j)
	   TEO(j) = TE(j)
	   TIO(j) = TI(j)
	   FPO(j) = FP(j)
	   do	jj=0,9
	      FJO(j,jj) = FJ(j,jj)
	   enddo
	enddo
	end
C======================================================================|
	subroutine	TRANSF(NO,FO,XO,N,FN,XN)
C-----------------------------------------------------------------------
C		The subroutine transmits a function FO(1:NO)
C		from an arbitrary grid XO(1:NO) to a functon F(1:N)
C		on an another arbitrary grid XN(1:N) by the method
C		of quadratic interpolation.
C	Input:	NO, FO(1:NO), XO(1:NO), N, XN(1:N)
C	Output:	FN(1:N)
C-----------------------------------------------------------------------
C       Note:	[XO(NO)-XO(1)]*[XN(N)-XN(1)] must be > 0
C-----------------------------------------------------------------------
	implicit none
	integer	NO,N,I,j
	double precision	XO(*),FO(*),XN(*),FN(*),YF1,YF2,YF3
	I = 1
	YF1 = FO(1)/((XO(1)-XO(2))*(XO(1)-XO(3)))
	YF2 = FO(2)/((XO(2)-XO(1))*(XO(2)-XO(3)))
	YF3 = FO(3)/((XO(3)-XO(1))*(XO(3)-XO(2)))
	do	5	j=1,N
		if(2.*XN(j).le.XO(I+1)+XO(I+2))	go to	4
 3		I = I+1
		if(I.gt.NO-2)	I=NO-2
	if(2.*XN(j).gt.XO(I+1)+XO(I+2).and.I.lt.NO-2)	go to	3
		YF1 = FO(I)/((XO(I)-XO(I+1))*(XO(I)-XO(I+2)))
		YF2 = FO(I+1)/((XO(I+1)-XO(I))*(XO(I+1)-XO(I+2)))
		YF3 = FO(I+2)/((XO(I+2)-XO(I))*(XO(I+2)-XO(I+1)))
 4		FN(j) = YF1*(XN(j)-XO(I+1))*(XN(j)-XO(I+2))
     .			+YF2*(XN(j)-XO(I))*(XN(j)-XO(I+2))
     .			+YF3*(XN(j)-XO(I))*(XN(j)-XO(I+1))
 5	continue
	end
C=======================================================================
	subroutine	CHEBFT(NIN,FIN,XIN,NCHCF,CHEBCF)
C-----------------------------------------------------------------------
C	The subroutine returns NCHCF coefficients of a Chebyshev
C	polynomial fit to the function FIN(1:NIN) given as a
C	function of any "radial" variable on the grid XIN(1:NIN)
C Example:
C	call CHEBFT(NA1,TE,FP,5,CHOUT)
C	out = PFITN(FP,CHOUT,5)
C
C-----------------------------------------------------------------------
      implicit none
      integer	NIN,NCHCF,NMAX,K,J
      double precision	FIN(NIN),XIN(NIN),CHEBCF(NCHCF),YA,YB,FAC,YD
      parameter	(NMAX=10)
      double precision	YC(NMAX),YF(NMAX),SUM,PI,PIOVN,BMA,BPA
      data	PI/3.141592654/
      if (NCHCF .gt. NMAX)	then
	 write(*,*)'>>> Chebyshev fit error: too high power'
	 call	a_stop
      endif
      YA = max(XIN(1),XIN(NIN))
      YB = min(XIN(1),XIN(NIN))
      PIOVN = PI/NCHCF

C The following two lines require that at the boundary
C     the fit concides with the original function
	YD = COS(0.5*PIOVN)
	YA = (2.*YA-YB*(1.-YD))/(1.+YD)

      BMA=0.5*(YB-YA)
      BPA=0.5*(YB+YA)
      do 11 K=1,NCHCF
        YC(K) = BPA+BMA*COS(PIOVN*(K-0.5))
 11   continue
	call	TRANSF(NIN,FIN,XIN,NCHCF,YF,YC)

      FAC=2./NCHCF
      do 13 J=1,NCHCF
        SUM=0.
        do 12 K=1,NCHCF
          SUM=SUM+YF(K)*COS(PIOVN*(K-0.5)*(J-1.))
 12     continue
        YC(J)=FAC*SUM
 13   continue

      call   CHEBPC(YC,CHEBCF,YF,NCHCF)

      FAC = 1./BMA
      do 14 J=2,NCHCF
        CHEBCF(J)=CHEBCF(J)*FAC
        FAC=FAC/BMA
 14   continue
      do 16 J=1,NCHCF-1
        do 15 K=NCHCF-1,J,-1
          CHEBCF(K)=CHEBCF(K)-BPA*CHEBCF(K+1)
 15     continue
 16   continue
      return
      end
C=======================================================================
      subroutine	CHEBPC(C,D,F,N)
C-----------------------------------------------------------------------
C	F(1:N) array for internal use
C
C   Input:	N	number of coefficients
C		C(1:N)	
C   Output:	D(1:N)	array of Chebyshev coefficients
C-----------------------------------------------------------------------
      implicit none
      integer  J,N,K
      double precision C(N),D(N),F(N),SV
      do 11 J=1,N
        D(J)=0.
        F(J)=0.
11    continue
      D(1)=C(N)
      do 13 J=N-1,2,-1
        do 12 K=N-J+1,2,-1
          SV=D(K)
          D(K)=2.*D(K-1)-F(K)
          F(K)=SV
12      continue
        SV=D(1)
        D(1)=-F(1)+C(J)
        F(1)=SV
13    continue
      do 14 J=N,2,-1
        D(J)=D(J-1)-F(J)
14    continue
      D(1)=-F(1)+0.5*C(1)
      return
      end
C======================================================================|
      double precision function PFITN(x,cp,N)
C-----------------------------------------------------------------------
C Nth order Polynomial FIT to a function f(a)
C CP are the N polynomial coefficients given
C PFITN = f(a) = \Sum {CP_j * a^(j-1)} ;   1<=j<=N
C-----------------------------------------------------------------------
      implicit none
      integer N,j
      double precision x,cp(N)
      PFITN = CP(N)
      do   j  = N-1,1,-1
	PFITN = PFITN*x+CP(j)
      enddo
      end
C======================================================================|
	subroutine	RUNN(A,B,DS,C,D,NO,N,GT,H,E,NN,VO,VN,G11)
C----------------------------------------------------------------------|
C Example call:
C      SNN(J)=SNNEU
C      YWA(J)=DN(J)
C      YWB(J)=-CN(J)		! CN=-VP*VRHH
C      NE(ND1)=...
C      NEO(ND1)=NE(ND1)
C      YWC(4)=1.
C      call RUNN(YWA,YWB,SNN,SN,YWD,NEO,ND,TAU,HRO,YWC,NE,VRO,VR,G11)
C----------------------------------------------------------------------|
C Exponential scheme
C----------------------------------------------------------------------|
C	The subroutine provides run of the equation
C	dn/dt=1/V'*d[V'<(\nabla\rho)^2>*(A*dn/dr+B*n)]/dr+C*n+D
C	Scheme:		h_j/GT*V'(j)*(NN(j)-NO(j))=
C	   =GT*G11(j+1/2)*A(j+1/2)/h(j+1/2)
C		*(NN(j+1)*f(j+1/2)-NN(j)*g(j+1/2))-
C	   -GT*G11(j-1/2)*A(j-1/2)/h(j-1/2)
C		*(NN(j)*f(j-1/2)-NN(j-1)*g(j-1/2))-
C	   +h_j*(V'(j)*C(j)*NN(j)+D(j))
C	Here	A,B,C,D		- arrays(1:N)
C		NOld,NNew 	- arrays(1:N+1)
C		VOld,VNew 	- arrays(1:N+1)
C	Input:	H	- radial step (m)
C		HB	- edge grid cell size (m)
C		N+1	- number of mesh points
C		GT	- time step (sec)
C		A(N),B(N),C(N),D(N),VO(N),VN(N),NO(N),NN(N+1),E(1:4)
C		E(1) = HROA
C		E(4) < 0 if (NEB isn't set) .and. (QNB .or. QNNB is set)
C			then E(2)=QNB or E(3)=QNNB
C	Not used: SN(N+1)
C		  SNN(N+1)
C	Output:	NN(N) - new quantity
C		A(N)  - coefficient at (j+1) for flux calculation
C		B(N)  - coefficient at  (j)  for flux calculation
C	Internal use:
C		E(*)
C----------------------------------------------------------------------|
	implicit none
	integer	N,j
	double precision
     1	       A(*),B(*),DS(*),C(*),D(*),NO(*),NN(*),VO(*),VN(*),E(*),
     2	       G11(*),H,HB,Q2,Q3,BC,HJ,GT,G0,G1,P0,P1,Q0,Q1,
     3	       Y0,Y1,YA,YB,YS,YJ,AJ,BJ,CJ,DJ
	HB = E(1)
	Q2 = E(2)
	Q3 = E(3)
	BC = E(4)
	G1 = 0.
	P1 = 0.
	Q1 = 0.
	Y1 = 0.
	YS = 0.
	HJ = H
	do	10	J=1,N
	   G0 = G1
	   P0 = P1
	   Q0 = Q1
	   G1 = GT*G11(j)
	   YA = A(j)
	   YB = B(j)
	   if (YA .lt. 0.d0)	then
              Y1 = j*H
              goto	99
           endif
C Power-law scheme: (1 line)
	   P1 = 0.5*(abs(YB)+YB)		! 0.5(|B|+B)
	   if (YA .eq. 0.d0)	goto	1
	   if (j .eq. N)	HJ = HB
C Exponential scheme: (7 lines)
	   YJ = HJ*YB/YA			! |\xi|
	   if (abs(YJ) .ge. 4.d1) goto	1	! Use (A/h)*f(\xi) = .5*(|B|+B)
	   if (abs(YJ) .ge. 1.d-5) then
	      P1 = YB/(1.-exp(-YJ))	   ! Use (A/h)*f(\xi) = B/(1-exp{-\xi})
	   else
	      P1 = YA/HJ*(1.+0.5*YJ)   	   ! Use (A/h)*f(\xi) = A/h/(1-\xi/2)
	   endif
C Power-law scheme: (6 lines)
C	   YJ = abs(HJ*YB/YA)			! \xi = h B/A
C	   if (YJ .ge. 1.d1)	goto	1
C	   Y1 = 1.d0-1.d-1*YJ			! (1 - 0.1|\xi|)
C	   Y2 = Y1*Y1
C	   YJ = Y2*Y2*Y1*YA/HJ			! (A/h)*(1 - 0.1|\xi|)^5
C	   P1 = P1+YJ				! (A/h)*f(\xi) is done
 1	   continue
	   Q1 = P1-YB
	   Y0 = Y1
	   Y1 = DS(j)/HJ
	   AJ = G1*(P1+Y1)
	   BJ = G1*(Q1+Y1)+G0*(P0+Y0)+H*(VN(J)-C(J)*VN(J)*GT)
	   CJ = G0*(Q0+Y0)
	   YJ = YS
	   YS = G1*Y1*(NO(j+1)-NO(j))
	   DJ = H*(D(J)*VN(J)*GT+NO(J)*VO(J))-YS+YJ
	   if(J .ne. 1)	then
	      BJ = BJ-CJ*E(J-1)
	      DJ = DJ+CJ*NN(J-1)
	   endif
	   E(J) = AJ/BJ
	   NN(J) = DJ/BJ
	   A(j) = P1
	   B(j) = Q1
 10	continue
	if (BC .lt. 0.)	NN(N+1) = (G11(N)*Q1*NN(N)-Q2)
     .			     /(G11(N)*(P1-E(N)*Q1)+Q3)
	do	J=N,1,-1
	   NN(J) = E(J)*NN(J+1)+NN(J)
	enddo
	return
 99     write(*,'(2A,F10.6,A,I4,A,F10.6)')
     &	  " >>> ERROR >>> Diffusion coefficient shouldn't be negative.",
     &	  "  RHO =",Y1,"  node =",j,"   D =",YA
	call	IFKEY(ichar(' '))
	end
C======================================================================|
	subroutine	RUNN1(A,B,C,D,G,NO,N,GT,H,E,NN,VO,VN,G11)
C----------------------------------------------------------------------|
C Linear scheme
C----------------------------------------------------------------------|
C this is a copy of old RUNN
C	The subroutine provides run of the equation
C	dn/dt=1/V'*d[V'<(\nabla\rho)^2>*(A*dn/dr+B*n)]/dr+C*n+D
C	Scheme:		H*H/GT*V'(j)*(NN(j)-NO(j))=
C		=V'(j+1/2)*A(j+1/2)*(NN(j+1)-NN(j))-
C		-V'(j-1/2)*A(j-1/2)*(NN(j)-NN(j-1))+
C		+0.5*H*V'(j+1/2)*B(j+1/2)*(NN(j+1)+NN(j))-
C		-0.5*H*V'(j-1/2)*B(j-1/2)*(NN(j)+NN(j-1))-
C		+H*H*(V'(j)*C(j)*NN(j)+D(j))
C	Here	A,B,C,D		- arrays(1:N)
C		NOld,NNew 	- arrays(1:N+1)
C		VOld,VNew 	- arrays(1:N+1)
C	Input:	H	- radial step (m)
C		HB	- edge grid cell size (m)
C		N+1	- number of mesh points
C		GT	- time step (sec)
C		A(N),B(N),C(N),D(N),VO(N),VN(N),NO(N)
C		NN(N+1),E(1),E(2),E(3),E(4)
C	Output:	NN(N)
C----------------------------------------------------------------------|
	implicit none
	integer	N,j,N1
	double precision
     1		A(*),B(*),C(*),D(*),G(*),NO(*),NN(*),VO(*),VN(*),E(4),
     2		G11(*),H,HB,GT,RA,RB,HH,H05,RA0,RB0,YG11,AJ,BJ,CJ,DJ
	HB = E(1)
	RA = 0.
	RB = 0.
	HH = H*H
	H05 = 0.5*H
	do	1	J=1,N
		RA0 = RA
		RB0 = RB
		YG11 = G11(J)*GT
		RA = A(J)*YG11
		if (j .eq. N)	then
			H05 = 0.5*HB
			HH  = H05*(H+HB)
			RA0 = RA0*HB/H
			RB0 = RB0*HB/H
		endif
		RB = B(J)*YG11*H05
		AJ = RA+RB
		BJ = RA+RA0-RB+RB0+HH*(VN(J)-C(J)*VN(J)*GT)
C		AJ = RA
C		BJ = RA+RA0-RB+RB0+HH*(VN(J)-C(J)*VN(J)*GT)
C		if (j .eq. N)	then
C		   BJ = BJ-RB
C		else
C		   AJ = AJ+RB
C		endif
		CJ = RA0-RB0
		DJ = HH*(D(J)*VN(J)*GT+NO(J)*VO(J))
		if(J .ne. 1)	then
			BJ = BJ-CJ*G(J-1)
			DJ = DJ+CJ*NN(J-1)
		endif
		G(J) = AJ/BJ
		NN(J) = DJ/BJ
 1	continue
	if (E(4) .lt. 0.)	NN(N+1) =
     .	  ((RA-RB)*NN(N)-HB*GT*E(2))/(RA+RB+G(N)*(RB-RA)+HB*GT*E(3))
C	RB = B(N)*G11(N)*GT*HB
C	if (E(4) .lt. 0.)	NN(N+1) = ((RA-RB)*NN(N)-GT*HB*E(2))
C     .				    /(RA+G(N)*(RB-RA)+GT*HB*E(3))
	do	2	J=N,1,-1
 2	NN(J) = G(J)*NN(J+1)+NN(J)
	end
C======================================================================|
	subroutine RUNT(A,B,C,P,G,NO,NN,TO,N,GT,H,Q,TN,VO,VN,G11)
C----------------------------------------------------------------------|
C Exponential scheme
C----------------------------------------------------------------------|
C	The subroutine provides run of the equation:
C		d(N*T)/dt=1/V'*d[V'*(A*dT/dr+B*T)]/dr+625.*(C*T+P)
C		V'(j)*(NN(j)*TN(j)-NO(j)*TO(j))*H*H/GT=
C			=V'(j+1/2)*A(j+1/2)*(TN(j+1)-TN(j))-
C			-V'(j-1/2)*A(j-1/2)*(TN(j)-TN(j-1))+
C			+0.5*V'(j+1/2)*B(j+1/2)*(TN(j+1)+TN(j))-
C			-0.5*V'(j-1/2)*B(j-1/2)*(TN(j)+TN(j-1))-
C			+V'(j)*H*H*625.*(C(j)*TN(j)+P(j))
C
C	Here	A,B,C,P		- arrays(1:N)
C		NOld,NNew 	- arrays(1:N+1)		(densities)
C		TOld,TNew 	- arrays(1:N+1)		(temperatures)
C		VOld,VNew 	- arrays(1:N+1)		(dV/drho)
C	Input:	H	- radial step (m)
C		N	- number of mesh points
C		GT	- time step (sec)
C		A(N),B(N),C(N),P(N),NO(N),NN(N),VO(N),VN(N),TO(N+1)
C		TN(N+1),Q(1),Q(2),Q(3),Q(4)
C	Output:	TN(N)
C----------------------------------------------------------------------|
	implicit none
	double precision
     1		A(*),B(*),C(*),P(*),G(*),NO(*),NN(*),TO(*),TN(*),Q(*),
     2		VO(*),VN(*),G11(*),H,HB,GT,GT23,Y625,AJ,BJ,CJ,DJ,
     3		H1,HJ,YA,YB,Y1,Y2,YJ,P0,P1,Q0,Q1,G0,G1
	integer	N,j
	HB = Q(1)
	GT23 = GT/1.5
	Y625 = 625.*GT23
	G1 = 0.
	P1 = 0.
	Q1 = 0.
	H1 = H
	HJ = H
	do	10	J=1,N
	   G0 = G1
	   P0 = P1
	   Q0 = Q1
	   if (j .eq. N)	then
C	      HJ = 0.5*(H+HB)
	      H1 = HB
	   endif
	   G1 = GT23*G11(j)
C	   if (DV(N+1) .gt. 0.)	then
C	      Y1 = 0.5*(NN(j)+NN(j+1))*DV(j)
C	      YA = A(j)+Y1
C	      YB = B(j)-Y1*log(TO(j+1)/TO(j))/H1
C	   else
	      YA = A(j)
	      YB = B(j)
C	   endif
	   if (YA .lt. 0.d0)	goto	99
	   P1 = 0.5*(abs(YB)+YB)
	   if (YA .eq. 0.d0)	goto	1
C Exponential scheme: (7 lines)
	   YJ = H1*YB/YA			! |\xi|
	   if (abs(YJ) .ge. 4.d1) goto	1	! Use (A/h)*f(\xi) = .5*(|B|+B)
	   if (abs(YJ) .ge. 1.d-5) then
	      P1 = YB/(1.-exp(-YJ))	   ! Use (A/h)*f(\xi) = B/(1-exp{-\xi})
	   else
	      P1 = YA/H1*(1.+0.5*YJ)   	   ! Use (A/h)*f(\xi) = A/h/(1-\xi/2)
	   endif
C Power-law scheme: (6 lines)
C	   YJ = abs(H1*YB/YA)			! |\xi|
C	   if (YJ .ge. 1.d1)	goto	1
C	   Y1 = 1.d0-.1d0*YJ
C	   Y2 = Y1*Y1
C	   YJ = Y2*Y2*Y1*YA/H1
C	   P1 = YJ+P1
 1	   continue
	   Q1 = P1-YB
	   AJ = G1*P1
	   BJ = G1*Q1+G0*P0+HJ*VN(J)*(NN(J)-Y625*C(J))
	   CJ = G0*Q0
	   YJ = (VO(J)/VN(J))**0.666667
	   DJ = HJ*(Y625*P(J)*VN(J)+NO(J)*TO(J)*VO(J)*YJ)
	   if(J .ne. 1)	then
	      BJ = BJ-CJ*G(J-1)
	      DJ = DJ+CJ*TN(J-1)
	   endif
	   G(J) = AJ/BJ
	   TN(J) = DJ/BJ
	   A(j) = P1
	   B(j) = Q1
 10	continue
C Here Q(2)=QjB, Q(3)=QjTB, where j should be replaced by E or I
	if (Q(4) .lt. 0.)	TN(N+1) = (G11(N)*Q1*TN(N)-625.*Q(2))
     .				     /(G11(N)*(P1-G(N)*Q1)+625.*Q(3))
C or equivalent
C	if (Q(4) .lt. 0.)
C     .	   TN(N+1) = (G1*Q1*TN(N)-Y625*Q(2))/(AJ-G(N)*G1*Q1+Y625*Q(3))
	do	J=N,1,-1
	   TN(J) = G(J)*TN(J+1)+TN(J)
	enddo
	do	J=1,N
	   Q(j) = -0.0016*G11(j)*(A(j)*TN(j+1)-B(j)*TN(j))
	enddo
	Q(N+1) = Q(N)
	do	J=1,N
	   P(j) = C(j)*TN(j)+P(j)
	enddo
	P(N+1) = P(N)
	do	J=1,N+1
	   TN(J) = MAX(TN(J),1.d-4)
	enddo
	return
 99	write(*,'(2A,1F10.6)')
     &	   " >>> ERROR >>> Heat conductivity shouldn't be negative.",
     &	   "  RHO =",j*H
	call	IFKEY(ichar(' '))
	end
C======================================================================|
	double precision function GETPEI(j)
	implicit none
	integer	 j
	double precision COULG,PEI
	include	'for/parameter.inc'
	include	'for/status.inc'
C	COULG=15.9-.5*log(NE(J))+log(TE(J))
C	PEI=0.00246*COULG*NE(J)*NI(J)*ZMAIN(J)*ZMAIN(J)/
C     .		(AMAIN(J)*TE(J)*sqrt(TE(J)))
	include	'fml/pei'
	GETPEI = PEI
	end
C======================================================================|
	subroutine RUNTT(A,B,C,D,NO,NN,TO,N,GT,H,HB,DV,DS,VO,VN,G11,W)
C----------------------------------------------------------------------|
C	call RUNTT (YWA,YWB,PET,YWC,NEO,NE,TEO,
C       		ND,TAU,HRO,QE(1),YWD,DSE,VRO,VR,G11,WORK1)
C----------------------------------------------------------------------|
C Exponential scheme
C----------------------------------------------------------------------|
C	The subroutine makes time step in the matrix equation:
C		d(N*T)/dt=1/V'*d[V'*(A*dT/dr+B*T)]/dr+625.*(C*T+D)
C
C	Here	A,B,C,D		- arrays(1:N)
C		NOld,NNew 	- arrays(1:N+1)		(densities)
C		VOld,VNew 	- arrays(1:N+1)		(dV/drho)
C		TO	 	- array (1:N+1)		(temperature)
C		W(*)		- work space (e.g., work1(*))
C	Input:	H	- radial step (m)
C		N	- number of mesh points
C		GT	- time step (sec)
C       	A(1:N)  - diffusivity n_e*\chi_e (or n_i*\chi_i)
C       	B(1:N)  - convective velocity e.g. 5/2*GNX(J)*SLAT(J)/G11(J)
C       	C(1:N)  - PET or PIT without equipartition
C       	D(1:N)  - PE  or PI
C       	NO(1:N) - old density (previous time step)
C       	NN(1:N) - new density (next time step)
C       	VO(1:N) - old V'
C       	VN(1:N) - new V'
C       	DV(1:N) - Aux. heat conductiviy compensated by advection
C       	DS(1:N) - Aux. heat conductiviy compensated by source
C       	DV(N+1) - Enable DV treatment if DV(N+1) is nonzero
C       	DS(N+1) - Enable DS treatment if DS(N+1) is nonzero
C       	TO(1:N+1) - old T_e (or T_i) 
C		HB 	  - edge cell size
C	Output:	DV(1:N) - contribution to rhs due to DV
C----------------------------------------------------------------------|
	implicit none
	double precision
     1		A(*),B(*),C(*),D(*),NO(*),NN(*),TO(*),DV(*),DS(*),VO(*),
     2		VN(*),G11(*),H,HB,GT,GT23,Y625,AJ,BJ,CJ,DJ,W0,W1,P0,P1,
     3		Q0,Q1,G0,G1,H1,HJ,YA,YB,Y0,Y1,Y2,YS,YJ,Y11,Y12,Y21,Y22
	integer	N,j,j1,jn,icall,N0
	double precision GETPEI,W(N,*)
	external	 GETPEI
	save	icall,N0
	data	icall/0/
        call add2loc("Subroutine RUNTT"//char(0))
	GT23 = GT/1.5
	Y625 = 625.*GT23
	HJ = H
	G1 = 0.
	P1 = 0.
	Q1 = 0.
	H1 = H
	CJ = 0.
	W1 = 0.
	Y1 = 0.
	YS = 0.
	do	2	J=1,N
	   G0 = G1
	   P0 = P1
	   Q0 = Q1
	   Y0 = Y1
	   if (j .eq. N)	H1 = HB
	   G1 = GT23*G11(j)
	   YA = A(j)
	   YB = B(j)
	   if (DV(N+1) .gt. 0.)	then
	      Y11 = DV(j)*0.5*(NN(j)+NN(j+1))
	      Y12 = Y11*log(TO(j)/TO(j+1))/H1
	      YA = YA+Y11
	      YB = YB+Y12
	      Y2 = TO(j)-TO(j+1)
	      if (abs(Y2) .lt. 1.d-6)	then
		 DV(j) = Y11*2./(TO(j)+TO(j+1))
	      else
		 DV(j) = Y12/Y2		   ! Contribution to the source
	      endif
	   endif
	   if (YA .lt. 0.d0)	goto	99
	   P1 = 0.5*(abs(YB)+YB)	   ! 0.5(|B|+B) (Used if vh/D is big)
	   if (YA .eq. 0.d0)	goto	1
	   YJ = H1*YB/YA		   ! |\xi| Peclet number
	   if (abs(YJ) .ge. 4.d1) goto	1  ! Use (A/h)*f(\xi) = .5*(|B|+B)
	   if (abs(YJ) .ge. 1.d-5) then
	      P1 = YB/(1.-exp(-YJ))	   ! Use (A/h)*f(\xi) = B/(1-exp{-\xi})
	   else
	      P1 = YA/H1*(1.+0.5*YJ)   	   ! Use (A/h)*f(\xi) = A/h/(1-\xi/2)
	   endif
 1	   continue
	   Q1 = P1-YB			   ! (A/h)*g(\xi) = (A/h)*f(\xi)-B
	   if (DS(N+1) .gt. 5.d-1)	then
	      Y1 = DS(j)*0.5*(NN(j)+NN(j+1))/H1
C	      A(j) = 1.6d-3*G11(j)*Y1		
	   endif

	   if (icall .eq. 0)	then
	      N0 = N
	      W(j,17) = P1*G11(j)
	      W(j,18) = Q1*G11(j)
	      W(j,10) = -GETPEI(J)
	   else
	      if (N0.ne.N)	goto	98
	      W(j,19) = P1*G11(j)
	      W(j,20) = Q1*G11(j)
	   endif
C	   AJ = G1*P1
	   AJ = G1*(P1+Y1)
C	   BJ = G1*Q1+G0*P0
	   BJ = G1*(Q1+Y1)+G0*(P0+Y0)
C	   CJ = G0*Q0
	   CJ = G0*(Q0+Y0)
	   BJ = BJ+HJ*VN(J)*(NN(J)-Y625*(W(j,10)+C(J)))
	   YJ = (VO(J)/VN(J))**0.666667
	   DJ = HJ*(Y625*D(J)*VN(J)+NO(J)*TO(J)*VO(J)*YJ)
	   YJ = YS
	   YS = G1*Y1*(TO(j+1)-TO(j))
	   DJ = DJ-YS+YJ
	   W(j,1+icall) = AJ			! A_e or A_i | P_k
	   W(j,3+icall) = BJ			! B_e or B_i | Q_k+P_{k-1}+C_k
	   W(j,5+icall) = CJ			! C_e or C_i | Q_{k-1}
	   W(j,7+icall) = DJ			! D_e or D_i | 
	   W(j,23+icall) = C(j)
	   W(j,21+icall) = D(j)+B(j)/(Y625*HJ*VN(j))	! PDE or PDI
 2	continue
	if (icall .eq. 0)	then
	   icall = 1
	   return
	else
	   icall = 0
	endif
C  2nd call:
C     Work array usage:		W(1:N,1<->24)
C	 Exchange with RUNTT:	W(1:N,1<->8)
C	 Exchange with NURTT:	W(1:N,10<->24)
C  W(1:N,1) - A_e	W(1:N,2) - A_i
C  W(1:N,3) - B_e	W(1:N,4) - B_i
C  W(1:N,5) - C_e	W(1:N,6) - C_i
C  W(1:N,7) - D_e	W(1:N,8) - D_i
C  W(1:N,9) - not used
C  W(1:N,10)  - Pe->i
C  W(1:N,11) - E_11	W(1:N,12) - E_12
C  W(1:N,13) - E_21	W(1:N,14) - E_22
C  W(1:N,15) - G_1	W(1:N,16) - G_2
C  W(1:N,17),	W(1:N,18) - Coefficients for flux evaluation
C  W(1:N,19),	W(1:N,20) -   (used in NURTT to compute QE, QI)
C  W(1:N,21),	W(1:N,22) - Coefficients for RHS evaluation
C  W(1:N,23),	W(1:N,24) -   (used in NURTT to compute PETOT, PITOT)
	do	3	j=1,N
	   Y1 = W(j,3)
	   Y2 = W(j,4)
	   YA = HJ*VN(j)*Y625*W(j,10)
	   YB = YA
	   if (j .gt. 1)	then
	      Y1 = Y1-W(j,5)*W(j-1,11)	! Direct = (B_k - Q_{k-1}*E_{k-1})
	      Y2 = Y2-W(j,6)*W(j-1,14)	!          (Y1 YA)
	      YA = YA-W(j,5)*W(j-1,12)	!          (YB Y2)
	      YB = YB-W(j,6)*W(j-1,13)
	   endif
	   YJ = 1./(Y1*Y2-YA*YB)
	   Y11 = Y2*YJ			!           ( Y2 -YA)      (Y11 Y12)
	   Y12 =-YA*YJ			! Inversed =         /det =
	   Y21 =-YB*YJ			!           (-YB  Y1)      (Y21 Y22)
	   Y22 = Y1*YJ
	   W(j,11) = Y11*W(j,1)		! 
	   W(j,12) = Y12*W(j,2)		!            (-YA  Y1) / det
	   W(j,13) = Y21*W(j,1)		!      (W) = ( 11 12 )
	   W(j,14) = Y22*W(j,2)		!          = ( 13 22 )
	   Y1 = W(j,7)
	   Y2 = W(j,8)
	   if (j .gt. 1)	then
	      Y1 = Y1+W(j,5)*W(j-1,15)
	      Y2 = Y2+W(j,6)*W(j-1,16)
	   endif
	   W(j,15) = Y1*Y11+Y2*Y12
	   W(j,16) = Y1*Y21+Y2*Y22
 3	continue
	return
 98	write(*,'(2A)')
     &	   ">>> ERROR >>> The same boundary is required",
     &     " for both TE and TI"
	call	IFKEY(ichar(' '))
 99	write(*,'(2A,1F10.6)')
     &	   " >>> ERROR >>> Heat conductivity shouldn't be negative.",
     &	   "  RHO =",j*H
	call	IFKEY(ichar(' '))
	end
C======================================================================|
	subroutine NURTT(T1,T2,Q1,Q2,P1,P2,N,W)
C----------------------------------------------------------------------|
C TE,TI or TI,TE in the same order as by calling RUNTT
C
C In:  Boundary conditions -
C	 Q1(4) > 0   ->   given	T1(N+1)
C	 Q2(4) > 0   ->   given	T2(N+1)
C	 Q1(4) < 0   ->   given	Q1(N+1) = Q1(2)+T{e,i}(N+1)*Q1(3)
C	 Q2(4) < 0   ->   given	Q2(N+1) = Q2(2)+T{e,i}(N+1)*Q2(3)
C
C Out: Fluxes
C        Q_j(1:N+1) = -0.0016*(p_j*T_{j,k+1}-q_j*T_{j,N})
C      RHSs
C        P_j(1:N+1) = 
C----------------------------------------------------------------------|
	implicit none
	double precision T1(*),T2(*),Q1(*),Q2(*),P1(*),P2(*)
	integer	N,j
	double precision W(N,*)
	double precision Y1,Y2,Y11,Y12,Y21,Y22,YD
	if (Q1(4).lt.0. .and. Q2(4).lt.0.)	then
C Both eqns use fluxes as boundary conditions:
	   Y11 = W(N,19)*W(N,11)-W(N,17)-625.*Q1(3)
	   Y22 = W(N,20)*W(N,14)-W(N,18)-625.*Q2(3)
	   Y12 = W(N,19)*W(N,12)
	   Y21 = W(N,20)*W(N,13)
	   YD = Y11*Y22-Y12*Y21
	   Y1 = 625.*Q1(2)-W(N,19)*W(N,15)
	   Y2 = 625.*Q2(2)-W(N,20)*W(N,16)
	   T1(N+1) = (Y11*Y1-Y12*Y2)/YD
	   T2(N+1) = (Y22*Y2-Y21*Y1)/YD
	elseif (Q1(4) .lt. 0.)	then
C Mixed boundary conditions: 1st eqn flux, 2nd eqn temperature
	   Y11 = W(N,19)*W(N,11)-W(N,17)-625.*Q1(3)
	   Y1  = 625.*Q1(2)-W(N,19)*(W(N,15)+W(N,12)*T2(N+1))
	   T1(N+1) = Y1/Y11
	elseif (Q2(4) .lt. 0.)	then
C Mixed boundary conditions: 2nd eqn flux, 1st eqn temperature
	   Y22 = W(N,20)*W(N,14)-W(N,18)-625.*Q2(3)
	   Y2  = 625.*Q2(2)-W(N,20)*(W(N,16)+W(N,13)*T1(N+1))
	   T2(N+1) = Y2/Y22
	endif
C Define new temperatures:
	do	J=N,1,-1
	   T1(J) = W(j,11)*T1(J+1)+W(j,12)*T2(J+1)+W(j,15)
	   T2(J) = W(j,13)*T1(J+1)+W(j,14)*T2(J+1)+W(j,16)
	enddo
C Define fluxes:
	do	J=1,N
	   Q1(j) = -0.0016*(W(j,17)*T1(j+1)-W(j,18)*T1(j))
	   Q2(j) = -0.0016*(W(j,19)*T2(j+1)-W(j,20)*T2(j))
	enddo
	Q1(N+1) = Q1(N)
	Q2(N+1) = Q2(N)
C Define RHSs:
	do j=1,N+1
	   P1(j) = W(j,21)+W(j,23)*T1(j)+W(j,10)*(T1(j)-T2(j))
	   P2(j) = W(j,22)+W(j,24)*T2(j)+W(j,10)*(T2(j)-T1(j))
	enddo
	P1(N+1) = P1(N)
	P2(N+1) = P2(N)
	do	J=1,N+1
	   T1(j) = max(T1(j),1.d-4)
	   T2(j) = max(T2(j),1.d-4)
	enddo
	end
C======================================================================|
	subroutine RUNTT1(A,B,C,D,NO,NN,TO,N,GT,H,HB,TN,VO,VN,G11,W)
C----------------------------------------------------------------------|
C Old version of RUNTT without the stiff-transport-correction DV
C----------------------------------------------------------------------|
C	call RUNTT1 (YWA,YWB,PET,PE,NEO,NE,
C       		TEO,ND,TAU,HRO,YWC,TE,VRO,VR,G11,WORK1)
C----------------------------------------------------------------------|
C Exponential scheme
C----------------------------------------------------------------------|
C	The subroutine makes time step in the matrix equation:
C		d(N*T)/dt=1/V'*d[V'*(A*dT/dr+B*T)]/dr+625.*(C*T+D)
C
C	Here	A,B,C,D		- arrays(1:N)
C		NOld,NNew 	- arrays(1:N+1)		(densities)
C		TOld,TNew 	- arrays(1:N+1)		(temperatures)
C		VOld,VNew 	- arrays(1:N+1)		(dV/drho)
C		W(*)		- work space (e.g., work1(*))
C	Input:	H	- radial step (m)
C		N	- number of mesh points
C		GT	- time step (sec)
C       	A(1:N)  - diffusivity n_e*\chi_e (or n_i*\chi_i)
C       	B(1:N)  - convective velocity e.g. 5/2*GNX(J)*SLAT(J)/G11(J)
C       	C(1:N)  - PET or PIT without equipartition
C       	D(1:N)  - PE  or PI
C       	NO(1:N) - old density (previous time step)
C       	NN(1:N) - new density (next time step)
C       	VO(1:N) - old V'
C       	VN(1:N) - new V'
C       	TO(1:N+1) - old T_e (or T_i) 
C       	TN(1:N+1) - not used
C		HB 	  - edge cell size
C----------------------------------------------------------------------|
	implicit none
	double precision
     1		A(*),B(*),C(*),D(*),NO(*),NN(*),TO(*),TN(*),
     2		VO(*),VN(*),G11(*),H,HB,GT,GT23,Y625,AJ,BJ,CJ,DJ,
     3		H1,HJ,YA,YB,Y1,Y2,YJ,P0,P1,Q0,Q1,G0,G1,Y11,Y12,Y21,Y22
	integer	N,j,j1,jn,icall,N0
	double precision GETPEI,W(N,*)
	external	 GETPEI
	save	icall,N0
	data	icall/0/
        call add2loc("Subroutine RUNTT1"//char(0))
	GT23 = GT/1.5
	Y625 = 625.*GT23
	HJ = H
	G1 = 0.
	P1 = 0.
	Q1 = 0.
	H1 = H
	do	2	J=1,N
	   G0 = G1
	   P0 = P1
	   Q0 = Q1
	   G1 = GT23*G11(j)
	   YA = A(j)
	   YB = B(j)
	   if (YA .lt. 0.d0)	goto	99
	   P1 = 0.5*(abs(YB)+YB)		! 0.5(|B|+B)
	   if (YA .eq. 0.d0)	goto	1
	   if (j .eq. N)	then
C	      HJ = 0.5*(H+HB)
	      H1 = HB
	   endif
C Exponential scheme: (7 lines)
	   YJ = H1*YB/YA		   ! |\xi|
	   if (abs(YJ) .ge. 4.d1) goto	1  ! Use (A/h)*f(\xi) = .5*(|B|+B)
	   if (abs(YJ) .ge. 1.d-5) then
	      P1 = YB/(1.-exp(-YJ))	   ! Use (A/h)*f(\xi) = B/(1-exp{-\xi})
	   else
	      P1 = YA/H1*(1.+0.5*YJ)   	   ! Use (A/h)*f(\xi) = A/h/(1-\xi/2)
	   endif
C Power-law scheme: (6 lines)
C	   YJ = abs(H1*YB/YA)			! |\xi|
C	   if (YJ .ge. 1.d1)	goto	1	! Use (A/h)*f(\xi) = .5*(|B|+B)
C	   Y1 = 1.d0-.1d0*YJ			! (1 - 0.1|\xi|)
C	   Y2 = Y1*Y1				! (1 - 0.1|\xi|)^2
C	   YJ = Y2*Y2*Y1*YA/H1			! (A/h)*(1 - 0.1|\xi|)^5
C	   P1 = P1+YJ				! (A/h)*f(\xi) is done
 1	continue
	   Q1 = P1-YB				! (A/h)*g(\xi) = (A/h)*f(\xi)-B
	   if (icall .eq. 0)	then
	      N0 = N
	      W(j,17) = P1*G11(j)
	      W(j,18) = Q1*G11(j)
	      W(j,10) = -GETPEI(J)
	   else
	      if (N0.ne.N)	goto	98
C	      A(j) = P1
C	      B(j) = Q1
	      W(j,19) = P1*G11(j)
	      W(j,20) = Q1*G11(j)
	   endif
	   AJ = G1*P1
	   BJ = G1*Q1+G0*P0+HJ*VN(J)*(NN(J)-Y625*(W(j,10)+C(J)))
	   CJ = G0*Q0
	   YJ = (VO(J)/VN(J))**0.666667
	   DJ = HJ*(Y625*D(J)*VN(J)+NO(J)*TO(J)*VO(J)*YJ)
	   W(j,1+icall) = AJ			! A_e or A_i | P_k
	   W(j,3+icall) = BJ			! B_e or B_i | Q_k+P_{k-1}+C_k
	   W(j,5+icall) = CJ			! C_e or C_i | Q_{k-1}
	   W(j,7+icall) = DJ			! D_e or D_i | 
	   W(j,21+icall) = D(j)
	   W(j,23+icall) = C(j)
 2	continue
	if (icall .eq. 0)	then
	   icall = 1
	   return
	else
	   icall = 0
	endif
C Work array usage:
C	Internally:		W(1:N,1<->8)
C	Exchange with NURTT:	W(1:N,10<->24)
C  W(1:N,1) - A_e	W(1:N,2) - A_i
C  W(1:N,3) - B_e	W(1:N,4) - B_i
C  W(1:N,5) - C_e	W(1:N,6) - C_i
C  W(1:N,7) - D_e	W(1:N,8) - D_i
C  W(1:N,9) - free
C  W(1:N,10) - Pe->i
C  W(1:N,11) - E_11	W(1:N,12) - E_12
C  W(1:N,13) - E_21	W(1:N,14) - E_22
C  W(1:N,15) - G_1	W(1:N,16) - G_2
C  W(1:N,17),	W(1:N,18) - Coefficients for
C  W(1:N,19),	W(1:N,20) - flux evaluation
C  W(1:N,21),	W(1:N,22) - Coefficients for
C  W(1:N,23),	W(1:N,24) - RHS evaluation
	do	3	j=1,N
	   Y1 = W(j,3)
	   Y2 = W(j,4)
	   YA = HJ*VN(j)*Y625*W(j,10)
	   YB = YA
	   if (j .gt. 1)	then
	      Y1 = Y1-W(j,5)*W(j-1,11)	! Direct = (B_k - Q_{k-1}*E_{k-1})
	      Y2 = Y2-W(j,6)*W(j-1,14)	!          (Y1 YA)
	      YA = YA-W(j,5)*W(j-1,12)	!          (YB Y2)
	      YB = YB-W(j,6)*W(j-1,13)
	   endif
	   YJ = 1./(Y1*Y2-YA*YB)
	   Y11 = Y2*YJ			!           ( Y2 -YA)      (Y11 Y12)
	   Y12 =-YA*YJ			! Inversed =         /det =
	   Y21 =-YB*YJ			!           (-YB  Y1)      (Y21 Y22)
	   Y22 = Y1*YJ
	   W(j,11) = Y11*W(j,1)		! 
	   W(j,12) = Y12*W(j,2)		!            (-YA  Y1) / det
	   W(j,13) = Y21*W(j,1)		!      (W) = ( 11 12 )
	   W(j,14) = Y22*W(j,2)		!          = ( 13 22 )
	   Y1 = W(j,7)
	   Y2 = W(j,8)
	   if (j .gt. 1)	then
	      Y1 = Y1+W(j,5)*W(j-1,15)
	      Y2 = Y2+W(j,6)*W(j-1,16)
	   endif
	   W(j,15) = Y1*Y11+Y2*Y12
	   W(j,16) = Y1*Y21+Y2*Y22
 3	continue
	return
 98	write(*,'(2A)')
     &	   ">>> ERROR >>> Sorry, this mode requires the same boundary",
     &     " for both TE and TI"
	call	IFKEY(ichar(' '))
 99	write(*,'(2A,1F10.6)')
     &	   " >>> ERROR >>> Heat conductivity shouldn't be negative.",
     &	   "  RHO =",j*H
	call	IFKEY(ichar(' '))
	end
C======================================================================|
	subroutine RUNT1(A,B,C,D,G,NO,NN,TO,N,GT,H,E,TN,VOL,VNE,G11)
C----------------------------------------------------------------------|
C Linear scheme
C----------------------------------------------------------------------|
C	The subroutine provides run of the equation:
C		d(N*T)/dt=1/V'*d[V'*(A*dT/dr+B*T)]/dr+625.*(C*T+D)
C		V'(j)*(NN(j)*TN(j)-NO(j)*TO(j))*H*H/GT=
C			=V'(j+1/2)*A(j+1/2)*(TN(j+1)-TN(j))-
C			-V'(j-1/2)*A(j-1/2)*(TN(j)-TN(j-1))+
C			+0.5*V'(j+1/2)*B(j+1/2)*(TN(j+1)+TN(j))-
C			-0.5*V'(j-1/2)*B(j-1/2)*(TN(j)+TN(j-1))-
C			+V'(j)*H*H*625.*(C(j)*TN(j)+D(j))
C
C	Here	A,B,C,D		- arrays(1:N)
C		NOld,NNew 	- arrays(1:N+1)		(densities)
C		TOld,TNew 	- arrays(1:N+1)		(temperatures)
C		VOLd,VNEw 	- arrays(1:N+1)		(dV/drho)
C	Input:	H	- radial step (m)
C		N	- number of mesh points
C		GT	- time step (sec)
C		A(N),B(N),C(N),D(N),NO(N),NN(N),VOL(N),VNE(N),TO(N+1)
C		TN(N+1),E(1),E(2),E(3),E(4)
C	Output:	TN(N)
C----------------------------------------------------------------------|
	implicit none
	double precision
     1		A(*),B(*),C(*),D(*),G(*),NO(*),NN(*),TO(*),TN(*),E(4),
     2		VOL(*),VNE(*),G11(*),H,HB,GT,GT23,Y625,HH,H05,YG11,
     3          RA,RA0,RB,RB0,AJ,BJ,CJ,DJ,VO,VN
	integer	N,j
	HB = E(1)
	GT23 = GT/1.5
	Y625 = 625.*GT23
	RA = 0.
	RB = 0.
	HH = H*H
	H05 = 0.5*H
	do	1	J=1,N
		VO = VOL(J)**1.666667
		VN = VNE(J)**0.666667
		RA0 = RA
		RB0 = RB
		YG11 = G11(J)*GT23
		RA = A(J)*YG11
		if (j .eq. N)	then
			H05 = 0.5*HB
			HH  = H05*(H+HB)
			RA0 = RA0*HB/H
			RB0 = RB0*HB/H
		endif
		RB = B(J)*YG11*H05
		AJ = RA+RB
		BJ = RA+RA0-RB+RB0+HH*(NN(J)-Y625*C(J))*VNE(J)
C		AJ = RA
C		BJ = RA+RA0-RB+RB0+HH*(VNE(J)-C(J)*VNE(J)*GT)
C		if (j .eq. N)	then
C		   BJ = BJ-RB
C		else
C		   AJ = AJ+RB
C		endif
		CJ = RA0-RB0
		DJ = HH*(Y625*D(J)*VNE(J)+NO(J)*TO(J)*VO/VN)
		if(J .ne. 1)	then
			BJ = BJ-CJ*G(J-1)
			DJ = DJ+CJ*TN(J-1)
		endif
		G(J) = AJ/BJ
		TN(J) = DJ/BJ
1	continue
C Here E(2)=QjB, E(3)=QjTB, where j should be replaced by E or I
	if (E(4) .lt. 0.)	TN(N+1) = ((RA-RB)*TN(N)-Y625*HB*E(2))
     .				    /(RA+RB+G(N)*(RB-RA)+Y625*HB*E(3))
C	RB = B(N)*G11(N)*GT23*HB
C	if (E(4) .lt. 0.)	TN(N+1) = ((RA-RB)*TN(N)-Y625*HB*E(2))
C     .				    /(RA+G(N)*(RB-RA)+Y625*HB*E(3))
	do	2	J=N,1,-1
2	TN(J) = G(J)*TN(J+1)+TN(J)
	do	3	J=1,N+1
3	TN(J) = MAX(TN(J),1.d-4)
	end
C======================================================================|
	subroutine	
     >		RUNF(AK,B,C,D,FO,N,GT,H,HB,F,FV)
C----------------------------------------------------------------------|
C	The subroutine provides inversion of the matrix equation for FP
C       The boundary condition is supplied in the following form:
C       	F(N-1)*Psi(N+1)+F(N)*Psi(N)=F(N+1)
C           where F(N-1), F(N) and F(N+1) are input parameters.
C	Scheme:
C	Input:	H	- radial step (m)
C		HB	- edge radial step (m)
C		N+1	- number of grid points
C		GT	- time step (sec)
C		AK(N)	- G22
C		B(N)	- conductivity
C		D(N)	- external (+bootstrap) current
C		FO(N)	- old poloidal flux
C		F(N-1),F(N),F(N+1)	- edge conditions
C	Output:	F(1:N+1)- new poloidal flux
C		C(1:N+1)- (1/rho)d{K*dF/d(rho)}/d(rho) ~ current density
C		B(1:N+1)- (1/rho)dF/d(rho)	~ rotational transform
C		D(1:N+1)- dF/dt	toroidal loop voltage
C----------------------------------------------------------------------|
	implicit none
	integer	N,j
	double precision
     1		AK(N+1),B(N+1),C(N+1),D(N+1),FV(N+1),F(N+1),FO(N+1),
     2		H,HB,GT,HH,AJ,BJ,CJ,DJ,RJ,RJHH,YHB
	HH = H*H
	AJ = 0.
	RJ = -0.5*H
	do	1	J=1,N
	   CJ = AJ
	   RJ = RJ+H
	   RJHH = RJ*HH
	   if (j .eq. N)	then
	      CJ = CJ*HB/H
	      RJHH = RJ*HB*H
	   endif
	   AJ = AK(J)
	   DJ = RJHH*B(J)/GT
	   BJ = AJ+CJ+DJ
	   DJ = DJ*FO(J)+RJHH*D(J)-AJ*(FV(J+1)-FV(J))
	   if(J .ne. 1)	then
	      BJ = BJ-CJ*C(J-1)
	      DJ = DJ+CJ*(D(J-1)-FV(J-1)+FV(J))
	   endif
	   C(J) = AJ/BJ
	   D(J) = DJ/BJ
 1	continue
	F(N+1) = (F(N+1)-F(N)*D(N))/(F(N-1)+F(N)*C(N))
	do	2	J=N,1,-1
	F(J) = C(J)*F(J+1)+D(J)
 2	continue
	do	j=1,N+1
	   D(j) = (F(j)-FO(j))/GT
C	   FO(j) = F(j)
	enddo
	end
C======================================================================|
	subroutine	
     >		RUNFOLD(AK,B,C,D,FO,N,GT,H,HB,F,FV)
C----------------------------------------------------------------------|
C	The subroutine provides inversion of the matrix equation for FP
C       The boundary condition is supplied in the following form:
C       	F(N-1)*Psi(N+1)+F(N)*Psi(N)=F(N+1)
C           where F(N-1), F(N) and F(N+1) are input parameters.
C	Scheme:
C	Input:	H	- radial step (m)
C		HB	- edge radial step (m)
C		N+1	- number of grid points
C		GT	- time step (sec)
C		AK(N)	- G22
C		B(N)	- conductivity
C		D(N)	- external (+bootstrap) current
C		FO(N)	- old poloidal flux
C		F(N-1),F(N),F(N+1)	- edge conditions
C	Output:	F(1:N+1)- new poloidal flux
C		C(1:N+1)- (1/rho)d{K*dF/d(rho)}/d(rho) ~ current density
C		B(1:N+1)- (1/rho)dF/d(rho)	~ rotational transform
C		D(1:N+1)- dF/dt	toroidal loop voltage
C----------------------------------------------------------------------|
	implicit none
	integer	N,j
	double precision
     1		AK(N+1),B(N+1),C(N+1),D(N+1),FV(N+1),F(N+1),FO(N+1),
     2		H,HB,GT,HH,AJ,BJ,CJ,DJ,RJ,RJHH,YHB,YDN,YDN1
	YDN = D(N)
	YDN1 = D(N-1)
	HH = H*H
	AJ = 0.
	RJ = -0.5*H
	do	1	J=1,N
	   CJ = AJ
	   RJ = RJ+H
	   RJHH = RJ*HH
	   if (j .eq. N)	then
	      CJ = CJ*HB/H
	      RJHH = RJ*0.5*HB*(H+HB)
	   endif
	   AJ = AK(J)
	   DJ = RJHH*B(J)/GT
	   BJ = AJ+CJ+DJ
	   DJ = DJ*FO(J)+RJHH*D(J)-AJ*(FV(J+1)-FV(J))
	   if(J .ne. 1)	then
	      BJ = BJ-CJ*C(J-1)
	      DJ = DJ+CJ*(D(J-1)-FV(J-1)+FV(J))
	   endif
	   C(J) = AJ/BJ
	   D(J) = DJ/BJ
 1	continue
	F(N+1) = (F(N+1)-F(N)*D(N))/(F(N-1)+F(N)*C(N))
	do	2	J=N,1,-1
	F(J) = C(J)*F(J+1)+D(J)
 2	continue
	goto	4
	AJ = 0.
	YHB = 0.5*HB/H
	do	3	J=1,N
		CJ = AJ
		if (j .lt. N)	then
		   AJ = (F(j+1)-F(j))/HH
		   B(j) = AJ/j
		   AJ = AK(j)*AJ
		   C(j) = (AJ-CJ)/H
		else
		   AJ = (F(j+1)-F(j))/H/HB
C		   B(j) = AJ/(j-0.5+YHB)
		   B(j) = AJ/j
		   AJ = AK(j)*AJ
		   C(j) = 2.*(AJ-CJ)/(H+HB)
C		   C(j) = C(N-1)
		endif
		C(j) = C(j)/(j-0.5)
 3	continue
	B(N+1) = B(N)*N/(N-0.5+HB/H)
C j_NA and j_NA1 are redifined in equftn
C	C(N+1) = C(N)+YDN+((C(N)+YDN)-(C(N-1)+YDN1))/H*HB
 4	continue
	do	j=1,N+1
	   D(j) = (F(j)-FO(j))/GT
C	   FO(j) = F(j)
	enddo
	end
C======================================================================|
	logical	function	ISTORY(NCH,CTYPE,ITIMES,TTOUT,TOUT)
C----------------------------------------------------------------------|
C Output to the file astory.dat
C----------------------------------------------------------------------|
C Input:
C	NCH - Logical units NCH and NCH+1 are used by the subroutine
C	CTYPE 
C	    = 'INITAL'
C	    = 'TIMOUT'
C	    = 'RADOUT'
C Output:
C	ISTORY  = .TRUE.  - Normal exit
C		= .FALSE. - Error. Blocks further use of the subroutine
C----------------------------------------------------------------------|
	implicit none
	character*(*) CTYPE
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
C	include	'for/timeoutput.inc'	! Can replace (ITIMES,TTOUT,TOUT)
	logical	EXITIM,EXIRAD
	integer	length,ITIMES
	integer	JSTMAX,JSRMAX,IERR,NCH,ITYPE,j,j1,JST,JSR,JCHT,JCHR
	parameter (JSTMAX = 3*NCONST,JSRMAX = NCONST)
	integer*2 IFALR(JSRMAX),IFALT(JSTMAX)
	character FILNAM*40,STRI*80
	double precision YARR(JSTMAX+NRW),TTOUT(ITIMES),TOUT(ITIMES,NRW)
	save	IFALT,IFALR,ITYPE,JST,JSR,JCHT,JCHR
	data	IFALT/JSTMAX*0/ IFALR/JSRMAX*0/
	data	JST/0/ JSR/0/ JCHT/0/ JCHR/0/
C----------------------------------------------------------------------|
C	write(*,*)CTYPE(1:6)
	if (CTYPE(1:6) .eq. 'TIMOUT')	goto	50
	if (CTYPE(1:6) .eq. 'RADOUT')	goto	60
	if (CTYPE(1:6) .ne. 'INITAL')	goto	94

Collect information:
	ITYPE = 0					! Formatted record
	inquire(file='tmp/astory.tim',exist=EXITIM)
	inquire(file='tmp/astory.rad',exist=EXIRAD)
	ISTORY = EXITIM .or. EXIRAD
C----------------------------------------------------------------------|
	do	j=1,NTOUT
	   if (MARKT(j) .eq. 1)	JCHT = JCHT+1
	enddo
	if (.not.EXITIM)	goto	6
	FILNAM = 'tmp/astory.tim'//char(0)
	call	OPENRD(NCH+1,FILNAM,ITYPE,IERR)
	if (IERR .ne. 0)	goto	2
 1	continue
	read(NCH+1,'(1A80)',END=2,ERR=2)STRI
	if (STRI(1:1) .eq. ' ')	goto	1
	j1 = index(STRI(1:),' ')
	if (j1 .gt. 7)	then
	   write(*,*)
	   write(*,*)'>>> Illegal record in the file "tmp/astory.tim"'
	   write(*,*)'    The line "', STRI(1:j1-1),'" is not processed' 
	   goto	1
	endif
	JST = JST+1
	if (JST .gt. JSTMAX)	goto	95
	IFALT(JST) = 0
C	write(*,*) STRI(1:j1-1)
C	write(*,*)(PRNAME(j),j=1,NPRNAM)
	do	j=1,NPRNAM
	   if (PRNAME(j) .eq. STRI(1:j1-1))	IFALT(JST) = j
C	   if (PRNAME(j) .eq. STRI(1:j1-1))Write(*,*)jst,j,PRNAME(j)
	enddo
C	write(*,*)(SRNAME(j),j=1,NINVAR)
	do	j=1,NINVAR
	   if (SRNAME(j) .eq. STRI(1:j1-1))	IFALT(JST) = j+NCONST
C	   if (SRNAME(j) .eq. STRI(1:j1-1))Write(*,*)jst,j,SRNAME(j)
	enddo
C	write(*,*)(CFNAME(j),j=1,NCFNAM)
	do	j=1,NCFNAM
	   if (CFNAME(j) .eq. STRI(1:j1-1))	IFALT(JST) = j+2*NCONST
C	   if (CFNAME(j) .eq. STRI(1:j1-1))Write(*,*)jst,j,CFNAME(j)
	enddo
	if (IFALT(JST) .eq. 0)	then
	   write(*,*)
	   write(*,*)'>>> Illegal request in the file "tmp/astory.tim"'
	   write(*,*)'    Variable "', STRI(1:j1-1),'" does not exist' 
	   JST = JST-1
	endif
	goto	1
 2	continue
	close(NCH+1)
C	if (j2 .ne. JSTMAX)	goto	95
C----------------------------------------------------------------------|
 6	continue
	do	j=1,NROUT
	   if (MARKR(j) .eq. 1)	JCHR = JCHR+1
	enddo
	if (.not.EXIRAD)	goto	9
	FILNAM = 'tmp/astory.rad'//char(0)
	call	OPENRD(NCH+1,FILNAM,ITYPE,IERR)
	if (IERR .ne. 0)	goto	8
 7	continue
	read(NCH+1,'(1A80)',END=8,ERR=8)STRI
	if (STRI(1:1) .eq. ' ')	goto	7
	j1 = index(STRI(1:),' ')
	if (j1 .gt. 7)	then
	   write(*,*)
	   write(*,*)'>>> Illegal record in the file "tmp/astory.rad"'
	   write(*,*)'    The line "', STRI(1:j1-1),'" is not processed' 
	   goto	7
	endif
	JSR = JSR+1			! Total number of rad. channels
	if (JSR .gt. JSRMAX)	goto	95
	IFALR(JSR) = 0
	do	j=1,NARNAM
	   if (ARNAME(j) .eq. STRI(1:j1-1))	IFALR(JSR) = j
	enddo
	if (IFALR(JSR) .eq. 0)	then
	   write(*,*)
	   write(*,*)'>>> Illegal request in the file "tmp/astory.rad"'
	   write(*,*)'    Array "', STRI(1:j1-1),'" does not exist' 
	   JSR = JSR-1
	endif
	goto	7
 8	continue
	close(NCH+1)
C	write(*,*)'Number of arrays to save',JSR
C     &,		',  Total number',NARNAM
C     &,		',  Number of channels to save',JCHR
C	do	j=1,JSR
C	   if (IFALR(j) .ne. 0) write(*,'(A)')ARNAME(IFALR(j))
C	enddo
C	write(*,*)jsr,'   "',FILNAM(1:length(FILNAM)),'"'
C	write(*,*)'ISTORY called for the 1st time ',TIME,JST+JCHT,JSR+JCHR
C----------------------------------------------------------------------|
 9	continue
	if (JST+JCHT+JSR+JCHR.eq.0)	goto	99

C Start writing
	write(STRI,'(80A1)')(' ',j=1,80)		! Empty string
	FILNAM = 'tmp/astory.dat'//char(0)
	call	OPENWT(NCH,FILNAM,ITYPE,IERR)
	if (IERR .eq. 2)	goto	91
C	if (length(SHOTID) .gt. 0)	then
C	   write(NCH,*)SHOTID
C	else
	   write(NCH,80)STRI
C	   write(*,80)STRI
C	endif
 80	format(1A80)
	write(NCH,'(1A32)')VERSION
	write(NCH,80)RUNID
	write(NCH,80)XLINE1(1:80)
	write(NCH,80)XLINE2(1:80)
	write(NCH,*)JST+JCHT,JSR+JCHR
C	write(*,'(1A32)')VERSION
C	write(*,80)RUNID
C	write(*,80)XLINE1(1:80)
C	write(*,80)XLINE2(1:80)
C	write(*,*)JST+JCHT,JSR+JCHR

	if (.not.EXITIM .or. JST.eq.0)	goto	12
	FILNAM = 'tmp/astory.tim'//char(0)
	call	OPENRD(NCH+1,FILNAM,ITYPE,IERR)
	if (IERR .ne. 0)		goto	12
	j1 = 0
 11	continue
	read(NCH+1,'(1A80)',END=12,ERR=12)STRI
	if (STRI(1:1) .eq. ' ')		goto	11
	if (index(STRI(1:),' ') .gt. 7)	goto	11
	if (IFALT(j1+1) .eq. 0)		goto	11
	j1 = j1+1
	write(NCH,'(A)')STRI(1:)
	goto	11
 12	close(NCH+1)
	if (JCHT .eq. 0)	goto	13
	do	j=1,NTOUT
	   if (MARKT(j) .ne. 0)	then
	      write(NCH,'(A)')NAMET(j)
	   endif
	enddo
 13	continue

	if (.not.EXIRAD .or. JSR.eq.0)	goto	15
	FILNAM = 'tmp/astory.rad'//char(0)
	call	OPENRD(NCH+1,FILNAM,ITYPE,IERR)
	if (IERR .ne. 0)		goto	15
	j1 = 0
 14	continue
	read(NCH+1,'(1A80)',END=15,ERR=15)STRI
	if (STRI(1:1) .eq. ' ')		goto	14
	if (index(STRI(1:),' ') .gt. 7)	goto	14
	if (IFALR(j1+1) .eq. 0)		goto	14
	j1 = j1+1
	write(NCH,'(A)')STRI(1:)
C	write(*,'(A)')STRI(1:)
	goto	14
 15	close(NCH+1)
	if (JCHR .eq. 0)	goto	17
	do	j=1,NROUT
	   if (MARKR(j) .ne. 0)	then
	      write(NCH,'(A)')NAMER(j)
	   endif
	enddo
 17	continue
	goto	51

 50	continue
	if (JST+JCHT .eq. 0)	goto	52
	FILNAM = 'tmp/astory.dat'//char(0)
	call	OPENAP(NCH,FILNAM,ITYPE,IERR)
	if (IERR .ne. 0)	goto	93
 51	continue
	if (JST+JCHT .eq. 0)	goto	52
C	write(*,*)'ISTORY called for a time record',TIME
	write(NCH,*)TTOUT(LTOUT-1),0
	j1 = 0
	if (JST .eq. 0)		goto	52
	do	j=1,JST
	   if (IFALT(j) .ne. 0)	then
	      j1 = j1+1
	      call	APPST(IFALT(j),YARR(j1))
	   endif
	enddo
	if (JCHT .eq. 0)	then
C	   write(*,*)J1,' =',JST
	   write(NCH,100)(YARR(j),j=1,JST)
	   goto	61
	endif
 52	continue
	if (JCHT .eq. 0)	goto	61
	do	j=1,NTOUT
	   if (MARKT(j) .ne. 0)	then
	      j1 = j1+1
	      YARR(j1) = TOUT(LTOUT-1,j)
	   endif
	enddo
C	write(*,*)J1,' =',JST+JCHT
	write(NCH,100)(YARR(j),j=1,JST+JCHT)
	if (CTYPE(1:6).eq.'INITAL' .and. JSR+JCHR.gt.0)  goto	61
	goto	77

 60	continue
	if (JSR+JCHR .eq. 0)	goto	77
	FILNAM = 'tmp/astory.dat'//char(0)
	call	OPENAP(NCH,FILNAM,ITYPE,IERR)
	if (IERR .ne. 0)	goto	92
 61	continue
	if (JSR+JCHR .eq. 0)	goto	77
C	write(*,*)'ISTORY called for a radial record',TIME
	write(NCH,*)TTOUT(LTOUT-1),NA1
	if (JSR .eq. 0)		goto	62
	do	j=1,JSR
	   if (IFALR(j) .ne. 0)	call	APPSR(NCH,NA1,IFALR(j))
	enddo
 62	continue
	if (JCHR .eq. 0)	goto	77
	do	j=1,NROUT
	   if (MARKR(j) .ne. 0)	write(NCH,100)(ROUT(j1,j),j1=1,NA1)
	enddo

 77	ISTORY = .TRUE.
	close(NCH)
	close(NCH+1)
	return
 91	j = length(FILNAM)
	write(*,*)'Cannot open storage file "',FILNAM(1:j),'"'
	goto	99
 92	j = length(FILNAM)
	write(*,*)'Cannot write radial record to file "',FILNAM(1:j),'"'
	goto	99
 93	j = length(FILNAM)
	write(*,*)'Cannot write time record to file "',FILNAM(1:j),'"'
	goto	99
 94	write(*,*)'>>> ISTORY ERROR >>> Wrong call'
	goto	99
 95	j = length(FILNAM)
	write(*,*)'>>> Illegal changes in the file "',FILNAM(1:j),'"'
	goto	99
 99	ISTORY = .FALSE.
	close(NCH)
	close(NCH+1)
 100	format(1p,5e14.6)
	end
C======================================================================|
	subroutine	APPST(JFALT,YARR)
C-----------------------------------------------------------------------
C This subroutine is created automatically by the write sequence 
C    which has to be placed in the file for/init.f after reading
C    the file for/const.inc, line ~ 250.
!	do	j=1,NPRNAM
!	   write(*,*)'	if (JFALT .eq.',j,')  YARR = ',PRNAME(j)
!	enddo
C Exception for DTEQ(4,20)
!	do	j=1,NINVAR
!	   write(*,*)'	if (JFALT .eq.',j,'+NCONST) YARR = ',SRNAME(j)
!	enddo
!	do	j=1,NCFNAM
!	   write(*,*)'	if (JFALT .eq.',j,'+2*NCONST) YARR = ',CFNAME(j)
!	enddo
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	integer*2 JFALT
	double precision	YARR
	if (JFALT .eq. 1)  YARR = AB
	if (JFALT .eq. 2)  YARR = ABC
	if (JFALT .eq. 3)  YARR = AIM1
	if (JFALT .eq. 4)  YARR = AIM2
	if (JFALT .eq. 5)  YARR = AIM3
	if (JFALT .eq. 6)  YARR = AMJ
	if (JFALT .eq. 7)  YARR = AWALL
	if (JFALT .eq. 8)  YARR = BTOR
	if (JFALT .eq. 9)  YARR = ELONG
	if (JFALT .eq. 10)  YARR = ELONM
	if (JFALT .eq. 11)  YARR = ENCL
	if (JFALT .eq. 12)  YARR = ENWM
	if (JFALT .eq. 13)  YARR = FECR
	if (JFALT .eq. 14)  YARR = FFW
	if (JFALT .eq. 15)  YARR = FICR
	if (JFALT .eq. 16)  YARR = FLH
	if (JFALT .eq. 17)  YARR = GN2E
	if (JFALT .eq. 18)  YARR = GN2I
	if (JFALT .eq. 19)  YARR = IPL
	if (JFALT .eq. 20)  YARR = LEXT
	if (JFALT .eq. 21)  YARR = NNCL
	if (JFALT .eq. 22)  YARR = NNWM
	if (JFALT .eq. 23)  YARR = QECR
	if (JFALT .eq. 24)  YARR = QFW
	if (JFALT .eq. 25)  YARR = QICR
	if (JFALT .eq. 26)  YARR = QLH
	if (JFALT .eq. 27)  YARR = QNBI
	if (JFALT .eq. 28)  YARR = RTOR
	if (JFALT .eq. 29)  YARR = SHIFT
	if (JFALT .eq. 30)  YARR = TRIAN
	if (JFALT .eq. 31)  YARR = TRICH
	if (JFALT .eq. 32)  YARR = UEXT
	if (JFALT .eq. 33)  YARR = UPDWN
	if (JFALT .eq. 34)  YARR = WNE
	if (JFALT .eq. 35)  YARR = WTE
	if (JFALT .eq. 36)  YARR = WTI
	if (JFALT .eq. 37)  YARR = ZMJ
	if (JFALT .eq. 38)  YARR = ZRD1
	if (JFALT .eq. 39)  YARR = ZRD10
	if (JFALT .eq. 40)  YARR = ZRD11
	if (JFALT .eq. 41)  YARR = ZRD12
	if (JFALT .eq. 42)  YARR = ZRD13
	if (JFALT .eq. 43)  YARR = ZRD14
	if (JFALT .eq. 44)  YARR = ZRD15
	if (JFALT .eq. 45)  YARR = ZRD16
	if (JFALT .eq. 46)  YARR = ZRD17
	if (JFALT .eq. 47)  YARR = ZRD18
	if (JFALT .eq. 48)  YARR = ZRD19
	if (JFALT .eq. 49)  YARR = ZRD2
	if (JFALT .eq. 50)  YARR = ZRD20
	if (JFALT .eq. 51)  YARR = ZRD21
	if (JFALT .eq. 52)  YARR = ZRD22
	if (JFALT .eq. 53)  YARR = ZRD23
	if (JFALT .eq. 54)  YARR = ZRD24
	if (JFALT .eq. 55)  YARR = ZRD25
	if (JFALT .eq. 56)  YARR = ZRD26
	if (JFALT .eq. 57)  YARR = ZRD27
	if (JFALT .eq. 58)  YARR = ZRD28
	if (JFALT .eq. 59)  YARR = ZRD29
	if (JFALT .eq. 60)  YARR = ZRD3
	if (JFALT .eq. 61)  YARR = ZRD30
	if (JFALT .eq. 62)  YARR = ZRD31
	if (JFALT .eq. 63)  YARR = ZRD32
	if (JFALT .eq. 64)  YARR = ZRD33
	if (JFALT .eq. 65)  YARR = ZRD34
	if (JFALT .eq. 66)  YARR = ZRD35
	if (JFALT .eq. 67)  YARR = ZRD36
	if (JFALT .eq. 68)  YARR = ZRD37
	if (JFALT .eq. 69)  YARR = ZRD38
	if (JFALT .eq. 70)  YARR = ZRD39
	if (JFALT .eq. 71)  YARR = ZRD4
	if (JFALT .eq. 72)  YARR = ZRD40
	if (JFALT .eq. 73)  YARR = ZRD41
	if (JFALT .eq. 74)  YARR = ZRD42
	if (JFALT .eq. 75)  YARR = ZRD43
	if (JFALT .eq. 76)  YARR = ZRD44
	if (JFALT .eq. 77)  YARR = ZRD45
	if (JFALT .eq. 78)  YARR = ZRD46
	if (JFALT .eq. 79)  YARR = ZRD47
	if (JFALT .eq. 80)  YARR = ZRD48
	if (JFALT .eq. 81)  YARR = ZRD5
	if (JFALT .eq. 82)  YARR = ZRD6
	if (JFALT .eq. 83)  YARR = ZRD7
	if (JFALT .eq. 84)  YARR = ZRD8
	if (JFALT .eq. 85)  YARR = ZRD9
	if (JFALT .eq. 1+NCONST) YARR = DROUT
	if (JFALT .eq. 2+NCONST) YARR = DTOUT
	if (JFALT .eq. 3+NCONST) YARR = DPOUT
	if (JFALT .eq. 4+NCONST) YARR = TIME
	if (JFALT .eq. 5+NCONST) YARR = TAUMIN
	if (JFALT .eq. 6+NCONST) YARR = TAUMAX
	if (JFALT .eq. 7+NCONST) YARR = TAUINC
	if (JFALT .eq. 8+NCONST) YARR = DELVAR
	if (JFALT .eq. 9+NCONST) YARR = ITEREX
	if (JFALT .eq. 10+NCONST) YARR = ITERIN
	if (JFALT .eq. 11+NCONST) YARR = TINIT
	if (JFALT .eq. 12+NCONST) YARR = TSCALE
	if (JFALT .eq. 13+NCONST) YARR = NB1
	if (JFALT .eq. 14+NCONST) YARR = NUF
	if (JFALT .eq. 15+NCONST) YARR = XOUT
	if (JFALT .eq. 16+NCONST) YARR = XINPUT
	if (JFALT .eq. 17+NCONST) YARR = NB2EQL
	if (JFALT .eq. 18+NCONST) YARR = NEQUIL
	if (JFALT .eq. 19+NCONST) YARR = NBND
	if (JFALT .eq. 20+NCONST) YARR = XFLAG
	if (JFALT .eq. 21+NCONST) YARR = DTEQL
	if (JFALT .eq. 22+NCONST) YARR = MEQUIL
	if (JFALT .eq. 23+NCONST) YARR = TPAUSE
	if (JFALT .eq. 24+NCONST) YARR = TEND
	if (JFALT .eq. 25+NCONST) YARR = TSTART
	if (JFALT .eq. 26+NCONST) YARR = TEND
	if (JFALT .eq. 27+NCONST) YARR = TAU
C	if (JFALT .eq. 28+NCONST) YARR = DTEQ
	if (JFALT .eq. 29+NCONST) YARR = HRO
	if (JFALT .eq. 30+NCONST) YARR = HROA
	if (JFALT .eq. 31+NCONST) YARR = ALBPL
	if (JFALT .eq. 32+NCONST) YARR = NNCX
	if (JFALT .eq. 33+NCONST) YARR = QETB
	if (JFALT .eq. 34+NCONST) YARR = QFF0B
	if (JFALT .eq. 35+NCONST) YARR = QFF1B
	if (JFALT .eq. 36+NCONST) YARR = QFF2B
	if (JFALT .eq. 37+NCONST) YARR = QFF3B
	if (JFALT .eq. 38+NCONST) YARR = QFF4B
	if (JFALT .eq. 39+NCONST) YARR = QFF5B
	if (JFALT .eq. 40+NCONST) YARR = QFF6B
	if (JFALT .eq. 41+NCONST) YARR = QFF7B
	if (JFALT .eq. 42+NCONST) YARR = QFF8B
	if (JFALT .eq. 43+NCONST) YARR = QFF9B
	if (JFALT .eq. 44+NCONST) YARR = QITB
	if (JFALT .eq. 45+NCONST) YARR = QNNB
	if (JFALT .eq. 46+NCONST) YARR = ROB
	if (JFALT .eq. 47+NCONST) YARR = ROWALL
	if (JFALT .eq. 48+NCONST) YARR = ROC
	if (JFALT .eq. 49+NCONST) YARR = RON
	if (JFALT .eq. 50+NCONST) YARR = ROE
	if (JFALT .eq. 51+NCONST) YARR = ROI
	if (JFALT .eq. 52+NCONST) YARR = RO0
	if (JFALT .eq. 53+NCONST) YARR = RO1
	if (JFALT .eq. 54+NCONST) YARR = RO2
	if (JFALT .eq. 55+NCONST) YARR = RO3
	if (JFALT .eq. 56+NCONST) YARR = RO4
	if (JFALT .eq. 57+NCONST) YARR = RO5
	if (JFALT .eq. 58+NCONST) YARR = RO6
	if (JFALT .eq. 59+NCONST) YARR = RO7
	if (JFALT .eq. 60+NCONST) YARR = RO8
	if (JFALT .eq. 61+NCONST) YARR = RO9
	if (JFALT .eq. 62+NCONST) YARR = VOLUME
	if (JFALT .eq. 63+NCONST) YARR = GP
	if (JFALT .eq. 64+NCONST) YARR = GP2
	if (JFALT .eq. 65+NCONST) YARR = NA
	if (JFALT .eq. 66+NCONST) YARR = NA1
	if (JFALT .eq. 67+NCONST) YARR = NA1N
	if (JFALT .eq. 68+NCONST) YARR = NA1E
	if (JFALT .eq. 69+NCONST) YARR = NA1I
	if (JFALT .eq. 70+NCONST) YARR = NA10
	if (JFALT .eq. 71+NCONST) YARR = NA11
	if (JFALT .eq. 72+NCONST) YARR = NA12
	if (JFALT .eq. 73+NCONST) YARR = NA13
	if (JFALT .eq. 74+NCONST) YARR = NA14
	if (JFALT .eq. 75+NCONST) YARR = NA15
	if (JFALT .eq. 76+NCONST) YARR = NA16
	if (JFALT .eq. 77+NCONST) YARR = NA17
	if (JFALT .eq. 78+NCONST) YARR = NA18
	if (JFALT .eq. 79+NCONST) YARR = NA19
	if (JFALT .eq. 80+NCONST) YARR = NAB
	if (JFALT .eq. 81+NCONST) YARR = QBEAM
	if (JFALT .eq. 1+2*NCONST) YARR = CF1
	if (JFALT .eq. 2+2*NCONST) YARR = CF2
	if (JFALT .eq. 3+2*NCONST) YARR = CF3
	if (JFALT .eq. 4+2*NCONST) YARR = CF4
	if (JFALT .eq. 5+2*NCONST) YARR = CF5
	if (JFALT .eq. 6+2*NCONST) YARR = CF6
	if (JFALT .eq. 7+2*NCONST) YARR = CF7
	if (JFALT .eq. 8+2*NCONST) YARR = CF8
	if (JFALT .eq. 9+2*NCONST) YARR = CF9
	if (JFALT .eq. 10+2*NCONST) YARR = CF10
	if (JFALT .eq. 11+2*NCONST) YARR = CF11
	if (JFALT .eq. 12+2*NCONST) YARR = CF12
	if (JFALT .eq. 13+2*NCONST) YARR = CF13
	if (JFALT .eq. 14+2*NCONST) YARR = CF14
	if (JFALT .eq. 15+2*NCONST) YARR = CF15
	if (JFALT .eq. 16+2*NCONST) YARR = CF16
	if (JFALT .eq. 17+2*NCONST) YARR = CV1
	if (JFALT .eq. 18+2*NCONST) YARR = CV2
	if (JFALT .eq. 19+2*NCONST) YARR = CV3
	if (JFALT .eq. 20+2*NCONST) YARR = CV4
	if (JFALT .eq. 21+2*NCONST) YARR = CV5
	if (JFALT .eq. 22+2*NCONST) YARR = CV6
	if (JFALT .eq. 23+2*NCONST) YARR = CV7
	if (JFALT .eq. 24+2*NCONST) YARR = CV8
	if (JFALT .eq. 25+2*NCONST) YARR = CV9
	if (JFALT .eq. 26+2*NCONST) YARR = CV10
	if (JFALT .eq. 27+2*NCONST) YARR = CV11
	if (JFALT .eq. 28+2*NCONST) YARR = CV12
	if (JFALT .eq. 29+2*NCONST) YARR = CV13
	if (JFALT .eq. 30+2*NCONST) YARR = CV14
	if (JFALT .eq. 31+2*NCONST) YARR = CV15
	if (JFALT .eq. 32+2*NCONST) YARR = CV16
	if (JFALT .eq. 33+2*NCONST) YARR = CHE1
	if (JFALT .eq. 34+2*NCONST) YARR = CHE2
	if (JFALT .eq. 35+2*NCONST) YARR = CHE3
	if (JFALT .eq. 36+2*NCONST) YARR = CHE4
	if (JFALT .eq. 37+2*NCONST) YARR = CHI1
	if (JFALT .eq. 38+2*NCONST) YARR = CHI2
	if (JFALT .eq. 39+2*NCONST) YARR = CHI3
	if (JFALT .eq. 40+2*NCONST) YARR = CHI4
	if (JFALT .eq. 41+2*NCONST) YARR = CNB1
	if (JFALT .eq. 42+2*NCONST) YARR = CNB2
	if (JFALT .eq. 43+2*NCONST) YARR = CNB3
	if (JFALT .eq. 44+2*NCONST) YARR = CNB4
	if (JFALT .eq. 45+2*NCONST) YARR = CNBI1
	if (JFALT .eq. 46+2*NCONST) YARR = CNBI2
	if (JFALT .eq. 47+2*NCONST) YARR = CNBI3
	if (JFALT .eq. 48+2*NCONST) YARR = CNBI4
	if (JFALT .eq. 49+2*NCONST) YARR = CCD1
	if (JFALT .eq. 50+2*NCONST) YARR = CCD2
	if (JFALT .eq. 51+2*NCONST) YARR = CCD3
	if (JFALT .eq. 52+2*NCONST) YARR = CCD4
	if (JFALT .eq. 53+2*NCONST) YARR = CRF1
	if (JFALT .eq. 54+2*NCONST) YARR = CRF2
	if (JFALT .eq. 55+2*NCONST) YARR = CRF3
	if (JFALT .eq. 56+2*NCONST) YARR = CRF4
	if (JFALT .eq. 57+2*NCONST) YARR = CNEUT1
	if (JFALT .eq. 58+2*NCONST) YARR = CNEUT2
	if (JFALT .eq. 59+2*NCONST) YARR = CNEUT3
	if (JFALT .eq. 60+2*NCONST) YARR = CNEUT4
	if (JFALT .eq. 61+2*NCONST) YARR = CPEL1
	if (JFALT .eq. 62+2*NCONST) YARR = CPEL2
	if (JFALT .eq. 63+2*NCONST) YARR = CPEL3
	if (JFALT .eq. 64+2*NCONST) YARR = CPEL4
	if (JFALT .eq. 65+2*NCONST) YARR = CBND1
	if (JFALT .eq. 66+2*NCONST) YARR = CBND2
	if (JFALT .eq. 67+2*NCONST) YARR = CBND3
	if (JFALT .eq. 68+2*NCONST) YARR = CBND4
	if (JFALT .eq. 69+2*NCONST) YARR = CFUS1
	if (JFALT .eq. 70+2*NCONST) YARR = CFUS2
	if (JFALT .eq. 71+2*NCONST) YARR = CFUS3
	if (JFALT .eq. 72+2*NCONST) YARR = CFUS4
	if (JFALT .eq. 73+2*NCONST) YARR = CIMP1
	if (JFALT .eq. 74+2*NCONST) YARR = CIMP2
	if (JFALT .eq. 75+2*NCONST) YARR = CIMP3
	if (JFALT .eq. 76+2*NCONST) YARR = CIMP4
	if (JFALT .eq. 77+2*NCONST) YARR = CMHD1
	if (JFALT .eq. 78+2*NCONST) YARR = CMHD2
	if (JFALT .eq. 79+2*NCONST) YARR = CMHD3
	if (JFALT .eq. 80+2*NCONST) YARR = CMHD4
	if (JFALT .eq. 81+2*NCONST) YARR = CRAD1
	if (JFALT .eq. 82+2*NCONST) YARR = CRAD2
	if (JFALT .eq. 83+2*NCONST) YARR = CRAD3
	if (JFALT .eq. 84+2*NCONST) YARR = CRAD4
	if (JFALT .eq. 85+2*NCONST) YARR = CSOL1
	if (JFALT .eq. 86+2*NCONST) YARR = CSOL2
	if (JFALT .eq. 87+2*NCONST) YARR = CSOL3
	if (JFALT .eq. 88+2*NCONST) YARR = CSOL4
	end
C======================================================================|
	subroutine	APPSR(JCH,JA1,NUMALR)
C-----------------------------------------------------------------------
C This subroutine is created automatically by the write sequence 
C    which has to be placed in the file for/init.f after reading
C    the file for/status.inc, line ~ 250.
!	do	j=1,NARNAM
!	   write(*,*)'     if (NUMALR .eq.',j,
!     >		') write(JCH,*)(',ARNAME(j),"(j),j=1,JA1)"	
!	enddo
C-----------------------------------------------------------------------
	implicit none
	integer*2 NUMALR
	integer	JCH,JA1,j
	include	'for/parameter.inc'
	include	'for/status.inc'
 100	format(1p,5e14.6)
      if (NUMALR .eq.   1) write(JCH,*)(AMAIN (j),j=1,JA1)
      if (NUMALR .eq.   2) write(JCH,*)(AMETR (j),j=1,JA1)
      if (NUMALR .eq.   3) write(JCH,*)(B0DB2 (j),j=1,JA1)
      if (NUMALR .eq.   4) write(JCH,*)(BDB0  (j),j=1,JA1)
      if (NUMALR .eq.   5) write(JCH,*)(BDB02 (j),j=1,JA1)
      if (NUMALR .eq.   6) write(JCH,*)(BMAXT (j),j=1,JA1)
      if (NUMALR .eq.   7) write(JCH,*)(BMINT (j),j=1,JA1)
      if (NUMALR .eq.   8) write(JCH,*)(CAR1  (j),j=1,JA1)
      if (NUMALR .eq.   9) write(JCH,*)(CAR10 (j),j=1,JA1)
      if (NUMALR .eq.  10) write(JCH,*)(CAR10X(j),j=1,JA1)
      if (NUMALR .eq.  11) write(JCH,*)(CAR11 (j),j=1,JA1)
      if (NUMALR .eq.  12) write(JCH,*)(CAR11X(j),j=1,JA1)
      if (NUMALR .eq.  13) write(JCH,*)(CAR12 (j),j=1,JA1)
      if (NUMALR .eq.  14) write(JCH,*)(CAR12X(j),j=1,JA1)
      if (NUMALR .eq.  15) write(JCH,*)(CAR13 (j),j=1,JA1)
      if (NUMALR .eq.  16) write(JCH,*)(CAR13X(j),j=1,JA1)
      if (NUMALR .eq.  17) write(JCH,*)(CAR14 (j),j=1,JA1)
      if (NUMALR .eq.  18) write(JCH,*)(CAR14X(j),j=1,JA1)
      if (NUMALR .eq.  19) write(JCH,*)(CAR15 (j),j=1,JA1)
      if (NUMALR .eq.  20) write(JCH,*)(CAR15X(j),j=1,JA1)
      if (NUMALR .eq.  21) write(JCH,*)(CAR16 (j),j=1,JA1)
      if (NUMALR .eq.  22) write(JCH,*)(CAR16X(j),j=1,JA1)
      if (NUMALR .eq.  23) write(JCH,*)(CAR17 (j),j=1,JA1)
      if (NUMALR .eq.  24) write(JCH,*)(CAR17X(j),j=1,JA1)
      if (NUMALR .eq.  25) write(JCH,*)(CAR18 (j),j=1,JA1)
      if (NUMALR .eq.  26) write(JCH,*)(CAR18X(j),j=1,JA1)
      if (NUMALR .eq.  27) write(JCH,*)(CAR19 (j),j=1,JA1)
      if (NUMALR .eq.  28) write(JCH,*)(CAR19X(j),j=1,JA1)
      if (NUMALR .eq.  29) write(JCH,*)(CAR1X (j),j=1,JA1)
      if (NUMALR .eq.  30) write(JCH,*)(CAR2  (j),j=1,JA1)
      if (NUMALR .eq.  31) write(JCH,*)(CAR20 (j),j=1,JA1)
      if (NUMALR .eq.  32) write(JCH,*)(CAR20X(j),j=1,JA1)
      if (NUMALR .eq.  33) write(JCH,*)(CAR21 (j),j=1,JA1)
      if (NUMALR .eq.  34) write(JCH,*)(CAR21X(j),j=1,JA1)
      if (NUMALR .eq.  35) write(JCH,*)(CAR22 (j),j=1,JA1)
      if (NUMALR .eq.  36) write(JCH,*)(CAR22X(j),j=1,JA1)
      if (NUMALR .eq.  37) write(JCH,*)(CAR23 (j),j=1,JA1)
      if (NUMALR .eq.  38) write(JCH,*)(CAR23X(j),j=1,JA1)
      if (NUMALR .eq.  39) write(JCH,*)(CAR24 (j),j=1,JA1)
      if (NUMALR .eq.  40) write(JCH,*)(CAR24X(j),j=1,JA1)
      if (NUMALR .eq.  41) write(JCH,*)(CAR25 (j),j=1,JA1)
      if (NUMALR .eq.  42) write(JCH,*)(CAR25X(j),j=1,JA1)
      if (NUMALR .eq.  43) write(JCH,*)(CAR26 (j),j=1,JA1)
      if (NUMALR .eq.  44) write(JCH,*)(CAR26X(j),j=1,JA1)
      if (NUMALR .eq.  45) write(JCH,*)(CAR27 (j),j=1,JA1)
      if (NUMALR .eq.  46) write(JCH,*)(CAR27X(j),j=1,JA1)
      if (NUMALR .eq.  47) write(JCH,*)(CAR28 (j),j=1,JA1)
      if (NUMALR .eq.  48) write(JCH,*)(CAR28X(j),j=1,JA1)
      if (NUMALR .eq.  49) write(JCH,*)(CAR29 (j),j=1,JA1)
      if (NUMALR .eq.  50) write(JCH,*)(CAR29X(j),j=1,JA1)
      if (NUMALR .eq.  51) write(JCH,*)(CAR2X (j),j=1,JA1)
      if (NUMALR .eq.  52) write(JCH,*)(CAR3  (j),j=1,JA1)
      if (NUMALR .eq.  53) write(JCH,*)(CAR30 (j),j=1,JA1)
      if (NUMALR .eq.  54) write(JCH,*)(CAR30X(j),j=1,JA1)
      if (NUMALR .eq.  55) write(JCH,*)(CAR31 (j),j=1,JA1)
      if (NUMALR .eq.  56) write(JCH,*)(CAR31X(j),j=1,JA1)
      if (NUMALR .eq.  57) write(JCH,*)(CAR32 (j),j=1,JA1)
      if (NUMALR .eq.  58) write(JCH,*)(CAR32X(j),j=1,JA1)
      if (NUMALR .eq.  59) write(JCH,*)(CAR3X (j),j=1,JA1)
      if (NUMALR .eq.  60) write(JCH,*)(CAR4  (j),j=1,JA1)
      if (NUMALR .eq.  61) write(JCH,*)(CAR4X (j),j=1,JA1)
      if (NUMALR .eq.  62) write(JCH,*)(CAR5  (j),j=1,JA1)
      if (NUMALR .eq.  63) write(JCH,*)(CAR5X (j),j=1,JA1)
      if (NUMALR .eq.  64) write(JCH,*)(CAR6  (j),j=1,JA1)
      if (NUMALR .eq.  65) write(JCH,*)(CAR6X (j),j=1,JA1)
      if (NUMALR .eq.  66) write(JCH,*)(CAR7  (j),j=1,JA1)
      if (NUMALR .eq.  67) write(JCH,*)(CAR7X (j),j=1,JA1)
      if (NUMALR .eq.  68) write(JCH,*)(CAR8  (j),j=1,JA1)
      if (NUMALR .eq.  69) write(JCH,*)(CAR8X (j),j=1,JA1)
      if (NUMALR .eq.  70) write(JCH,*)(CAR9  (j),j=1,JA1)
      if (NUMALR .eq.  71) write(JCH,*)(CAR9X (j),j=1,JA1)
      if (NUMALR .eq.  72) write(JCH,*)(CU    (j),j=1,JA1)
      if (NUMALR .eq.  73) write(JCH,*)(CUBM  (j),j=1,JA1)
      if (NUMALR .eq.  74) write(JCH,*)(CUBS  (j),j=1,JA1)
      if (NUMALR .eq.  75) write(JCH,*)(CUECR (j),j=1,JA1)
      if (NUMALR .eq.  76) write(JCH,*)(CUFI  (j),j=1,JA1)
      if (NUMALR .eq.  77) write(JCH,*)(CUFW  (j),j=1,JA1)
      if (NUMALR .eq.  78) write(JCH,*)(CUICR (j),j=1,JA1)
      if (NUMALR .eq.  79) write(JCH,*)(CULH  (j),j=1,JA1)
      if (NUMALR .eq.  80) write(JCH,*)(CUTOR (j),j=1,JA1)
      if (NUMALR .eq.  81) write(JCH,*)(CUX   (j),j=1,JA1)
      if (NUMALR .eq.  82) write(JCH,*)(CV    (j),j=1,JA1)
      if (NUMALR .eq.  83) write(JCH,*)(DIMP1 (j),j=1,JA1)
      if (NUMALR .eq.  84) write(JCH,*)(DIMP2 (j),j=1,JA1)
      if (NUMALR .eq.  85) write(JCH,*)(DIMP3 (j),j=1,JA1)
      if (NUMALR .eq.  86) write(JCH,*)(DRODA (j),j=1,JA1)
      if (NUMALR .eq.  87) write(JCH,*)(DRODAX(j),j=1,JA1)
      if (NUMALR .eq.  88) write(JCH,*)(ELON  (j),j=1,JA1)
      if (NUMALR .eq.  89) write(JCH,*)(ELX   (j),j=1,JA1)
      if (NUMALR .eq.  90) write(JCH,*)(EQFF  (j),j=1,JA1)
      if (NUMALR .eq.  91) write(JCH,*)(EQPF  (j),j=1,JA1)
      if (NUMALR .eq.  92) write(JCH,*)(ER    (j),j=1,JA1)
      if (NUMALR .eq.  93) write(JCH,*)(F0    (j),j=1,JA1)
      if (NUMALR .eq.  94) write(JCH,*)(F0O   (j),j=1,JA1)
      if (NUMALR .eq.  95) write(JCH,*)(F0X   (j),j=1,JA1)
      if (NUMALR .eq.  96) write(JCH,*)(F1    (j),j=1,JA1)
      if (NUMALR .eq.  97) write(JCH,*)(F1O   (j),j=1,JA1)
      if (NUMALR .eq.  98) write(JCH,*)(F1X   (j),j=1,JA1)
      if (NUMALR .eq.  99) write(JCH,*)(F2    (j),j=1,JA1)
      if (NUMALR .eq. 100) write(JCH,*)(F2O   (j),j=1,JA1)
      if (NUMALR .eq. 101) write(JCH,*)(F2X   (j),j=1,JA1)
      if (NUMALR .eq. 102) write(JCH,*)(F3    (j),j=1,JA1)
      if (NUMALR .eq. 103) write(JCH,*)(F3O   (j),j=1,JA1)
      if (NUMALR .eq. 104) write(JCH,*)(F3X   (j),j=1,JA1)
      if (NUMALR .eq. 105) write(JCH,*)(F4    (j),j=1,JA1)
      if (NUMALR .eq. 106) write(JCH,*)(F4O   (j),j=1,JA1)
      if (NUMALR .eq. 107) write(JCH,*)(F4X   (j),j=1,JA1)
      if (NUMALR .eq. 108) write(JCH,*)(F5    (j),j=1,JA1)
      if (NUMALR .eq. 109) write(JCH,*)(F5O   (j),j=1,JA1)
      if (NUMALR .eq. 110) write(JCH,*)(F5X   (j),j=1,JA1)
      if (NUMALR .eq. 111) write(JCH,*)(F6    (j),j=1,JA1)
      if (NUMALR .eq. 112) write(JCH,*)(F6O   (j),j=1,JA1)
      if (NUMALR .eq. 113) write(JCH,*)(F6X   (j),j=1,JA1)
      if (NUMALR .eq. 114) write(JCH,*)(F7    (j),j=1,JA1)
      if (NUMALR .eq. 115) write(JCH,*)(F7O   (j),j=1,JA1)
      if (NUMALR .eq. 116) write(JCH,*)(F7X   (j),j=1,JA1)
      if (NUMALR .eq. 117) write(JCH,*)(F8    (j),j=1,JA1)
      if (NUMALR .eq. 118) write(JCH,*)(F8O   (j),j=1,JA1)
      if (NUMALR .eq. 119) write(JCH,*)(F8X   (j),j=1,JA1)
      if (NUMALR .eq. 120) write(JCH,*)(F9    (j),j=1,JA1)
      if (NUMALR .eq. 121) write(JCH,*)(F9O   (j),j=1,JA1)
      if (NUMALR .eq. 122) write(JCH,*)(F9X   (j),j=1,JA1)
      if (NUMALR .eq. 123) write(JCH,*)(FOFB  (j),j=1,JA1)
      if (NUMALR .eq. 124) write(JCH,*)(FP    (j),j=1,JA1)
      if (NUMALR .eq. 125) write(JCH,*)(FPO   (j),j=1,JA1)
      if (NUMALR .eq. 126) write(JCH,*)(FV    (j),j=1,JA1)
      if (NUMALR .eq. 127) write(JCH,*)(G11   (j),j=1,JA1)
      if (NUMALR .eq. 128) write(JCH,*)(G11X  (j),j=1,JA1)
      if (NUMALR .eq. 129) write(JCH,*)(G22   (j),j=1,JA1)
      if (NUMALR .eq. 130) write(JCH,*)(G22X  (j),j=1,JA1)
      if (NUMALR .eq. 131) write(JCH,*)(G33   (j),j=1,JA1)
      if (NUMALR .eq. 132) write(JCH,*)(G33X  (j),j=1,JA1)
      if (NUMALR .eq. 133) write(JCH,*)(GN    (j),j=1,JA1)
      if (NUMALR .eq. 134) write(JCH,*)(GNX   (j),j=1,JA1)
      if (NUMALR .eq. 135) write(JCH,*)(GRADRO(j),j=1,JA1)
      if (NUMALR .eq. 136) write(JCH,*)(IPOL  (j),j=1,JA1)
      if (NUMALR .eq. 137) write(JCH,*)(IPOLX (j),j=1,JA1)
      if (NUMALR .eq. 138) write(JCH,*)(MU    (j),j=1,JA1)
      if (NUMALR .eq. 139) write(JCH,*)(MUX   (j),j=1,JA1)
      if (NUMALR .eq. 140) write(JCH,*)(MV    (j),j=1,JA1)
      if (NUMALR .eq. 141) write(JCH,*)(MVX   (j),j=1,JA1)
      if (NUMALR .eq. 142) write(JCH,*)(NALF  (j),j=1,JA1)
      if (NUMALR .eq. 143) write(JCH,*)(NDEUT (j),j=1,JA1)
      if (NUMALR .eq. 144) write(JCH,*)(NE    (j),j=1,JA1)
      if (NUMALR .eq. 145) write(JCH,*)(NEO   (j),j=1,JA1)
      if (NUMALR .eq. 146) write(JCH,*)(NEX   (j),j=1,JA1)
      if (NUMALR .eq. 147) write(JCH,*)(NHE3  (j),j=1,JA1)
      if (NUMALR .eq. 148) write(JCH,*)(NHYDR (j),j=1,JA1)
      if (NUMALR .eq. 149) write(JCH,*)(NI    (j),j=1,JA1)
      if (NUMALR .eq. 150) write(JCH,*)(NIBM  (j),j=1,JA1)
      if (NUMALR .eq. 151) write(JCH,*)(NIO   (j),j=1,JA1)
      if (NUMALR .eq. 152) write(JCH,*)(NIX   (j),j=1,JA1)
      if (NUMALR .eq. 153) write(JCH,*)(NIZ1  (j),j=1,JA1)
      if (NUMALR .eq. 154) write(JCH,*)(NIZ2  (j),j=1,JA1)
      if (NUMALR .eq. 155) write(JCH,*)(NIZ3  (j),j=1,JA1)
      if (NUMALR .eq. 156) write(JCH,*)(NN    (j),j=1,JA1)
      if (NUMALR .eq. 157) write(JCH,*)(NNBM1 (j),j=1,JA1)
      if (NUMALR .eq. 158) write(JCH,*)(NNBM2 (j),j=1,JA1)
      if (NUMALR .eq. 159) write(JCH,*)(NNBM3 (j),j=1,JA1)
      if (NUMALR .eq. 160) write(JCH,*)(NTRIT (j),j=1,JA1)
      if (NUMALR .eq. 161) write(JCH,*)(PBEAM (j),j=1,JA1)
      if (NUMALR .eq. 162) write(JCH,*)(PBLON (j),j=1,JA1)
      if (NUMALR .eq. 163) write(JCH,*)(PBOL1 (j),j=1,JA1)
      if (NUMALR .eq. 164) write(JCH,*)(PBOL2 (j),j=1,JA1)
      if (NUMALR .eq. 165) write(JCH,*)(PBOL3 (j),j=1,JA1)
      if (NUMALR .eq. 166) write(JCH,*)(PBPER (j),j=1,JA1)
      if (NUMALR .eq. 167) write(JCH,*)(PEBM  (j),j=1,JA1)
      if (NUMALR .eq. 168) write(JCH,*)(PEECR (j),j=1,JA1)
      if (NUMALR .eq. 169) write(JCH,*)(PEFW  (j),j=1,JA1)
      if (NUMALR .eq. 170) write(JCH,*)(PEICR (j),j=1,JA1)
      if (NUMALR .eq. 171) write(JCH,*)(PELH  (j),j=1,JA1)
      if (NUMALR .eq. 172) write(JCH,*)(PELON (j),j=1,JA1)
      if (NUMALR .eq. 173) write(JCH,*)(PEPER (j),j=1,JA1)
      if (NUMALR .eq. 174) write(JCH,*)(PETOT (j),j=1,JA1)
      if (NUMALR .eq. 175) write(JCH,*)(PEX   (j),j=1,JA1)
      if (NUMALR .eq. 176) write(JCH,*)(PIBM  (j),j=1,JA1)
      if (NUMALR .eq. 177) write(JCH,*)(PIFW  (j),j=1,JA1)
      if (NUMALR .eq. 178) write(JCH,*)(PIICR (j),j=1,JA1)
      if (NUMALR .eq. 179) write(JCH,*)(PITOT (j),j=1,JA1)
      if (NUMALR .eq. 180) write(JCH,*)(PIX   (j),j=1,JA1)
      if (NUMALR .eq. 181) write(JCH,*)(PRAD  (j),j=1,JA1)
      if (NUMALR .eq. 182) write(JCH,*)(PRADX (j),j=1,JA1)
      if (NUMALR .eq. 183) write(JCH,*)(PSXR1 (j),j=1,JA1)
      if (NUMALR .eq. 184) write(JCH,*)(PSXR2 (j),j=1,JA1)
      if (NUMALR .eq. 185) write(JCH,*)(PSXR3 (j),j=1,JA1)
      if (NUMALR .eq. 186) write(JCH,*)(QE    (j),j=1,JA1)
      if (NUMALR .eq. 187) write(JCH,*)(QF0   (j),j=1,JA1)
      if (NUMALR .eq. 188) write(JCH,*)(QF1   (j),j=1,JA1)
      if (NUMALR .eq. 189) write(JCH,*)(QF2   (j),j=1,JA1)
      if (NUMALR .eq. 190) write(JCH,*)(QF3   (j),j=1,JA1)
      if (NUMALR .eq. 191) write(JCH,*)(QF4   (j),j=1,JA1)
      if (NUMALR .eq. 192) write(JCH,*)(QF5   (j),j=1,JA1)
      if (NUMALR .eq. 193) write(JCH,*)(QF6   (j),j=1,JA1)
      if (NUMALR .eq. 194) write(JCH,*)(QF7   (j),j=1,JA1)
      if (NUMALR .eq. 195) write(JCH,*)(QF8   (j),j=1,JA1)
      if (NUMALR .eq. 196) write(JCH,*)(QF9   (j),j=1,JA1)
      if (NUMALR .eq. 197) write(JCH,*)(QI    (j),j=1,JA1)
      if (NUMALR .eq. 198) write(JCH,*)(QN    (j),j=1,JA1)
      if (NUMALR .eq. 199) write(JCH,*)(RHO   (j),j=1,JA1)
      if (NUMALR .eq. 200) write(JCH,*)(SCUBM (j),j=1,JA1)
      if (NUMALR .eq. 201) write(JCH,*)(SF0TOT(j),j=1,JA1)
      if (NUMALR .eq. 202) write(JCH,*)(SF1TOT(j),j=1,JA1)
      if (NUMALR .eq. 203) write(JCH,*)(SF2TOT(j),j=1,JA1)
      if (NUMALR .eq. 204) write(JCH,*)(SF3TOT(j),j=1,JA1)
      if (NUMALR .eq. 205) write(JCH,*)(SF4TOT(j),j=1,JA1)
      if (NUMALR .eq. 206) write(JCH,*)(SF5TOT(j),j=1,JA1)
      if (NUMALR .eq. 207) write(JCH,*)(SF6TOT(j),j=1,JA1)
      if (NUMALR .eq. 208) write(JCH,*)(SF7TOT(j),j=1,JA1)
      if (NUMALR .eq. 209) write(JCH,*)(SF8TOT(j),j=1,JA1)
      if (NUMALR .eq. 210) write(JCH,*)(SF9TOT(j),j=1,JA1)
      if (NUMALR .eq. 211) write(JCH,*)(SHEAR (j),j=1,JA1)
      if (NUMALR .eq. 212) write(JCH,*)(SHIF  (j),j=1,JA1)
      if (NUMALR .eq. 213) write(JCH,*)(SHIV  (j),j=1,JA1)
      if (NUMALR .eq. 214) write(JCH,*)(SHX   (j),j=1,JA1)
      if (NUMALR .eq. 215) write(JCH,*)(SLAT  (j),j=1,JA1)
      if (NUMALR .eq. 216) write(JCH,*)(SLATX (j),j=1,JA1)
      if (NUMALR .eq. 217) write(JCH,*)(SNEBM (j),j=1,JA1)
      if (NUMALR .eq. 218) write(JCH,*)(SNIBM1(j),j=1,JA1)
      if (NUMALR .eq. 219) write(JCH,*)(SNIBM2(j),j=1,JA1)
      if (NUMALR .eq. 220) write(JCH,*)(SNIBM3(j),j=1,JA1)
      if (NUMALR .eq. 221) write(JCH,*)(SNNBM (j),j=1,JA1)
      if (NUMALR .eq. 222) write(JCH,*)(SNTOT (j),j=1,JA1)
      if (NUMALR .eq. 223) write(JCH,*)(SNX   (j),j=1,JA1)
      if (NUMALR .eq. 224) write(JCH,*)(SQEPS (j),j=1,JA1)
      if (NUMALR .eq. 225) write(JCH,*)(TE    (j),j=1,JA1)
      if (NUMALR .eq. 226) write(JCH,*)(TEO   (j),j=1,JA1)
      if (NUMALR .eq. 227) write(JCH,*)(TEX   (j),j=1,JA1)
      if (NUMALR .eq. 228) write(JCH,*)(TI    (j),j=1,JA1)
      if (NUMALR .eq. 229) write(JCH,*)(TIO   (j),j=1,JA1)
      if (NUMALR .eq. 230) write(JCH,*)(TIX   (j),j=1,JA1)
      if (NUMALR .eq. 231) write(JCH,*)(TN    (j),j=1,JA1)
      if (NUMALR .eq. 232) write(JCH,*)(TRIA  (j),j=1,JA1)
      if (NUMALR .eq. 233) write(JCH,*)(TRX   (j),j=1,JA1)
      if (NUMALR .eq. 234) write(JCH,*)(ULON  (j),j=1,JA1)
      if (NUMALR .eq. 235) write(JCH,*)(UPL   (j),j=1,JA1)
      if (NUMALR .eq. 236) write(JCH,*)(VIMP1 (j),j=1,JA1)
      if (NUMALR .eq. 237) write(JCH,*)(VIMP2 (j),j=1,JA1)
      if (NUMALR .eq. 238) write(JCH,*)(VIMP3 (j),j=1,JA1)
      if (NUMALR .eq. 239) write(JCH,*)(VOLUM (j),j=1,JA1)
      if (NUMALR .eq. 240) write(JCH,*)(VP    (j),j=1,JA1)
      if (NUMALR .eq. 241) write(JCH,*)(VPFP  (j),j=1,JA1)
      if (NUMALR .eq. 242) write(JCH,*)(VPOL  (j),j=1,JA1)
      if (NUMALR .eq. 243) write(JCH,*)(VPOLX (j),j=1,JA1)
      if (NUMALR .eq. 244) write(JCH,*)(VR    (j),j=1,JA1)
      if (NUMALR .eq. 245) write(JCH,*)(VRO   (j),j=1,JA1)
      if (NUMALR .eq. 246) write(JCH,*)(VRS   (j),j=1,JA1)
      if (NUMALR .eq. 247) write(JCH,*)(VRX   (j),j=1,JA1)
      if (NUMALR .eq. 248) write(JCH,*)(VTOR  (j),j=1,JA1)
      if (NUMALR .eq. 249) write(JCH,*)(VTORX (j),j=1,JA1)
      if (NUMALR .eq. 250) write(JCH,*)(ZEF   (j),j=1,JA1)
      if (NUMALR .eq. 251) write(JCH,*)(ZEF1  (j),j=1,JA1)
      if (NUMALR .eq. 252) write(JCH,*)(ZEF2  (j),j=1,JA1)
      if (NUMALR .eq. 253) write(JCH,*)(ZEF3  (j),j=1,JA1)
      if (NUMALR .eq. 254) write(JCH,*)(ZEFX  (j),j=1,JA1)
      if (NUMALR .eq. 255) write(JCH,*)(ZIM1  (j),j=1,JA1)
      if (NUMALR .eq. 256) write(JCH,*)(ZIM2  (j),j=1,JA1)
      if (NUMALR .eq. 257) write(JCH,*)(ZIM3  (j),j=1,JA1)
      if (NUMALR .eq. 258) write(JCH,*)(ZMAIN (j),j=1,JA1)
	return
	end
C======================================================================|
	subroutine	SETVAR_NOTUSED
C Obsolete. Since July 2007 the subroutine is not used.
	implicit none
	character str*2
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer j,js,jj,length
	double precision	YV,YF,YMU,YN,YNE,YNI,YTE,YTI,YZF
	call	markloc("SETVAR"//char(0))
	YZF = 0.
	YNE = 0.
	YNI = 0.
	YTE = 0.
	YTI = 0.
	YV = GP2*BTOR*HRO
	do	j=1,NA1
	   SQEPS(j) = SQRT(AMETR(j)/(RTOR+SHIF(j)))
! VP = c<E_||/h>/B_p	in Hinton & Hazeltine Rev.Mod.Phys.48(1976) p.297
	   if (MU(j) .le. 0.d0)	then
	      write(*,'(/A,F10.6,A)')'>>> ERROR >>>  Time =',TIME,' sec'
	      write(*,'(A,F10.6)')
     >	'               The rotational transform is less or equal zero '
     >,			'at rho_N =',RHO(j)/ROC
	      if (TIME .le. TSTART+TAU/2.)	write(*,'(A/2A,1H"/)')
     >	'               Most probably it has not been properly defined'
     >,	'               by the data file "',RDNAME(1:length(RDNAME))
	      call	a_stop
	   endif
	   VP(j) = ULON(j)/(YV*j*MU(j))
	   YZF = max(YZF,ZEF(j))
	   YNE = max(YNE,NE(j))
	   YNI = max(YNI,NI(j))
	   YTE = max(YTE,TE(j))
	   YTI = max(YTI,TI(j))
	enddo
	str = '3M'
	if (min(YZF,YNE,YNI,YTE,YTI) .gt. 0.)	then
	   goto	11
	elseif (YZF .le. 0.)	then
	   str = "ZF"
	elseif (YNE .le. 0.)	then
	   str = "NE"
	elseif (YNI .le. 0.)	then
	   str = "NI"
	elseif (YTE .le. 0.)	then
	   str = "TE"
	elseif (YTI .le. 0.)	then
	   str = "TI"
	endif
	write(*,'(/A,F6.3)')'>>> ERROR >>>  Time =',TIME
	if (str .ne. '3M' .and. str .ne. 'ZF') write(*,'(3A/)')
     >	'               The variable "',str,'" is less or equal zero'
	if (str .eq. '3M') write(*,'(A/)')
     >	'               The 3M equilibrium solver does not converge'
	if (str .eq. 'ZF') write(*,'(A/)')
     >	'               The variable "ZEF" is less or equal zero'
	call	a_stop

 11	continue
	if (NA1 .eq. NB1)	goto	12

	do	j=NA1+1,NB1
Common /STATUS/: TE -> CV	! excluding FV, MV, AMETR, RHO
	   do	jj=1,66
		STAARR(j,jj)=STAARR(NA1,jj)
	   enddo
Common /A_IONS/: NN -> PSXIM3
	   do	jj=1,31
		EXTARR(j,jj)=EXTARR(NA1,jj)
	   enddo
Common /A_CARS/: CAR1 -> CAR32
	   do	jj=1,16
		CAR(j,jj)=CAR(NA1,jj)
	   enddo
	   YV = exp((ABC-AMETR(j))/WNE)
	   NE(j) = NE(NA1)*YV
	   NI(j) = NI(NA1)*YV
	   NALF(j) = NALF(NA1)*YV
	   NHE3(j) = NHE3(NA1)*YV
	   NHYDR(j) = NHYDR(NA1)*YV
	   NDEUT(j) = NDEUT(NA1)*YV
	   NTRIT(j) = NTRIT(NA1)*YV
	   TE(j) = TE(NA1)*exp((ABC-AMETR(j))/WTE)
	   TI(j) = TI(NA1)*exp((ABC-AMETR(j))/WTI)
C Should be changed when a SOL equilibrium is implemented
	   MU(j) = MU(NA)*G22(NA)*ROC/(HRO*j)**2
	   FP(j) = FP(NA1)
	   FV(j) = FV(NA1)
	enddo
 12	continue
	end
C======================================================================|
