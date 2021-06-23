C ./neutex /home/grp/TAstra/.tsk/Intel/neuts.exe 14638 772309392 1
C----------------------------------------------------------------------|
C Make as
C $> cc -c xpr/neut.c -o xpr/neut.o
C $> cc -c xpr/ipc.c -o xpr/ipc.o
C $> f95f xpr/neut.f xpr/neut.o xpr/ipc.o -o xpr/neut
C----------------------------------------------------------------------|
	program   main
	implicit  none
	integer	  j,iargc,getpid,mampid,mamkey,eignr
	character*132 STRING,eigpath,mampath
	if (iargc() .ne. 4)	then
	   write(*,*)"Error"
	   call a_stop
	endif
	call	getarg(0,STRING)
	eigpath = STRING(1:len_trim(STRING))//char(0)
	call	getarg(1,STRING)
	mampath = STRING(1:len_trim(STRING))//char(0)
	call	getarg(2,STRING)
	read(STRING,*)mampid
	call	getarg(3,STRING)
	read(STRING,*)mamkey
	call	getarg(4,STRING)
	read(STRING,*)eignr
	write(*,*)
	call sbp2shm(eigpath, mampath, mampid, mamkey, eignr)
C	do j=0,iargc()
C	   write(*,*)
C	   call	getarg(j,STRING)
C	   write(*,*)j,'"',STRING(1:len_trim(STRING)),'"'
C	   write(*,*)len(STRING),len_trim(STRING)
C	enddo
	end
C----------------------------------------------------------------------|
C ./neutex /home/grp/MAstra/.tsk/Intel/neuts.exe 17527 1996962867 1 &
C======================================================================|
	subroutine A_NEUTEX(
     >		A_ABC,
     >		A_AMJ,
     >		A_NA1,
     >		A_NAB,
     >		A_NB1,
     >		A_NNCX,
     >		A_NNCL,
     >		A_NNWM,
     >		A_ENCL,
     >		A_ENWM,
     >		A_TE,
     >		A_TI,
     >		A_NE,
     >		A_NI,
     >		A_AMAIN,
     >		A_METR,
     >		A_SNNBM,
     >		A_NN,
     >		A_TN,
     >		A_ALBPL
     >			)
C-----------------------------------------------------------22.01.97---|
C	Input:	ABC,AMJ,NA1,NAB,NB1,NNCX
C		AMETR(j),AMAIN(j),TE(j),TI(j),NE(j),NI(j)
C		ENCL,ENWM or wall neutral distribution
C		NNCL,NNWM
C	Warning:	ENCL > 0.5;	ENWM =/= 0 if NNWM =/= 0
C	Output:	NN,	TN,    ALBPL
C---------------------------------------------CHANGED BY POLEVOY-------|
	implicit none
	integer A_NA1, A_NAB, A_NB1, A_NNCX
	double precision A_ABC,A_AMJ,A_ENCL,A_ENWM,A_NNCL,A_NNWM,A_ALBPL
	double precision
     >		A_TE(*),
     >		A_TI(*),
     >		A_NE(*),
     >		A_NI(*),
     >		A_AMAIN(*),
     >		A_METR(*),
     >		A_SNNBM(*),
     >		A_NN(*),
     >		A_TN(*)
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	J
	call A_NEUTEX_(
     >		A_ABC,
     >		A_AMJ,
     >		A_NA1,
     >		A_NAB,
     >		A_NB1,
     >		A_NNCX,
     >		A_NNCL,
     >		A_NNWM,
     >		A_ENCL,
     >		A_ENWM,
     >		A_TE,
     >		A_TI,
     >		A_NE,
     >		A_NI,
     >		A_AMAIN,
     >		A_METR,
     >		A_SNNBM,
     >		A_NN,
     >		A_TN,
     >		A_ALBPL
     >			)
	return
	end
C----------------------------------------------------------------------|
	subroutine A_NEUTEX_(
     >		A_ABC,
     >		A_AMJ,
     >		A_NA1,
     >		A_NAB,
     >		A_NB1,
     >		A_NNCX,
     >		A_NNCL,
     >		A_NNWM,
     >		A_ENCL,
     >		A_ENWM,
     >		A_TE,
     >		A_TI,
     >		A_NE,
     >		A_NI,
     >		A_AMAIN,
     >		A_METR,
     >		A_SNNBM,
     >		A_NN,
     >		A_TN,
     >		A_ALBPL
     >			)
C-----------------------------------------------------------22.01.97---|
C	Input:	ABC,AMJ,NA1,NAB,NB1,NNCX
C		AMETR(j),AMAIN(j),TE(j),TI(j),NE(j),NI(j)
C		ENCL,ENWM or wall neutral distribution
C		NNCL,NNWM
C	Warning:	ENCL > 0.5;	ENWM =/= 0 if NNWM =/= 0
C	Output:	NN,	TN,    ALBPL
C---------------------------------------------CHANGED BY POLEVOY-------|
	implicit none
	integer A_NA1, A_NAB, A_NB1, A_NNCX
	double precision A_ABC,A_AMJ,A_ENCL,A_ENWM,A_NNCL,A_NNWM,A_ALBPL
	double precision
     >		A_TE(*),
     >		A_TI(*),
     >		A_NE(*),
     >		A_NI(*),
     >		A_AMAIN(*),
     >		A_METR(*),
     >		A_SNNBM(*),
     >		A_NN(*),
     >		A_TN(*)
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	J
C----------------------------------------------------------------------|
	ABC = A_ABC			! Not used
	AMJ = A_AMJ			! Wall neutral mass
	NA1 = A_NA1
	NAB = A_NAB
	NB1 = A_NB1
	NNCX = A_NNCX
	NNCL = A_NNCL
	NNWM = A_NNWM
	ENCL = A_ENCL
	ENWM = A_ENWM
	do j=1,NB1
	   TE(j) = A_TE(j)
	   TI(j) = A_TI(j)
	   NE(j) = A_NE(j)
	   NI(j) = A_NI(j)
	   AMETR(j) = A_METR(j)
	   AMAIN(j) = A_AMAIN(j)	! Main ion mass
	   SNNBM(j) = A_SNNBM(j)
	enddo
	call	NEUTEX
	A_ALBPL = ALBPL
	do j=1,NAB
	   A_NN(j) = NN(j)
	   A_TN(j) = TN(j)
	enddo
	end
C----------------------------------------------------------------------|
	subroutine	NEUTEX
C-------------------------------------- Pereverzev ---------07.09.07---|
C	Input:	ABC, NA1, NB1, NNCX
C		AMETR(j), AMAIN(j), TE(j), TI(j), NE(j), NI(j), SNNBM(j)
C		ENCL, ENWM or wall neutral distribution
C		NNCL, NNWM
C	Warning:	ENCL > 0.5;	ENWM =/= 0 if NNWM =/= 0
C	Output:	  NN,	TN,   ALBPL
C----------------------------------------------------------------------|
	implicit none
C----------------------------------------------------------------------|
	include	'for/parameter.inc'
C	integer NRD
C	parameter(NRD=501)
C----------------------------------------------------------------------|
	include 'for/const.inc'
C	double precision	! Only these parameters are used:
C     1		ABC, AMJ, TIME, NNCL, NNWM, ENCL, ENWM, NNCX, ALBPL
C	integer NA1, NAB, NB1
C----------------------------------------------------------------------|
	include 'for/status.inc'
C	double precision	! Only these arrays are used:
C     &		TE(NRD),   TI(NRD),   NE(NRD),  NI(NRD),  AMAIN(NRD),
C     1		AMETR(NRD),SNNBM(NRD),NN(NRD),  TN(NRD)
C	double precision	WORK1(NRD,2*NRD+7)
C----------------------------------------------------------------------|
	integer j,jj,JI,JNN,JNA,jcall
	double precision
     1		AKDLT,SVIE,SVII,SVREC,SVCX,YAB,
     2		YV1,YV2,YN0,YZ1,YZ2,YHA,YH1,YG1,YH2,YG2,YY1,YY2,YE1,YE2,
     3		YF(NRD),YC(NRD),YS(NRD),YA(NRD),YR(NRD),YX(NRD),YV(NRD),
     4		YN(NRD),YT(NRD),YK(NRD,NRD)
	equivalence (WORK1(1,1),YF(1)),   (WORK1(1,2),YC(1)),
     1		    (WORK1(1,3),YS(1)),   (WORK1(1,4),YA(1)),
     2		    (WORK1(1,5),YR(1)),	  (WORK1(1,6),YX(1)),
     3		    (WORK1(1,7),YN(1)),	  (WORK1(1,8),YT(1)),
     4		    (WORK1(1,9),YV(1)),	  (WORK1(1,10),YK(1,1))
	data jcall/0/
	save jcall
C----------------------------------------------------------------------|
C----------------------------------------------------------------------|
	JNA = NA1			! Entry NEUTEX (ignore SOL)
	goto	2
C	entry	NEUTAB
	JNA = NAB			! Entry NEUTAB (include SOL)
 2	continue
C-------------------------------       Set the internal grid
C	JNN = NRD
	JNN = NB1
	YAB = AMETR(JNA)		! ABC(NEUTEX) or AB(NEUTAB)
	YHA = YAB/(JNN-1.)
	do	j=1,JNN
	   YA(j) = (j-1.)/(JNN-1.)
	enddo
	do	j=1,JNA
	   YR(j) = AMETR(j)/AMETR(JNA)
	   YS(j) = 2.52d5*SQRT(max(TI(j),1.d-3)/AMAIN(j)) ! v_Ti/sqrt(3)
	enddo
	call	SMOOTH(1.d-3,JNA,YS,YR,JNN,YV,YA)
C-------------------------------	Some checks
	if (ENCL .lt. 0.0001)	then
	   write(*,'(A,F6.3,A)')
     1	      ">>> NEUT >>> Too low energy of incoming neutrals ENCL ="
     2	      ,1000.*ENCL,"eV"
	   write(*,*) "             Setting ENCL = 2 eV"
	   ENCL = 0.002
	endif
	if (ENWM .lt. 0.002)	then
	   write(*,'(A,F6.3,A)')
     1	      ">>> NEUT >>> Too low energy of incoming neutrals ENWM ="
     2	      ,1000.*ENWM,"eV"
	   write(*,*) "             Setting ENWM = 2 eV"
	   ENWM = 0.002
	endif
C-------------------------------	Atomic reaction rates:
	do	j=1,JNA
	   include 'fml/svrec'
	   YC(j) = SVREC*NE(j)*NI(j)+SNNBM(j)
	   include 'fml/svcx'
	   YN(j) = SVCX*NI(j)
	   include 'fml/svie'
	   include 'fml/svii'
	   YX(j) = YN(j)+SVIE*NE(j)+SVII*NI(j)		! YX=s(a)
	enddo
	call	SMOOTH(1.d-3,JNA,YX,YR,JNN,YS,YA)	! YS=s(x)
	YZ1 = 0.5*YHA
	YX(1) = 0.
	do	j=1,JNN-1
	   YX(j+1) = YX(j)+(YS(j)+YS(j+1))*YZ1		! YX(j)=X(j)
	enddo
C	call	SMOOTH(1.d-3,JNN,YS,YA,JNA,CAR(1,29),YR)
C	call	SMOOTH(1.d-3,JNN,YX,YA,JNA,CAR(1,30),YR)
C-------------------------------	Compute kernel:
	do	jj=1,JNN				! jj -> x
	   YH2 = 0.
	   YG2 = 0.
	   YZ2 = abs(YX(1)-YX(jj))/YV(1)		! h(x_jj,\xi_1)
	   YE2 = AKDLT(YZ2)
	   do	 j=1,JNN				! j  -> \xi
	      YH1 = YH2
	      YG1 = YG2
	      YZ1 = YZ2
	      YE1 = YE2
	      if (j .lt. JNN)	then
		 YZ2 = abs(YX(j+1)-YX(jj))/YV(j+1)	! g(x_jj,\xi_{j+1})
		 YE2 = AKDLT(YZ2)			! exp{-g}
		 YH2 = (YE2-YE1)/(YZ2-YZ1)
		 YG2 = YHA/(YZ2-YZ1)			! 1/g'_{i+1/2}
	      else
		 YG2 = 0.
	      endif
	      YY1 = YG1*(YE1+YH1)			! H_j 1st term
	      YY2 = YG2*(YE1+YH2)			! G_j 1st term
	      YK(jj,j) = YY2-YY1
	   enddo
	enddo
	do	jj=1,JNN
	   YH2 = 0.
	   YG2 = 0.
	   YZ2 = (YX(1)+YX(jj))/YV(1)			! h(x_jj,\xi_1)
	   YE2 = AKDLT(YZ2)
	   do	 j=1,JNN
	      YH1 = YH2
	      YG1 = YG2
	      YZ1 = YZ2
	      YE1 = YE2
	      if (j .lt. JNN)	then
		 YZ2 = (YX(j+1)+YX(jj))/YV(j+1)		! h(x_jj,\xi_{j+1})
		 YE2 = AKDLT(YZ2)			! exp{-h}
	         YH2 = (YE2-YE1)/(YZ2-YZ1)		! d[e^(-h)]/dh
	         YG2 = YHA/(YZ2-YZ1)			! 1/h'_{j+1/2}(x_jj)
	      else
		 YG2 = 0.
	      endif
	      YY1 = YG1*(YE1+YH1)			! H_j 2nd term
	      YY2 = YG2*(YE1+YH2)			! G_j 2nd term
	      YK(jj,j) = YK(jj,j)+YY2-YY1		! G_j(x_jj)-H_j(x_jj)
	   enddo
	enddo			! YK(jj,j)=K(x_jj,\xi_j) - intermediate result
C	jj = 1
C	do	j=1,20
C	   call	SMOOTH(1.d-3,JNN,YK(1,jj),YA,JNA,CAR(1,j),YR)
C	   jj = jj+2
C	enddo
C---------------------------------	Zero guess N_0(x)+R(x)
	YN0 = NNCL+NNWM+1.d-11
	YY1 = (NNCL+1.d-11)/YN0
	YY2 = NNWM/YN0
	YV1  = 4.37d5*SQRT(ENCL/AMJ)
	YV2  = 4.37d5*SQRT(ENWM/AMJ)
	call	SMOOTH(1.d-3,JNA,TI,YR,JNN,YT,YA)	! T_i -> Int. grid
	call	SMOOTH(1.d-3,JNA,YN,YR,JNN,YF,YA)	! YF=svcx*n_i
	call	SMOOTH(1.d-3,JNA,YC,YR,JNN,YS,YA)	! YS =svrec*n_e*n_i
 	do	j=1,JNN
	   YF(j) = 0.5*YF(j)/YV(j)			! F(j) factor in kernel
	   YS(j) = 0.5*YS(j)/YV(j)			! Factor in R(j)
	enddo
	do	j=1,JNN					! x_j
	   YZ1 = (YX(JNN)+YX(j))/YV1
	   YZ2 = (YX(JNN)-YX(j))/YV1
	   YH1 = YY1*(AKDLT(YZ1)+AKDLT(YZ2))		! K(x,a,v_1)
	   NN(j) = YH1				! 1st species, 0th generation
	   TN(j) = ENCL*YH1
	   if (NNWM .le. 0.)	goto	9
	   YZ1 = (YX(JNN)+YX(j))/YV2
	   YZ2 = (YX(JNN)-YX(j))/YV2
	   YH2 = YY2*(AKDLT(YZ1)+AKDLT(YZ2))		! K(x,a,v_2)
	   NN(j) = NN(j)+YH2			! 2nd species, 0th generation
	   TN(j) = TN(j)+ENWM*YH2
 9	   continue
	   YG1 = 0.
	   YG2 = 0.
	   do	jj=1,JNN
	      YZ1 = YS(jj)*YK(j,jj)		! YS=[<sv^rec>*n_e*n_i/(2v_Ti)]
	      YG1 = YG1+YZ1
	      YG2 = YG2+YZ1*YT(jj)
	   enddo
	   NN(j) = NN(j)+YG1/YN0		! NN_j=N_0(x_j)
	   TN(j) = TN(j)+YG2/YN0		! TN_j=T_0(x_j)
	enddo
C---------------------------------	Auxiliary output
C	call	SMOOTH(1.d-5,JNN,NN,YA,JNA,CAR(1,31),YR)
C	call	SMOOTH(1.d-5,JNN,TN,YA,JNA,CAR(1,32),YR)
C---------------------------------	Finalize the kernel \tilde{K} = F*K
 	do	j=1,JNN
	   do	jj=1,JNN
	      YK(jj,j) = YF(j)*YK(jj,j)		! \tilde{K} = F*K
	   enddo
	   YS(j)  = NN(j)			! YS=NN_(0)
C	   YS(j)  = 1.				! Check convergence
	enddo
C---------------------------------	Kernel is ready
C---------------------------------	Iterations:
	do	10	JI=1,NNCX
	   do	j=1,JNN
	      YN(j)  = YS(j)			! YN=NN_{i}
	   enddo
C	   if (ji .le. 20)
C     >		call	SMOOTH(1.d-5,JNN,YN,YA,JNA,CAR(1,JI),YR)
	   YY2 = 0.
	   do	j=1,JNN
	      YY1 = 0.
	      do       jj=1,JNN
		 YY1 = YY1+YN(jj)*YK(j,jj)	! YN=NN_(i-1)
	      enddo
	      YS(j) = YY1			! YS=NN_(i)
	      YY2 = max(YY2,YY1)
	   enddo
	   YZ1 = YZ2
	   YZ2 = YY2
	   do	j=1,JNN
	      NN(j)  = NN(j)+YS(j)		! NN(1:) at the internal grid
	   enddo
	   if (JI .eq. 1)	YN0 = YZ2
	   if (JI .eq. 1)	YZ1 = YZ2
	   if (abs(YZ2/YN0) .lt. 1.d-12)	goto	11
C	   write(*,'(1I5,1P,4E12.3)')JI,YZ2/YN0,YZ2/YZ1,YZ1,YZ2
 10	continue
	if (abs(YZ2/YN0) .gt. 1.d-4)	goto	12
 11	continue
C	write(*,100)JI,'    Accuracy =',YZ2/YZ1,YZ2/YN0
 100	format('Number of iterations = ',I4,A,1P,5E12.3)
C---------------------------------	End of iterations
	do	j=1,JNN
	   YS(j) = YT(j)*NN(j)			! T_i*N
	enddo
	do	j=1,JNN
	   YY1	= 0.
	   do	jj=1,JNN
	      YY1  = YY1+YS(jj)*YK(j,jj)
	   enddo
	   YX(j) = (TN(j)+YY1)/NN(j)
C	   TN(j) = max(0.d0,YT(j))
	   YN(j) = NN(j)
	enddo

C albpl: [d/l]	Plasma albedo
C		Pereverzev	15-05-95
C (Neutral_outflux)/(Neutral_influx) =
C	= sqrt{(NN-N1-N2)(NN*TN-N1*E1-N2*E2)}/(N1*sqrt(E1)+N2*sqrt(E2))
	ALBPL=sqrt((NNCL+NNWM)*(YN(JNN)-1)*
     .		((NNCL+NNWM)*YN(JNN)*YX(JNN)-NNCL*ENCL-NNWM*ENWM))/
     .		(NNCL*sqrt(ENCL)+NNWM*sqrt(ENWM))
	call	SMOOTH(1.d-5,JNN,YX,YA,JNA,TN,YR)
	call	SMOOTH(1.d-5,JNN,YN,YA,JNA,NN,YR)
C	write(*,'(1P,6E12.5)')(ABC*YA(j),j=JNN-5,JNN)
C	write(*,'(1P,6E12.5)')(YN(j),j=JNN-5,JNN)
C	stop
	return
 12	write(*,101)TIME
	if (jcall .eq. -1)	return
	write(*,*)"                        Increase NNCX can help"
	jcall = -1
 101	format(' >>> NEUTEX  Warning >>> Poor iteration convergence',
     >         ' @ TIME =',1P,4E12.3)
	end
C======================================================================|
	double precision function AKDLT(X)
	double precision X
	if (X .gt. 50.)	then
	   AKDLT = 0.
	   return
	else
	   AKDLT = exp(-X)
	endif
	end
C======================================================================|
C  Subroutine minimizes the value of functional
C  INTEGRAL(alfa*P(x)*(dU/dx)**2+(U-F)**2)*dx,
C  where FO(NO) is a function, given on the grid XOld(NOld)
C	P(x) is equal to unit now
C	ALFA=alfa<<0.01*XO(NO)**2 is regularizator
C	NO - number of old grid points
C	N=<NRD - number of new grid points
C	0<=XO(NO) - old grid |     both grids are arbitrary
C	0<=XN(N)  - new grid |     but XO(NO)=XN(N)
C	FO(NO) - origin function, given on the grid XO(NO)
C	FN(N) - smoothed function on grid XN(N)
C  The result is function FN(XN), given on the new grid
C	with additional conditions:
C	dFN/dx(x=0)=0 - cylindrical case and
C	FN(XN(N))=FO(XO(NO))
C----------------------------------------------------------------------|
	subroutine	SMOOTH(ALFA,NO,FO,XO,N,FN,XN)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	integer	NO,N,J,I
	double precision	ALFA,XO(*),FO(*),XN(*),FN(*),P(NRD)
	double precision	YF,YX,YP,YQ,YD,FJ
	if (N .gt. NRD .or. NO .le. 0)	then
		write(*,*)' >>> SMOOTH: array is out of limits'
		call	a_stop
	endif
	if (NO .eq. 1)	then
	   do	j=1,N
		FN(j) = FO(1)
	   enddo
	   return
	endif
	if (NO .eq. 2)	then
	  do	j=1,N
	   FN(j)=(FO(2)*(XN(j)-XO(1))-FO(1)*(XN(j)-XO(2)))/(XO(2)-XO(1))
	  enddo
	  return
	endif
	if (N .lt. 2)	then
		write(*,*)' >>> SMOOTH: no output grid is provided'
		call	a_stop
	endif
	if (abs(XO(NO)-XN(N)) .gt. XN(N)/N)	then
	    write(*,*)'>>> SMOOTH: grids are not aligned'
	    write(*,'(1A23,I4,F8.4)')'     Old grid size/edge',NO,XO(NO)
	    write(*,'(1A23,I4,F8.4)')'     New grid size/edge',N,XN(N)
	    call	a_stop
	endif
	do	1	j=2,N
	   YP = (XN(j)-XN(j-1))
	   if (YP .le. 0.d0)	then
	write(*,*)'>>> SMOOTH: new grid is not increasing monotonically'
	      write(*,'(A,I4,A,F8.4)')'Node ',j-1,'   Value',XN(j-1)
	      write(*,'(A,I4,A,F8.4)')'Node ',j,  '   Value',XN(j)
	      call	a_stop
	   endif
	   P(j)	=ALFA/YP/XO(NO)**2
 1	continue
	P(1)	=0.
	FN(1)	=0.
	I	=1
	YF	=(FO(2)-FO(1))/(XO(2)-XO(1))
	YX	=2./(XN(2)+XN(1))
	YP	=0.
	YQ	=0.
	do	5	j=1,N-1
		if(XO(I) .gt. XN(j))	GO TO 4
 3		I	=I+1
		if(I .gt. NO)	I=NO
		if(I .ne. NO .and. XO(I) .lt. XN(j))	GOTO	3
		YF	=(FO(I)-FO(I-1))/(XO(I)-XO(I-1))
 4		FJ	=FO(I)+YF*(XN(j)-XO(I))
		YD=1.+YX*(YP+P(j+1))
		P(j)	=YX*P(j+1)/YD
		FN(j)	=(FJ+YX*YQ)/YD
		if (j .eq. N-1)	goto	5
		YX	=2./(XN(j+2)-XN(j))
		YP	=(1.-P(j))*P(j+1)
		YQ	=FN(j)*P(j+1)
 5	continue
	FN(N)	=FO(NO)
	do	6	j=N-1,1,-1
		FN(j)	=P(j)*FN(j+1)+FN(j)
 6	continue
	end
C======================================================================|
