C======================================================================|
	subroutine	A2G
C----------------------------------------------------------------------|
	implicit none
	integer	JEQUIL,ierr,NP,NY,j,i,length
	parameter(NP=21,NY=32)
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision	C(NRD),BA(NRD),BB(NRD),BC(NRD),BD(NRD),
     .		GBD(NRD),GL(NRD),GSD(NRD),DGBD(NRD),DGL(NRD),DGSD(NRD),
     .		XTR(NRD),XEQ(NRD),AGR(NP),GR(NRD)
	character*80	action*8,timems*3,asfile,file	! XDR file name
	integer	ll,ml,nl,sym,info,ierror
	double precision PRD(NP),PCD(NP),ETA(NP),THETA(NY),PHI(1)
C----------------------------------------------------------------------|
	double precision BNT(NP,NY),BNP(NP,NY),LABEL(NP),HZZ(NP,NY)
	double precision HZDR(NP,NY),HZDZ(NP,NY),HZR(NP,NY)
	double precision YB1,YTH2,YG,ALFA,YAI,YRI,YTHJ,YSJ,YCJ,YAILI,YY
	double precision GSDOA,GROA,YBNT,YS2J,YSQG
	data	nl/1/
	JEQUIL = NEQUIL
C	write(*,*)NEQUIL,JEQUIL
	if (NP .ne. JEQUIL)	then
	write(*,*)'>>> Astra_to_Garbo: cannot write XDR file.'
     &		,'Change parameter Nequil to ',NP
		return
	endif
	YB1	=(RTOR+SHIFT)/RTOR
	do	J=1,NA1
		XTR(J)	=AMETR(J)/AMETR(NA1)
		BD(J)	=YB1*EQPF(J)
		YTH2	=RHO(j)*G22(j)*(MU(j)/RTOR)**2
		YG	=(1.+YTH2)*G33(j)
		C(J)	=YB1*CU(J)/IPOL(J)/YG+BD(J)*(1.-YB1*YB1/YG)
C dI/d\psi in ASTRA notations:	-EQFF/(5.*BTOR*IPOL(j))
C Multiply by 2\pi/\mu_0 = 5.E6    (HPZ)
		BC(j)	=-1.E6*EQFF(J)/(BTOR*IPOL(j))
	enddo
* From the transport grid in "a" to the equidistant grid in "a"
	do	J=1,JEQUIL
		XEQ(J)	=(J-1.)/(JEQUIL-1.)
	enddo
	ALFA	=.001
	CALL	SMOOTH(ALFA,NA1,C,XTR,JEQUIL,BA,XEQ)
	CALL	SMOOTH(ALFA,NA1,BD,XTR,JEQUIL,BB,XEQ)
	CALL	SMOOTH(ALFA,NA1,BC,XTR,JEQUIL,C,XEQ)
	CALL	SMOOTH(ALFA,NA1,CC,XTR,JEQUIL,BD,XEQ)

	do	j=1,JEQUIL
		PCD(j) = C(j)
		ETA(j) = 1.E-6/BD(j)
C \prti{p}{\psi} in [J/(V s m^3)]=[J/(T m^5)]=[Pascal/(V s)]
		PRD(j) = -BB(j)/(1.E-6*GP2*(RTOR+SHIFT))
	enddo
C-------------------------------------------------------------------------
	CALL	BFIELD(BA,BB,RTOR+SHIFT,ABC,ELONG,TRIAN*ABC,JEQUIL,
     .		GR,GBD,GL,GSD,DGBD,DGL,DGSD,BA,BB,BC,C,BD,
     .		BTOR*RTOR/(RTOR+SHIFT),IPL)
C-------------------------------------------------------------------------
C	Output:	GR  - rho
C		GBD - shift
C		GL  - elongation
C		GSD - triangularity	= \delta^{Zakh}	= a*\delta^{Astra}
C		DGBD- \Delta'
C		DGL - \lambda'
C		DGSD- (a\delta^{Astra})'
C		BA  - <g^{33}>		= <1/r^2>=G33/RTOR**2
C		BB  - IPOL*RTOR*BTOR	= I
C		BC  - d\rho/da
C		C   - (dV/da)/(4\pi^2)
C		BD  - V(a)
C-------------------------------------------------------------------------
	do	j=1,NY
	    THETA(j) = GP2*(j-1.)/NY
	enddo
	do	i=1,JEQUIL
	    YAI = AMETR(NA1)*XEQ(i)
	    YRI = RTOR+GBD(i)
	    do	j=1,NY
		YTHJ = THETA(j)
		YSJ = sin(YTHJ)
		YCJ = cos(YTHJ)
		HZR(i,j) = YRI+YAI*YCJ-GSD(i)*YSJ**2
		HZZ(i,j) = YAI*GL(i)*YSJ
		HZDR(i,j) = DGBD(i)+YCJ
		HZDZ(i,j) = (YAI*DGL(i)+GL(i))*YSJ
	    enddo
	enddo
	do	i=1,NA-1
	    AGR(i) = i*HRO/GR(JEQUIL)
	enddo
	AGR(NA) = 0.5*(1.+XTR(NA-1))
	AGR(NA1)= 1.
	CALL	SMOOTH(ALFA,NA1,MU,AGR,JEQUIL,BA,GR)
	do	i=1,JEQUIL
	    YAI = AMETR(NA1)*XEQ(i)
	    LABEL(i) = YAI
	    YAILI = YAI*DGL(i)
	    if(i.ne.1)	then
		GSDOA = GSD(i)/YAI
		GROA = GR(i)/YAI
	    else
		GSDOA = DGSD(1)
		GROA = BC(1)
	    endif
	    YY = 2.*DGL(i)*GSD(i)+2.*GL(i)*GSDOA-GL(i)*DGSD(i)
	    YBNT = BTOR*BA(i)*GROA*BC(i)
	    do	j=1,NY
		YTHJ = THETA(j)
		YCJ = cos(YTHJ)
		YS2J = sin(YTHJ)**2
		BNP(i,j) = BB(i)/HZR(i,j)**2
		YSQG = GL(i)*(1.+DGBD(i)*YCJ)+(YAILI+YY*YCJ)*YS2J
		BNT(i,j) = YBNT/(YSQG*HZR(i,j))
	    enddo
	enddo
	ll = JEQUIL
	ml = NY
	sym = 11
	PHI(1) = 0.
	asfile = 'ToGarbo.inp'
	call OPENWT(8,asfile,0,ierr)
	if (ierr .eq. 2)write(*,*)'Cannot open file "',asfile
	write(8,*)'Array: ll*ml R-coordinate values'
	write(8,'(1P7E15.7)')HZR
	write(8,*)'Array: ll*ml z-coordinate values'
	write(8,'(1P7E15.7)')HZZ
	write(8,*)'Array: ll*ml dR/da - values'
	write(8,'(1P7E15.7)')HZDR
	write(8,*)'Array: ll*ml dz/da - values'
	write(8,'(1P7E15.7)')HZDZ
	write(8,*)'Array: ll*ml BDotNablaTheta values'
	write(8,'(1P7E15.7)')BNT
	write(8,*)'Array: ll*ml BDotNablaPhi values'
	write(8,'(1P7E15.7)')BNP
	write(8,*)'Array: ll flux surface label values'
	write(8,'(1P7E15.7)')LABEL
	write(8,*)'Number of label values'
	write(8,'(I8)')NP
	write(8,*)'Array: ml equidistant theta values'
	write(8,'(1P7E15.7)')THETA
	write(8,*)'Number of angle theta values (power of 2)'
	write(8,'(I8)')NY
	write(8,*)'Array: ll pressure derivative values'
	write(8,'(1P7E15.7)')PRD
	write(8,*)'Array: ll poloidal current derivative values'
	write(8,'(1P7E15.7)')PCD
	write(8,*)'Array: ll resistivity values in Ohm*m'
	write(8,'(1P7E15.7)')ETA
	write(8,*)'Symmetry flag'
	write(8,'(3I8)')sym
	close(8)
	write(timems,'(1I3.3)')int(time*1000)
	write(*,*)'File "',asfile(1:length(asfile)),'" is created.',
     &			'  Error code =',ierr
	file = 'Astra/aug6905.'//timems
	action = 'write'
	call	rwGarboXDR(file,action, 
     &  HZR, HZZ,    ! Arrays: ll*ml*1 R,z-coordinate values
     &  HZDR,HZDZ,   ! Arrays: ll*ml*1 dR/da, dz/da values
     &  BNT, BNP,    ! ll*ml*1 Scalar products (B.Nabla(Theta)),(B.Nabla(Phi))
     &  LABEL,ll,    ! Array(ll) flux surface label values
     &  THETA,ml,    ! Array(ml) equidistant theta values
     &  PHI,  nl,    ! Array(1) phi value
     &  PRD,PCD,     ! Arrays(ll) pressure and poloidal current derivatives
     &  ETA,         ! Array: ll resistivity values
     &  sym,info,ierror)        ! integer symmetry,info and error flags
	write(*,*)'XDR file "',file(1:length(file)),'" is created.'
     &		,'  Error code =',ierror
	END
C======================================================================|
      subroutine BFIELD(BA,BB,BR00,SA0,GL0,GD30,NA1,GR,GBD,GL,GSD
     .		,GRA,SQGRA,GRAR,AVR2,AI0,dgrda,avsqg,Vol,B0T,PLCUR)
C Input:
C     &	(BA
C     &	,BB	! j_zeta = BA*(R_0/r) + BB*(r/R_0-R_0/r)
C     &	,BR00	! R_0+\Delta_edge		! RTOR+SHIFT
C     &	,SA0	! a_edge			! ABC
C     &	,GL0	! \lambda_edge			! ELONG
C     &	,GD30	! \delta_edge			! TRIAN*ABC
C     &	,NA1	! radial grid point No.		! JEQUIL
C Output:
C     &	,GR	! \Rho(a)			-> 
C     &	,GBD	! \Delta(a)			-> SHIF
C     &	,GL	! \lambda(a)			-> ELON
C     &	,GSD	! \delta(a)			-> TRIA
C     &	,GRA	! \Delta'
C     &	,SQGRA	! \lambda'
C     &	,GRAR	! \delta'
C     &	,AVR2	! <g^{33}>			-> <1/r^2>=G33/RTOR**2
C     &	,AI0	! I				-> IPOL*RTOR*BTOR
C     &	,dgrda	! \prti{\rho}{a}
C     &	,avsqg	! \prti{V}{a}/(4\pi^2)
C     &	,Vol	! V(a)
C Input:
C     &	,B0T	! B_0 at the magnetic axis	! BTOR*RTOR/(RTOR+SHIFT)
C     &	,PLCUR	! Total plasma current		! IPL
C     &	)
C - BR00,SA0 - MAJOR & MINOR RADII /METER/
C - GLO,GD3O - ELONGATION AND TRIANGULARITY /BOUNDARY VALUE/
C - GSD = TRIANGULARITY
C GRA,GRAR,SQGRA change meaning here to		!!!!!!!!!!!!!
C GRA	\Delta' =WDSD1(I)*WSA(I)
C SQGRA	\lambda'=WDGL(I)*WGL(I)*WSAA(I)
C GRAR	\delta' =WDSD3(I)*WSA(I)

	implicit none
	include 'for/emeq.inc'
      double precision BA(NP),BB(NP),GR(NP),GBD(NP),GL(NP),GSD(NP)
     &		,GRA(NP),SQGRA(NP),GRAR(NP),AVR2(NP),AI0(NP)
     &		,dgrda(NP),avsqg(NP),Vol(NP)
      double precision AOLD,GLOLD,G3DOLD,cgp,SA0,GL0,GD30,BR00
     &		,B0T,PLCUR,D0,GR2,fi,fj,s,sqg
      integer NA1,NAOLD,NA,NT,NITER,I,j
      save AOLD,GLOLD,G3DOLD,NAOLD
      external EQC1

      data AOLD/0./GLOLD/0./G3DOLD/0./NAOLD/1/,cgp/3.14159265359/

*************************************************
C - TO CALCULATE EQUILIBRIUM SURFACE CHARACTERISTICS B0(T) AND Ipl(MA)
C - ARE REQUIRED, TAKE THEM FROM ASTRA COMMON BLOCKS
      WBBS0=B0T
      WBJ0=0.2*PLCUR
**************************************************

      NA=NA1-1
      NT=12*MAX(1,NA1/8)
      IF(NAOLD.EQ.NA.AND.ABS(AOLD-SA0).LT.1.E-6.AND.ABS(GLOLD-GL0).LT.
     *	.1E-6.AND.ABS(G3DOLD-GD30).LT.1E-6)GOTO 4
      CALL EQGB3(BR00,SA0,GL0,GD30,NA)
      AOLD=SA0
      GLOLD=GL0
      G3DOLD=GD30
      NAOLD=NA
      WBR00=BR00
 4    continue
      NITER=20
      DO  I=1,NA1
         WSJP(I)=BA(I)
         WSP(I)=BB(I)
      enddo
      CALL EQAB3(NA,NT,NITER,0.1,EQC1)
      CALL EQPPAB(NA)
      D0=WSD1(NA1)*WSAA(NA1)
      GR(1)=0.
      Vol(1)=0.
      GR2=0.
      fi=0.
      s=4.*cgp*cgp
      DO I=1,NA1
         Vol(i)=s*WSAA(i)*(WBR0*wsl0(i)+WSAA(i)*wsl1(i))
         j=i-1
         fj=fi
         AI0(I)=WBF(I)*WBR0
         sqg=wbg0(i)*wgl(i)*wbr0
         avsqg(i)=WSA(i)*sqg
         GRA(I)=WDSD1(I)*WSA(I)
         SQGRA(I)=WDGL(I)*WGL(I)*WSAA(I)
         GRAR(I)=WDSD3(I)*WSA(I)
         AVR2(I)=wbg33(i)*wgl(i)/(wbr0*sqg)
         GL(I)=WGL(I)
         GBD(I)=D0-WSD1(I)*WSAA(I)
         GSD(I)=WSD3(I)*WSAA(i)
         fi=wbg33(i)*wgl(i)*AI0(I)/(wbr0*WBBS0)
         IF(I.GT.1) then
	    GR2=GR2+fj*WSCJ1(j)+fi*WSCI1(j)
            GR(I)=SQRT(GR2*2.)
            dgrda(i)=fi*WSA(i)/gr(i)
         else
            dgrda(i)=sqrt(fi)
         endif
      enddo
 30   format(5(f10.4))
      END
C======================================================================|
