      SUBROUTINE INITAL
C
C ************ ************ ************ ************ ************
C  Ray tracing code RAYIC
C  Initializatyion routine
C ************ ************ ************ ************ ************
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
         ISTOP = 0
C
C  Consistency and charge neutrality
C
      IF(NUMPHI.GT.NPHRAY)  THEN
         ISTOP = 1
         WRITE(6,9910)
 9910 FORMAT(' Arrays underdimensioned, NUMPHI > NPHRAY')
         RETURN
      END IF
C
      IF(NTHIN.GT.NTHRAY)  THEN
         ISTOP = 1
         WRITE(6,9920)
 9920 FORMAT(' Arrays underdimensioned, NTHIN > NTHRAY')
         RETURN
      END IF
C
      IF(NPROF.GT.NPSIPT)  THEN
         ISTOP = 1
         WRITE(6,9940)
 9940 FORMAT(' Arrays underdimensioned, NPROF > NPSIPT')
         RETURN
      END IF
C
C  Constants
C
         PI = 3.14159265359
C
C  Density, temperature and MHD profiles
C
      IF(NSPEC.EQ.1)  THEN
         MAINSP = 1
         MINSP = 1
      ELSE
C
         IF(MAINSP.GT.NSPEC.OR.MAINSP.LE.0)  THEN
            ISTOP = 1
            WRITE(6,9950)
 9950 FORMAT(' MAINSP not correctly specified')
            RETURN
         END IF
C
         IF(IPROEQ.NE.0)  THEN
            DO 10  I=1,NSPEC
               ZNI = ACONC(I)*PNE(1)/PNI(1,MAINSP)
               ZTI = TEMPIC(I)/PTI(1,MAINSP)
            DO 10  N=1,NPROF
               PNI(N,I) = ZNI*PNI(N,MAINSP)
   10          PTI(N,I) = ZTI*PTI(N,MAINSP)
         END IF
C
         DO 30  J=1,NPROF
            ZDEN = PNE(J)
         DO 20 I=1,NSPEC
            IF(I.NE.MAINSP)  ZDEN = ZDEN - AZI(I)*PNI(J,I)
   20    CONTINUE
C
            IF(ZDEN.LT.0.)   THEN
               ISTOP = 1
               WRITE(6,9960)
 9960 FORMAT(' Charge neutrality cannot be satisfied')
               RETURN
            END IF
C
   30       PNI(J,MAINSP) = ZDEN/AZI(MAINSP)
C
      END IF
C
         MINSP = MAINSP
      DO 40  I=1,NSPEC
         ATZI(I) = ATM(I)/AZI(I)
         ATOPI(I) = AZI(I)*AZI(I)/ATM(I)
         ACONC(I) = PNI(1,I)/PNE(1)
      IF(ACONC(I).LT.ACONC(MINSP).AND.ATZI(I).NE.ATZI(MINSP)) MINSP=I
      IF(ACONC(I).GT.ACONC(MINSP).AND.ATZI(I).EQ.ATZI(MINSP)) MINSP=I
   40 CONTINUE
C
C   Frequency
C
         WLENGT = 2.997925E8/FREQCY
         UKZERO = 2.*PI/WLENGT
         ULENGT = 1./UKZERO
C
         IF(MAXHRM.LT.2)  MAXHRM=2
         IF(MAXHRM.GT.5)  MAXHRM=5
C
C  Parameters for the MHD configuration
C
         DPSI = 1./(NPROF-1)
      DO 60  I=1,NPROF
   60    PPSI(I) = FLOAT(I-1)*DPSI
C
            CALL METRIN
C
C  Splines for the plasma density and temperature profiles
C
            CALL PROFIN
C
      IF(ISTOP.NE.0)  RETURN
C
C  FACILITY FOR PLOTTING THE MHD CONFIGURATION ALONE
C
      IF(IPLOPR.EQ.-1)  THEN
C
         NCOPSI = 5
         NPTHET = 101
         IPRINT = 2
C
C            CALL PLOPSI(NCOPSI,NPTHET)
            CALL OUTRUN(1)
C
         WRITE(6,9970)
 9970 FORMAT(/,' Only the MHD configuration plotted')
         ISTOP = 1
         RETURN
C
      END IF
C
C  Cyclotron and two-ion hybrid resonances
C
            CALL CYTWIN
C
C  Antenna evaluation
C
      IF(NUMPHI.EQ.1 .OR. NTHIN.LE.2)  NANTS = 0
      IF(JOUTA.EQ.-1)  NANTS = 1
C
      IF(NANTS.NE.0)  THEN
C
            CALL ANTIC
C
C  Facility to evaluate only the antenna    
C
         IF(JOUTA.EQ.-1)  THEN
            WRITE(6,9980)
 9980 FORMAT(/,' Only the antenna evaluated')
            ISTOP = 1
            RETURN
         END IF
C
      END IF
C
      IF(NANTS.EQ.0)  THEN
C
C  Antenna evaluation skipped or failed
C
         PSASYM = 1.
         APWEDG = 0.
C
C  Model toroidal wavenumber power spectrum
C
         IF(NUMPHI.EQ.1)  THEN
            NPHI(1) = NPHI1
            PWCPL(1) = 1.
         ELSE
C
            DO 110  N=1,NUMPHI
               NPHI(N) = NPHI1 + (N-1)*JUMPHI
               ZARG = FLOAT(NPHI(N))*WIDTH*UKZERO/2.
               PWCPL(N) = 1.
            IF(ZARG.NE.0)   PWCPL(N) = (SIN(ZARG)/ZARG)**2
  110       CONTINUE
C
            IASYM = 1
            IF(NPHI(1).NE.0)  IASYM = 0
            IF(DPHASE.NE.0.)  IASYM = 0
            IF(THANTN.NE.0.AND.THANTN.NE.180.)  IASYM = 0
            IF(JPOLE.EQ.2.AND.EPHASE.NE.0.)  IASYM = 0
               ZW = 2.
            IF(IASYM.EQ.0)  ZW = 1.
C
            ZPWNOR = PWCPL(1)
            DO  120  N=2,NUMPHI
               PWCPL(N) = ZW*PWCPL(N)
  120          ZPWNOR = ZPWNOR + PWCPL(N)
            DO 130  N=1,NUMPHI
  130          PWCPL(N) = POWER*PWCPL(N)/ZPWNOR
C
         END IF
C
C  Model azimutal distribution of the power flux
C
         ZTHAPR = 180.*HEIGTH/(PI*RPLASM)
C
         IF(NTHIN.EQ.1)  THEN
C
            THSTRT(1) = THANTN
            DO 140  L=1,NUMPHI
  140          POYNT(L,1) = 1.
C
         ELSE
C
            ZPWNOR = 0.
            ZDTH = 2.*ZTHAPR/FLOAT(NTHIN-1)
         DO 150  N=1,NTHIN
            THSTRT(N) = THANTN - ZTHAPR + FLOAT(N-1)*ZDTH
            ZCOSN = COS(PI*(THSTRT(N)-THANTN)*RPLASM/180.)
            ZPWNOR = ZPWNOR + ZCOSN*ZCOSN
         DO 150  L=1,NUMPHI
  150       POYNT(L,N) = ZCOSN*ZCOSN
         DO 160  N=1,NTHIN
         DO 160  L=1,NUMPHI
  160       POYNT(L,N) = POYNT(L,N)/ZPWNOR
C
         END IF
C
      END IF
C
C  Initialization of the output tables
C
         PWTOTE = 0.
         PWIBWS = 0.
      DO 310 I=1,NSPEC  
  310    PWTOTI(I) = 0.
C
      DO 320 N=1,NPROF
         PWE(N) = 0.
      DO 320 I=1,NSPEC
  320    PWI(N,I) = 0.
C
      DO 330  N=1,NUMPHI
         PWPARE(N) = 0.
         PWPARB(N) = 0.
      DO 330  I=1,NSPEC
  330    PWPARI(N,I) = 0.
      DO 340  N=1,NUMPHI
  340    PWTRA(N) = 0.
         PWEDGE = APWEDG
C
         PWCONF = 1.E-5*POWER/FLOAT(NUMPHI*NTHIN)
C
         XIJUMP = 1.75
C
C  Initialization of the smoothing facility
C
         ITWIHR = 0
      IF(ISMOOT.NE.0)  THEN
         DXIEQ = (XMAX - XMIN)/FLOAT(NPROF-1)
         UX = UKZERO*XMAX
         UZ = 0.
         PSI = 1.
         IGUESS = 1
         IXIEQ(1) = NPROF
         ZDUX = DXIEQ*UKZERO
         ACCUR = 1.E-7
      DO 350  I=2,NPROF
         UX = UX - ZDUX
            CALL INVERT(IGUESS)
  350    IXIEQ(I) = IFIX(PSI/DPSI)+1
C
      END IF
C
C  Parameters for the Runge-Kutta integration
C
         ACCUR = 1.E-5
         IND = 1
         NDIM = 5
C
         PRECDY(5) = 0.04
         ZPREC = (1. - PRECDY(5))/4.
      DO 410 I=1,4
  410    PRECDY(I) = ZPREC
         DSTEPH = 0.1
C
         PSIMAX = PSASYM + (1.-PSASYM)/2.
C
C  Control of the graphical output
C
         ZDEYKP = 2.*PI*DEYKPR
      IF(NUMPHI.GT.1)  ZDEYKP = 2000.*PI
      IF(IPRINT.EQ.3)  IGRAPH = 0
      IF(IGRAPH.EQ.0) RETURN
      IF(NUMPHI.GT.1)  IGRAPH = -1
C
      RETURN
      END

      SUBROUTINE CYTWIN
C
C ************ ************ ************ ************ ************
C   This subroutine determines the position
C   of cyclotron and ion-ion resonances.
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
C ************ ************ ************ ************ ************
      DATA     ZCHIH   /6.5594E-08 /
C ************ ************ ************ ************ ************
C
         ZXMAX = UKZERO*XMAX
         ZXMIN = UKZERO*XMIN
C
C  Position of the cyclotron resonances
C
      DO  10  I=1,NSPEC
         ZBRES = ZCHIH*FREQCY*ATZI(I)
      DO  10  L=1,2
         ZX = UKZERO*RAXMAG*(FLOAT(L)*BAXIS/ZBRES-1.)
         UXCYCL(L,I) = ZX
         NCYCLR(L,I) = 0
      IF(ZX.LT.ZXMAX .AND. ZX.GT.ZXMIN)  NCYCLR(L,I) = 1
   10 CONTINUE
C
C  Position of ion-ion resonances and cutoffs
C
         NCOF = 0
         NRES = 0
      IF(NSPEC.EQ.1)  RETURN
C
         ZCONC = 1.
      DO 20  I=1,NSPEC
         IF(ACONC(I).LT.ZCONC .AND. NCYCLR(1,I).NE.0)  ZCONC = ACONC(I)
   20 CONTINUE
C
         NIX = 100
   30    ZDX = (ZXMAX-ZXMIN)/FLOAT(NIX)
C
      IF(ZDX.LE.0.25*ZCONC)  GO TO 40
C
      IF(NIX.GE.500)  GOTO 40
         NIX = NIX+100
      GO TO 30
C
   40    ACCUR = 1.E-5
         PKPAR = 0.
         IGUESS = 1
         ZLOLD = 0.
         ZSOLD = 0.
         THETA = 0.
         PSI = 1.
C
            CALL COORDS
C
            CALL PROFIX
C
            CALL DISPIC(1)
C
      DO 100  IX=1,NIX
C
      IF(NRES.GE.NSPEC-1.AND.NCOF.GE.NSPEC-1)  GO TO 110
C
         UX = UX - ZDX
         ZLWOLD = ZLOLD
         ZSWOLD = ZSOLD
         ZLOLD = AL
         ZSOLD = AS
C
            CALL INVERT(IGUESS)
C
            CALL PROFIX
C
            CALL DISPIC(1)
C
      IF(ZLOLD*AL.LE.0. .AND.
     >   ABS(ZLOLD).LT.ABS(ZLWOLD))  THEN
C
C  FINE SEARCH FOR A CUT-OFF
C
         ZL = AL
         ZS = AS
         ZUX = UX
         ZZDX = AL*ZDX/(AL-ZLOLD)
         IFALSE = 0
C
      DO  50  ITRIN=1,10
         ZXOLD = UX
         ZLOLD = AL
         UX = UX + ZZDX
      IF(UX.GT.ZXMAX.OR.UX.LT.ZXMIN)  GO TO  60
C
            CALL INVERT(IGUESS)
C
            CALL PROFIX
C
            CALL DISPIC(1)
C
         ZZDX = -AL*(UX-ZXOLD)/(AL-ZLOLD)
      IF(ABS(AL/ZLOLD).GT.2.)  IFALSE=IFALSE+1
      IF(IFALSE.GT.5)  GO TO  60
C
      IF(ABS(ZZDX)/URPLAS.LT.ACCUR)  THEN
         NCOF = NCOF+1
         UXCOF(NCOF) = UX + ZZDX
         GO TO  60
      END IF
C
   50 CONTINUE
   60    AL = ZL
         AS = ZS
         UX = ZUX
      END IF
C
      IF(AS*ZSOLD.LE.0. .AND.
     >   ABS(ZSOLD).LT.ABS(ZSWOLD))  THEN
C
C  FINE SEARCH FOR A RESONANCE
C
         ZS = AS
         ZL = AL
         IFALSE = 0
         ZUX = UX
         ZZDX = AS*ZDX/(AS-ZSOLD)
C
      DO  70  ITRIN=1,10
         ZXOLD = UX
         ZSOLD = AS
         UX = UX + ZZDX
      IF(UX.GT.ZXMAX.OR.UX.LT.ZXMIN)  GO TO 80
C
            CALL INVERT(IGUESS)
C
            CALL PROFIX
C
            CALL DISPIC(1)
C
         ZZDX = -AS*(UX-ZXOLD)/(AS-ZSOLD)
      IF(ABS(AS/ZSOLD).GT.2.)  IFALSE = IFALSE+1
      IF(IFALSE.GT.5)  GO TO 80
C
      IF(ABS(ZZDX)/URPLAS.LT.ACCUR)  THEN
         NRES = NRES+1
         UXRES(NRES) = UX
         GOTO 80
      END IF
C
   70 CONTINUE
C
   80    AS = ZS
         AL = ZL
         UX = ZUX
C
      END IF
C
  100 CONTINUE
C
  110 RETURN
      END

      SUBROUTINE PROFIN
C
C ************ ************ ************ ************ ************
C
C      THIS SUBROUTINE INITIALIZES THE DENSITY, TEMPERATURE,
C      AND EVALUATES THE COEFFICIENT FOR SPLINE INTERPOLATION.
C
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
C ************ ************ ************ ************ ************
C  Transfer of the profiles for density & temperatures
C     to the tables for the cubic interpolation
C
         DENC = PNE(1)
         TEMPEC = PTE(1)
      DO 10  J=1,NPROF
         TNE(J,1) = PNE(J)/DENC
   10    TTE(J,1) = PTE(J)/TEMPEC
      DO 30  I=1,NSPEC
         DENIC(I) = PNI(1,I)
         TEMPIC(I) = PTI(1,I)
      DO 20  J=1,NPROF
         TNI(J,1,I) = PNI(J,I)/DENIC(I)
         TTIX(J,1,I) = PTI(J,I)/TEMPIC(I)
   20    TTIZ(J,1,I) = PTI(J,I)/TEMPIC(I)
   30 CONTINUE
C
C  Coefficients for the cubic spline interpolation
C
            CALL CUBSPL(DPSI,TNE,NPROF)
            CALL CUBSPL(DPSI,TTE,NPROF)
C
            CALL CUBSPL(DPSI,TNI(1,1,1),NPROF)
            CALL CUBSPL(DPSI,TTIX(1,1,1),NPROF)
            CALL CUBSPL(DPSI,TTIZ(1,1,1),NPROF)
C
      IF(NSPEC.EQ.1 .OR. IPROEQ.EQ.1) RETURN
C
      DO 110  I=2,NSPEC
            CALL CUBSPL(DPSI,TNI(1,1,I),NPROF)
            CALL CUBSPL(DPSI,TTIX(1,1,I),NPROF)
            CALL CUBSPL(DPSI,TTIZ(1,1,I),NPROF)
  110 CONTINUE
C
      RETURN
      END

      SUBROUTINE METRIN
C
C   ****************************************************************
C   *                                                              *
C   *      THIS MODULUS INITIALIZES THE PARAMETRIC REPRE-          *
C   *           SENTATION OF  THE PLASMA EQUILIBRIUM CON-          *
C   *           FIGURATION.                                        *
C   *                                                              *
C   ****************************************************************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
      DIMENSION ZS(3),  ZV(3), ZL(3), ZQ(3)
      DIMENSION ZSURF(NPSIPT),   ZLAMB(NPSIPT),   ZQFAC(NPSIPT),
     +          ZCURR(NPSIPT)
C
C ************ ************ ************ ************ ************
C
C   MHD equilibrium parameters
C
         URPLAS = RPLASM*UKZERO
         RAXMAG = RMAJ(1)
         RSHIFT = RAXMAG - RTORUS
         POLASP(1) = ELONG(1)
         POLASP(2) = ELONG(NMHD)
         DTRIAN = TRIANG(NMHD)
         URTOR  = UKZERO*RAXMAG
         URSHFT = UKZERO*RSHIFT
         BAXIS  = BZERO*RTORUS/RAXMAG
C
         DPSMHD = 1./FLOAT(NMHD-1)
      DO  10  N=1,NMHD
         PSIMHD(N) = FLOAT(N-1)*DPSMHD
         BX0(N,1) = UKZERO*(RMAJ(N) - RMAJ(1))
         BX2(N,1) = TRIANG(N)
   10    BZ0(N,1) = URPLAS*PSIMHD(N)*ELONG(N)
C
            CALL CUBSPL(DPSMHD,BX0,NMHD)
C
            CALL CUBSPL(DPSMHD,BX2,NMHD)
C
            CALL CUBSPL(DPSMHD,BZ0,NMHD)
C
C   Plasma shape and dimensions
C
         PSI = 1.
         THETA = 0.
            CALL COORDS
         XMAX = UX/UKZERO
C
         THETA = PI
            CALL COORDS
         XMIN = UX/UKZERO
C
C  Specific volume of the magnetic surfaces,
C  Poloidal magnetic field and safety factor.
C
         ZDTH = PI/200.
         ZDTHX = 2.*ZDTH/3.
         ZK2 = UKZERO*UKZERO
         ZK3 = UKZERO*ZK2
      DO 130  J=2,NPROF
         PSI = PPSI(J)
         ZSINT = 0.
         ZVINT = 0.
         ZLINT = 0.
         ZQINT = 0.
         THETA = 0.
            CALL COORDS
         ZS(3) = AJ
         ZV(3) = URHS*AJ
         ZL(3) = ANTAU*ANTAU/(AJ*URHS)
         ZQ(3) = AJ/URHS
      DO 120  N=1,100
         ZS(1) = ZS(3)
         ZV(1) = ZV(3)
         ZL(1) = ZL(3)
         ZQ(1) = ZQ(3)
      DO 110  L=2,3
         THETA = THETA + ZDTH
            CALL COORDS
         ZS(L) = AJ
         ZV(L) = URHS*AJ
         ZL(L) = ANTAU*ANTAU/(AJ*URHS)
  110    ZQ(L) = AJ/URHS
         ZSINT = ZSINT + ZS(1) + 4.*ZS(2) + ZS(3)
         ZVINT = ZVINT + ZV(1) + 4.*ZV(2) + ZV(3)
         ZLINT = ZLINT + ZL(1) + 4.*ZL(2) + ZL(3)
  120    ZQINT = ZQINT + ZQ(1) + 4.*ZQ(2) + ZQ(3)
         ZSURF(J) = ZSINT*ZDTHX*DPSI/ZK2
         VOLUME(J) = 2.*PI*ZVINT*DPSI*ZDTHX/ZK3
         ZLAMB(J) = ZLINT*URTOR*ZDTHX/(2.*PI)
  130    ZQFAC(J) = ZQINT*ZDTHX/(2.*PI*URPLAS)
         VOLUME(1) = VOLUME(2)/3.
C
         ZCURR(1) = 0.
      DO 140  J=2,NPROF
  140    ZCURR(J) = ZCURR(J-1) + 0.5*(ZSURF(J-1)*PJAVG(J-1)+
     +        ZSURF(J)*PJAVG(J))
      DO 150  J=2,NPROF
         ZCURR(J) = TCURR*ZCURR(J)/ZCURR(NPROF)
         TQR(J,1) = 2.E-7*ZCURR(J)/(RPLASM*BAXIS*ZLAMB(J))
  150    PQR(J) = ZQFAC(J)/TQR(J,1)
         TQR(1,1) = 0.
         PQR(1) = PQR(2)
C
            CALL CUBSPL(DPSI,TQR,NPROF)
C
C ************ ************ ************ ************ ************
C
      RETURN
      END

      SUBROUTINE INIWFR
C
C ************ ************ ************ ************ ************
C
C  Initial wavefront: evaluation of the wavevector
C      The function PSIWFR must be provided by the user:
C      PSIWFR(THETA,ZDPSDT) = PSI(THETA) on the first wavefront;
C      ZDPSDT = d(PSI)/d(THETA)
C
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
         ZPSI   = PSASYM
         IPASYM = 0
C
   10    IPASYM = IPASYM+1
         ISTART = 1
C
      DO 30 JTH=1,NTHIN
C
         INZWFR(JTH) = 0
         THETA  = PI*THSTRT(JTH)/180.
         PSI    = PSIWFR(ZPSI,THETA,ZDPSDT)
C
            CALL COORDS
            CALL PROFIX
C
         PKPHI = ANPHI/URHS
         PKPAR = PKPHI*ACOSTQ
         PKPAR2 = PKPAR*PKPAR
C
            CALL DISPIC(1)
C
         ZPKP2 = -ANR*ANL/ANS
C
      IF(ZPKP2 .LT. 25.)  THEN
         ISTART = 0
         IF(IPASYM.GE.1)
     >      WRITE(6,9980)  PSI,THSTRT(JTH),ZPKP2
 9980 FORMAT(5X,'Starting point too close to the R-cutoff',/,
     +   5X,'Psi =',F10.4,' Theta =',F8.2,' PKP2 =',F10.2)
         GOTO 40
      END IF
C
         ZSINZ = -AJ*ZDPSDT
         ZCOSZ = ANTAU*ANTAU+AG*ZDPSDT
         ZNORMZ = SQRT(ZSINZ*ZSINZ+ZCOSZ*ZCOSZ)
         ZSINZ = ZSINZ/ZNORMZ
         ZCOSZ = ZCOSZ/ZNORMZ
C
         PKP2 = ZPKP2
         ZA = ZCOSZ*ZCOSZ + ZSINZ*ZSINZ*ACOSTQ*ACOSTQ
         ZB = ZSINZ*ACOSTQ*ASINTQ*PKPHI
C
C  Evaluation of the initial value of KPSI & KTHETA
C
         NITER = 0
C
   20    NITER = NITER+1
C
      IF(NITER.GT.20)  THEN
         ISTART = 0
         IF(IPASYM.GE.1)
     >      WRITE(6,9990) PKP2,ZDKP2,DISPH
 9990 FORMAT(5X,'Search for the initial Nperp not convergent',/,
     +   5X,'Last values: KP2 =',1P,E13.4,' DKP2 =',E13.4,
     +      '  DISPH =',E13.4,0P)
         GOTO 40
      END IF
C
         ZC = PKPHI*PKPHI*ASINTQ*ASINTQ - PKP2
         ZDISCR = ZB*ZB - ZA*ZC
C
      IF(ZDISCR.LE.0.)  THEN
         ISTART = 0
         IF(IPASYM.GE.1)
     >      WRITE(6,9970)  PSI,THSTRT(JTH),PKP2,ZPKP2
 9970 FORMAT(5X,'Starting point too close to the R-cutoff',/,
     +   5X,'Psi =',F10.4,' Theta =',F8.2,' PKP2 =',1P,E13.4,
     +      ' QF2 =',E13.4,0P)
         GOTO 40
      END IF
C
         ZKAPPA = (-ZB + SQRT(ZDISCR))/ZA
         PKPSI = -ZKAPPA*ZCOSZ
         PKTAU = -ZKAPPA*ZSINZ
         AKPSI = (AJ*PKPSI + AG*PKTAU)/ANTAU
         AKTH = ANTAU*PKTAU
C
            CALL WAVECT
            CALL DISPIC(1)
C
         ZDKP2 = DISPH/DHKP2
         PKP2 = PKP2 - ZDKP2
C
      IF(ABS(ZDKP2/PKP2).GT.0.1*ACCUR) GO TO 20
C
      IF(ISTART.EQ.0)  GOTO 40
C
         ZPSINI(JTH) = PSI
         ZTHINI(JTH) = THETA
         ZKPSIN(JTH) = AKPSI
         ZKTHIN(JTH) = AKTH
         INZWFR(JTH) = 1
C
   30 CONTINUE
C
   40 IF(ISTART.EQ.0 .AND. IPASYM.LT.5) THEN
         ZPSI = ZPSI-0.05
         GOTO 10
      END IF
C
      RETURN
      END
