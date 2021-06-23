         SUBROUTINE ANTIC
C
C =================================================================
C  Simplified antenna model - Master routine
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
C ===============================================
C Initialisation and output
C ===============================================
C
            CALL INIANT
C
      IF(IBREAK.NE.0)  GOTO 9000
C
C ===============================================
C Loop over NPHI
C ===============================================
C
      DO 100  NZ=1,NUMPHI
C
C ===============================================
C Evaluation of the surface impedance
C ===============================================
C
            CALL CYSFIC
C
      IF(IBREAK.NE.0)  GO TO 9000
C
C ===============================================
C Evaluation of the function ZB(NY,NZ)
C ===============================================
C
            CALL ZBFAC2
C
C ===============================================
C Power coupled to the partial wave NPHI
C ===============================================
C
            CALL POYNTX(1)
C
  100 CONTINUE
C
C ===============================================
C Normalisation of the power flux
C ===============================================
C
            CALL POYNTX(2)
C
C ===============================================
C Output of the antenna package
C ===============================================
C
            CALL ANTOUT
C
C ===============================================
      RETURN
C ===============================================
C
 9000    NANTS = 0
C
C ===============================================
      RETURN
C ===============================================
C
      END

         SUBROUTINE INIANT
C
C =================================================================
C  Simplified antenna model - Initialisation and output.
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
      include 'COMHPCSD.F'
C
C =================================================================
      DATA ZCOPE2 /8.0615E+07/,      ZCOPH2 /4.3905E+04 /,
     +     ZCHE   /3.5724E-11/,      ZCHH   /6.5594E-08 /
C =================================================================
      COFUN(X) = COS(X) + CIM*SIN(X)
      ZMOD2(C) = REAL(C)*REAL(C)+AIMAG(C)*AIMAG(C)
C =================================================================
C
         IBREAK = 0
         IOUTA = JOUTA
         IND = -1
         NDIM = 4
C
         PI = 3.14152965359
         CRE = CMPLX(1.,0.)
         CIM = CMPLX(0.,1.)
         CZERO = CMPLX(0.,0.)
C
C ===============================================
C  Parameters for the evaluation of profiles
C ===============================================
C
         DUAW = DISTAW/ULENGT
         DUPW = (DISTAP+DISTAW)/ULENGT
C
         ACOPE2 = ZCOPE2/(FREQCY*FREQCY)
         ACOHE = ZCHE*FREQCY
         ACOPH2 = ZCOPH2/(FREQCY*FREQCY)
         ACHIH = ZCHH*FREQCY
C
C ===============================================
C  Plasma edge
C ===============================================
C
         THETA = PI*THANTN/180.
         PSI = 1.
C
            CALL COORDS
C
            CALL PROFIX
C
         DENSX = DENS
         BTOTX = BTOT
C
            CALL DISPAN
C
         DENEDG = DENS
         BZEDGE = BTOT
         ZEX = UX
         ZEZ = UZ
C
C ===============================================
C  Asymptotic position
C ===============================================
C
         PSI = 0.9
C
   10    CONTINUE
            CALL COORDS
C
            CALL PROFIX
C
         DENSX = DENS
         BTOTX = BTOT
C
            CALL DISPAN
C
         ZQFSQ = VR*VL/VS
      IF(ZQFSQ.LE.400.)  GOTO 20
         ZDX = UX - ZEX
         ZDZ = UZ - ZEZ
         ULENGR = SQRT(ZDX*ZDX + ZDZ*ZDZ)
         ZCONF = SQRT(ZQFSQ)*ULENGR
      IF(ZCONF.GT.5.)  GO TO 30
   20    PSI = PSI - 0.05
      IF(PSI.LT.0.1)  GOTO 9100
      GOTO 10
C
   30    PSASYM = PSI
         PSIINT = 1. - PSASYM
         DNASYM = DENS
         BZASYM = BTOT
         ARASYM = VR
         ALASYM = VL
         ASASYM = VS
         ADASYM = VD
C
C ===============================================
C Values of NPHI and evaluation of JZ(NPHI)
C ===============================================
C
         MUMPHI = NUMPHI
         NUMPHI = 0
         ZPHDIF = PI*DPHASE/180.
         ZNANTS = FLOAT(NANTS)
         ZDRHF = 2.*PI/ZNANTS
      DO 40  N=1,MUMPHI
         NPH = NPHI1 + (N-1)*JUMPHI
         ZNZ = FLOAT(NPH)*ULENGT/RTORUS
         VNZ2 = ZNZ*ZNZ
         ZQFSQ = -(VNZ2-VR)*(VNZ2-VL)/(VNZ2-VS)
C
      IF(ZQFSQ.LT.25.)  GO TO 40
C
         NUMPHI = N
         NPHI(N) = NPH
         ANZ(N) = ZNZ
C
         ZPERFC = 1.
         ZARG = 0.5*(ZPHDIF + FLOAT(NPH)*ZDRHF)
         ZSNARG = SIN(ZARG)
      IF(ABS(ZSNARG).GT.1.E-7)
     +   ZPERFC = SIN(ZNANTS*ZARG)/(ZNANTS*ZSNARG)
C
         CJZ(N) = ZPERFC*CFUNJZ(NPH)
C
   40 CONTINUE
C
      IF(NUMPHI.EQ.0)  GOTO 9200
         IASYM = 1
      IF(NPHI(1).NE.0)  IASYM = 0
      IF(DPHASE.NE.0.)  IASYM = 0
      IF(THANTN.NE.0.AND.THANTN.NE.180.)  IASYM = 0
      IF(JPOLE.EQ.2.AND.EPHASE.NE.0.)  IASYM = 0
C
C ===============================================
C Representative values of NY and JY(NY)
C ===============================================
C
         ZRA = (RPLASM+DISTAP)/ULENGT
         ZNYMAX = 2.*PI*(RPLASM+DISTAP)/HEIGTH
         NYMAX = IFIX(ZNYMAX) + 1
         JUMPNY = (2*NYMAX + 1)/NYDIM + 1
         NUMNY  = (2*NYMAX + 1)/JUMPNY
         IF(2*(NUMNY/2).EQ.NUMNY)  NUMNY = NUMNY+1
C
      DO  60  N=1,NUMNY
         NNY = -NYMAX + (N-1)*JUMPNY
         ANY(N) = FLOAT(NNY)/ZRA
         ZNY = ANY(N)
   60    CJY(N) = CFUNJY(ZNY)
C
C ===============================================
C Points at which PX has to be evaluated
C ===============================================
C
      IF(NTHIN.EQ.1) THEN
         ZYMIN = 0.
         ZDY = 0.
      ELSE
         ZYMIN = HEIGTH
         ZDY = 2.*ZYMIN/FLOAT(NTHIN-1)
      END IF
      DO  70  J=1,NTHIN
         YY(J) = -ZYMIN + FLOAT(J-1)*ZDY
         THSTRT(J) = THANTN + 180.*YY(J)/(PI*(RPLASM+DISTAP))
      DO  70  N=1,NUMNY
         ZARG = ANY(N)*YY(J)/ULENGT
   70    CEXNYY(N,J) = CMPLX(COS(ZARG),SIN(ZARG))
C
C ===============================================
C  Output of the ANTIC package
C ===============================================
C
      IF(IOUTA.GT.0)  THEN
      WRITE(6,1000)
 1000 FORMAT(/,45X,'DETAILS OF THE ANTENNA EVALUATION',/)
      WRITE(6,1010)
 1010 FORMAT(/,25X,'Fourier spectrum of the antenna current',
     +     1X,'in the toroidal direction'/)
      WRITE(6,1020)  (NPHI(N),ZMOD2(CJZ(N)),N=1,NUMPHI)
 1020 FORMAT(3(5X,'Nphi =',I4,' Jz sqrd =',1P,E13.3))
C
      WRITE(6,1030)
 1030 FORMAT(/,25X,'Fourier spectrum of the antenna current',
     +     1X,'in the poloidal direction'/)
      WRITE(6,1040)  (ANY(N),ZMOD2(CJY(N)),N=1,NUMNY)
 1040 FORMAT(3(5X,'Ny =',1PE13.3,' Jy sqrd =',E13.3))
      END IF
C ===============================================
      RETURN
C ===============================================
 9100 WRITE(6,9110)
 9110 FORMAT(/,5X,'Antenna evaluation failed,',
     +   ' central density too loo to be asymptotic')
      GOTO 9999
C ===============================================
 9200 WRITE(6,9210)
 9210 FORMAT(/,5X,'Antenna evaluation failed,',
     +   ' no mode is asymptotic in the required range')
C ===============================================
 9999    IBREAK = 1
C ===============================================
      RETURN
C ===============================================
      END

         SUBROUTINE DISPAN
C
C =================================================================
C  Simplified antenna model
C  Elements of the cold plasma dielectric tensor
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
      include 'COMHPCSD.F'
C
C =================================================================
C
         OPE2 = ACOPE2*DENSX
         OHE = ACOHE/BTOTX
         OPCYE2 = OPE2*OHE*OHE
C
         VR = 1. + OPCYE2
         VL = 1. + OPCYE2
C
         ZHIH = ACHIH/BTOTX
         ZOPH2 = ACOPH2*DENSX
      DO 10  I=1,NSPEC
         OPI2(I) = ZOPH2*ATOPI(I)*ACONC(I)
         OHI(I) = ZHIH*ATZI(I)
         OPCYI2(I) = OPI2(I)*OHI(I)*OHI(I)
         VR = VR + OPCYI2(I)/(OHI(I)+1.)
   10    VL = VL - OPCYI2(I)/(OHI(I)-1.)
C
         VS = 0.5*(VR + VL)
         VD = 0.5*(VR - VL)
C
C ===============================================
C
      RETURN
      END

         SUBROUTINE ANTOUT
C
C =================================================================
C  Simplified antenna model - output.
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
C =================================================================
C
      WRITE(6,2000)
 2000 FORMAT(/,45X,'ANTENNA PARAMETERS')
C
      WRITE(6,2010)  NANTS,DPHASE
 2010 FORMAT(  20X,'Number of antennas              =',I10,/,
     +         20X,'Phase diff. between antennas    =',F10.1)
      WRITE(6,2020)  JPOLE
 2020 FORMAT(20X,' Number of elements in each antenna',I10)
C
      IF(JALIM.EQ.1)  THEN
         WRITE(6,2030)
 2030 FORMAT(20X,'Feeders up and down, short central')
      ELSE
         WRITE(6,2040)
 2040 FORMAT(20X,'Shorts up and down, feeder central')
      END IF
C
      WRITE(6,2110) DISTAP,DISTAW,THANTN,WIDTH,HEIGTH
 2110 FORMAT(/,20X,'Distance antenna-plasma         =',F10.3,' m',/
     +         20X,'Distance antenna-wall           =',F10.3,' m',/
     +         20X,'Position of the antenna         =',F10.2,' deg.',/
     +         20X,'Width of the antenna conductor  =',F10.3,' m',/
     +         20X,'Length of the antenna conductor =',F10.3,' m')
      IF(JPOLE.GE.2)
     >   WRITE(6,2120) WGAP,EPHASE
 2120 FORMAT(  20X,'Gap between conductors          =',F10.3,' m',/
     +         20X,'Phase diff. between conductors  =',F10.1)
         ZLENGR = ULENGR*ULENGT
      WRITE(6,2140)  ZLENGR,DNASYM,DENEDG,BZEDGE
 2140 FORMAT(/,20X,'Length of the near field region =',F10.3,' m',/
     +         20X,'Asymptotic point density        =',E14.4,' cm-3',/
     +         20X,'Edge density                    =',E14.4,' cm-3',/
     +         20X,'Edge magnetic field             =',F10.3,' T')
      WRITE(6,2150) PWEDGE
 2150 FORMAT(/,20X,'Power absorbed in the near field region =',
     +   E13.4,' MW',/)
C
      IF(IPRINT.EQ.0.OR.PWEDGE.EQ.0.)  RETURN
C
      WRITE(6,2170)  (NPHI(N),PWNEAR(N),N=1,NUMPHI)
 2170 FORMAT(2(5X,'Nphi =',I3,' Abs. power =',E13.4,5X))
C
C ===============================================
C
      RETURN
      END

      SUBROUTINE CYSFIC
C
C =================================================================
C  Simplified antenna model - Surface impedance of the plasma
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
      include 'COMHPCSD.F'
C
C =================================================================
      EXTERNAL FUANTY,OUTANT
C =================================================================
C
         VNZ2 = ANZ(NZ)*ANZ(NZ)
C
         ZNR = VNZ2 - ARASYM
         ZNL = VNZ2 - ALASYM
         ZNS = VNZ2 - ASASYM
         ZQFSQ = -ZNR*ZNL/ZNS
         ANUCOL = 1.E-3
         ZSTEP = PSIINT/10.
C
      DO 100  NY=1,NUMNY
         VANY = ANY(NY)
         VANY2 = VANY*VANY
         IREANT = 0
         ADFLUX(NY) = 0.
C
         ZQ2 = -VANY2 + ZQFSQ
C
      IF(ZQ2)  10,20,30
C
   10    CZYAS = CIM*(ZNR*ZNL/(ZNS*SQRT(-ZQ2)-VANY*ADASYM))
         CY(NY) = CZERO
      GO TO 100
C
   20    CZYAS = -CIM*(ZNR*ZNL/(VANY*ADASYM))
         CY(NY) = CZERO
      GO TO 100
C
   30    CZYAS = -ZNR*ZNL/(ZNS*SQRT(ZQ2)*CRE-VANY*ADASYM*CIM)
C
   40    Y(1) = 1.
         Y(2) = 0.
         Y(3) = REAL(CZYAS)
         Y(4) = AIMAG(CZYAS)
C
         PRMT(1) = 0.
         PRMT(2) = PSIINT
         PRMT(3) = ZSTEP
         PRMT(4) = 1.E-3
      DO 50  MM=1,4
   50    DY(MM) = 0.25
C
            CALL HPCSD(FUANTY,OUTANT)
C
      IF(IBREAK.NE.0)  THEN
C
      IF(ANUCOL.GT.0.05)  GOTO 9000
C
         IBREAK = 0
         ANUCOL = 2.*ANUCOL
         ZSTEP = ZSTEP/2.
         GOTO 40
      END IF
C
         CZEY = CMPLX(Y(1),Y(2))
         CZBZ = CMPLX(Y(3),Y(4))
         CY(NY) = CZBZ/CZEY
C
      IF(IREANT.NE.0)
     >   ADFLUX(NY) = REAL(CONJG(CZEY)*CZBZ) - REAL(CZYAS)
C
  100    CYASIM(NY) = CZYAS
C
      IF(IOUTA.GT.0)  THEN
         WRITE(6,1010)  NZ,ANZ(NZ),ZQFSQ
 1010 FORMAT(1H0,/,20X,'MODE N. =',I3,5X,'N// =',F7.3,
     +   5X,'N_perp sqrd (asympt.) =',1P,E13.4,/)
         WRITE(6,1020)  (ANY(NY),CYASIM(NY),CY(NY),NY=1,NUMNY)
 1020 FORMAT(5X,'Ny =',1PE13.4,'   CYAS =',2E13.4,'   CY =',2E13.4)
      END IF
C
C ===============================================
      RETURN
C ===============================================
C
 9000    IBREAK = 1
      WRITE(6,9999) JHLF,HLAST,VNS,ANUCOL
 9999 FORMAT(1H0,5X,'Antenna evaluation failed, JHLF =',I4,
     +   '  last step =',E13.4,/,11X,'Nzsq - S =',E13.4,' Nu_coll =',
     +   E13.4)
C
C ===============================================
C
      RETURN
      END

      SUBROUTINE FUANTY(X)
C
C =================================================================
C  Simplified antenna model
C  Diff. equations for the evaluation of the surface impedance
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
      include 'COMHPCSD.F'
C
C =================================================================
      DIMENSION ZDY(4)
C =================================================================
C
         ZX = X/PSIINT
         ZX1 = 1.-ZX
C
         DENSX = ZX*DENEDG + ZX1*DNASYM
         BTOTX = ZX*BZEDGE + ZX1*BZASYM
C
            CALL DISPAN
C
         VNR = VNZ2 - VR
         VNL = VNZ2 - VL
         VNS = VNZ2 - VS
C
         ZANID = VANY*VD
         ZANYS = VANY2 + VNS
         ZNRL = VNR*VNL
C
         ZDY(1) = ZANID*Y(1) - ZANYS*Y(4)
         ZDY(2) = ZANID*Y(2) + ZANYS*Y(3)
         ZDY(3) = -ZANID*Y(3) + ZNRL*Y(2)
         ZDY(4) = -ZANID*Y(4) - ZNRL*Y(1)
C
         ZEPS = ANUCOL*OPCYI2(MAINSP)
      IF(ABS(VNS).LT.10.*ZEPS)  GOTO 20
C
         ZMULT = -ULENGR/VNS
      DO 10  I=1,4
   10    DY(I) = ZMULT*ZDY(I)
C
      RETURN
C
   20    IREANT = 1
         ZMULT = -ULENGR/(VNS*VNS + ZEPS*ZEPS)
         DY(1) = ZMULT*(VNS*ZDY(1) - ZEPS*ZDY(2))
         DY(2) = ZMULT*(VNS*ZDY(2) + ZEPS*ZDY(1))
         DY(3) = ZMULT*(VNS*ZDY(3) - ZEPS*ZDY(4))
         DY(4) = ZMULT*(VNS*ZDY(4) + ZEPS*ZDY(3))
C
      RETURN
C
      END

         SUBROUTINE OUTANT(X)
C
C =================================================================
C  Simplified antenna model
C  Output of the integration routine
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
      include 'COMHPCSD.F'
C
C =================================================================
C
      IF(JHLF.LT.11) RETURN
         PRMT(5) = 1
         IBREAK = 1
C
C ===============================================
C
      RETURN
      END

      SUBROUTINE ZBFAC2
C
C =================================================================
C  Simplified antenna model
C  Form factor due to the wall-plasma guide, ZB(NY,NY)
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
C =================================================================
C
         ZANU2 = VNZ2 - 1.
C
      DO 100  NY=1,NUMNY
         ZANUX2 = VNZ2 + ANY(NY)*ANY(NY) - 1.
C  NX SQRD = 1
      IF(ZANUX2.EQ.0.) THEN
         CZEDG(NY) = CZERO
      ELSE
C  NX SQRD > 1
         IF(ZANUX2.GT.0.)  THEN
            ZANUX = SQRT(ZANUX2)
            ZSAW = -ZANUX*SINH(ZANUX*DUAW)
            ZANUW = ZANUX*DUPW
            CDIV = ZANU2*COSH(ZANUW)-CIM*CY(NY)*ZANUX*SINH(ZANUW)
C  NX SQRD < 1
         ELSE
            ZANUX = SQRT(-ZANUX2)
            ZSAW = ZANUX*SIN(ZANUX*DUAW)
            ZANUW = ZANUX*DUPW
            CDIV = ZANU2*COS(ZANUW)+CIM*CY(NY)*ZANUX*SIN(ZANUW)
         END IF
         CZEDG(NY) = CJZ(NZ)*CJY(NY)*ZSAW/CDIV
      END IF
  100 CONTINUE
C
C ===============================================
      IF(IOUTA.NE.0)  THEN
         WRITE(6,1000)
 1000 FORMAT(/,20X,'Surface factor at the plasma edge',/)
         WRITE(6,1010)  (ANY(NY),CZEDG(NY),NY=1,NUMNY)
 1010 FORMAT(2(5X,'Ny =',F8.3,' CZ =',1P,2E13.3,0P))
      END IF
C ===============================================
C
      RETURN
      END

      SUBROUTINE POYNTX(IENTRY)
C
C =================================================================
C  Simplified antenna model
C  Power coupled to the partial wave NPHI
C  Power flux distribution along the antenna P(NZ,Y)
C =================================================================
C
      IMPLICIT COMPLEX (C)
C
      include 'COMMON.F'
C
      include 'COMANTIC.F'
C
C =================================================================
      ZMOD2(C) = REAL(C)*REAL(C)+AIMAG(C)*AIMAG(C)
C =================================================================
C
      GOTO (100,200) IENTRY
C
  100    PWCPL(NZ) = 0.
         PWNEAR(NZ) = 0.
      DO 110  L=1,NUMNY
         ZB2 = ZMOD2(CZEDG(L))
         PWCPL(NZ) = PWCPL(NZ) + ZB2*REAL(CY(L))
  110    PWNEAR(NZ) = PWNEAR(NZ) + ZB2*ADFLUX(L)
C
      DO 130 J=1,NTHIN
         CZEY = 0.
         CZBZ = 0.
      DO 120  N=1,NUMNY
         CTERM = CZEDG(N)*CEXNYY(N,J)
         CZEY = CZEY + CTERM
  120    CZBZ = CZBZ + CTERM*CYASIM(N)
  130    POYNT(NZ,J) = REAL(CZEY*CONJG(CZBZ))
C
         ZNORM = 0.
      DO 140  J=1,NTHIN
  140    ZNORM = ZNORM + POYNT(NZ,J)
C
      IF(ZNORM.EQ.0.)  RETURN
C
      DO 150  J=1,NTHIN
  150    POYNT(NZ,J) = POYNT(NZ,J)/ZNORM
C
      IF(IOUTA.EQ.0)  RETURN
C
      WRITE(6,1010) NPHI(NZ)
 1010 FORMAT(/,20X,'Theta power flux distribution for Nphi =',I3)
      WRITE(6,1020) (YY(J),POYNT(NZ,J),J=1,NTHIN)
 1020 FORMAT(4(3X,'Y=',1P,E11.3,' Px=',E11.3,0P))
C
      RETURN
C
C =================================================================
C  Normalisation of the NPHI power spectrum
C =================================================================
C
  200    ZW = 2.
      IF(IASYM.EQ.0)  ZW = 1.
C
         ZPWNOR = PWCPL(1) + PWNEAR(1)
      DO  210  N=2,NUMPHI
         PWCPL(N) = ZW*PWCPL(N)
         PWNEAR(N) = ZW*PWNEAR(N)
  210    ZPWNOR = ZPWNOR + PWCPL(N) + PWNEAR(N)
C
      IF(ZPWNOR.EQ.0.)  GOTO 9000
C
         PWEDGE = 0.
         ZPWNOR = POWER/ZPWNOR
      DO 220  N=1,NUMPHI
         PWNEAR(N) = PWNEAR(N)*ZPWNOR
         PWEDGE = PWEDGE + PWNEAR(N)
  220    PWCPL(N) = PWCPL(N)*ZPWNOR
C
C ===============================================
C  Output
C ===============================================
C
      IF(IOUTA.EQ.0)  RETURN
      WRITE(6,1000)
 1000 FORMAT(/,40X,'NORMALIZED POWER SPECTRUM',/)
      WRITE(6,1030)  (NPHI(N),ANZ(N),PWCPL(N),PWNEAR(N),N=1,NUMPHI)
 1030 FORMAT(5X,'Nphi =',I4,5X,'N par. =',F10.4,5X,
     +      'Power: radiated =',E14.4,' absorbed =',E14.4)
C
C ===============================================
      RETURN
C ===============================================
C
 9000    IBREAK = 1
      WRITE(6,9990)
 9990 FORMAT(/,5X,'Antenna evaluation failed,',
     +   ' no power coupled to the spectral range')
C
C ===============================================
C
      RETURN
      END
