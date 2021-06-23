       SUBROUTINE DISPIC(IENTRY)
C
C ************ ************ ************ ************ ************
C  Ion Cyclotron Ray Tracing code RAYIC
C  Last modified: January 29, 1996
C ************ ************ ************ ************ ************
C
C  Dispersion relation, derivatives, and absorption for IC waves
C
C      IENTRY =  0 : Parameters for the evaluation of GAMMA
C                1 : dispersion relation H and dH/dkperp
C                2   same and all derivatives of H
C                    If (IPWDP.NE.0) GAMMA is reevaluated.
C
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
      DIMENSION PZR(NSPION),      PZLW(NSPION),
     +          PZRW(NSPION),     PZTW(NSPION)
C
      DIMENSION PZLDEL(NSPION),   PZLWDL(NSPION)
C
      DIMENSION OCYDIV(2,NSPION), ZQL(2,NSPION),    ZQR(2,NSPION),
     +          ZQLDOM(2,NSPION), ZQRDOM(2,NSPION), ZQLDKZ(2,NSPION),
     +          ZQALFA(NSPION)
C
C ----------------------------------------------------------------
C
      DATA ZCOPE2 /8.0615E+07/,      ZCOPH2 /4.3905E+04 /,
     +     ZCHE   /3.5724E-11/,      ZCHH   /6.5594E-08 /,
     +     ZCXE   /15.9844   /,      ZCXI   /684.945    /,
     +     ZCVTE2 /3.91389E-3/,      ZCVTH2 /2.13152E-6 /
C
      DATA ZSQRTH /0.7071607  /,     ZDSQPI /1.7724538509/,
     +     ZMIN   /1.E-20/
C
C ************ ************ ************ ************ ************
C
         IPWDPE = 0
         IPWDPI = 0
         GAMMA = 0.
         GAME(1) = 0.
         GAME(2) = 0.
         GAME(3) = 0.
      DO 10  I=1,NSPEC
         GAMI(I) = 0.
   10    IRES(I) = 0
C
C ************ ************ ************ ************ ************
C  Dimensionless parameters
C ************ ************ ************ ************ ************
C
         ZKPAR = PKPAR
      IF(ABS(ZKPAR).LT.ZMIN)  ZKPAR = ZMIN
C
         OPE2   = ZCOPE2*DENS/(FREQCY*FREQCY)
         OHE    = ZCHE*FREQCY/BTOT
         OPCYE2 = OPE2*OHE*OHE
         OMCYE  = 1./OHE
         ZVTHE2 = ZCVTE2*TEMPE
         BETAE  = OPCYE2*ZVTHE2
C
         XE     = ZCXE/(ZKPAR*SQRT(TEMPE))
      IF(ABS(XE).LE.10.)  THEN
         IPWDPE = 1
            CALL ZETA(XE,ZFE,ZDE,ZWE)
         ZEXPTT = ZDSQPI*ABS(XE)*EXP(-XE*XE)
         ZEXPE  = 2.*XE*XE*ZEXPTT
      ELSE
         ZFE    = 1.
         ZDE    = 1.
         ZWE    = 0.
         ZEXPE  = 0.
         ZEXPTT = 0.
      END IF
C
         ZOPH2  = ZCOPH2/(FREQCY*FREQCY)
         ZXIH   = ZCXI/ZKPAR
         ZHIH   = ZCHH*FREQCY/BTOT
      DO 20  I=1,NSPEC
         OPI2(I)   = ZOPH2*ATOPI(I)*DENSI(I)
         OHI(I)    = ZHIH*ATZI(I)
         OPCYI2(I) = OPI2(I)*OHI(I)*OHI(I)
         OMCYI(I)  = 1./OHI(I)
         OXI(I)    = ZXIH*SQRT(ATM(I)/TEMPIZ(I))
         ZVTHI2(I) = ZCVTH2*TEMPIX(I)/ATM(I)
         BETAI(I)  = OPCYI2(I)*ZVTHI2(I)
         ZQALFA(I) = 0.5*BETAI(I)
      DO 20  IH=1,MAXHRM
         OCYDIV(IH,I) = 1.-FLOAT(IH)*OMCYI(I)
   20    XI(IH,I)     = OXI(I)*OCYDIV(IH,I)
C
      DO 30  I=1,NSPEC
      DO 30  IH=1,MAXHRM       
      IF(ABS(XI(IH,I)).LE.10.)  THEN
         IPWDPI = 1
         IRES(I)  = IH
         ZXI = XI(IH,I)
         ZEXPI(IH,I) = ZDSQPI*ABS(OXI(I))*EXP(-ZXI*ZXI)
      ELSE
         ZEXPI(IH,I) = 0.
      END IF
   30 CONTINUE
C
         IPWDP = IPWDPI + IPWDPE
C
      IF(IENTRY.EQ.0)  RETURN
C
C ****************************************************************
C  Elements of the dielectric tensor
C ****************************************************************
C
         ZQRE   = OHE/(OHE - 1.)
         ZQLE   = OHE/(OHE + 1.)
C
      DO 40  I=1,NSPEC
      DO 40  IH=1,2
C
         ZDIV   = OHI(I)/(OHI(I) + FLOAT(IH))
         ZQR(IH,I) = ZDIV
         ZQRDOM(IH,I) = -FLOAT(IH)*OMCYI(I)*ZDIV*ZDIV
C
      IF(ABS(OCYDIV(IH,I)) .GT. 1.E-12)  THEN
         ZDIV   = 1./OCYDIV(IH,I)
      ELSE
         ZDIV   = 1.E12
      END IF

         ZXI    = XI(IH,I)
      IF(ABS(ZXI) .LE. 10.)  THEN
            CALL ZETA(ZXI,ZF1,ZD1,ZW1)
         ZQL(IH,I)    = ZDIV*ZF1
         ZQLDOM(IH,I) = FLOAT(IH)*OMCYI(I)*ZDIV*ZDIV*ZD1
         ZQLDKZ(IH,I) = -ZDIV*ZW1
      ELSE
         ZQL(IH,I)    = ZDIV
         ZQLDOM(IH,I) = FLOAT(IH)*OMCYI(I)*ZDIV*ZDIV
         ZQLDKZ(IH,I) = 0.
      END IF
C
   40 CONTINUE
C
         PZRE   = -OPE2*ZQRE
         PZLE   = -OPE2*ZQLE
         PZTWE  = 0.5*BETAE*ZFE
C
         AR     = 1. + PZRE
         AL     = 1. + PZLE
         ARWARM = 0.
         ALWARM = 0.
         ATWARM = PZTWE
C
      DO 50  I=1,NSPEC
         PZR(I) = -OPI2(I)*ZQR(1,I)
         AR     = AR + PZR(I)
         PZL(I) = -OPI2(I)*ZQL(1,I)
         AL     = AL + PZL(I)
         PZRW(I) = ZQALFA(I)*ZQR(2,I)
         PZLW(I) = ZQALFA(I)*ZQL(2,I)
         PZTW(I) = ZQALFA(I)
         ARWARM = ARWARM + PZRW(I)
         ALWARM = ALWARM + PZLW(I)
   50    ATWARM = ATWARM + PZTW(I)
C
         AS     = 0.5*(AR + AL)
         ANL    = PKPAR2 - AL
         ANR    = PKPAR2 - AR
         ANS    = PKPAR2 - AS
C
C ***************************************************************
C  Dispersion relation
C ***************************************************************
C
         ZCWL  = ATWARM + ALWARM + 0.5
         ZCWR  = ATWARM + ARWARM + 0.5
C
         HA = 0.5*(ALWARM + ARWARM)
         HB = ANR*ZCWL + ANL*ZCWR
         HC = ANR*ANL
C
         DISPH = (HA*PKP2 + HB)*PKP2 + HC
C
C  k-perp derivative
C
         DHKP2 = 2.*HA*PKP2 + HB
C
      IF(IENTRY .EQ. 1)  RETURN
C
C ***************************************************************
C  Derivatives
C ***************************************************************
C
C  k-parallel derivative
C
         ZLDEL  = 0.
         ZLWDEL = 0.
      IF(PKPAR.NE.0.)  THEN
         DO 60  I=1,NSPEC
            PZLDEL(I) = -OPI2(I)*ZQLDKZ(1,I)
            PZLWDL(I) = ZQALFA(I)*ZQLDKZ(2,I)
            ZLDEL  = ZLDEL + PZLDEL(I)
   60       ZLWDEL = ZLWDEL + PZLWDL(I)
         ZTWDEL = -0.5*BETAE*ZWE
      ELSE
         DO 70  I=1,NSPEC
            PZLDEL(I) = 0.
   70       PZLWDL(I) = 0.
         ZTWDEL = 0.
      END IF
C
         ZAKZ = 0.5*ZLWDEL/ZKPAR
         ZBKZ = (ANR*ZLWDEL + 2*ANS*ZTWDEL -
     +              ZLDEL*ZCWR)/ZKPAR + 2.*PKPAR*(ZCWR + ZCWL)
         ZCKZ = -ANR*ZLDEL/ZKPAR + 4.*PKPAR*ANS
         DHKZ = (ZAKZ*PKP2 + ZBKZ)*PKP2 + ZCKZ
C
C  Magnetic field derivative
C
         ZROM   = -OPE2*OMCYE*ZQRE*ZQRE
         ZLOM   = OPE2*OMCYE*ZQLE*ZQLE
         ZRWOM  = 0.
         ZLWOM  = 0.
C
      DO 80  I=1,NSPEC
         ZROM   = ZROM - OPI2(I)*ZQRDOM(1,I)
         ZLOM   = ZLOM - OPI2(I)*ZQLDOM(1,I)
         ZRWOM  = ZRWOM + ZQALFA(I)*ZQRDOM(2,I)
   80    ZLWOM  = ZLWOM + ZQALFA(I)*ZQLDOM(2,I)
C
         ZRWOM  = -2.*ARWARM + ZRWOM
         ZLWOM  = -2.*ALWARM + ZLWOM
         ZTWOM  = -2.*ATWARM
C
         ZAOM   = 0.5*(ZRWOM+ZLWOM)
         ZBOM   = -(ZROM*ZCWL + ZLOM*ZCWR) +
     +              ANR*ZLWOM + ANL*ZRWOM + 2.*ANS*ZTWOM
         ZCOM   = -(ANR*ZLOM + ANL*ZROM)
         DHOM   = (ZAOM*PKP2 + ZBOM)*PKP2 + ZCOM
C
C  Density and temperature derivatives
C
         ZBN    = -(PZLE*ZCWR + PZRE*ZCWL)
         ZCN    = -(PZLE*ANR + PZRE*ANL)
         DHNE   = ZBN*PKP2 + ZCN
C
         ZTWTE  = PZTWE + ZTWDEL
         DHTE = 2*ANS*ZTWTE*PKP2
C
      DO 90  I=1,NSPEC
C
         ZATX    = 0.5*(PZLW(I) + PZRW(I))
         ZBTX    = ANR*PZLW(I) + ANL*PZRW(I) + 2.*ANS*PZTW(I)
         DHTX(I) = (ZATX*PKP2 + ZBTX)*PKP2
C
         ZBN     = -(PZL(I)*ZCWR + PZR(I)*ZCWL) + ZBTX
         ZCN     = -(PZL(I)*ANR + PZR(I)*ANL)
         DHNI(I) = (ZATX*PKP2 + ZBN)*PKP2 + ZCN
C
         ZATZ    = 0.5*PZLWDL(I)
         ZBTZ    = ANR*PZLWDL(I) - PZLDEL(I)*ZCWR
         ZCTZ    = -ANR*PZLDEL(I)
         DHTZ(I) = -0.5*((ZATZ*PKP2 + ZBTZ)*PKP2 + ZCTZ)
C
   90 CONTINUE
C
      IF(IPWDP .EQ. 0)  RETURN
C
C ****************************************************************
C  Coefficients of the power transport equation.
C ****************************************************************
C
C  Polarization unit vectors (rotating frame)
C  E_z is evaluated as if x_e >> 1. No other choice is meaningful
C     within the eikonal approximation.
C
         EPOL(1) = AR - PKP2*ARWARM - PKPAR2
         EPOL(2) = -(AL - PKP2*ALWARM - PKPAR2)
         EPOLX = ZSQRTH*(EPOL(1) + EPOL(2))
         EPOLY = ZSQRTH*(EPOL(1) - EPOL(2))
         EPOL(3) = PKP*PKPAR*(EPOLX + BETAE*EPOLY/OHE)/(PKP2 + OPE2)
C
          ZENORM = 0.
       DO 110  I=1,3
  110    ZENORM = ZENORM + EPOL(I)*EPOL(I)
         ZENORM = SQRT(ZENORM)
      DO 120  I=1,3
  120    EPOL(I) = EPOL(I)/ZENORM
         EPOLX = EPOLX/ZENORM
         EPOLY = EPOLY/ZENORM
C
C  GAMMAS (Normalized absorption constants)
C
            CALL SGAMMA
C
C ************ ************ ************ ************ ************
C
         RETURN
         END

         SUBROUTINE SGAMMA
C
C ************ ************ ************ ************ ************
C  Ion Cyclotron Ray Tracing code RAYIC
C  Normalized absorption constants
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
C  Poynting vector projection along the ray
C
         EPEPSQ = EPOL(1)*EPOL(1)
         EPEMSQ = EPOL(2)*EPOL(2)
         EPEZSQ = EPOL(3)*EPOL(3)
         EPEYSQ = EPOLY*EPOLY
C
         EPEPM  = EPOL(1)*EPOL(2)
         EPEYZ  = EPOLY*EPOL(3)
C
         POYNTN = PKP*(ARWARM*EPEMSQ + ALWARM*EPEPSQ +
     +              (1. + 2.*ATWARM)*EPEYSQ)
C
      IF(IPWDPE.NE.0)  THEN
         ZLBE = 0.5*PKP2*ZVTHE2*OHE*OHE
         GAME(1) = PKP2*EPEYSQ*BETAE*ZEXPTT
         GAME(2) = OPE2*EPEZSQ*(1. - ZLBE)*ZEXPE
         GAME(3) = -PKP*PKPAR*EPEYZ*BETAE*ZEXPE/OHE
      END IF
C
         GAME(4) = GAME(1) + GAME(2) + GAME(3)
         GAMMA  = GAME(4)
C
      IF(IPWDPI.NE.0)  THEN
         DO 140  I=1,NSPEC
         IF(IRES(I).NE.0)  THEN
            IH = IRES(I)
            ZLAMB = 0.5*PKP2*ZVTHI2(I)*OHI(I)*OHI(I)
C
            CALL IBESSN(ZLAMB,ZIBN,ZIBNM1,ZIBNP1,IH)
C
            ZIH = FLOAT(IH)
            ZIBESP = (ZIH-ZLAMB)*ZIBNM1 + ZLAMB*ZIBN
            ZIBESM = -(ZIH+ZLAMB)*ZIBNP1 + ZLAMB*ZIBN
            ZIBESQ = ZLAMB*(ZIBNM1+ZIBNP1 - 2.*ZIBN)
            GAMI(I) = OPI2(I)*ZEXPI(IH,I)*(EPEPSQ*ZIBESP +
     +           EPEMSQ*ZIBESM + EPEPM*ZIBESQ)
            GAMMA  = GAMMA + GAMI(I)
         END IF
  140    CONTINUE
      END IF
C
C ************ ************ ************ ************ ************
C
      RETURN
      END

      SUBROUTINE ABSORB(IENTRY)
C
C ************ ************ ************ ************ ************
C  Ion Cyclotron Ray Tracing code RAYIC
C  Power absorbed along the ray
C  ABSORB is called by CYCRES, OUTPUT, TIHRES.
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
      DIMENSION       ZWGT(6,5),      ZGAMI(NSPION)
C
C ************ ************ ************ ************ ************
C
      DATA ZWGT /0.576116885,  0.211941558,  0.000000000,
     +           0.0        ,  0.0        ,  0.0        ,
     +           0.303641225,  0.236476024,  0.111703364,
     +           0.0        ,  0.0        ,  0.0        ,
     +           0.207995415,  0.186122475, 0.133362581,
     +           0.076517237,  0.0        ,  0.0        ,
     +           0.158434610,  0.148835542,  0.123388998,
     +           0.0902733194, 0.058284836,  0.0        ,
     +           0.128015356,  0.122995802,  0.109087491,
     +           0.089313283,  0.067501527,  0.047094218/
C
C ************ ************ ************ ************ ************
C
         ZPWABS = POWOLD - PWX
C
C  Check that PSI has been advanced. This should also
C  eliminate power creation.
C
      IF(ZPWABS.LE.1.E-12)  RETURN
C
         POWOLD = PWX
C
C  Electrons
C
         ZGAME  = ZPWABS*GAME(4)/GAMMA
         PWTOTE = PWTOTE + ZGAME
C
C  Bernstein waves
C
      IF(ITWIHR.EQ.1)  PWIBWS = PWIBWS + PWMODC
C
C  Ions
C
         ZRES = 0.
      DO 10  I=1,NSPEC
         ZGAMI(I) = ZPWABS*GAMI(I)/GAMMA
         ZRES = ZRES + ZGAMI(I)
   10    PWTOTI(I) = PWTOTI(I) + ZGAMI(I)
C
C  Power deposition profiles
C
         L = LEFT
      IF(HH.GT.0.5*DPSI)  L = LEFT+1
C
      IF(ISMOOT.EQ.0 .OR. ITWIHR+KRES.EQ.0)  THEN
C
C  Power deposition without smoothing
C
         PWE(L) = PWE(L) + ZGAME/VOLUME(L)
         DO 20  I=1,NSPEC
   20       PWI(L,I) = PWI(L,I) + ZGAMI(I)/VOLUME(L)
C
         GOTO 500
C
      END IF
C
C  Smoothing option
C
      IF(ITWIHR.EQ.0)
     >   PWE(L) = PWE(L) + ZGAME/VOLUME(L)
C
         ZDZ = 6.28/(UKZERO*PKP)
         N = NPROF*INT(ZDZ/RPLASM) + 1
         KX = IFIX((XMAX-XU)/DXIEQ) + 1
         KPSMIN = IXIEQ(KX)
      IF(N.GT.5)  N=5
         N1 = N+1
C
      IF(ITWIHR.NE.0)
     >   PWE(L) = PWE(L) + ZWGT(1,N)*ZGAME/VOLUME(L)
C
      DO 110  I=1,NSPEC
  110    PWI(L,I) = PWI(L,I) + ZWGT(1,N)*ZGAMI(I)/VOLUME(L)
C
      DO 130  JP=2,N1
         ISIGN = -1
      DO 130  IS=1,2
         ISIGN = -ISIGN
         J = L + ISIGN*(JP-1)
         IF(J.LT.KPSMIN)  J = 2*KPSMIN-J
C
         IF(J.LE.NPROF)  THEN
            IF(ITWIHR.NE.0)
     >         PWE(J) = PWE(J) + ZWGT(JP,N)*ZGAME/VOLUME(J)
               DO 120  I=1,NSPEC
  120             PWI(J,I) = PWI(J,I) + ZWGT(JP,N)*ZGAMI(I)/VOLUME(J)
         ELSE
            PWEDGE = PWEDGE + ZWGT(JP,N)*(ZGAME + ZRES)
         END IF
C
  130 CONTINUE
C
C  Check wether a resonance is too close
C
  500 IF(IENTRY.EQ.0)  RETURN
C
         KRES = 0
      IF(ZRES.EQ.0.)  RETURN
C
      DO 540  I=1,NSPEC
         J = IRES(I)
      IF(J.EQ.0.OR.J.GT.2) GO TO 540
      IF(ABS(XI(J,I)).GE.XIJUMP)  GO TO 540
      IF(AKXINI*XI(J,I)*PKPAR.GT.0.)  GO TO 540
C
      GO TO (510,520), J
C
  510    ZCONF = OPI2(I)*OXI(I)/AL
      GO TO 530
C
  520    ZCONF = BETAI(I)*OXI(I)/ALWARM
C
  530 IF(ABS(ZCONF).LT.0.01)  GO TO 540
C
         KRES = I
      GO TO 550
C
  540 CONTINUE
C
      RETURN
C
  550 IF(PRMT(5).EQ.0.)  PRMT(5) = -1.
C
      RETURN
      END
