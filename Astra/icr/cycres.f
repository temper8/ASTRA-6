      SUBROUTINE CYCRES(X)
C
C ************ ************ ************ ************ ************
C    Ion cyclotron ray-tracing code RAYIC
C    Last modified: January 26, 1996
C ************ ************ ************ ************ ************
C    Reinitialization of the ray after a break
C    at an ion cyclotron resonance.
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
      DIMENSION  ZPSI(51),     ZTHETA(51)
      DIMENSION  ZPOLY1(3),    ZPOLY2(3)
C
      DATA   ZSQRTH /0.7071607  /,    ZCXI   /684.945    /
C
C ************ ************ ************ ************ ************
C
C  Meaning of ZPARM5:
C
C   0 - Reinitialization successful
C
C   1 - k// changes sign within the layer or
C       Search for new k// not convergent.
C
C   2 - Search for new NX not converged; or
C       PKP2 too small; or
C       k// jump too large; or
C       group velocity vanishing at exit point;
C
C  If ZPARM5 = 2 the subroutine TIHRES (ion-ion resonance)
C  is called, after which ZPARM5 is set to 1 again.
C
         ZPARM5 = 0.
C
         JJRES = IRES(KRES)
         XU = UX/UKZERO
         ZU = UZ/UKZERO
      WRITE(6,1010)  NPHI(INPHI),THSTRT(ITH),JJRES,KRES,EIKON,XU,ZU,PWX,
     +                     PKPAR,PKP,XI(JJRES,KRES)
 1010 FORMAT(/,' Nphi =',I3,' The ray from theta =',F7.2,/,
     +  5X,'crosses the',I2,
     +               ' cycl. harm. of sp.',I2,' at phase',F8.2,/,
     +  5X,'Entry point: X =',F8.3,' Z =',F8.3,' power =',1P,E12.3,/,
     +  5X,'N par =',E12.3,' N perp =',E12.3,' x(i) =',E12.3,0P)
C
C  Store the old values of coordinates and wavevector
C
         ZUX     = UX
         ZUZ     = UZ
         ZKPAR1 = PKPAR
         ZPKP1  = PKP
         POYNT1 = POYNTN
         ZDIVY1 = ZDIVY
         ZNUMY1 = ZNUMY
         ZKPAR2 = PKPAR
         ZPKP2  = PKP
         POYNT2 = POYNTN
         ZDIVY2 = ZDIVY
         ZNUMY2 = ZNUMY
      DO 10  I=1,3
         ZPOLY1(I) = EPOL(I)
   10    ZPOLY2(I) = EPOL(I)
         ZDPNEW = DY(5)
C
C  Direction of the ray at the break
C
         ZQX = DXPSI*DY(2) + DXTH*DY(4)
         ZQZ = DZPSI*DY(2) + DZTH*DY(4)
C
C  Required jump in the horizontal direction
C
         ZOMKVT = ZCXI*SQRT(ATM(KRES)/(TEMPIZ(KRES)*PKPAR2))
         ZADQX = ZQX/ABS(ZQX)
         ZJX = 2.*XIJUMP*URHS*ZADQX/ZOMKVT
         ZJZ = ZJX*ZQZ/ZQX
C ====== test begs
C     write(6,7001)  ZQX,ZQZ,ZJX,ZJZ
C7001 format(' Ray direction at break:',1P,2E13,4,/,
C    +   ' Required jump in x, z:',2E13.4)
C ====== test ends
C
C  Test feasibility of reinitialization
C
      IF(ABS(ZQX).LE.0.1*ABS(ZQZ))  THEN
         WRITE(6,9500) ZQX,ZQZ
 9500 FORMAT(10X,'Ray almost vertically incident, DX =',1P,E13.4,/,
     +    10X,'DZ =',E13.4,0P,' reinitialization might fail')
      END IF
C
      IF(ABS(ZJX).GT.0.30*URPLAS)  THEN
         ZMX = ZJX/UKZERO
         ZKMX = PKP*ABS(ZJX)
         WRITE(6,9510) ZMX,ZKMX
 9510 FORMAT(10X,'Warning: cyclotron layer too large, DX =',1P,E13.4,/,
     +    10X,'estimated optical thickness =',E13.4,0P)
      END IF
C
         ITRY = 0
C
C  Step to advance the ray
C
         NEXT = 50
C
   30    ZDX = ZJX/FLOAT(NEXT)
         ZDZ = ZJZ/FLOAT(NEXT)
         ZDS = ZDX/ZQX
C
C  Next point along the ray
C
         IGUESS = 1
         IEXT = 0
C
   40    IEXT = IEXT+1
         UX = UX + ZDX
         UZ = UZ + ZDZ
C
C  Values od PSI & THETA at the new point
C
            CALL INVERT(IGUESS)
C ====== test begs
C        XU = UX/UKZERO
C        ZU = UZ/UKZERO
C     WRITE(6,7002)  IEXT,XU,ZU,PSI,THETA
C7002 FORMAT(' Step n.',i3,' X, z =',1P,2E13.4,' Psi,th =',2E13.4,0P)
C ====== test ends
C
      IF(IGUESS.EQ.0)  THEN
         WRITE(6,9520) XU,ZU,IEXT
 9520 FORMAT(10X,'Search for Psi not convergent at X =',1P,E13.4,/,
     +   10X,'Z =',E13.4,0P,' after',I3,' steps')
         PRMT(5) = 1.
         RETURN
      END IF
C
      IF(PSI.GE.1.)  THEN
         WRITE(6,9530)
 9530 FORMAT(10X,'The ray is leaving the plasma')
         ZPARM5 = 1.
         GO TO 200
      END IF
C
         ZPSI(IEXT) = PSI
         ZTHETA(IEXT) = THETA
C
      IF(IEXT.LT.NEXT) GOTO 40
C
C  Reevaluation of the wavevector at the new point
C
            CALL PROFIX
C
         ZCTAU = DZTH/ANTAU
         ZSTAU = -DXTH/ANTAU
         PKPHI = ANPHI/URHS
         ZPKZ2 = PKZ*PKZ
C
      DO 130  KTER=1,20
C
C  New K PERP with the old estimate of KPAR
C
         ZKPAR = PKPAR
C
            CALL DISPIC(1)
C
         ZDISCR = HB*HB - 4.*HA*HC
C
      IF(ZDISCR.LE.0.)  THEN
         WRITE(6,9540)
 9540 FORMAT(10X,'New N PERP is complex')
         ZPARM5 = 2.
         GO TO 200
      END IF
C
         ZQF2 = -ANR*ANL/ANS
         ZKP1 = 0.5*(-HB + SQRT(ZDISCR))/HA
         ZKP2 = 0.5*(-HB - SQRT(ZDISCR))/HA
         PKP2 = ZKP1
      IF(ABS(ZKP2-ZQF2).LT.ABS(ZKP1-ZQF2))  PKP2 = ZKP2
C
C  New KX (keeping KZ constant)
C
         ZA = 1. - ZSTAU*ZSTAU*ASINTQ*ASINTQ
         ZK = PKZ*ZCTAU*ASINTQ + PKPHI*ACOSTQ
         ZB = ZK*ZSTAU*ASINTQ
         ZC = -PKP2 + ZPKZ2 + PKPHI*PKPHI - ZK*ZK
         ZDISCR = ZB*ZB - ZA*ZC
C
      IF(ZDISCR.LE.0.)  THEN
         WRITE(6,9550) PKP2
 9550 FORMAT(10X,'New point too close or beyond a cutoff, PKP2 =',
     +         1P,E13.4,0P)
         ZPARM5 = 2.
         GO TO 200
      END IF
C
         PKX = (-ZB + ZADQX*SQRT(ZDISCR))/ZA
C
         PKP = SQRT(PKP2)
      DO 110  ITER=1,5
         PKPSI = PKX*ZCTAU + PKZ*ZSTAU
         PKTAU = -PKX*ZSTAU + PKZ*ZCTAU
         PKHETA = PKTAU*ACOSTQ - PKPHI*ASINTQ
         ZFUN = PKPSI*PKPSI + PKHETA*PKHETA - PKP2
         ZDIF = 2.*(PKPSI*ZCTAU - PKHETA*ZSTAU*ACOSTQ)
         ZDPKX = ZFUN/ZDIF
      IF(ABS(ZDPKX).LE.ACCUR*(ABS(PKX)+1.))  GO TO 120
  110    PKX = PKX - ZDPKX
C
C  Procedure not convergent after 5 iterations
C
         WRITE(6,9560)
 9560 FORMAT(10X,'Search for NX not convergent')
         ZPARM5 = 2.
         GO TO 200
C
C  New K PARALLEL estimate
C
  120    PKPAR = PKTAU*ASINTQ + PKPHI*ACOSTQ
         PKPAR2 = PKPAR*PKPAR
      IF(ABS(PKPAR-ZKPAR).LE.ACCUR*(ABS(PKPAR)+1.))  GO TO 140
  130 CONTINUE
C
C  Procedure not convergent after 20 iterations
C
         WRITE(6,9570)
 9570 FORMAT(10X,'Search for new k paralell not convergent')
         ZPARM5 = 1.
         GO TO 200
C
C  Check that the restarting point is 'asymptotic'
C
  140 IF(PKPAR*ZKPAR1.LE.0.)  THEN
         WRITE(6,9580) ZKPAR1,PKPAR
 9580 FORMAT(10X,'N parallel changes sign within the IC layer',/,
     +    10X,'old =',F10.3,', new =',F10.3)
         ZPARM5 = 1.
         GO TO 200
      END IF
C
      IF(PKP2.LT.PKPMIN*PKPMIN)  THEN
         WRITE(6,9550) PKP2
         ZPARM5 = 2.
         GO TO 200
      END IF
C
         ZXIC = ABS(XI(JJRES,KRES))/XIJUMP
C
C  New point still too close to resonance (too large change in k//)
C
      IF(ZXIC.LT.0.8 .OR. ZXIC.GT.1.25)  THEN
         IF(ITRY.EQ.2)  THEN
            WRITE(6,9590) ZKPAR1,PKPAR,XI(JJRES,KRES)
 9590 FORMAT(10X,'Warning: N par jump across the layer is too large',/,
     +    10X,'old =',F10.3,', new =',F10.3,/,
     +    10X,'xn at the exit point =',F8.3)
            IF(ZXIC.GT.2.0.OR.ZXIC.LT.0.5)  ZPARM5 = 2.
            GO TO 200
         END IF
C
C  Try again with a different end point
C
         ITRY = ITRY + 1
         UX = ZUX
         UZ = ZUZ
         ZJX = 2.*ZJX/(1. + ZXIC)
         ZJZ = 2.*ZJZ/(1. + ZXIC)
         GO TO 30
      END IF
C
C  Reinitialize the dependent variables
C
  200    AKPSI = (AJ*PKPSI + AG*PKTAU)/ANTAU
         AKTH = ANTAU*PKTAU
         Y(1) = AKPSI
         Y(2) = PSI
         Y(3) = AKTH
         Y(4) = THETA
C
      IF(ZPARM5.NE.0)  GOTO 300
C
C  Test the direction of the group velocity at the new point
C
            CALL DISPIC(2)
C
      IF(POYNTN*POYNT1.LE.0.)  THEN
         WRITE(6,9600)
 9600 FORMAT(10X,'Group velocity not OK at the exit point')
         ZPARM5 = 2.
         GO TO 300
      END IF
C
C  Reinitialization successful
C
         XU = UX/UKZERO
         ZU = UZ/UKZERO
         ZEIKON = (X + ZJX/ZQX)/(2.*PI)
      WRITE(6,1020)  XU,ZU,ZEIKON,PKPAR,PKP,XI(JJRES,KRES)
 1020 FORMAT(5X,'Exit point: X =',F8.3,' Z =',F8.3,
     +    ' at phase',1P,E13.4,/,
     +  5X,'N par =',E13.4,' N perp =',E13.4,' x(i) =',E13.4,0P)
C
         ZKPAR2 = PKPAR
         ZPKP2 = PKP
         ZPOYNT2 = POYNTN
         ZDIVY2 = ZDIVY
         ZNUMY2 = ZNUMY
         ZSIGN = SIGN(1.,EPOL(2)*ZPOLY1(2))
      DO 210  I=1,3
  210    ZPOLY2(I) = ZSIGN*EPOL(I)
C
C  Absorption in the resonance layer
C
  300    NEXT = IEXT
C
      DO 340  IEXT=1,NEXT
C
         PSI = ZPSI(IEXT)
         THETA = ZTHETA(IEXT)
C
         ZH2 = FLOAT(IEXT)/FLOAT(NEXT)
         ZH1 = 1.-ZH2
         PKPAR = ZH1*ZKPAR1 + ZH2*ZKPAR2
         PKP = ZH1*ZPKP1 + ZH2*ZPKP2
         PKP2 = PKP*PKP
C
      DO 310  I=1,3
  310    EPOL(I) = ZH1*ZPOLY1(I) + ZH2*ZPOLY2(I)
         ZDIVY  = ZH1*ZDIVY1 + ZH2*ZDIVY2
         ZNUMY  = ZH1*ZNUMY1 + ZH2*ZNUMY2
         POYNTN = ZH1*POYNT1 + ZH2*POYNT2
C
         ZENORM = 0.
       DO 320  I=1,3
  320    ZENORM = ZENORM + EPOL(I)*EPOL(I)
         ZENORM = SQRT(ZENORM)
      DO 330  I=1,3
  330    EPOL(I) = EPOL(I)/ZENORM
         EPOLX = ZSQRTH*(EPOL(1) + EPOL(2))
         EPOLY = ZSQRTH*(EPOL(1) - EPOL(2))
C
            CALL COORDS
            CALL PROFIX
            CALL DISPIC(0)
            CALL SGAMMA
C
      IF(IPWDP .NE. 0)  THEN
C
         ZDPOLD = ZDPNEW
         ZDPNEW = -PKP*ZNUMY*GAMMA/POYNTN
C ====== test begs
         ZPWX = PWX
         PWX = EXP(0.5*ZDS*(ZDPNEW+ZDPOLD))*PWX
      IF(PWX .GT. ZPWX)  THEN
         write(6,7003)  IEXT,PKP,ZNUMY,GAMMA,POYNTN
 7003 format(' Power is increasing at point',I3,/,
     + ' kperp =',F10.4,' numy =',1P,E13.4,' gamma =',E13.4,
     + ' Poyntn =',E13.4,0P)
      END IF
C ====== test ends
C         PWX = EXP(0.5*ZDS*(ZDPNEW+ZDPOLD))*PWX
C
            CALL ABSORB(0)
C
      END IF
C
C  New value of the phase - graphical output
C
         X = X + ZDS
         EIKON = X/(2.*PI)
C
      IF(ZPARM5.NE.2 .AND. X.GT.EIKREF)  THEN
C
         EIKREF = EIKREF + ZDEYKP
         EIKON = X/(2.*PI)
         XU = UX/UKZERO
         ZU = UZ/UKZERO
C
         IF(IGRAPH.GT.0)  THEN
            ISXZ = ISXZ+1
               CALL OUTGRA(2)
         END IF
C
      END IF
C
  340 CONTINUE
C
      IF(ZPARM5.GT.1.)  THEN
C
         NTYPE = 0
      IF(JJRES.EQ.1 .AND. KRES.NE.MAINSP)   NTYPE=-1
C
C  Modified july 16, 1993 to allow second harm. heating of minority
C
C *old* IF(JJRES.EQ.2 .AND. KRES.EQ.MAINSP)   NTYPE=-2
      IF(JJRES.EQ.2)   NTYPE=-2
            CALL TIHRES(NTYPE)
         ZPARM5 = 1.
      END IF
C
         PRMT(5) = ZPARM5
         EIKON = X/(2.*PI)
         Y(5) = ALOG(PWX/PWINIT)
C
      IF(PRMT(5).NE.0.)
     >      CALL OUTPUT(X)
C
      RETURN
      END

      SUBROUTINE TIHRES(NTYPE)
C
C ************ ************ ************ ************ ************
C
C  THIS MODULE EVALUATES THE WAVE BEHAVIOUR NEAR RESONANCES
C       (TWO-ION HYBRID, CYCLOTRON HARMONIC)
C
C ************ ************ ************ ************ ************
C
C  CAUSE OF CALL TO THE TIHRES SUBROUTINE:
C         NTYPE = -2 - 1ST HARMONIC OF MAIN SPECIES
C         NTYPE = -1 - FAILURE NEAR AN IC RESONANCE
C         NTYPE =  0 - OTHER FAILURE
C         NTYPE =  1 - HORIZONTAL REFLECTION
C         NTYPE =  2 - CONFLUENCE WITH THE ACOUSTIC WAVE
C
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
      DIMENSION    ZABCYC(8)
C ************ ************ ************ ************ ************
C
C   Decide how to proceed
C
      IF(NSPEC.EQ.1)  THEN
C
C   Only one species
C
         J = IRES(1)
C
         IF(J.LE.0)  THEN
C
C           No resonance, impossible with confluence
C                         STOP in all other cases.
C
            IF(NTYPE.EQ.2)  GO TO 9000
            RETURN
C
         END IF
C
         IF(J.EQ.1)  THEN
C
C           Fundamental, impossible with confluence or horiz. reflection
C
            IF(NTYPE.NE.0)  GO TO 9000
            RETURN
C
         END IF
C
C        2d harmonic
C
         IF(J.EQ.2)  THEN
C
            IF(AKXINI.EQ.0.)  GOTO 9000
C
C           2d harmonic from the low field side
C                    impossible with a confluence
C
            IF(AKXINI.LT.0. .AND. NTYPE.EQ.2)  GO TO 9000
C
C           2d harmonic from the high field side,
C                    impossible with horiz. reflection 
C
            IF(AKXINI.GT.0. .AND. NTYPE.EQ.1)  GO TO 9000
C
C           2d harmonic, cases which can be treated 
C
            MHARM = 1
            IHARM(1) = 1
            ZBHARM = BETAI(1)
            ZSCREN = UX - UXCYCL(2,1)
            IMIN = 0
            ZCNMIN = 0.
            PKZ2 = PKZ*PKZ
            ZQF2 = -ANR*ANL/ANS
            ZQFX2 = ZQF2 - PKZ2
            ZNS = ANS
C
            GO TO 200
C
        END IF
C
      END IF
C
C ********* ********* ********* ********* *********
C  Several species, test for ion-ion resonance
C
      IF(AKXINI.EQ.0.)  GOTO 9030
C
C     Ray from the LFS - incompatible with confluence
C
      IF(AKXINI.LT.0. .AND. 
     >      (NTYPE.EQ.2.OR.NCOF.EQ.0))  GO TO 9000
C
C     Ray from the HFS - Confluence only near an ion-ion resonance
C
      IF(AKXINI.GT.0. .AND. NRES.EQ.0)  GO TO 9000
C
C  Identification of the minority species
C
         ZPLMIN = 0.05
      IF(NTYPE.EQ.-1)  ZPLMIN = 0.1
         ZCNTOT = 0.
         ZCNMIN = 0.
         IMIN = 0
      DO 120  I=1,NSPEC
C
         ZCNTOT = ZCNTOT + OPCYI2(I)
         ZNU = 0.
         DO 110  J=1,NSPEC
            IF(ATZI(J).EQ.ATZI(I))  
     >            ZNU = ZNU + DENSI(J)/DENS
  110    CONTINUE
C
         IF(AKXINI.GT.0.)  THEN
            ZXTEOR = 0.5*AZI(I)*ZNU*(ATZI(I)/ATZI(MAINSP) -
     +        ATZI(MAINSP)/ATZI(I))
         ELSE
            ZXTEOR = AZI(I)*ZNU*(ATZI(I)/ATZI(MAINSP)-1.)
         END IF
C
         ZCONF = OHI(I)/(ZXTEOR+1.)
         IF(ABS(ZCONF-1.).LE.ZPLMIN)  THEN
            IMIN = IMIN+1
            MINOR(IMIN) = I
            ZCNMIN = ZCNMIN + OPCYI2(I)
         END IF
C
  120 CONTINUE
C
C   Check of consistency
C
      IF(IMIN.EQ.0)  GO TO 9030
      IF(ZCNMIN/ZCNTOT.GT.0.333)  GO TO 9010
C
C  Check that Doppler broadening does not wash out the siongularity
C
         IDOPPL = 0
         ZDL = 0.
      DO 130  J=1,IMIN
         I = MINOR(J)
         ZSCREN = UX - UXCYCL(1,I)
         ZDL = ZDL + PZL(I)
         ZUXCYC = UXCYCL(1,I)
         ZXITIH = OXI(I)*ZXTEOR
      IF(ABS(ZXITIH).LT.XIJUMP)  IDOPPL=1
  130 CONTINUE
C
C  Coincidences with second harmonics
C
         MHARM = 0
         ZBHARM = 0.
      DO 140  I=1,NSPEC
      IF(ABS(XI(2,I)).GT.5. .OR.
     >      ABS(ZUXCYC-UXCYCL(2,I)).LE.0.01*URHS) THEN
         MHARM = MHARM+1
         IHARM(MHARM) = I
         ZBHARM = ZBHARM + BETAI(I)
      END IF
  140 CONTINUE
C
C   The ray cannot be continued after a ion-ion resonance
C       or a first harmonic resonance with evanescence
C
      IF(IMIN+MHARM.EQ.0)  RETURN
         PRMT(5) = 1.
C
C   Approximate dispersion relation
C
         ZL = AL-ZDL
         ZS = 0.5*(AR+ZL)
         ZNL = PKPAR2 - ZL
         ZNS = PKPAR2 - ZS
         ZQF2 = -ANR*ZNL/ZNS
C
         ZPKZ2 = PKZ*PKZ
         ZQFX2 = ZQF2 - ZPKZ2
         ZQFXC2 = -2.*ANR-ZPKZ2
C
      IF(ZQFXC2.LE.0.)  GO TO 9020
         ZNXCYC = SQRT(ZQFXC2)
C
C  Investigation of the ion-ion layer
C
  200 IF(ZQFX2.LE.0.) GO TO 9020
C
         ZQFX = SQRT(ZQFX2)
C
         ZFAPOL = ANR/ZNS
         ZFAOPT = PI*ZFAPOL*ZFAPOL/2.
C
C  Optical thickness
C
         ZDELTA = 0.5*ZCNMIN*URHS
         ZSIGMA = 0.25*ZBHARM*URHS
C
         HTOPT1 = ZFAOPT*ZDELTA/ZQFX
         HTOPT2 = ZFAOPT*ZSIGMA*ZQF2/ZQFX
         ZETOPT = EXP(-2.*(HTOPT2 + HTOPT1))
C
C  Provenience of the ray 
C
      IF(AKXINI.LT.0.)  THEN
C
C        Ray from the low field side
C
         REFTWH = (1. - ZETOPT)**2
         TRATWH = ZETOPT
C
      ELSE
C
C        Ray from the high field side
C
         REFTWH = 0.
         TRATWH = ZETOPT
C
      END IF
C
         ABSTWH = 1. - TRATWH - REFTWH
C
C  Wave subject to cyclotron damping
C
      IF(NTYPE.LT.0)  ZSCREN = -ZSCREN
C
         ZAMPLI = 0.
      IF(AKXINI*ZSCREN.LT.0.)  THEN
         ZAMPLI = TRATWH
      ELSE
         IF(IMIN+MHARM.NE.0.) ZAMPLI = REFTWH
      END IF
C
C  Power balance
C
         ZABTOT = 0.
      DO 210  I=1,NSPEC
         UPWABI(I) = 0.
  210    ZABCYC(I) = 0.
C
      IF(IMIN.NE.0) THEN
C
C   Absorption by the minority species
C
         ZOPXI = 0.
         DO 220  J=1,IMIN
            I = MINOR(J)
  220       ZOPXI = ZOPXI + OPI2(I)*ABS(OXI(I))
            ZCYFAC = PI*ANR*ANR*URHS/(ZNXCYC*ZOPXI*ZOPXI)
         DO 230  J=1,IMIN
            I = MINOR(J)
               CALL INTHRM(I,ZOPXI,ZASUM)
            ZABCYC(I) = ZCYFAC*OPI2(I)*ZASUM
  230       ZABTOT = ZABTOT + ZABCYC(I)
C
      END IF
C
      IF(MHARM.NE.0)  THEN
C
C    Absorption at the first harmonic
C
         DO 240  J=1,MHARM
            I = IHARM(J)
               CALL INTHRM(I,ZOPXI,ZASUM)
            ZABCYC(I) = ZCYFAC*ABS(ANR)*BETAI(I)*ZASUM
  240       ZABTOT = ZABTOT + ZABCYC(I)
C
      END IF
C
C  Mode converted power
C
         PWMODC = ABSTWH*PWX
C
C  Power to the electrons
C
         ZABSE = 1.
C
      IF(IPWDP.NE.0 .AND. MHARM.EQ.0)  
     >   ZABSE = GAME(4)/GAMMA
C
         UPWABE = ZABSE*PWMODC
         ZPWTIH = UPWABE
C
C  Power to the ions
C
      IF(ZABTOT.NE.0.)  THEN
         ZPWCYI = (1.-EXP(-ZABTOT))*ZAMPLI*PWX/ZABTOT
         DO 310  I=1,NSPEC
  310    UPWABI(I) = ZABCYC(I)*ZPWCYI
      END IF
C
      IF(IPWDP.NE.0 .AND. MHARM.EQ.0)  THEN
         DO 320  I=1,NSPEC
  320    UPWABI(I) = UPWABI(I) + GAMI(I)*PWMODC/GAMMA
      END IF
C
      DO 330  I=1,NSPEC
  330    ZPWTIH = ZPWTIH + UPWABI(I)
C
C   Redefining the GAMMAS
C
      IF(ZPWTIH.EQ.0.)  GO TO 9030
C
         GAME(4) = UPWABE/ZPWTIH
      DO 340  I=1,NSPEC
  340    GAMI(I) = UPWABI(I)/ZPWTIH
         GAMMA = 1.
C
         PWX = PWX - ZPWTIH
         Y(5) = ALOG(PWX/PWINIT)
C
         ITWIHR = 1
C
            CALL ABSORB(0)
C
         ITWIHR = 0
C
C   OUTPUT STATEMENTS
C
            CALL OUTRUN(3)
C
      RETURN
C
C   DIAGNOSTIC STATEMENTS
C
 9000 WRITE(6,9900)
 9900 FORMAT(10X,'Incompatible with further analysis')
         PRMT(5) = 1.
      RETURN
C
 9010 WRITE(6,9910)
 9910 FORMAT(10X,'Minority concentration too large, all power',
     +   ' reflected')
         PRMT(5) = 1.
      RETURN
C
 9020 WRITE (6,9920) ZQF2,ZQFX2,ZQFXC2
 9920 FORMAT(10X,'QF2 =',E13.4,' QFX2 =',E13.4,' QFX2(0) =',
     +     E13.4,' Absorption calculation unsuccessful')
         PRMT(5) = 1.
      RETURN
C
 9030 WRITE(6,9930)
 9930 FORMAT(10X,'Apparently unrelated to an ion-ion resonance')
         PRMT(5) = 1.
      RETURN
C
      END

      SUBROUTINE INTHRM(ISPEC,ZOPXI,ZASUM)
C
C ************ ************ ************ ************ ************
C  THIS MODULE EVALUATES THE INTEGRAL FOR ABSORPTION
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
C ************ ************ ************ ************ ************
C
      DIMENSION ZS(502)
C
      DIMENSION    ZNU(NSPION),   ZAKI(NSPION)
C
      DATA    ZSQPI /1.7724538509/,   ZDXI /0.01/
C
C ************ ************ ************ ************ ************
C
         ZASUM = 1.
C
      IF(IMIN.EQ.0)  RETURN
C
      DO 10  J=1,IMIN
         I = MINOR(J)
         ZNU(I) = OPI2(I)*ABS(OXI(I))/ZOPXI
   10    ZAKI(I) = SQRT(ATM(I)*TEMPIZ(ISPEC)/(ATM(ISPEC)*TEMPIZ(I)))
C
         ZS(1) = 1./(ZSQPI*ZSQPI)
C
      DO 30  K=2,502
         ZXI = FLOAT(K-1)*ZDXI
         ZDIVR = 0.
         ZDIVI = 0.
      DO 20  J=1,IMIN
         I = MINOR(J)
         ZX = ZAKI(I)*ZXI
            CALL ZETA(ZX,ZXR,DZR,DUMMY)
         ZR = -ZXR/ZX
         ZI = ZSQPI*EXP(-ZX*ZX)
         ZDIVR = ZDIVR + ZNU(I)*ZR
   20    ZDIVI = ZDIVI + ZNU(I)*ZI
   30    ZS(K) = EXP(-ZXI*ZXI)/(ZDIVR*ZDIVR+ZDIVI*ZDIVI)
C
         ZSUM = 0.
      DO 40 I=1,500,2
   40    ZSUM = ZSUM + ZS(I) + 4.*ZS(I+1) + ZS(I+2)
C
         ZASUM = (4.*ZDXI/(3.*ZSQPI))*ZSUM
C
      RETURN
      END
