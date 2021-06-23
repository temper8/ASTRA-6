      SUBROUTINE COORDS
C
C   ****************************************************************
C   *                                                              *
C   *      THIS MODULUS EVALUATES X,Z, THEIR DERIVATIVES           *
C   *           WITH RESPECT TO PSI AND THETA, AND THE             *
C   *           METRIC ELEMENTS FOR USE IN THE ION CYCLOTRON       *
C   *           RAY TRACING PROGRAM                                *
C   *                                                              *
C   *      Version with Pereverzev representation of the equilib.  *
C   *                                                              *
C   ****************************************************************
C
      include 'COMMON.F'
C
C   ****************************************************************
C
      DIMENSION ZBX0(3), ZBX2(3), ZBZ0(3)
C
C   ****************************************************************
C
      IF(PSI.LT.ACCUR)  PSI=ACCUR
C
         ACOSTH = COS(THETA)
         ASINTH = SIN(THETA)
C
C   ****************************************************************
C   EVALUATIONS OF RSHIFT, RHO, DTRIAN, ELONG AND THEIR DERIVATIVES
C   ****************************************************************
C
         LMHD = IFIX(PSI/DPSMHD) + 1
         ZH = PSI - PSIMHD(LMHD)
C
      DO 10  I=1,3
         ID = I-1
         ZBX0(I) = VALSPL(ZH,BX0,LMHD,ID,NMHD)
         ZBX2(I) = VALSPL(ZH,BX2,LMHD,ID,NMHD)
   10    ZBZ0(I) = VALSPL(ZH,BZ0,LMHD,ID,NMHD)
C
C   ****************************************************************
C   EVALUATION OF XU(PSI,THETA), ZU(PSI,THETA), AND DERIVATIVES
C   ****************************************************************
C
         ASINSQ = ASINTH*ASINTH
C
         UX = ZBX0(1) + URPLAS*PSI*(ACOSTH - ZBX2(1)*ASINSQ)
         UZ = ZBZ0(1)*ASINTH
C
         URHS = URTOR + UX
C
         DXPSI = ZBX0(2) + URPLAS*(ACOSTH -
     +              (ZBX2(1) + PSI*ZBX2(2))*ASINSQ)
         DXTH = -URPLAS*PSI*ASINTH*(1. + 2.*ZBX2(1)*ACOSTH)
C
         DZPSI = ZBZ0(2)*ASINTH
         DZTH = ZBZ0(1)*ACOSTH
C
         ZDDXPP = ZBX0(3) - URPLAS*(2.*ZBX2(2) + PSI*ZBX2(3))*ASINSQ
         ZDDXPT = -URPLAS*ASINTH*(1. +
     +               2.*(ZBX2(1) + PSI*ZBX2(2))*ACOSTH)
         ZDDXTT = -URPLAS*(PSI*ACOSTH + 
     +               2*ZBX2(1)*(1. - 2.*ASINSQ))
C
         ZDDZPP = ZBZ0(3)*ASINTH
         ZDDZPT = ZBZ0(2)*ACOSTH
         ZDDZTT = -UZ
C
C   ****************************************************************
C   EVALUATION OF THE JACOBIAN, G, NTAU, AND DERIVATIVES
C   ****************************************************************
C
         ANTAU2 = DXTH*DXTH+DZTH*DZTH
         ANTAU = SQRT(ANTAU2)
         AJ = DXPSI*DZTH - DZPSI*DXTH
         AG = DXPSI*DXTH + DZPSI*DZTH
C
C  LOG. DERIV. OF ANTAU & AJ, ORD. DERIV. OF AG
C
         DNTAUP = (DXTH*ZDDXPT+DZTH*ZDDZPT)/ANTAU2
         DNTAUT = (DXTH*ZDDXTT+DZTH*ZDDZTT)/ANTAU2
         DJP = (ZDDXPP*DZTH+ZDDZPT*DXPSI-ZDDXPT*DZPSI-ZDDZPP*DXTH)/
     +         AJ
         DJT = (ZDDXPT*DZTH+ZDDZTT*DXPSI-ZDDXTT*DZPSI-ZDDZPT*DXTH)/
     +         AJ
         DGP = ZDDXPP*DXTH+ZDDXPT*DXPSI+ZDDZPP*DZTH+ZDDZPT*DZPSI
         DGT = ZDDXPT*DXTH+ZDDXTT*DXPSI+ZDDZPT*DZTH+ZDDZTT*DZPSI
C
      IF(IPRINT.GE.4)  THEN
         WRITE(6,7010)  PSI,THETA,UX,UZ,URHS,ASINTH,ACOSTH
 7010 FORMAT(/,' Metric elements at psi =',F7.3,' theta =',F7.3,/,
     +   ' Ux     =',1P,E13.4,' Uz     =',E13.4,' URhs   =',E13.4,/,
     +   ' Asinth =',E13.4,' Acosth =',E13.4,0P)
         WRITE(6,7020)  DXPSI,DXTH,DZPSI,DZTH,ANTAU2,AJ,AG
 7020 FORMAT(' dXdPsi =',1P,E13.4,' dXdTh  =',E13.4,/,
     +   ' dZdPsi =',E13.4,' dZdTh  =',E13.4,/,
     +   ' Ntau2  =',E13.4,' J      =',E13.4,' G      =',E13.4,0P)
      END IF
C
      RETURN
      END

      SUBROUTINE WAVECT
C
C   ****************************************************************
C   *                                                              *
C   *      THIS MODULUS EVALUATES THE PHISICAL COMPONENTS          *
C   *           OF THE WAVEVECTOR IN THE PSI, THETA, PHI           *
C   *           FRAME AND IN THE LOCAL STIX FRAME.                 *
C   *                                                              *
C   ****************************************************************
C
      include 'COMMON.F'
C
C   ****************************************************************
C
         PKPSI = (ANTAU*AKPSI - AG*AKTH/ANTAU)/AJ
         PKTAU = AKTH/ANTAU
         PKPHI = ANPHI/URHS
C
C  STIX COMPONENTS
C
         PKHETA = PKTAU*ACOSTQ - PKPHI*ASINTQ
         PKPAR = PKTAU*ASINTQ + PKPHI*ACOSTQ
         PKPAR2 = PKPAR*PKPAR
         PKP2 = PKPSI*PKPSI + PKHETA*PKHETA
         PKP = SQRT(PKP2)
C
C  CARTESIAN COMPONENTS
C
         PKX = (AKPSI*DZTH - AKTH*DZPSI)/AJ
         PKZ = (-AKPSI*DXTH + AKTH*DXPSI)/AJ
C
C  DERIVATIVES
C
         DKPSIP = (DNTAUP-DJP)*PKPSI + (2.*AG*DNTAUP-DGP)*PKTAU/AJ
         DKPSIT = (DNTAUT-DJT)*PKPSI + (2.*AG*DNTAUT-DGT)*PKTAU/AJ
         DKHEP = -DNTAUP*ACOSTQ*PKTAU+PKPHI*DXPSI*ASINTQ/URHS
     +      -DQPSI*PKPAR
         DKHET = -DNTAUT*ACOSTQ*PKTAU+PKPHI*DXTH*ASINTQ/URHS-DQTH*PKPAR
         DKP2P = 2.*(PKPSI*DKPSIP+PKHETA*DKHEP)
         DKP2T = 2.*(PKPSI*DKPSIT+PKHETA*DKHET)
         DKPZP = -DNTAUP*ASINTQ*PKTAU - PKPHI*DXPSI*ACOSTQ/URHS +
     +      DQPSI*PKHETA
         DKPZT = -DNTAUT*ASINTQ*PKTAU - PKPHI*DXTH*ACOSTQ/URHS +
     +      DQTH*PKHETA
C
      IF(IPRINT.GE.4)  THEN
         WRITE(6,7010) AKPSI,AKTH,ANPHI
 7010 FORMAT(' ak_psi  =',1P,E13.4,' ak_th   =',E13.4,
     +       ' ak_phi  =',E13.4,0P)
         WRITE(6,7020) PKPSI,PKTAU,PKPHI,PKHETA,PKPAR,PKP
 7020 FORMAT(' k_psi  =',1P,E13.4,' k_tau  =',E13.4,
     +       ' k_phi  =',E13.4,/,
     +       ' k_eta  =',E13.4,' k_par  =',E13.4,' k_perp =',E13.4,0P)
      END IF
C
      RETURN
      END

      SUBROUTINE INVERT(IGUESS)
C
C ************ ************ ************ ************ ************
C
C  This modulus evaluates PSI amd THETA given X and Z
C  (inversion of SUBROUTINE COORDS)
C
C  INPUT:
C      IGUESS = 0 - No initial guess for PSI and THETA available
C      IGUESS = 1 - Initial guess for PSI and THETA available
C
C  OUTPUT
C      IGUESS = 0 - Search not convergent
C      IGUESS = 1 - Search successful
C
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
C ************ ************ ************ ************ ************
C
C  Avoid the origin
C
      IF((ABS(UX)+ABS(UZ)).GT.1.E-12)  GOTO 10
         PSI = 1.E-12
         THETA = ATAN2(UZ,UX)
      RETURN
C
C  Store the goal values of X, Z
C
   10    ZXU = UX
         ZZU = UZ
C
         ITRY = 0
      IF(IGUESS.EQ.1)  GOTO 100
C
C  Initial guess (when not already available)
C
   20    ITRY = 1
         THETA = ATAN2(ZZU,ZXU)
         ZPX = (URPLAS - URSHFT)*COS(THETA)
         ZPZ = URPLAS*POLASP(1)*SIN(THETA)
         PSI = SQRT((ZXU*ZXU+ZZU*ZZU)/(ZPX*ZPX+ZPZ*ZPZ))
C
C  Iteration
C
  100 DO 110  ITER=1,20
C
      IF(PSI.GT.1.001.OR.PSI.LT.0.)  GO TO 120
C
            CALL COORDS
C
         ZDX = ZXU - UX
         ZDZ = ZZU - UZ
         ZDPSI = (DZTH*ZDX - DXTH*ZDZ)/AJ
         ZDTH = (-DZPSI*ZDX + DXPSI*ZDZ)/AJ
         PSI = PSI + ZDPSI
         THETA = THETA + ZDTH
C
      IF((ABS(ZDX)+ABS(ZDZ))/URPLAS.LE.ACCUR)  GOTO 130
C
  110 CONTINUE
C
C   Search not convergent
C
  120 IF(ITRY.EQ.0)  GO TO 20
         IGUESS = 0
      IF(IPRINT.LE.1)  RETURN
      WRITE(6,9000) ZDPSI,ZDTH,ITER,ZXU,ZZU,PSI,ZDX,ZDZ
 9000 FORMAT(' Search for PSI, THETA not convergent:',/,
     +   ' DPSI, ZDTH =',1P,2E13.4,0P,' after',I3,' iterations',/,
     +   ' UX,UZ,PSI,ZDX,ZDZ =',1P,5E13.4,0P)
C
      RETURN
C
  130    IGUESS = 1
C
      RETURN
      END
