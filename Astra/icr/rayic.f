      SUBROUTINE RAYIC
C
C ************ ************ ************ ************ ************
C  RAY-TRACING CODE RAYIC - Monitor subroutine
C  Last modified February 9, 1996
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
C ************ ************ ************ ************ ************
C   Plasma parameters and initialization
C
            CALL INITAL
C
      IF(ISTOP.NE.0)  RETURN
C
C ************ ************ ************ ************
C   Printout of the run parameters
C
            CALL OUTRUN(1)
C
      IF(IGRAPH.NE.0)
     >      CALL OUTGRA(1)
C
C ************ ************ ************ ************
C   Ray tracing subroutine
C
            CALL RAYPAT
C
C ************ ************ ************ ************
C   Power deposition profiles and plots
C
            CALL OUTRUN(2)
C
      IF(IGRAPH.NE.0)
     >      CALL OUTGRA(3)
C
C ************ ************ ************ ************
C
      RETURN
      END

      SUBROUTINE RAYPAT
C
C ************ ************ ************ ************ ************
C  RAY-TRACING CODE RAYIC - Ray Tracing implementation.
C  Last modified February 9, 1996
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
      DIMENSION   ZPWOLI(NSPION)
C
      EXTERNAL DERHIC,OUTPUT
C
C ************ ************ ************ ************ ************
C   Loop over the initial value of N-parallel
C ************ ************ ************ ************ ************
C
      DO 120  JNPHI=1,NUMPHI
C
      IF(IPRINT.NE.0)
     >   WRITE(6,1010) NPHI(JNPHI)
 1010 FORMAT(/,20X,'------- Nphi =',I3,' -------')
C
         INPHI  = JNPHI
         ANPHI  = FLOAT(NPHI(INPHI))
         INDGR  = INPHI
         ZPWOLE = PWTOTE
      DO 10  I=1,NSPEC
   10    ZPWOLI(I) = PWTOTI(I)
C
            CALL INIWFR
C
C ************ ************ ************ ************ ************
C   Loop over the initial poloidal angle 
C ************ ************ ************ ************ ************
C
      DO 100  JTH=1,NTHIN
C
      IF(IPRINT.NE.0)
     >   WRITE(6,1020)  THSTRT(JTH)
 1020 FORMAT(/,'================== ',
     +         'Ray from theta =',F7.1,' =====================')
C
         ITH = JTH
      IF(IGRAPH.EQ.1)  INDGR=ITH
         ISXZ = 0
C
C ************ ************ ************ ************ ************
C  Initialization for HPCSD (integrating subroutine)
C ************ ************ ************ ************ ************
C
         PWINIT = PWCPL(INPHI)*POYNT(INPHI,ITH)
C
      IF(PWINIT.LT.PWCONF)  THEN
         WRITE(6,9010) PWINIT
 9010 FORMAT(/,20X,'No power coupled to this ray, PWINIT =',E13.4)
         GOTO 100
      END IF
C
      IF(INZWFR(ITH) .EQ. 0)  THEN
         WRITE(6,9020)
 9020 FORMAT(/,20X,'Ray could not be initialized')
         GOTO 100
      END IF
C
         PSI    = ZPSINI(ITH)
         THETA  = ZTHINI(ITH)
         AKPSI  = ZKPSIN(ITH)
         AKTH   = ZKTHIN(ITH)
         PWX    = PWINIT
C
         Y(1)   = AKPSI
         Y(2)   = PSI
         Y(3)   = AKTH
         Y(4)   = THETA
         Y(5)   = 0.
C
         EIKON  = 0.
         X      = 0.
         EIKREF = ZDEYKP
         POWOLD = PWINIT
         PKPMIN = 4.
         HLAST  = 0.5*DSTEPH
         NHLAST = 0
         NUREFL = 0
         KEIK   = 1
C
   20       CALL DERHIC(X)
C
         DKDRIN = DY(2)/ABS(DY(2))
         AKRINI = Y(1)/ABS(Y(1))
         ZDXINI = DXPSI*DY(2) + DXTH*DY(4)
         AKXINI = ZDXINI/ABS(ZDXINI)
C
         NEIK    = KEIK
         PRMT(4) = ACCUR
         PRMT(5) = 0.
C
C ************ ************ ************ ************ ************
C  Loop over equiphase surfaces
C ************ ************ ************ ************ ************
C
      DO 50  JEIK=NEIK,1000
         KEIK    = JEIK
         PRMT(1) = X
         PRMT(2) = EIKREF
         PRMT(3) = 2.*HLAST
      DO 30  NZK=1,5
   30    DY(NZK) = PRECDY(NZK)
C
C ************ ************ ************ ************ ************
C  Call the integration subroutine
C ************ ************ ************ ************ ************
C
            CALL HPCSD(DERHIC,OUTPUT)
C
         X      = 2.*PI*EIKON
C
         IF(PRMT(5) .NE. 0.)  THEN
C
C  Integration failed
C
             IF(PRMT(5).EQ.1.)  GO TO 60
C
C  Cyclotron or ion-ion resonance to be crossed
C
            IF(PRMT(5).EQ.-1.)  THEN
C
               CALL CYCRES(X)
C
               IF(PRMT(5).NE.0.)  GO TO 60
               GO TO 50
C
            END IF
C
            IF(PRMT(5).GT.-2.)  GO TO 60
C
C  Integration can be resumed
C
            PKPSI  = -PKPSI
            Y(1)   = (AJ*PKPSI + AG*PKTAU)/ANTAU
            IF(PRMT(5).EQ.-2.)  PSIMAX = PSI
            IF(PRMT(5).EQ.-3.)  PKPMIN = PKP
C
            GO TO 20
C
         END IF
C
         EIKREF = EIKREF + ZDEYKP
C
   50    CONTINUE
C
C  Integration must be stopped
C
   60    PWTRA(INPHI) = PWTRA(INPHI) + PWX
C
  100 CONTINUE
C
         PWPARE(INPHI) = PWTOTE - ZPWOLE
      DO 110  I=1,NSPEC
  110    PWPARI(INPHI,I) = PWTOTI(I) - ZPWOLI(I)
C
  120 CONTINUE
C
C ************ ************ ************ ************ ************
C
      RETURN
      END

      SUBROUTINE DERHIC(X)
C
C ************ ************ ************ ************ ************
C  RAY TRACING CODE RAYIC
C  RHS of the ray tracing equations
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
C  Ray through the origin?
C
      IF(Y(2).LT.0.)  THEN
         Y(2) = ABS(Y(2))
         IF(DY(4).GE.0.)  THEN
            Y(4) = Y(4) + PI
         ELSE
            Y(4) = Y(4) - PI
         END IF
      END IF
C
C  Transfer to physical variables
C
   10    PSI    = Y(2)
         THETA  = Y(4)
         AKPSI  = Y(1)
         AKTH   = Y(3)
C
      IF(PSI.GT.1.)   PSI=1.
C
C  Metric elements
C
            CALL COORDS
C
C  Profile factors
C
            CALL PROFIX
C
C  Components of the wavevector
C
            CALL WAVECT
C
C  Dispersion relation and power transport
C
            CALL DISPIC(2)
C
C  Derivative of AKPSI
C
         DY(1) = -(DHNE*DERDEN + DHTE*DERTEM)
      DO 20  I=1,NSPEC
   20    DY(1) = DY(1) - (DHTX(I)*DERTIX(I) + DHTZ(I)*DERTIZ(I) + 
     +                    DHNI(I)*DERDNI(I))
         DY(1) = DY(1) - (DHOM*DOMPSI + DKP2P*DHKP2 + DKPZP*DHKZ)
C
C  Derivative of PSI
C
         DY(2) = 2.*PKPSI*DHKP2*ANTAU/AJ
C
C  Derivative of AKTH
C
         DY(3) = -(DHOM*DOMTH + DKP2T*DHKP2 + DKPZT*DHKZ)
C
C  Derivative of THETA
C
         DY(4) = (2.*(-AG*PKPSI/AJ + ACOSTQ*PKHETA)*DHKP2 +
     +         ASINTQ*DHKZ)/ANTAU
C
C  Phase as independen variable
C
         ZDIVY = Y(1)*DY(2) + Y(3)*DY(4)
      DO 30  I=1,4
   30    DY(I) = DY(I)/ZDIVY
C
C  Derivative of the power flux along the ray
C
      IF(IPWDP.EQ.0)  THEN
         DY(5) = 0.
         ZNUMY = 0.
      ELSE
         ANPSI2 = (AJ*AJ + AG*AG)/ANTAU2
         ZNUMY = ANPSI2*DY(2)*DY(2) + 2.*AG*DY(2)*DY(4) +
     +               ANTAU2*DY(4)*DY(4)
         DY(5) = -PKP*ZNUMY*GAMMA/POYNTN
      END IF
C
      RETURN
      END
