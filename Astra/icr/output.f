      SUBROUTINE OUTPUT(X)
C
C ************ ************ ************ ************ ************
C
C  Output subroutine for the Ray Tracing code RAYIC
C  Tests for failures, absorption, reflection, etc.
C  Prints the current results of the integration
C  Stores the results for the graphical output.
C
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
C  Phase and position in meters
C
         EIKON = X/(2.*PI)
         XU = UX/UKZERO
         ZU = UZ/UKZERO
         PWX = EXP(Y(5))*PWINIT
C
      IF(IPRINT.EQ.-11)  GO TO 300
      IF(PRMT(5).NE.0.)  GO TO 200
C
C  Accuracy check
C
         ZCONF = AMAX1(ABS(HC),ABS(HB*PKP2))
         ZERR = ABS(DISPH/ZCONF)
      IF(HLAST.LT.5.E-3)  NHLAST = NHLAST+1
      IF(NHLAST.GT.150)  JHLF = 11
C
      IF(JHLF.GT.10 .OR. ZERR.GT.0.1)  THEN
         PRMT(5) = 1.
         IF(JHLF.GT.11)  THEN
            IF(IPRINT.NE.-1)   WRITE(6,9999) JHLF
 9999 FORMAT(/,' Integration failed, JHLF =',I3)
            GO TO 100
         END IF
         IF(IPRINT.NE.-1)  WRITE(6,2010)  EIKON,XU,ZU,PWX
 2010 FORMAT(/,' Eikonal approximation fails at phase',F8.3,/,
     +   5X,'X =',F8.3,' Z =',F8.3,' power =',1P,E12.2,0P)
            CALL TIHRES(0)
         GO TO 100
      END IF
C
      IF(X.GT.(PRMT(2)+0.001*HLAST)) RETURN
C
C  Check the sign of the phase velocity (now superfluous)
C
      IF(IPWDP.NE.0 .AND. DY(5).GT.0.)  THEN
         PRMT(5) = 1.
         IF(IPRINT.NE.-1)  WRITE(6,2020)  EIKON,XU,ZU,PWX
 2020 FORMAT(/,' Group velocity changes sign at phase',F8.3,
     +   ' X =',F8.3,' Z =',F8.3,' power =',1P,E12.3,0P)
            CALL TIHRES(0)
         GO TO 100
      END IF
C
C  Ray leaving the plasma
C
      IF(Y(2).GT.PSIMAX .AND. DY(2).GT.0.)  THEN
         PRMT(5) = -2.
         NUREFL = NUREFL+1
         IF(NUREFL.GE.MAXREF)  PRMT(5) = 1.
         IF(IPRINT.NE.-1)  WRITE(6,2030) EIKON,XU,ZU,PWX
 2030 FORMAT(/,' The ray is leaving the plasma, phase',F8.3,
     +   ' X =',F8.3,' Z =',F8.3,' power =',1P,E12.3,0P)
         GO TO 100
      END IF
C
C  Radial reflection
C
      IF(DY(2)/DKDRIN.LE.0.)  THEN
         DKDRIN = -DKDRIN
C
C  Reflection from a magnetic surface
C
         IF(PKPSI/AKRINI.LE.0.) THEN
            AKRINI = -AKRINI
            GO TO 100
         END IF
C
C  Confluence of two waves
C
         IF(IPRINT.NE.-1) WRITE(6,2040) EIKON,XU,ZU,PWX
 2040 FORMAT(/,' A confluence is detected at phase',F8.3,
     +   ' X =',F8.3,' Z =',F8.3,' power =',1P,E12.3,0P)
            CALL TIHRES(2)
         GO TO 100
      END IF
C
C  Horizontal reflection
C
         ZDX = DXPSI*DY(2) + DXTH*DY(4)
      IF(ZDX/AKXINI.LE.0.)  THEN
         NUREFL = NUREFL+1
         IF(NUREFL.GE.MAXREF)  PRMT(5) = 1.
         IF(IPRINT.NE.-1)  WRITE(6,2050) EIKON,XU,ZU,PWX
 2050 FORMAT(/,' Reflection from a vertical plane at phase',
     +   F8.3,/,5X,'X =',F8.3,' Z =',F8.3,' power =',1P,E12.3,0P)
            CALL TIHRES(1)
         IF(PRMT(5).NE.1.)  AKXINI = -AKXINI
         GO TO 100
      END IF
C
C  Cut off approached
C
      IF(PKP.LT.PKPMIN) THEN
         IF(IPRINT.NE.-1)  WRITE(6,2060) EIKON,XU,ZU,PWX
 2060 FORMAT(/,' The ray approaches a cut-off at phase',F8.3,/,
     +   5X,'X =',F8.3,' Z =',F8.3,' power =',1P,E12.3,0P)
            CALL TIHRES(1)
         IF(PRMT(5).NE.1.)  THEN
            AKXINI = -AKXINI
            PRMT(5) = -3.
            NUREFL = NUREFL+1
            IF(NUREFL.GE.MAXREF)  PRMT(5) = 1.
         END IF
      END IF
C
  100 CONTINUE
C
      IF(IPWDP.NE.0) THEN
C
C  Power absorbed since the last step 
C
            CALL ABSORB(1)
C
C  Is the power completely absorbed?
C
         IF(PWX.LE.PWCONF)  THEN
            PRMT(5) = 1.
            IF(IPRINT.EQ.-1)  WRITE(6,2070)  EIKON,XU,ZU,PWX
 2070 FORMAT(/,' Power completely absorbed at phase',F8.3,
     +   ' X =',F8.3,' Z =',F8.3,' POWER =',1P,E12.3,0P)
         END IF
      END IF
C
C  Partial printout
C
      IF(IPRINT.LE.2 .AND. PRMT(5).EQ.0. .AND. X.NE.0. .AND.
     >   ABS(X-PRMT(2)).GT.0.001)  RETURN
C
  200 CONTINUE
C
      IF(IPRINT.NE.0)
     >       CALL OUTRUN(4)
C
C  Storage for the fraphical output
C
      IF(PRMT(5).NE.0.) RETURN
C
  300 CONTINUE
C
      IF(IGRAPH.LE.0)  RETURN
         ISXZ = ISXZ+1
C
C            CALL OUTGRA(2)
C
      RETURN
      END

      SUBROUTINE OUTGRA(IENTRY)
C
C ************ ************ ************ ************ ************
C INTERFACE FOR THE OUTPUT OF GRAPHS
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
C ************ ************ ************ ************ ************
C
      DIMENSION STABX(NSEQMX,2*NTHRAY), STABY(NSEQMX,7*NTHRAY)
      EQUIVALENCE (STABX(1,1),SEIKON(1,1)),(STABY(1,1),SKP(1,1))
C
C ************ ************ ************ ************ ************
C
      GO TO (1000,2000,3000), IENTRY
C
C ********** ********** ********** ********** ********** **********
C
 1000 CONTINUE
C
C  Initialization of the plot storage
C
         NPS    = NSEQMX
         NGR    = NTHIN
C
      DO 10  N=1,NGR
         ISIMPL(N) = 0
      DO 10 K=1,77
   10    STABY(N,K) = 0.
C
         KGRAPH = IGRAPH
C
      RETURN
C
C ************ ************ ************ ************ ************
C
 2000 CONTINUE
C
C   Storage for the graphical output
C
         ISIMPL(INDGR) = ISXZ
C
C  Plot tables full?
C
      IF(ISXZ.GE.NPS)  PRMT(5) = 1.
C
         SEIKON(ISXZ,INDGR) = EIKON
         SPSI(ISXZ,INDGR) = PSI
         SXU(ISXZ,INDGR) = XU
         SZU(ISXZ,INDGR) = ZU
         SKP(ISXZ,INDGR) = PKP
         SKPAR(ISXZ,INDGR) = PKPAR
C
         SKR(ISXZ,INDGR) = PKPSI
         SKTH(ISXZ,INDGR) = Y(3)
         SKX(ISXZ,INDGR) = PKX
         SKZ(ISXZ,INDGR) = PKZ
         SPOWER(ISXZ,INDGR) = PWX
C
      RETURN
C
C ************ ************ ************ ************ ************
C
 3000 CONTINUE
C
C ************ ************ ************ ************
C  OUTPUT OF GRAPHS
C ************ ************ ************ ************
C
C      IF(IPLOPR.EQ.1)  CALL GRPROF
C
C ************ ************ ************ ************
C  DRAWING THE TOKAMAK CROSSECTION
C ************ ************ ************ ************
C
C      IF(IGRAPH.EQ.1)  THEN
C
C         NCOPSI = 5
C         NPTHET = 201
C            CALL PLOPSI(NCOPSI,NPTHET)
C
C            CALL GCROSS
C
C            CALL GTITLE
C
C         IF(KGRAPH.EQ.-999)  THEN
C            WRITE(6,9990)
C 9990 FORMAT(1H0,20X,'Further plots aborted for lack of points')
C            RETURN
C         END IF
C
C         DO 70  NXVARB=1,2
C            CALL GSCALE
C         DO 70  NYVARB=1,7
C            IF(IPLOXY(NXVARB,NYVARB).NE.0)  CALL GKAPPA
C   70    CONTINUE
C
C      END IF
C
C            CALL GRAPOW
C
      RETURN
      END
