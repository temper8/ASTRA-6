      SUBROUTINE OUTRUN(IENTRY)
C
C ************ ************ ************ ************ ************
C
C  Ion Cyclotron Ray Tracing code RAYIC
C  Output of the run parameters
C
C ************ ************ ************ ************ ************
C
      include 'COMMON.F'
C
      include 'COMHPCSD.F'
C
      DIMENSION QPN(NSPION), QPWE(NPSIPT), QPWI(NPSIPT,NSPION),
     +          ZGAME(3), ZGAMI(NSPION)
C
C ************ ************ ************ ************ ************
C
      GO TO (1000,2000,3000,4000), IENTRY
C
C ************ ************ ************ ************ ************
C
 1000 CONTINUE
C
         ZTCURR = 1.E-6*TCURR
      WRITE(6,5010)  RTORUS,RAXMAG,RPLASM,XMIN,XMAX,BZERO,BAXIS,
     +   ZTCURR,PQR(NPROF),PQR(1)
 5010 FORMAT(/,10X,'Major radius =                    ',F9.3,' m',/,
     +       10X,'Radius of the magnetic axis =     ',F9.3,' m',/,
     +       10X,'Minor radius =                    ',F9.3,' m',/,
     +       10X,'Inner plasma edge at X =          ',F9.3,' m',/,
     +       10X,'Outer plasma edge at X =          ',F9.3,' m',/,
     +       10X,'Toroidal magnetic field, center   ',F9.2,' Tesla',/,
     +       10X,'Toroidal magnetic field, axis =   ',F9.2,' Tesla',/,
     +       10X,'Toroidal current =                ',F9.2,' MA',/,
     +       10X,'Safety factor at the edge =       ',F9.2,/,
     +       10X,'Safety factor on the axis =       ',F9.2)
C
      WRITE(6,5020)  RSHIFT,POLASP,DTRIAN
 5020 FORMAT(10X,'Shift of the magnetic axis =      ',F9.3,' m',/,
     +       10X,'Ellipticity at the axis =         ',F9.3,/,
     +       10X,'Ellipticity of the outer surface =',F9.3,/,
     +       10X,'Triangularity at outer surface =  ',F9.3)
C
      WRITE(6,5030)  DENC,TEMPEC
 5030 FORMAT(/,10X,'Central density =             ',1P,E13.3,' cm-3',/,
     +         10X,'Central electron temperature =',0P,F13.3,' keV')
      WRITE(6,5040) MAINSP
 5040 FORMAT(/,10X,'Majority ion species is species n.',I3)
      DO 110 J=1,NSPEC
         ZCONC = AZI(J)*DENIC(J)/DENC
  110 WRITE(6,5050) J,ATM(J),AZI(J),ZCONC,TEMPIC(J)
 5050 FORMAT(10X,'Sp. ',I2,' A =',F5.1,
     +   ' Z =',F5.1,' Z*n(j)/n(e) =',F7.3,' Ti =',F7.3,' keV')
C
C  FREQUENCY, WAVELENGTH
C
         ZFREMH = 1.D-6*FREQCY
         WLENGT = 2.*PI/UKZERO
      WRITE(6,5060)  ZFREMH,WLENGT
 5060 FORMAT(/,10X,'Applied frequency =           ',F13.3,' Mhz',/,
     +         10X,'Wavelength in vacuum =        ',F13.3,' m')
C
      IF(NUMPHI.EQ.1)  THEN
         WRITE(6,5070)  NPHI1,NTHIN
 5070 FORMAT(10X,'Toroidal wavenumber Nphi =             ',I4,/,
     +       10X,'Rays in the poloidal cross-section     ',I4)
      ELSE
         NPHI2 = NPHI1 + (NUMPHI-1)*JUMPHI
         WRITE(6,5080)  NPHI1,NPHI2,NTHIN
 5080 FORMAT(10X,'Toroidal wavenumber Nphi from     ',I4,' to',I4,/,
     +       10X,'Rays in the poloidal cross-section     ',I4)
      END IF
C
      WRITE(6,5090)
 5090 FORMAT(/,10X,'Singular layers:')
      DO 130  I=1,NSPEC
         ZUX1 = UXCYCL(1,I)/UKZERO
         ZUX2 = UXCYCL(2,I)/UKZERO
         RUX1 = RTORUS + RSHIFT + ZUX1
         RUX2 = RTORUS + RSHIFT + ZUX2
  130 WRITE(6,5100)  I,ZUX1,RUX1,ZUX2,RUX2
 5100 FORMAT(10X,'Sp.',I3,' Fundam. res. at X    =',F9.4,
     +             ' (R =',F9.4,')',/,
     +    16X,' First harm. res at X =',F9.4,' (R =',F9.4,')')
C
      IF(NCOF.GT.0)  THEN
         DO 140  I=1,NCOF
            ZUXCR = UXCOF(I)/UKZERO
            RUXCR = RTORUS + RSHIFT + ZUXCR
  140    WRITE(6,5110)  ZUXCR,RUXCR
 5110 FORMAT(10X,'Cutoff at X                 =',F9.4,
     +           ' (R =',F9.4,')')
      END IF
C
      IF(NRES.GT.0)  THEN
         DO 150  I=1,NRES
            ZUXCR = UXRES(I)/UKZERO
            RUXCR = RTORUS + RSHIFT + ZUXCR
  150    WRITE(6,5120)  ZUXCR,RUXCR
 5120 FORMAT(10X,'Hybrid resonance at X       =',F9.4,
     +           ' (R =',F9.4,')')
      END IF
C
C ************ ************ ************ ************ ************
C
      IF(IPRINT.GT.1)  THEN
C
         WRITE(6,5210)
 5210 FORMAT(/,10X,'COEFFICIENTS OF THE MHD CONFIGURATION:'/,
     +   6X,'Psi',12X,'Rmaj(PSI) (m)',8X,'triang(PSI)',
     +   8X,'Elong(PSI)')
         WRITE(6,5220)  (PSIMHD(N),RMAJ(N),TRIANG(N),
     +                                     ELONG(N),N=1,NMHD)
 5220 FORMAT(4E16.4)
C
         WRITE(6,5230)
 5230 FORMAT(/,10X,'DENSITY AND TEMPERATURE PROFILES',/)
         WRITE(6,5240) (PPSI(L),TQR(L,1),PQR(L),PNE(L),PTE(L),
     +       L=1,NPROF)
 5240 FORMAT(' Psi =',1P,E11.3,' Bpol =',E11.3,' q =',
     +   E11.3,' Ne =',E11.3,0P,' Te =',F7.3)
C
         JS = NSPEC
  160    WRITE(6,5250)
 5250 FORMAT(/)
         WRITE(6,5260) (PPSI(L),PNI(L,JS),PTI(L,JS),L=1,NPROF)
 5260 FORMAT(' PSI =',1P,E12.3,' Ni  =',E12.3,0P,' Ti  =',F7.3)
C
      END IF
C
C ************ ************ ************ ************ ************
C
      RETURN
C
C ************ ************ ************ ************ ************
C
 2000 CONTINUE
C
C ************ ************ ************ ************ ************
C
C   PRINT-OUT OF THE POWER DEPOSITION PROFILES AND POWER SPECTRA
C
C ************ ************ ************ ************ ************
C
         QPWE(1) = 0.
      DO 2010  I=1,NSPEC
 2010    QPWI(1,I) = 0.
C
      DO 2020  IP=2,NPROF
         QPWE(IP) = QPWE(IP-1) + PWE(IP)*VOLUME(IP)
      DO 2020 I=1,NSPEC
 2020    QPWI(IP,I) = QPWI(IP-1,I) + PWI(IP,I)*VOLUME(IP)
C
      IF(PWTOTE.NE.0.)  THEN
         QPNE = PWTOTE/QPWE(NPROF)
      ELSE
         QPNE = 0.
      END IF
      DO 2030  I=1,NSPION
      IF(PWTOTI(I).NE.0.)  THEN
         QPN(I) = PWTOTI(I)/QPWI(NPROF,I)
      ELSE
         QPN(I) = 0.
      END IF
 2030 CONTINUE
C
      DO 2040  IP=1,NPROF
         QPWE(IP) = QPNE*QPWE(IP)
      DO 2040  I=1,NSPEC
 2040    QPWI(IP,I) = QPN(I)*QPWI(IP,I)
C
      WRITE(6,5500)  (I,I=1,NSPEC)
 5500 FORMAT(/,20X,'POWER DEPOSITION PROFILES',/,
     + 5X,'PSI',6X,'ELECTRONS',4X,'IONS:',I3,7I10,/)
      DO 210  N=1,NPROF
  210 WRITE(6,5510) PPSI(N),PWE(N),(PWI(N,I),I=1,NSPEC)
 5510 FORMAT(2X,F8.4,1P,9E12.2)
C
      WRITE(6,5505)  (I,I=1,NSPEC)
 5505 FORMAT(/,20X,'Integrated profiles',/,
     + 5X,'PSI',6X,'ELECTRONS',4X,'IONS:',I3,7I10,/)
      DO 215  N=1,NPROF
  215 WRITE(6,5510) PPSI(N),QPWE(N),(QPWI(N,I),I=1,NSPEC)
C
C  RENORM. TO UNITY THE POSITIVE PART OF THE SPECTRUM ALONE
C      IF IASYM=1 (SYMMETRIC SPECTRUM)
C
      IF(IASYM.NE.0)  THEN    
         ZNORM = 1.+PWCPL(1)
         PWCPL(1) = 2.*PWCPL(1)
         PWTRA(1) = 2.*PWTRA(1)
         PWPARE(1) = 2.*PWPARE(1)
      DO 220  I=1,NSPEC
  220    PWPARI(1,I) = 2.*PWPARI(1,I)
      DO 230  N=1,NUMPHI
         PWCPL(N) = PWCPL(N)/ZNORM
         PWTRA(N) = PWTRA(N)/ZNORM
         PWPARE(N) = PWPARE(N)/ZNORM
      DO 230  I=1,NSPEC
  230    PWPARI(N,I) = PWPARI(N,I)/ZNORM
      END IF   
C  
      IF(NUMPHI.GT.1)  THEN     
         WRITE(6,5520)
 5520 FORMAT(/,10X,'SPECTRE OF THE NON-ABSORBED POWER')
         DO 240  N=1,NUMPHI
            ZPWABS = PWCPL(N) - PWTRA(N)
  240    WRITE(6,5530) NPHI(N),PWCPL(N),ZPWABS,PWTRA(N)
 5530 FORMAT(2X,'NPHI =',I3,2X,'POWER: RADIATED =',1P,E13.3,/,
     +    13X,'ABSORBED =',E13.3,2X,'NON ABSORBED =',E13.3,2X,'MW')
C
         WRITE(6,5540) (I,I=1,NSPEC)
 5540 FORMAT(/,10X,'ABSORPTIONS FOR EACH TOROIDAL MODE',/,
     +   13X,'ELECTRONS',3X,'IBW   ',3X,'IONS: ',I3,7I12)
         DO 250  N=1,NUMPHI
  250    WRITE(6,5550)  NPHI(N),PWPARE(N),PWPARB(N),
     +                 (PWPARI(N,I),I=1,3)
 5550 FORMAT(' NPHI =',I3,1P5E12.3,0P)
         IF(NSPEC.GT.3)  THEN
            WRITE(6,5555)  (PWPARI(N,I),I=4,NSPEC)
 5555 FORMAT(10X,1P,5E12.3,0P)
         END IF
      END IF
C
         ZPWNAB = POWER - PWEDGE - PWTOTE
      DO 260  I=1,NSPEC
  260    ZPWNAB = ZPWNAB - PWTOTI(I)
         ZPWABS = POWER - PWEDGE - ZPWNAB
C
      IF(ZPWABS.EQ.0.)  ZPWABS = POWER
C
         ZZGAME = 100.*PWTOTE/ZPWABS
         ZGAMBW = 100.*PWIBWS/ZPWABS
      DO 270  I=1,NSPEC
  270    ZGAMI(I) = 100.*PWTOTI(I)/ZPWABS
C
      WRITE(6,5560) POWER,ZPWABS,PWEDGE,ZPWNAB
 5560 FORMAT(/,25X,'TOTAL POWER BALANCE',/,
     +        10X,'Launched power             =',1P,E13.4,' MW',/,
     +        10X,'Absorbed by the plasma     =',E13.4,' MW',/,
     +        10X,'Absorbed near the antenna  =',E13.4,' MW',/,
     +        10X,'Non absorbed               =',E13.4,' MW')
      WRITE(6,5570) PWTOTE,ZZGAME,PWIBWS,ZGAMBW
 5570 FORMAT(/,10X,'Power to the electrons    =',1P,E13.4,' MW',
     +       0P,5X,'(',F10.2,' %)',/,
     +  10X,'Power to Bernst. waves    =',1P,E13.4,' MW',
     +       0P,5X,'(*',F9.2,' %)')
      WRITE(6,5580) (I,PWTOTI(I),ZGAMI(I),I=1,NSPEC)
 5580 FORMAT(  10X,'        to species ',I6,' =',1P,E13.3,' MW',
     +       0P,5X,'(',F10.2,' %)')
      WRITE(6,5590)
 5590 FORMAT(/,' *: The power to IB waves is also attributed',1X,
     +         'to electrons and ions',/,
     +  '    depending on the local absorption rates.')
C
C ************ ************ ************ ************ ************
C
      RETURN
C
C ************ ************ ************ ************ ************
C
 3000 CONTINUE
C
C ************ ************ ************ ************ ************
C
C  OUTPUT FROM TWO-ION HYBRID LAYERS
C
C ************ ************ ************ ************ ************
C
      IF(IMIN.LE.0)  THEN
         IF(MHARM.EQ.0)  RETURN
         WRITE(6,6000) (IHARM(J),J=1,MHARM)
 6000 FORMAT(10X,'First IC harmonic, species:',2I3)
      ELSE
         WRITE(6,6010)  (MINOR(J),J=1,IMIN)
 6010 FORMAT(10X,'Two-ion hybrid resonance, minority species:',2I3)
      END IF
C
      IF(IPRINT.LT.1)  RETURN
C
      IF(IDOPPL.NE.0)  WRITE(6,6020)
 6020 FORMAT(10X,'Warning: Doppler broadening too large')
C
      WRITE(6,6030)  HTOPT1,HTOPT2,REFTWH,TRATWH,ABSTWH
      WRITE(6,6040)  UPWABE
      WRITE(6,6050)  (UPWABI(I),I=1,NSPEC)
 6030 FORMAT(10X,'Optical thickness =',1P,2E12.3,0P,/,
     +       10X,'R =',F8.3,2X,'T =',F8.3,2X,'A =',F8.3)
 6040 FORMAT(10X,'Absorption by the electrons',1P,E12.3)
 6050 FORMAT(10X,'Absorption by the ions     ',
     +    1P,4E12.3,/50X,4E12.3)
C
      RETURN
C
C ************ ************ ************ ************ ************
C
C  Partial output during integration
C
C ************ ************ ************ ************ ************
C
 4000 CONTINUE
C
         THG = THETA*180./PI
      WRITE(6,4010)  EIKON,PSI,THG,XU,ZU,PWX
 4010 FORMAT(/,10X,'Information at Phase =',F7.3,/,
     +   ' Psi =',F7.3,' Th =',F7.2,' X =',F7.3,' Z =',F7.3,
     +         ' Pw =',1P,E11.3,0P)
      IF(IPRINT.GT.1)
     >   WRITE(6,4020)  DENS,TEMPE,TEMPIX(MAINSP),
     >                  TEMPIZ(MAINSP),BTOR,BPOL
 4020 FORMAT(' Ne =',1P,E11.3,' Te =',0P,F7.3,
     +       ' Ti =',2F7.3,' Btor =',F7.3,' Bpol =',F7.3)
      WRITE(6,4030)  PKPAR,PKP,PKPSI,PKX,PKZ
 4030 FORMAT(' Npar =',F7.3,' Nperp =',F7.3,' Npsi =',F8.3,
     +       ' NX =',F8.3,' NZ =',F8.3)
C
      IF(IPWDP.NE.0 .AND. IPRINT.NE.0) THEN
         IF(IPWDPE.NE.0)  THEN
            DO 410  I=1,3
  410          ZGAME(I) = GAME(I)/GAMMA
            WRITE(6,4040)  XE,(ZGAME(I),I=1,3)
 4040 FORMAT(' xe =',1P,E11.3,' Power(e): MP =',E11.3,
     +           ' ELD =',E11.3,' Mxd =',E11.3,0P)
         END IF
         DO 420  I=1,NSPEC
         IF(IRES(I).NE.0)  THEN
            IH = IRES(I)
            ZGAMI(I) = GAMI(I)/GAMMA
            WRITE(6,4050)  I,IRES(I),XI(IH,I),ZGAMI(I)
 4050 FORMAT(' Species',I2,': xi(',I1,') =',1P,E11.3,
     +            ' Power(i) =',E11.3,0P)
         END IF
  420    CONTINUE
      END IF
C
C  Detailed output
C
      IF(JHLF.GT.11 .OR. IPRINT.GT.1) THEN
         ZQFSQ = -ANR*ANL/ANS
         WRITE(6,4060) AR,AL,AS,ZQFSQ
 4060 FORMAT(' R   =',1P,E11.3,' L   =',E11.3,
     +       ' S   =',E11.3,' Nsq F  =',E11.3,0P)
         WRITE(6,4070) ARWARM,ALWARM,ATWARM,HA
 4070 FORMAT(' Rho =',1P,E11.3,' Lbd =',E11.3,
     +       ' Tau =',E11.3,' Sig =',E11.3)
         WRITE(6,4080) (Y(L),L=1,NDIM)
 4080 FORMAT(' Y (KPSI,PSI,KTH,TH,PW)',1P,5E12.3,0P)
         WRITE(6,4090) (DY(L),L=1,NDIM)
 4090 FORMAT(' DY(KPSI,PSI,KTH,TH,PW)',1P,5E12.3,0P)
         IF(IPWDP.NE.0)
     >       WRITE(6,4100)  (EPOL(I),I=1,3),EPOLX,EPOLY
 4100 FORMAT(' EPOL  =',1P,5E12.3,0P)
      END IF
C
      IF(IPRINT.GE.1)
     >   WRITE(6,4110)  DISPH,ZERR,HLAST
 4110 FORMAT(' DISPH  =',1P,E11.3,' Err  =',E11.3,
     +  ' Step   =',E11.3,0P)
C
      RETURN
      END
