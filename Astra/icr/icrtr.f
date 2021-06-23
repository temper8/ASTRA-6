C ************ ************ ************ ************ ************
C *                                                              *
C *     ION CYCLOTRON RAY TRACING CODE RAYIC   -   15/04/93      *
C *                                                              *
C ************ ************ ************ ************ ************
	subroutine	icrtr	(
C Input:
     &		ARRAY,NRAD,ASHIF,AELON,ATRIA,APNE,APTE,ACUR,
     &		NIONSP,APNI,APTI,IPLOT,
C Output:
C     &		APTE,APTI
     &			ierror	)
	integer	NRAD,NIONSP,IPLOT,ierror
	real	ARRAY(*),ASHIF(NRAD),AELON(NRAD),ATRIA(NRAD)
	real	APNE(NRAD),APTE(NRAD),ACUR(NRAD)
	real	APNI(NRAD,NIONSP),APTI(NRAD,NIONSP)
      include 'COMMON.F'
	if (NRAD.gt.51 .or. NRAD.lt.11)	then
	    ierror = 1
	    return
	endif
C	write(*,'(5I5)')NRAD,NIONSP,IPLOT,ierror
C	write(*,'(10F7.2)')(ARRAY(j),j=1,10)
C	write(*,'(10F7.2)')(ASHIF(j),j=1,NRAD)
C	write(*,'(10F7.2)')(AELON(j),j=1,NRAD)
C	write(*,'(10F7.2)')(ATRIA(j),j=1,NRAD)
C	write(*,'(10F7.2)')(APNE(j),j=1,NRAD)
C	write(*,'(10F7.2)')(APTE(j),j=1,NRAD)
C	write(*,'(10F7.2)')(ACUR(j),j=1,NRAD)
C	write(*,'(10F7.2)')(APTI(j),j=1,NRAD)
C	write(*,'(10F7.2)')((APNI(j,i),j=1,NRAD),i=1,NIONSP)
C =================================================================
C  MHD CONFIGURATION AND PLASMA PARAMETERS
C =================================================================
         NEQUIL = 1
         RTORUS = ARRAY(1)
         RPLASM = ARRAY(2)
         BZERO = ARRAY(3)
         TCURR = ARRAY(4)*1.E6
         NMHD = NRAD
         ZDPSI = RPLASM/FLOAT(NMHD-1)
      DO 31  J=1,NMHD
C ************ ************ ************ ************
         RMAJ(J) = RTORUS + ASHIF(J)
         ELONG(J) = AELON(J)
   31    TRIANG(J) = ATRIA(J)
C
C ************ ************ ************ ************
C  N.of ion species (up to 8 species allowed) presently 3 !!!
         NSPEC = NIONSP
C ************ ************ ************ ************
C  Reference species (charge neutrality and output)
         MAINSP = 1
C ************ ************ ************ ************
C  mass and charge (A.U.)
         ATM(1) = ARRAY(5)
         AZI(1) = ARRAY(6)
         ATM(2) = ARRAY(7)
         AZI(2) = ARRAY(8)
C         ATM(3) = ARRAY(9)
C         AZI(3) = ARRAY(10)
C ************ ************ ************ ************
C  Concentrations (n(i)/n(e); charge neutrality is
C     implemented by reevaluating ACONC(MAINSP)
         ACONC(2) = 0.05
C        ACONC(3) = 0.005
C ************ ************ ************ ************
C  Density and temperature profiles
C     NPROF = N. of points in the radial profiles
C             accepted if 11 <= NPROF <= 51
C
         NPROF = NRAD
         ACONC(MAINSP) = 1./AZI(MAINSP)
      DO 61  I=1,NSPEC
      IF(I.EQ.MAINSP)  GO TO 61
         ACONC(MAINSP) = ACONC(MAINSP) - AZI(I)*ACONC(I)/AZI(MAINSP)
   61    continue
      DO 71  J=1,NPROF
C ************ ************ ************ ************
C  Normalized "radial" variable
         PPSI(J) = FLOAT(J-1)/FLOAT(NPROF-1)
         PNE(J) = 1.E13*APNE(J)
         PTE(J) = APTE(J)
         PJAVG(J) = ACUR(J)
      DO 71  I=1,NSPEC
C         PNI(J,I) =  ACONC(I)*PNE(J)
         PNI(J,I) =  APNI(J,I)
C  In this version the ion temperature profiles are all equal
   71    PTI(J,I) = APTI(J,I)
C ************ ************ ************ ************
C  Central ion temperatures (keV)
      DO 72  I=1,NSPEC
 72         TEMPIC(I) = PTI(1,I)
C
C ************ ************ ************ ************
C  Density and temperature profiles
C     IPROEQ = 0 - Separate profiles for each species
C     IPROEQ = 1 - Same profile for all species
         IPROEQ = 0





C =================================================================
C  WAVE FREQUENCY, POWER, AND ANTENNA PARAMETERS
C =================================================================
C
C  Alternative procedures for the frequency choice are not applicable
C  Applied frequency (hz)
         FREQCY = ARRAY(21)
C     IRESP = 0 - frequency and magnetic field given;
         IRESP = 0
C ************ ************ ************ ************
C  Higher harmonics contributiong to the absorption
C    (Cannot exceed 5, default value = 2, first harmonic)
         MAXHRM = ARRAY(22)+.001
C ************ ************ ************ ************
C  Total power coupled (MW)
         POWER = ARRAY(23)
C ************ ************ ************ ************
C  NANTS : Number of antennas
C     if NANTS = 0 the antenna evaluation is skipped.
         NANTS = ARRAY(24)+.001
C ************ ************ ************ ************
C  DPHASE: phase difference between antennas
         DPHASE = ARRAY(25)
C ************ ************ ************ ************
C  JPOLE : Symmetry of the antenna:
C     JPOLE = 1 - monopole antenna
C     JPOLE = 2 - dipole antenna
         JPOLE = ARRAY(26)+.001
C ************ ************ ************ ************
C  JALIM : feeders position:
C     JALIM = 1 - Feeders up and down, equatorial short
C     JALIM = 2 - Shorts up and down, equatorial feeder
         JALIM = ARRAY(27)+.001
C ************ ************ ************ ************
C  DISTAP: distance antenna plasma
C  DISTAW: distance antenna wall;
C  WIDTH : z-width of the antenna
C  HEIGTH: half y-height of an antenna element
C  WGAP  : gap between conductors in the dipole case
C          (all lengths in meters)
         DISTAP = ARRAY(28)
         DISTAW = ARRAY(29)
         WIDTH = ARRAY(30)
         WGAP = ARRAY(31)
         HEIGTH = ARRAY(32)
C ************ ************ ************ ************
C  Poloidal position of the antenna middle point (degs)
C     THANTN = 0 - outer equatorial plane
         THANTN = ARRAY(33)
C ************ ************ ************ ************
C  Phase between conductors in the dipole case
C  Phase between feeders in the monopole JALIM=1 case
         EPHASE = ARRAY(34)
C ************ ************ ************ ************
C  Effective propagation constant along the antenna
         ANTKY = ARRAY(35)
C ************ ************ ************ ************
C  Output of the antenna package:
C     JOUTA = -1 - Only antenna evaluation
C     JOUTA =  0 - No details of the antenna evaluation
C     JOUTA =  1 - Output of details of the antenna evaluation
         JOUTA = ARRAY(36)+.001
C
C =================================================================
C  CONTROL OF RAY-TRACING AND OUTPUT
C =================================================================
C
C  Evaluation of the antenna power spectre:
C     NUMPHI - N. of values of NPHI
C     NPHI1  - first value of NPHI
C     JUMPHI - step of NPHI (ignored IF NUMPHI=1)
         NUMPHI = ARRAY(37)+.001
         NPHI1  = ARRAY(38)+.001
         JUMPHI = ARRAY(39)+.001
C ************ ************ ************ ************
C  NTHIN - N. of rays in the poloidal cross-section
C          (should not exceed NTHRAY)
         NTHIN = ARRAY(40)+.001
C ************ ************ ************ ************
C  DEYKPR - Step of EIKON (wave phase) between output
C           (only if NUMPHI = 1)
         DEYKPR = ARRAY(41)
C ************ ************ ************ ************
C  MAXREF = N. of horizontal reflections allowed
         MAXREF = ARRAY(42)+.001
C ************ ************ ************ ************
C  ISMOOT = 0  Power deposition at the closest mesh point
C  ISMOOT = 1  power deposition smoothed over one wavelength
C              (to be checked)
         ISMOOT = 1
C ************ ************ ************ ************
C  Details in the output:
C     IPRINT = 0 - output reduced;
C     IPRINT = 1 - normal output;
C     IPRINT = 2 - detailed output.
C     IPRINT = 3 - output at each call of OUTPUT(X)
C                  for tests only (no plots)
         IPRINT = ARRAY(43)+.001
C ************ ************ ************ ************
C  Graphical output: density and temperature profiles
C     IPLOPR = 1  -  profiles plotted;
C     IPLOPR = 0  -  profiles not plotted;
C     IPLOPR = -1 -  plots of the magnetic surfaces only
C                   (no ray tracing, for tests)
         IPLOPR = 0
C ************ ************ ************ ************
C  Graphical output: Ray tracing
C     Plots of ray tracing is suppressed by IGRAPH = 0
         IGRAPH = IPLOT
C ************ ************ ************ ************
C  List of plots to be made:
      DO 10 NX=1,2
      DO 10 NY=1,7
   10    IPLOXY(NX,NY) = 0
      IF(IGRAPH.EQ.0)  GO TO 20
C ===================================================================
C  NPOLXY(NXVARB,NYVARB) = 1 produces a plot
C                        according to the following table:
C ====================================================================
C  NXVARB=1  Indep. var. SEIKON  (wave phase in the poloidal plane)
C         2              PSI (magnetyic surface label)
C ====================================================================
C  NYVARB=1  Dep. var.   SKP
C         2              SKPAR
C         3              SKR
C         4              SKTH (M-THETA)
C         5              SKX
C         6              SKZ
C         7              SPOWER
C ===================================================================
         IPLOXY(2,1) = 1
         IPLOXY(2,2) = 1
C        IPLOXY(2,5) = 1
         IPLOXY(2,6) = 1
         IPLOXY(2,7) = 1
C        IPLOXY(1,1) = 1
C        IPLOXY(1,2) = 1
C        IPLOXY(1,7) = 1
   20 CONTINUE
C
C ************ ************ ************ ************
C  Title & date
C
         ZFREQ = 1.E-6*FREQCY
      WRITE (6,5000) ZFREQ,NPHI1
 5000 FORMAT(/,20X,'ION CYCLOTRON RAY-TRACING CODE RAYIC',/,
     +       20X,'           April 18 1994          ',/,
     +       20X,'     f =',F6.2,' Mhz - n-phi =',I4)
C
C ************ ************ ************ ************
C  Call to the ray-tracing code
C
            CALL RAYIC
      do   j=1,NPROF
         APNE(j) = PWE(j)
      do   i=1,NSPEC
         APNI(j,i) = PWI(j,i)
      enddo
      enddo
C
C ************ ************ ************ ************ ************
C
      return
      END
