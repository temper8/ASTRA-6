C----------------------------------------------------------------------|
	subroutine	ICRAY(PLOT)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
C ************ ************ ************ ************
C  IONSP - No. of ion species (can be extended up to 8 species)
C  NICR  - No. of radial grid points
	integer		NICR,IONSP,ierror,j,i,IPLOT,IPOW
	parameter	(NICR=51, IONSP=2)
	double precision	A(NRD),B(NRD),C(NRD),YNE(NICR)
	double precision	YTI(NICR,IONSP),YNI(NICR,IONSP)
	double precision	ARRAY(50),YTE(NICR),D(NRD)
	double precision	XTR(NRD),XPSI(NICR)
	double precision	ALFA,PLOT,VINT,POWSC
C	call	icrtr(	ARRAY,NICR,A,B,C,YNE,YTE,D,
C     &			IONSP,YNI,YTI,IPLOT,ierror)
C ************ ************ ************ ************
C  Major and plasma radius (m)
	ARRAY(1) = RTOR
	ARRAY(2) = ABC
C ************ ************ ************ ************
C  Toroidal magnetic field at the geometrical center (Tesla)
	ARRAY(3) = BTOR
C ************ ************ ************ ************
C  Toroidal current (A)
	ARRAY(4) = IPL
C ************ ************ ************ ************
C Grid size for ICRH module 
C		(same for the configuration and plasma parameters)
	do  J = 1,NICR
	    XPSI(J) = (J-1.)/(NICR-1.)
	enddo
	do  J = 1,NA1
	    XTR(J) = AMETR(J)/AMETR(NA1)
	enddo
	ALFA =.001
C =================================================================
C  MHD CONFIGURATION AND PROFILES
C =================================================================
C
C  Coefficients of the MHD configuration
C     NICR = N. of points in the radial profiles
C            of RMAJ(PSI), RRHO(PSI), TRIANG(PSI), ELLIPT(PSI)
C            accepted if 11 <= NICR <= 51
C
C ************ ************ ************ ************
C Shafranov shift, elongation and triangularity
C
	CALL	SMOOTH(ALFA,NA1,SHIF,XTR,NICR,A,XPSI)
	CALL	SMOOTH(ALFA,NA1,ELON,XTR,NICR,B,XPSI)
	CALL	SMOOTH(ALFA,NA1,TRIA,XTR,NICR,C,XPSI)
C ************ ************ ************ ************
C n_e, n_i, T_e, T_i, j
C
	CALL	SMOOTH(ALFA,NA1,NE,XTR,NICR,YNE,XPSI)
	CALL	SMOOTH(ALFA,NA1,TE,XTR,NICR,YTE,XPSI)
	CALL	SMOOTH(ALFA,NA1,CU,XTR,NICR,D,XPSI)
C	CALL	SMOOTH(ALFA,NA1,NI,XTR,NICR,YNI(1,1),XPSI)
	CALL	SMOOTH(ALFA,NA1,NHYDR,XTR,NICR,YNI(1,2),XPSI)
C	CALL	SMOOTH(ALFA,NA1,NHE3, XTR,NICR,YNI(1,3),XPSI)
	CALL	SMOOTH(ALFA,NA1,TI,XTR,NICR,YTI(1,1),XPSI)
	CALL	SMOOTH(ALFA,NA1,TI,XTR,NICR,YTI(1,2),XPSI)
C	CALL	SMOOTH(ALFA,NA1,TI,XTR,NICR,YTI(1,3),XPSI)
C ************ ************ ************ ************
C  mass and charge (A.U.)
C ************ ************ ************ ************
C  Reference (majority) species 
C			is adjusted to implement charge neutrality
	ARRAY(5) = AMJ
	ARRAY(6) = ZMJ
C ************ ************ ************ ************
C  mass and charge for minority ion species
C		(density distributions can be taken from the arrays:
C		 NHYDR, NDEUT, NHE3, NALF, NIM1, NIM2, NIM3)
C ************ ************ ************ ************
C  mass and charge for 1st minority
	ARRAY(7) = 1.
	ARRAY(8) = 1.
C ************ ************ ************ ************
C  mass and charge for 2nd minority (not used if IONSP < 3)
	ARRAY(9) = 2.
	ARRAY(10) = 1.
C ************ ************ ************ ************
C  Further minorities (not used if IONSP < 4)
C  	ARRAY(j) elements 11<= j <= 20 
C  	are reserved for mass and charge of further minorities
C ************ ************ ************ ************
C  Majoriry species is determined to fulfil quasineutrality
	do  i=1,IONSP
	    do  j=1,NICR
		if (i .eq. 1)  then
		   YNI(j,1) = YNE(j)
		else
		   YNI(j,1) = YNI(j,1)-ARRAY(6+2*i)*YNI(j,i)
		endif
	    enddo
	enddo

C =================================================================
C  WAVE FREQUENCY, POWER, AND ANTENNA PARAMETERS
C =================================================================
C
C  Applied frequency (Hz)
	ARRAY(21) = 1.E6*FICR
C ************ ************ ************ ************
C  Higher harmonics contributiong to the absorption
C    (Cannot exceed 5, default value = 2, first harmonic)
        ARRAY(22) = 2.
C ************ ************ ************ ************
C  Total power coupled (MW)
	ARRAY(23) = QICR
C ************ ************ ************ ************
C  Number of antennas. If 0 the antenna evaluation is skipped.
	ARRAY(24) = 1.
C ************ ************ ************ ************
C  Phase difference between antennas
	ARRAY(25) = 0.
C ************ ************ ************ ************
C  Symmetry of the antenna:
C     1 - monopole antenna
C     2 - dipole antenna
	ARRAY(26) = 2.
C ************ ************ ************ ************
C  Feeders position:
C     1 - Feeders up and down, equatorial short
C     2 - Shorts up and down, equatorial feeder
	ARRAY(27) = 1.
C ************ ************ ************ ************
C  Distance antenna plasma [m]
	ARRAY(28) = 0.07
C  Distance antenna wall [m]
	ARRAY(29) = 0.14
C  z-width of the antenna [m]
	ARRAY(30) = 0.18
C  half y-height of an antenna element [m]
	ARRAY(31) = 0.20
C  gap between conductors in the dipole case [m]
	ARRAY(32) = 0.50
C ************ ************ ************ ************
C  Poloidal position of the antenna middle point (degs)
C     0 - outer equatorial plane
	ARRAY(33) = 0.
C ************ ************ ************ ************
C  Phase between conductors in the dipole case
C  Phase between feeders in the monopole JALIM=1 case
	ARRAY(34) = 180.
C ************ ************ ************ ************
C  Effective propagation constant along the antenna
	ARRAY(35) = 0.
C ************ ************ ************ ************
C  Output of the antenna package:
C     -1 - Only antenna evaluation
C      0 - No details of the antenna evaluation
C      1 - Output of details of the antenna evaluation
	ARRAY(36) = 1
C
C =================================================================
C  CONTROL OF RAY-TRACING AND OUTPUT
C =================================================================
C
C  Evaluation of the antenna power spectre:
C     No. of values of NPHI
	ARRAY(37) = 10.
C     1st value of NPHI
	ARRAY(38) = 0.
C     step of NPHI (ignored IF ARRAY(37)=1)
	ARRAY(39) = 4.
C ************ ************ ************ ************
C  No. of rays in the poloidal cross-section
C          (should not exceed NTHRAY in COMMON.F)
	ARRAY(40) = 20.
C ************ ************ ************ ************
C  Step of EIKON (wave phase) between output
C           (only if ARRAY(37) = 1)
	ARRAY(41) = 0.5
C ************ ************ ************ ************
C  No. of horizontal reflections allowed
	ARRAY(42) = 2.
C ************ ************ ************ ************
C  Details in the output:
C     0 - output reduced;
C     1 - normal output;
C     2 - detailed output.
C     3 - output at each call of OUTPUT(X)
C                  for tests only (no plots)
	ARRAY(43) = 0.
C ************ ************ ************ ************
C  Graphical output: Ray tracing
C     Plots of ray tracing is suppressed by IPLOT = 0
	IPLOT = PLOT
	ierror = 0

C ************ ************ ************ ************
	write(*,'(">>> ICRH Entry  >>>  ",$)')
C	open(unit=6,file='/dev/null',STATUS='OLD')
	open(unit=6,file='tmp/raytr.out',err=91,STATUS='UNKNOWN')
C =================================================================
	call	icrtr(	ARRAY,NICR,A,B,C,YNE,YTE,D,
     &			IONSP,YNI,YTI,IPLOT,ierror)
C =================================================================
	close(6)
	open(unit=6,file='/dev/tty',STATUS='OLD')
C ************ ************ ************ ************
	if (ierror .eq. 0)write(*,'("ICRH >>> Normal exit",$)')
	if (ierror .ne. 0)write(*,*)'ICRH >>> Error #',ierror
C ************ ************ ************ ************
C Power deposition to electrons
	CALL	SMOOTH(ALFA,NICR,YNE,XPSI,NA1,PEICR,XTR)
C ************ ************ ************ ************
C Power deposition to 1st ion species
	CALL	SMOOTH(ALFA,NICR,YNI(1,1),XPSI,NA1,A,XTR)
	do	j=1,NA1
	    PIICR(j) = A(j)
	enddo
C ************ ************ ************ ************
C Power deposition to 2nd ion species
	CALL	SMOOTH(ALFA,NICR,YNI(1,2),XPSI,NA1,A,XTR)
	do	j=1,NA1
	    PIICR(j) = PIICR(j)+A(j)
	enddo
C ************ ************ ************ ************
C Power deposition to 3rd ion species
C	CALL	SMOOTH(ALFA,NICR,YNI(1,3),XPSI,NA1,A,XTR)
C	do	j=1,NA1
C	    PIICR(j) = PIICR(j)+A(j)
C	enddo
C ************ ************ ************ ************
C Note !!! All power deposition profiles are scaled so that
C	   the total absorbed power is equal to QICR
	POWSC = (VINT(PEICR,ROC)+VINT(PIICR,ROC))/QICR
	if (POWSC .lt. 0.6)	then
		IPOW = 100*POWSC
		write(*,'(" -> First pass absorption ",1I2,"%")')IPOW
	else
		write(*,*)
	endif
	do	j=1,NA1
	    PEICR(j) = PEICR(j)/POWSC
	    PIICR(j) = PIICR(j)/POWSC
	enddo

	return
 91	open(unit=6,file='/dev/tty',STATUS='OLD')
	write(*,*)'>>> ICRH >>> File open error'
	end
