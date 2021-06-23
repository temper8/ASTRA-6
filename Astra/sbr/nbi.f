C======================================================================|
	subroutine	NBI
C----------------------------------------------------------------------|
C Interface to the Neutral Beam Injection package by A.R.Polevoi
C						(Edition 19-APR-2000)
C 
C This interface is compatible with the Astra version 5.3 and later
C						(Pereverzev 27-09-01)
C----------------------------------------------------------------------|
C The subroutine NBI is called from ASTRA model and provides
C	the multisource neutral beam power, momentum and driven current
C	deposition including ripple losses.
C	Optionally, the time dependent Fokker-Planck solver (NBIBCE.f) 
C       can be used
C----------------------------------------------------------------------|
C Input:
C	The input is split in three parts. (See also description at the end)
C	  (1)	The tokamak parameters:
C		(i)  The NBI subroutine input parameters. These parameters: 
C		     CNB1,CNB2,CNB3,CNB4,CNBI1,CNBI2,CNBI3,CNBI4 either
C		     have to be set explicitly in this subroutine (see
C		     an example below) or by the data file "exp/DATA_FILE"
C		     or in the calling model. 
C		(ii) Parameters to user's function RIPRAD (sbr/nbuser.f).
C	  (2)	Beam box parameters:
C		(i)  The input via the file FILENA (the file name is
C		     constructed as "exp/DATA_FILE.nbi"),
C		(ii) User's functions NBFRY, NBFHZ (see file sbr/nbuser.f).
C	  (3)	Input plasma profiles as TE, TI, NE and so on are 
C		transferred through common blocks.
C----------------------------------------------------------------------|
C Output:	QNBI [MW] NBI power
C		PBEAM,PEBM,PIBM,NIBM,CUFI,CUBM,PBLON,PBPER,
C		SCUBM,SNEBM,SNNBM,NNBM1,2,3	for MAIN
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/outcmn.inc'
	integer		jj,length
C----------------------------------------------------------------------|
	if (NBFILE(1:1).eq.'*')	goto	99
	jj = length(NBFILE)
C	write(*,*)NBFILE(1:jj)
C	CNB1  = 8
C	CNB2  = 1   ? Explicit form of CX losses ?
	call NBINJ(
     +		NBFILE(1:jj),BTOR,RTOR,ABC,AB,ROC,SHIFT,UPDWN,
     +		HRO,TAU,NA1,NB1,AIM1,AIM2,AIM3,AMJ,ZMJ,NNCL,NNWM,QNBI,
C Any allowed control parameters can be used below
     +		CNB1,	! No. of NB sources			     (8)
     +		CNB2,	! bulk-ion CX losses:  1/0 on/off	     (1)
     +		CNB3,	! NBI space grid: N1=(NA1-1)/CNB3+1	     (2)
     +		CNB4,	! 1(2) steady st. (time dep.) FP solver	     (1)
C	     the next 4 parameters are effective if CNB4=2 only
     +		CNBI1,	! fast-ion CX losses due to cold neutrals
     +		CNBI2,	! FP solver time step: TAUNBI=CNBI2*TAU	
     +		CNBI3,  ! fast-ion CX losses due to NBI  neutrals
     +		CNBI4)	! FP solver velocity grid: IV1=80/CNBI4+1
	NBFLAG = 1	! Used for footprint drawing
	return
 99	write(*,*)'>>> NBI Error >>> Configuration file not found'
	stop
	end
C======================================================================|
C		Control parameters: (as described by A.R.Polevoi)
C CNB1	A number of PINIs with different geometry.		(Def=1)
C	When CNB1 is negative or changes its value during a run then 
C	a pop-up menu appears for interactive setting PINI parameters
C
C CNB2	Control of CX bulk-ion heat losses due to NBI		(Def=1)
C	CNB2 > 0 (.and. CNB4.eq.1) explicit treatment of CX heat losses
C	that is formula PBICX is included in PIBM. This corresponds 
C	to the explicit calculations of the heat sink from the bulk 
C	ion component at the moment when NBI is called.
C	= 0 No losses from bulk ions due to NBI charge exchange
C	This setting supposes that the implicit form is used:
C				CNB2=0;	PIT=...-.0024*SNNBM
C CNB3	=1,2... NBI X-mesh points number N1=(NA1-1)/CNBI3+1
C	NOTE: if the fast ion gyroradius Rg(EBEAM, ABEAM)
C	exceeds the Astra space step, i.e. Rg > ABC/NA1, 
C	then the internal NBI solver space step will inrcease 
C	automatically: CNB3 = Rg*NA1/ABC: NA1/CNB3,max = 5
C
C		Fast ion solver control:
C	CNBI1,2,3,4 are active for CNB4=2 only
C CNB4	=1 (2) steady (time dependent) Fokker-Planck Solver	(Def=1)
C
C CNBI1	=1 losses (=0 No losses) from fast ions due to CX	(Def=1)
C	with the cold neutrals
C

C CNBI2	=TAUNBI/TAU ratio of the NBI time-step to TAU	(Def=1)
C	WARNING: TAUNBI does not coincide automatically with the 
C	interval between sequental calls of NBI (that saves
C	the calculation time at the NBI-steady state phase) 
C CNBI3	=1 losses (=0 No losses) from fast ions due to CX
C	with the NBI neutrals itself
C	The same CX cross-section as for cold neutrals is used.
C
C
C CNBI4	= free parameter. 	            (Def=1)
C (At present number of mesh points in V space is fixed IV1=161)
C
C----------------------------------------------------------------------|
C Unlike the previous versions a number of parameters have to be
C	defined in the input file exp/data_file_name.nbi (FILENA)
C       This file provides a description of each beam line and 
C       includes N groups, where N is the total number of beam lines.
C
C Each group consists of 21 records. 
C The first record is the ordinal group number in a separate line.
C The rest of a group consisits of 20 character*12 fields.
C Every field can be either a real number or a name of variable.
C The allowed variables are: 
C   "ZRDn" (or "ZRDnX"), where n stands for an integer number 1<=n<=48
C   "Astra_Constant" (see the full list in "for/const.inc" 
C                     or just press "C" in the run mode).
C 
C Significance of input parameters in each group:
C
C Ordinal_beam_number
C QBEAM  [MW] Beam power
C CONTR       Counter injection fraction (0 co, 1 counter)
C ABEAM [m_p] Mass of beam ions in the proton mass
C ZBEAM       Beam ion charge in the proton charge units
C EBEAM [keV] Beam energy
C DBM1        EBEAM power fraction
C DBM2        EBEAM/2 power fraction
C DBM3        EBEAM/3 power fraction
C Orb_av      Type of averaging over ion orbits
C		  =0 No averaging (deposition at the birth point)
C		  =1 Averaging with a finite orbit width
C	    	  =2 Averaging with zero orbit width
C Penc_num    Number of pencils in the horizontal plane
C HBEAM   [m] Beam footprint center height
C RBMAX   [m] Beam footprint maximum radius
C RBMIN   [m] Beam footprint minimum radius
C tg(A)       A is the angle between the beam and the midplane
C Aspect      Footprint aspect ratio: Beam_height=Aspect*(RBMAX-RBMIN)
C Cver1       These parameters describe exponential (or any other)
C Cver2         beam power distribution across the beam cross-section
C Chor1         as described by the user functions NBFRY & NBFHZ
C Chor2         (see file sbr/nbuser.f)
C
C Parameter meaning:
C	QBEAM [MW] NBI power 					(Def=0)
C	CONTR=0-1 fraction of contr-inj. power			(Def=1)
C	ABEAM=1,2,3 Hydrogen,Deuterium,Tritium NBI only		(Def=1)
C	EBEAM [keV] the main component energy			(Def=1)
C	DBM1=0-1 ~ full energy component EBEAM power fraction 	(Def=1)
C	DBM2=0-1 ~ half energy component EBEAM/2 power fraction (Def=0)
C	DBM3=0-1 ~third energy component EBEAM/2 power fraction (Def=0) 
C Warning: DBM1 > 0 must be used;    SUM(DBMi).ne.1. is allowed
C	Power is renormalized withinNBI: Power(i)~DBMi/SUM(DBMi)	
C	HBEAM [m] NBI footprint center height in respect to Z=0	(Def=0)
C	RBMAX [m] max. major radius of the NBI footprint	(Def=0)
C	RBMIN [m] min. major radius of the NBI footprint	(Def=0)
C
Comments: the footprint position is interpreted as:
C
C	for tangential NBI: (RBMAX+RBMIN)/2 > RTOR - AB
C	as NBI cross-section by the meridianal plane perpendicular 
C	to the vertical plane of the central NBI pencil beam
C	with tangential radius of (RBMAX+RBMIN)/2 
C
C	for perpendicular NBI: (RBMAX+RBMIN)/2 < RTOR - AB 
C	as NBI cross-section with the vertical plane perpendicular to 
C	the central NBI pencil beam (RBMAX+RBMIN)/2 vertical plane,
C	placed at R = RTOR (the corrections, connected with 
C	Shafranov's shift are calculated into the NBI block
C
C	Penc_num  Each PINI is approxomated by Nh x Nr pencils	(Def=1)
C		where Nr = Penc_num (horizontal direction),
C		Nh (vertical direction) is automatically 
C		calculated as a number of magn. surfaces
C		in the vertical aperture of each PINI in 
C		the footprint crossection.
C	
C	tg(A)	tg(angle between central pencil and midplane)	(Def=1)
C		is used for perpendicular NBI only		
C	Aspect =(HBmax-HBmin)/(RBMAX-RBMIN): HBEAM=(HBmax+HBmin)/2(Def=1)
C
C User's functions
C
C	RIPRAD(ZUPDWN,J) [m] - Ripple loss cone boundary (major 
C		radius) for each magnetic surface: ZUPDWN [m]- 
C		shift in respect to the midplane of ripple 
C		simmetry, J - surface index
C	RIPRAD -depends on tokamak mag. field coils' and 
C		plasma configurations  (Def: No Ripple losses) (Def=999)
C======================================================================|
C Parameter mapping: 
C      External name (Astra_5.3) <-> Internal name (NBI and old Astra)
C			CNB1	 <->	CBM1
C			CNB2	 <->	CBM2
C			CNB3	 <->	CBMI3
C			CNB4	 <->	CBMI1
C			CNBI1	 <->	CBM4
C			CNBI2	 <->	CBMI2
C			CNBI3	 <->	CBM3
C			CNBI4	 <->	CBMI4
C======================================================================|
