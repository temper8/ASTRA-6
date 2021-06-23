C 	ASTRA constant name list
C		Each group below: 
C			Variables, Constans, Internal program parameters
C			cannot include more than NCONST=256 components
C			model compiler allows extension up to 256
C Variables
C AB            [m] - Limiter position in the mid-plane
C ABC           [m] - Calculation boundary in the mid-plane
C AIM1          [ ] - Impurity 1 ions atom mass/Mp
C AIM2          [ ] - Impurity 2 ions atom mass/Mp
C AIM3          [ ] - Impurity 3 ions atom mass/Mp
C AMJ           [ ] - (Main ions atom mass)/Mp
C AWALL         [m] - Wall position in horisontal axis
C BTOR          [T] - Toroidal field at vacuum chamber center
C ELONG         [ ] - Elongation of boundary surface
C ELONM         [ ] - Chamber elongation
C ENCL        [keV] - Edge cold neutrals energy
C ENWM        [keV] - Edge warm neutrals energy
C FECR        [GHz] - EC frequency
C FFW         [GHz] - FW frequency
C FICR        [GHz] - IC frequency
C FLH         [GHz] - LH frequency
C GN2E          [ ] - Factor in convective electron heat flux
C GN2I          [ ] - Factor in convective ion heat flux
C IPL          [MA] - Plasma current
C LEXT        [mkH] - External inductance
C NNCL  [10^19/m^3] - Edge cold neutrals density
C NNWM  [10^19/m^3] - Edge warm neutrals density
C QECR         [MW] - Total EC power input
C QFW          [MW] - Total FW power input
C QICR         [MW] - Total IC power input
C QLH          [MW] - Total LH power input
C QNBI         [MW] - Total beam power
C RTOR          [m] - Major radius (vacuum chamber center)
C SHIFT         [m] - Edge magnetic surface shift
C TRIAN         [ ] - Triangularity of the boundary surface
C TRICH         [ ] - Chamber triangularity
C UEXT          [V] - External voltage
C UPDWN         [m] - Plasma up- or down-shift
C WNE           [m] - Width of NE exponential decay in SOL (j > NA1)
C WTE           [m] - Width of TE exponential decay in SOL (j > NA1)
C WTI           [m] - Width of TI exponential decay in SOL (j > NA1)
C ZMJ           [ ] - Main-ion charge
C ZRD1          [ ] - Reserved variable for reading from "exp" file
C ZRD10         [ ] - Reserved variable for reading from "exp" file
C ZRD11         [ ] - Reserved variable for reading from "exp" file
C ZRD12         [ ] - Reserved variable for reading from "exp" file
C ZRD13         [ ] - Reserved variable for reading from "exp" file
C ZRD14         [ ] - Reserved variable for reading from "exp" file
C ZRD15         [ ] - Reserved variable for reading from "exp" file
C ZRD16         [ ] - Reserved variable for reading from "exp" file
C ZRD17         [ ] - Reserved variable for reading from "exp" file
C ZRD18         [ ] - Reserved variable for reading from "exp" file
C ZRD19         [ ] - Reserved variable for reading from "exp" file
C ZRD2          [ ] - Reserved variable for reading from "exp" file
C ZRD20         [ ] - Reserved variable for reading from "exp" file
C ZRD21         [ ] - Reserved variable for reading from "exp" file
C ZRD22         [ ] - Reserved variable for reading from "exp" file
C ZRD23         [ ] - Reserved variable for reading from "exp" file
C ZRD24         [ ] - Reserved variable for reading from "exp" file
C ZRD25         [ ] - Reserved variable for reading from "exp" file
C ZRD26         [ ] - Reserved variable for reading from "exp" file
C ZRD27         [ ] - Reserved variable for reading from "exp" file
C ZRD28         [ ] - Reserved variable for reading from "exp" file
C ZRD29         [ ] - Reserved variable for reading from "exp" file
C ZRD3          [ ] - Reserved variable for reading from "exp" file
C ZRD30         [ ] - Reserved variable for reading from "exp" file
C ZRD31         [ ] - Reserved variable for reading from "exp" file
C ZRD32         [ ] - Reserved variable for reading from "exp" file
C ZRD33         [ ] - Reserved variable for reading from "exp" file
C ZRD34         [ ] - Reserved variable for reading from "exp" file
C ZRD35         [ ] - Reserved variable for reading from "exp" file
C ZRD36         [ ] - Reserved variable for reading from "exp" file
C ZRD37         [ ] - Reserved variable for reading from "exp" file
C ZRD38         [ ] - Reserved variable for reading from "exp" file
C ZRD39         [ ] - Reserved variable for reading from "exp" file
C ZRD4          [ ] - Reserved variable for reading from "exp" file
C ZRD40         [ ] - Reserved variable for reading from "exp" file
C ZRD41         [ ] - Reserved variable for reading from "exp" file
C ZRD42         [ ] - Reserved variable for reading from "exp" file
C ZRD43         [ ] - Reserved variable for reading from "exp" file
C ZRD44         [ ] - Reserved variable for reading from "exp" file
C ZRD45         [ ] - Reserved variable for reading from "exp" file
C ZRD46         [ ] - Reserved variable for reading from "exp" file
C ZRD47         [ ] - Reserved variable for reading from "exp" file
C ZRD48         [ ] - Reserved variable for reading from "exp" file
C ZRD5          [ ] - Reserved variable for reading from "exp" file
C ZRD6          [ ] - Reserved variable for reading from "exp" file
C ZRD7          [ ] - Reserved variable for reading from "exp" file
C ZRD8          [ ] - Reserved variable for reading from "exp" file
C ZRD9          [ ] - Reserved variable for reading from "exp" file
C End variables		  - 85 entries (must be <=NCONST=256)
     	double precision
     1		AB,    ABC,   AIM1,  AIM2,  AIM3,  AMJ,   AWALL, BTOR,
     2		ELONG, ELONM, ENCL,  ENWM,  FECR,  FFW,   FICR,  FLH,
     3		GN2E,  GN2I,  IPL,   LEXT,  NNCL,  NNWM,  QECR,  QFW,
     4		QICR,  QLH,   QNBI,  RTOR,  SHIFT, TRIAN, TRICH, UEXT,
     5		UPDWN, WNE,   WTE,   WTI,   ZMJ,
     6		ZRD1,  ZRD10, ZRD11, ZRD12, ZRD13, ZRD14, ZRD15, ZRD16,
     8		ZRD17, ZRD18, ZRD19, ZRD2,  ZRD20, ZRD21, ZRD22, ZRD23,
     9		ZRD24, ZRD25, ZRD26, ZRD27, ZRD28, ZRD29, ZRD3,  ZRD30,
     &		ZRD31, ZRD32, ZRD33, ZRD34, ZRD35, ZRD36, ZRD37, ZRD38,
     1		ZRD39, ZRD4,  ZRD40, ZRD41, ZRD42, ZRD43, ZRD44, ZRD45,
     2		ZRD46, ZRD47, ZRD48, ZRD5,  ZRD6,  ZRD7,  ZRD8,  ZRD9,
     3		RESA1(NCONST-85)
     	double precision
     1		ABX,   ABCX,  AIM1X, AIM2X, AIM3X, AMJX,  AWALLX,BTORX,
     2		ELONGX,ELONMX,ENCLX, ENWMX, FECRX, FFWX,  FICRX, FLHX,
     3		GN2EX, GN2IX, IPLX,  LEXTX, NNCLX, NNWMX, QECRX, QFWX,
     4		QICRX, QLHX,  QNBIX, RTORX, SHIFTX,TRIANX,TRICHX,UEXTX,
     5		UPDWNX,WNEX,  WTEX,  WTIX,  ZMJX,
     6		ZRD1X, ZRD10X,ZRD11X,ZRD12X,ZRD13X,ZRD14X,ZRD15X,ZRD16X,
     8		ZRD17X,ZRD18X,ZRD19X,ZRD2X, ZRD20X,ZRD21X,ZRD22X,ZRD23X,
     9		ZRD24X,ZRD25X,ZRD26X,ZRD27X,ZRD28X,ZRD29X,ZRD3X, ZRD30X,
     &		ZRD31X,ZRD32X,ZRD33X,ZRD34X,ZRD35X,ZRD36X,ZRD37X,ZRD38X,
     1		ZRD39X,ZRD4X, ZRD40X,ZRD41X,ZRD42X,ZRD43X,ZRD44X,ZRD45X,
     2		ZRD46X,ZRD47X,ZRD48X,ZRD5X, ZRD6X, ZRD7X, ZRD8X, ZRD9X,
     3		RESA1X(NCONST-85)
     	common /A_VARIABLES/
     1		AB,    ABC,   AIM1,  AIM2,  AIM3,  AMJ,   AWALL, BTOR,
     2		ELONG, ELONM, ENCL,  ENWM,  FECR,  FFW,   FICR,  FLH,
     3		GN2E,  GN2I,  IPL,   LEXT,  NNCL,  NNWM,  QECR,  QFW,
     4		QICR,  QLH,   QNBI,  RTOR,  SHIFT, TRIAN, TRICH, UEXT,
     5		UPDWN, WNE,   WTE,   WTI,   ZMJ,
     6		ZRD1,  ZRD10, ZRD11, ZRD12, ZRD13, ZRD14, ZRD15, ZRD16,
     8		ZRD17, ZRD18, ZRD19, ZRD2,  ZRD20, ZRD21, ZRD22, ZRD23,
     9		ZRD24, ZRD25, ZRD26, ZRD27, ZRD28, ZRD29, ZRD3,  ZRD30,
     &		ZRD31, ZRD32, ZRD33, ZRD34, ZRD35, ZRD36, ZRD37, ZRD38,
     1		ZRD39, ZRD4,  ZRD40, ZRD41, ZRD42, ZRD43, ZRD44, ZRD45,
     2		ZRD46, ZRD47, ZRD48, ZRD5,  ZRD6,  ZRD7,  ZRD8,  ZRD9,
     3		RESA1
     	common /X_VARIABLES/
     1		ABX,   ABCX,  AIM1X, AIM2X, AIM3X, AMJX,  AWALLX,BTORX,
     2		ELONGX,ELONMX,ENCLX, ENWMX, FECRX, FFWX,  FICRX, FLHX,
     3		GN2EX, GN2IX, IPLX,  LEXTX, NNCLX, NNWMX, QECRX, QFWX,
     4		QICRX, QLHX,  QNBIX, RTORX, SHIFTX,TRIANX,TRICHX,UEXTX,
     5		UPDWNX,WNEX,  WTEX,  WTIX,  ZMJX,
     6		ZRD1X, ZRD10X,ZRD11X,ZRD12X,ZRD13X,ZRD14X,ZRD15X,ZRD16X,
     8		ZRD17X,ZRD18X,ZRD19X,ZRD2X, ZRD20X,ZRD21X,ZRD22X,ZRD23X,
     9		ZRD24X,ZRD25X,ZRD26X,ZRD27X,ZRD28X,ZRD29X,ZRD3X, ZRD30X,
     &		ZRD31X,ZRD32X,ZRD33X,ZRD34X,ZRD35X,ZRD36X,ZRD37X,ZRD38X,
     1		ZRD39X,ZRD4X, ZRD40X,ZRD41X,ZRD42X,ZRD43X,ZRD44X,ZRD45X,
     2		ZRD46X,ZRD47X,ZRD48X,ZRD5X, ZRD6X, ZRD7X, ZRD8X, ZRD9X,
     3		RESA1X
	double precision      DEVAR(NCONST),	 DEVARX(NCONST)
     	equivalence   	     (DEVAR(1),AB),	(DEVARX(1),ABX)
C Constants:
C CF1           [ ] - Model control constant (=1 default)
C CF2           [ ] - Model control constant (=1 default)
C CF3           [ ] - Model control constant (=1 default)
C CF4           [ ] - Model control constant (=1 default)
C CF5           [ ] - Model control constant (=1 default)
C CF6           [ ] - Model control constant (=1 default)
C CF7           [ ] - Model control constant (=1 default)
C CF8           [ ] - Model control constant (=1 default)
C CF9           [ ] - Model control constant (=1 default)
C CF10          [ ] - Model control constant (=1 default)
C CF11          [ ] - Model control constant (=1 default)
C CF12          [ ] - Model control constant (=1 default)
C CF13          [ ] - Model control constant (=1 default)
C CF14          [ ] - Model control constant (=1 default)
C CF15          [ ] - Model control constant (=1 default)
C CF16          [ ] - Model control constant (=1 default)
C CV1           [ ] - Model control constant (=0 default)
C CV2           [ ] - Model control constant (=0 default)
C CV3           [ ] - Model control constant (=0 default)
C CV4           [ ] - Model control constant (=0 default)
C CV5           [ ] - Model control constant (=0 default)
C CV6           [ ] - Model control constant (=0 default)
C CV7           [ ] - Model control constant (=0 default)
C CV8           [ ] - Model control constant (=0 default)
C CV9           [ ] - Model control constant (=0 default)
C CV10          [ ] - Model control constant (=0 default)
C CV11          [ ] - Model control constant (=0 default)
C CV12          [ ] - Model control constant (=0 default)
C CV13          [ ] - Model control constant (=0 default)
C CV14          [ ] - Model control constant (=0 default)
C CV15          [ ] - Model control constant (=0 default)
C CV16          [ ] - Model control constant (=0 default)
C CHE1          [ ]
C CHE2          [ ]
C CHE3          [ ]
C CHE4          [ ]
C CHI1          [ ]
C CHI2          [ ]
C CHI3          [ ]
C CHI4          [ ]
C CNB1          [ ] - CNB1=1,...,82    number of 'pencils'
C CNB2          [ ] - No. of internal grid points N1=(NA1-1)/CNB2+1
C CNB3          [ ]
C CNB4          [ ]
C CNBI1         [ ] - NBI FP solver control
C CNBI2         [ ] - not used
C CNBI3         [ ] - not used
C CNBI4         [ ] - 
C CCD1          [ ]
C CCD2          [ ]
C CCD3          [ ]
C CCD4          [ ]
C CRF1          [ ]
C CRF2          [ ]
C CRF3          [ ]
C CRF4          [ ]
C CNEUT1        [ ] - Neutrals control constant
C CNEUT2        [ ] - Neutrals control constant
C CNEUT3        [ ] - Neutrals control constant
C CNEUT4        [ ] - Neutrals control constant
C CPEL1         [ ] - Pellet control constant
C CPEL2         [ ] - Pellet control constant
C CPEL3         [ ] - Pellet control constant
C CPEL4         [ ] - Pellet control constant
C CBND1         [ ] - Boundary control constant
C CBND2         [ ] - Boundary control constant
C CBND3         [ ] - Boundary control constant
C CBND4         [ ] - Boundary control constant
C CFUS1         [ ] - Fusion control constant
C CFUS2         [ ] - Fusion control constant
C CFUS3         [ ] - Fusion control constant
C CFUS4         [ ] - Fusion control constant
C CIMP1         [ ] - Impurity control constant
C CIMP2         [ ] - Impurity control constant
C CIMP3         [ ] - Impurity control constant
C CIMP4         [ ] - Impurity control constant
C CMHD1         [ ] - MHD modes control constant
C CMHD2         [ ] - MHD modes control constant
C CMHD3         [ ] - MHD modes control constant
C CMHD4         [ ] - MHD modes control constant
C CRAD1         [ ] - Radiation control constant
C CRAD2         [ ] - Radiation control constant
C CRAD3         [ ] - Radiation control constant
C CRAD4         [ ] - Radiation control constant
C CSOL1         [ ] - Boundary control constant
C CSOL2         [ ] - Boundary control constant
C CSOL3         [ ] - Boundary control constant
C CSOL4         [ ] - Boundary control constant
C End constants	    - 88 entries (must be <=NCONST=256)
	double precision
     1		CF1,   CF2,   CF3,   CF4,   CF5,   CF6,   CF7,   CF8,
     2		CF9,   CF10,  CF11,  CF12,  CF13,  CF14,  CF15,  CF16,
     3		CV1,   CV2,   CV3,   CV4,   CV5,   CV6,   CV7,   CV8,
     4		CV9,   CV10,  CV11,  CV12,  CV13,  CV14,  CV15,  CV16,
     5		CHE1,  CHE2,  CHE3,  CHE4,  CHI1,  CHI2,  CHI3,  CHI4,
     6		CNB1,  CNB2,  CNB3,  CNB4,  CNBI1, CNBI2, CNBI3, CNBI4,
     7		CCD1,  CCD2,  CCD3,  CCD4,  CRF1,  CRF2,  CRF3,  CRF4,
     8		CNEUT1,CNEUT2,CNEUT3,CNEUT4,CPEL1, CPEL2, CPEL3, CPEL4,
     9		CBND1, CBND2, CBND3, CBND4, CFUS1, CFUS2, CFUS3, CFUS4,
     &		CIMP1, CIMP2, CIMP3, CIMP4, CMHD1, CMHD2, CMHD3, CMHD4,
     1		CRAD1, CRAD2, CRAD3, CRAD4, CSOL1, CSOL2, CSOL3, CSOL4,
     2		RESCON(NCONST-88)
C Old list: Versions 5.2 and earlier
	common /A_PARAMETERS/
C     1		CF1,   CF2,   CF3,   CF4,   CF5,   CF6,   CF7,   CF8,
C     2		CF9,   CF10,  CF11,  CF12,  CF13,  CF14,  CF15,  CF16,
C     3		CV1,   CV2,   CV3,   CV4,   CV5,   CV6,   CV7,   CV8,
C     4		CV9,   CV10,  CV11,  CV12,  CV13,  CV14,  CV15,  CV16,
C     5		CBND1, CBND2, CBND3, CBND4, CBM1,  CBM2,  CBM3,  CBM4,
C     6		CBMH1, CBMH2 ,CBMH3, CBMH4, CBMR1, CBMR2, CBMR3, CBMR4,
C     7		CBMS1, CBMS2 ,CBMS3, CBMS4, CBMI1, CBMI2, CBMI3, CBMI4,
C     8		CNEUT1,CNEUT2,CNEUT3,CNEUT4,CCD1,  CCD2,  CCD3,  CCD4,
C     9		CHE1,  CHE2,  CHE3,  CHE4,  CIMP1, CIMP2, CIMP3, CIMP4,
C     &		CRAD1, CRAD2, CRAD3, CRAD4, CMHD1, CMHD2, CMHD3, CMHD4,
C     1		CPEL1, CPEL2, CPEL3, CPEL4, CFUS1, CFUS2, CFUS3, CFUS4,
C     2		RESCON
C New list: Version 5.3
C     1		CA1,   CA2,   CA3,   CA4,   CA5,   CA6,   CA7,   CA8,
C     1		CB1,   CB2,   CB3,   CB4,   CB5,   CB6,   CB7,   CB8,
C     1		CG1,   CG2,   CG3,   CG4,   CG5,   CG6,   CG7,   CG8,
     1		CF1,   CF2,   CF3,   CF4,   CF5,   CF6,   CF7,   CF8,
     2		CF9,   CF10,  CF11,  CF12,  CF13,  CF14,  CF15,  CF16,
     3		CV1,   CV2,   CV3,   CV4,   CV5,   CV6,   CV7,   CV8,
     4		CV9,   CV10,  CV11,  CV12,  CV13,  CV14,  CV15,  CV16,
     5		CHE1,  CHE2,  CHE3,  CHE4,  CHI1,  CHI2,  CHI3,  CHI4,
     6		CNB1,  CNB2,  CNB3,  CNB4,  CNBI1, CNBI2, CNBI3, CNBI4,
     7		CCD1,  CCD2,  CCD3,  CCD4,  CRF1,  CRF2,  CRF3,  CRF4,
     8		CNEUT1,CNEUT2,CNEUT3,CNEUT4,CPEL1, CPEL2, CPEL3, CPEL4,
     9		CBND1, CBND2, CBND3, CBND4, CFUS1, CFUS2, CFUS3, CFUS4,
     &		CIMP1, CIMP2, CIMP3, CIMP4, CMHD1, CMHD2, CMHD3, CMHD4,
     1		CRAD1, CRAD2, CRAD3, CRAD4, CSOL1, CSOL2, CSOL3, CSOL4,
     2		RESCON
	double precision  CONSTF(NCONST)
     	equivalence	 (CONSTF(1),CF1)
C Internal program parameters
C DROUT       [s] - Radial profiles output interval
C DTOUT       [s] - Time variables output interval
C DPOUT       [s] - Profiles/time output interval
C TIME        [s] - Current time, t
C TAUMIN      [s] - Minimum time step
C TAUMAX      [s] - Maximum time step
C TAUINC      [ ] - Maximum time step increment (=TAUold/TAUnew)
C DELVAR      [ ] - Maximum relative change of variables at time step
C ITEREX      [ ] - No. of iterations in the external loop
C ITERIN      [ ] - Not used
C TINIT       [s] - Time axis initial time
C TSCALE      [s] - Time axis scale 
C NB1         [ ] - No. of points of the main transport grid
C NUF         [ ] - No. of radial points to be written in a U-file
C XOUT        [ ] - Output control: X-axis for radial curves (C Xaxis )
C XINPUT      [ ] - Default input mode for a radial variable (C Xdeflt)
C NB2EQL      [ ] - [0,1] contribution of beam pressure to equilibrium
C NEQUIL      [ ] - Eq. solver grid [<= none, 1 freeze eq., >= 42 - ESC)
C NBND        [ ] - No. of points on the boundary
C XFLAG       [ ] - X-point index
C DTEQL       [s] - Time interval to recalculate equilibrium
C MEQUIL      [ ] - Poloidal M (for ESC), grid size (for SPIDER)
C TPAUSE      [s] - Time to set "pause" mode
C TEND        [s] - End time, t_end, t_start <= t <= t_end
C TSTART      [s] - Start time, t_start <= t <= t_end
C TAU         [s] - Time step
C DTEQ        [s] - Subroutine J call interval  DTEQ(1,J)
C HRO         [m] - Radial grid size
C HROA        [m] - Radial grid edge size
C ALBPL       [ ] - Albedo of plasma (Neutral outflux)/(Neutral influx)
C NNCX        [ ] - No. of iterations in the neutral solver
C QETB   [MJ/keV] - Coefficient for Electron_heat_outflux = QETB*TEB
C QFF0B     [1/s] - Coefficient for F0_outflux = QFF0B*F0B
C QFF1B     [1/s] - Coefficient for F1_outflux = QFF1B*F1B
C QFF2B     [1/s] - Coefficient for F2_outflux = QFF2B*F2B
C QFF3B     [1/s] - Coefficient for F3_outflux = QFF3B*F3B
C QFF4B     [1/s] - Coefficient for F3_outflux = QFF4B*F4B
C QFF5B     [1/s] - Coefficient for F3_outflux = QFF5B*F5B
C QFF6B     [1/s] - Coefficient for F3_outflux = QFF6B*F6B
C QFF7B     [1/s] - Coefficient for F3_outflux = QFF7B*F7B
C QFF8B     [1/s] - Coefficient for F3_outflux = QFF8B*F8B
C QFF9B     [1/s] - Coefficient for F3_outflux = QFF9B*F9B
C QITB   [MJ/keV] - Coefficient for Ion_heat_outflux = QITB*TIB
C QNNB    [m^3/s] - Coefficient for Density_outflux = QNNB*NEB
C ROB         [m] - Limiter position (in "rho")
C ROWALL      [m] - Wall position    (in "rho")
C ROC         [m] - Limiter/separatrix radius (in the variable "rho")
C RON         [m] - Boundary position for density equation
C ROE         [m] - Boundary position for electron temperature equation
C ROI         [m] - Boundary position for ion temprature equation
C RO0         [m] - Boundary position for F0 equation
C RO1         [m] - Boundary position for F1 equation
C RO2         [m] - Boundary position for F2 equation
C RO3         [m] - Boundary position for F3 equation
C RO4         [m] - Boundary position for F4 equation
C RO5         [m] - Boundary position for F5 equation
C RO6         [m] - Boundary position for F6 equation
C RO7         [m] - Boundary position for F7 equation
C RO8         [m] - Boundary position for F8 equation
C RO9         [m] - Boundary position for F9 equation
C VOLUME    [m^3] - Total volume
C GP          [ ] = 3.1415926
C GP2         [ ] = 2*GP
C NA          [ ] = NA1-1
C NA1         [ ] - Limiter/separatrix position: ROC=RHO(NA1)
C NA1N        [ ] - Eq:NE boundary position:	 RON=RHO(NA1N)
C NA1E        [ ] - Eq:TE boundary position:	 ROE=RHO(NA1E)
C NA1I        [ ] - Eq:TI boundary position:	 ROI=RHO(NA1I)
C NA10        [ ] - Eq:F0 boundary position:	 RO1=RHO(NA11)
C NA11        [ ] - Eq:F1 boundary position:	 RO1=RHO(NA11)
C NA12        [ ] - Eq:F2 boundary position:	 RO2=RHO(NA12)
C NA13        [ ] - Eq:F3 boundary position:	 RO3=RHO(NA13)
C NA14        [ ] - Eq:F4 boundary position:	 RO4=RHO(NA14)
C NA15        [ ] - Eq:F5 boundary position:	 RO5=RHO(NA15)
C NA16        [ ] - Eq:F6 boundary position:	 RO6=RHO(NA16)
C NA17        [ ] - Eq:F7 boundary position:	 RO7=RHO(NA17)
C NA18        [ ] - Eq:F8 boundary position:	 RO8=RHO(NA18)
C NA19        [ ] - Eq:F9 boundary position:	 RO9=RHO(NA19)
C NAB         [ ] - AB=AMETR(NAB),		 ABC=AMETR(NA1)
C NSTEPS      [ ] - Total number of time steps
C QBEAM      [MW] - Total beam power (alias QNBI)
C CUT             - ASTRA function: CUT(x,y) = max(-x,min(x,y))
C FIXVAL          - ASTRA function: returns V(t<=t0) fixval(V,t0)
C FJUMP           - ASTRA function: unity step(t1)
C FRAMP           - ASTRA function: linear ramp-up(t1,t2)
C FRMIN           - ASTRA function: minimal value of array
C FRMAX           - ASTRA function: maximal value of array
C FTAV            - ASTRA function: Floating Time AVerage (see for/stdfun.f)
C FTAV2           - ASTRA function: 10^2 Time Step Average (see for/stdfun.f)
C FTAV3           - ASTRA function: 10^3 Time Step Average (see for/stdfun.f)
C FTBOX           - ASTRA function: FTBOX(t1,t2)=[1 if t1<t<t2][0 otherwise]
C FTMIN           - ASTRA function: FTMIN(Y)=min_t(y), where y=y(t[,r])
C FTMAX           - ASTRA function: FTMAX(Y)=max_t(y)
C AFVAL           - ASTRA function: a-position of array=value
C RFVAL           - ASTRA function: rho-position of array=value
C AFVEX           - ASTRA function: Outermost a-position of array=value
C RFVEX           - ASTRA function: Outermost rho-position of array=value
C AFVIN           - ASTRA function: Innermost a-position of array=value
C RFVIN           - ASTRA function: Innermost rho-position of array=value
C RFMIN           - ASTRA function: position of minimal array value
C RFMAX           - ASTRA function: position of maximal array value
C AFR             - ASTRA function: a=a(rho)
C AFX             - ASTRA function: a=a(x), x=\rho/ROC
C ATR             - ASTRA function: array at given \rho
C ATX             - ASTRA function: array at given x=\rho/ROC
C RFA             - ASTRA function: \rho=\rho(a)
C RFAN            - ASTRA function: \rho=\rho(a_N), a_N=a/ABC
C XFA             - ASTRA function: x=\rho(a)
C XFAN            - ASTRA function: x=\rho(a_N), a_N=a/ABC
C NODE            - ASTRA function: node(\rho)
C TIMAVG          - ASTRA function: time average
C TIMDER          - ASTRA function: time derivative
C TIMINT          - ASTRA function: time integral
C RECR            - ASTRA function: RECR(z,n) rho_tor of resonance surface
C STEP            - ASTRA Heaviside function = 0 if arg < 0, = 1 otherwise
C SIN             - FORTRAN intrinsic function
C COS             - FORTRAN intrinsic function
C TAN             - FORTRAN intrinsic function
C EXP             - FORTRAN intrinsic function
C ABS             - FORTRAN intrinsic function
C INT             - FORTRAN intrinsic function
C MIN             - FORTRAN intrinsic function
C MAX             - FORTRAN intrinsic function
C SQRT            - FORTRAN intrinsic function
C LOG             - FORTRAN intrinsic function
C LOG10           - FORTRAN intrinsic function
C ANINT           - FORTRAN intrinsic function: Nearest whole number
C TANH            - FORTRAN intrinsic function: Hyperbolic tangent
C REAL            - FORTRAN intrinsic function: real(N)
C DBLE            - FORTRAN intrinsic function: dble(N)
C SIGN            - FORTRAN intrinsic function: SIGN(a1,a2)=|a1|*signum(a2)
C End internal program parameters, 131 entries (must be <=NCONST=256)
	double precision
     1		DROUT,	DTOUT,	DPOUT,	TIME,  
     2		TAUMIN,	TAUMAX,	TAUINC,	DELVAR,
     3		ITEREX,	ITERIN,	TINIT,	TSCALE,
     4		NB1R,	NUFR,	XOUT,	XINPUT,
     5		NB2EQL,	NEQUIL,	NBNDR,	XFLAGR,
     6		DTEQL,  MEQUIL, TPAUSE,	TEND,   TEQ
	common /A_OUTPUT/
     1		DROUT,	DTOUT,	DPOUT,	TIME,
     2		TAUMIN,	TAUMAX,	TAUINC,	DELVAR,
     3		ITEREX,	ITERIN,	TINIT,	TSCALE,
     4		NB1R,	NUFR,	XOUT,	XINPUT,
     5		NB2EQL,	NEQUIL,	NBNDR,	XFLAGR,
     6		DTEQL,  MEQUIL, TPAUSE,	TEND,
     7		DTEQ(4,NSBMX),  TEQ(NSBMX)
	double precision	DELOUT(24+4*NSBMX),MESHEQ,	DTEQ
     	equivalence	(DELOUT(1),DROUT),     (MESHEQ,NEQUIL)
	double precision
     &			HRO,  HROA, VOLUME,ROB, ROC, ROWALL, GP,  GP2,
     1			ALBPL,TAU,  QETB, QITB, QNNB,QFF0B,QFF1B,QFF2B,
     2			QFF3B,QFF4B,QFF5B,QFF6B,QFF7B,QFF8B,QFF9B,
     3			RON,  ROE,  ROI,  RO0,  RO1, RO2,  RO3,
     4			RO4,  RO5,  RO6,  RO7,  RO8, RO9,
     5			TIMEQL,TSTART,QBEAM
     	equivalence	(QNBI,QBEAM)
	integer		NA,   NA1,  NAB,   NB1,  NUF,  NNCX, KEY, 
     1			NBND, XFLAG,IFESC, LEQ(NEQNS),
     2			NA1N, NA1E, NA1I,  NA10, NA11, NA12,NA13,
     3			NA14, NA15, NA16,  NA17, NA18, NA19,NSTEPS,
     4			AS_PID,ES_PID
	common /A_GRIDS/
     1			HRO,  HROA, VOLUME,ROB,  ROC,  GP,  GP2,
     2			TAU,  QETB, QITB,  QNNB, QFF0B,QFF1B,QFF2B,
     3			QFF3B,QFF4B,QFF5B, QFF6B,QFF7B,QFF8B,QFF9B,
     4			RON,  ROE,  ROI,   RO0,  RO1,  RO2,  RO3,
     5			RO4,  RO5,  RO6,   RO7,  RO8,  RO9,  ALBPL,
     6			ROWALL,     TIMEQL,TSTART,
     7			NA,   NA1,  NAB,   NB1,  NUF,  NNCX, KEY,
     8			NBND, XFLAG,IFESC,
     9			NA1N, NA1E, NA1I,  NA10, NA11, NA12, NA13,
     &			NA14, NA15, NA16,  NA17, NA18, NA19, NSTEPS,
     4			LEQ,  AS_PID,ES_PID
	integer		ARXUSE(NARRX)
	common /A_XUSAGE/ ARXUSE
	character*132  	XLINE1,XLINE2
	common /A_XCOMMENT/ XLINE1,XLINE2
C======================================================================|
