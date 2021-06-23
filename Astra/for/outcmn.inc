C NDTNAM - Number of names in the "Time step" table
C NCFNAM - Number of constant names
C NARNAM - Number of global arrays (described in the file "for/status.inc")
C NPRNAM - Number of global variables (between 'Variab' and 'End va')
C NSRNAM - Number of internal names (between 'Intern' and 'End in')
C NINVAR - Number of int. variables (between 'Intern' and 'ASTRA function:')
C NROUT  - Number of radial output profiles
C NTOUT  - Number of time output channels
C NSBR   - Number of subroutines
C LTOUT  - Current time output index
C IPOUT  - Current radial/time output index
C NSCR   - (0:9) Screen number in the current mode
C MOD10  - Output mode number
C NWIND1 - Windows for 8 boxes on screen
C NWIND2 - Windows for 2 boxes on screen
C NWIND3 - Windows for time output mode
C NWINDX - Windows for experimental data output
C SCALER - Radial profiles output scales
C SCALET - Time variables output scales
C ROUT   - Radial profiles output array
C TASK   - Task type
C RDNAME - Data file name
C EQNAME - Model file name
C RSNAME - PostView file name
C DTNAME - Time step name list
C CFNAME - Constant name list
C PRNAME - Variable name list
C SRNAME - Internal variable name list
C ARNAME - Array name list
C NAMER  - Radial profiles output names
C NAMET  - Time variables output names
C RUNID  - Run time, date, model and exp data
C NBND   - No. of boundary points (descr. in const.inc) 
C NBNT   - No. of time points for boundary
C BNDARR - array(25000) for boundary evolution (2*NBND+1)*NBNT <= 25000
C----------------------------------------------------------------------|
	integer
     1		NDTNAM,	NCFNAM,	NPRNAM,	NSRNAM,	NARNAM,
     2		NTOUT,  NROUT,	NSBR,	NSBP,	NBNT,
     3		LTOUT,	IPOUT,	MOD10,	NGR,	NXOUT,	
     4		NINVAR,	NWIND1,	NWIND2,	NWIND3,	NWINDX,
     5		IFDFVX,	IFDFAX,	NGRIDX,	NTYPEX, SIGNSB,
     6		KTO,	KOGDA,	GDEX,	GDEY,   IFSBX,
     7		NPTM,   NSCR,	MARKT,	MARKR
	common /A_EXCHAI/
     1		NDTNAM, NCFNAM, NPRNAM, NSRNAM,	NARNAM,
     2		NTOUT,  NROUT,  NSBR,   NSBP,	NBNT,	LTOUT,
     3		IPOUT,	MOD10,	NGR,	NXOUT,  NINVAR,	NSCR(9),
     4		SIGNSB(NSBMX),  IFSBX(NSBMX),	IFDFVX(NCONST),
     5		NWIND1(NRW),  NWIND2(NRW),  NWIND3(NRW), NWINDX(NRW),
     6		MARKT(NRW),   MARKR(NRW),   KTO(NTARR),  NGRIDX(NTARR),
     7		NTYPEX(NTARR),IFDFAX(NARRX),KOGDA(NARRX),NPTM(NARRX)
	double precision
     1		SCALET,	SCALER,	OSHIFT,	OSHIFR,
     2		ROUT,	TIMEX,	TOUTX,	XAXES,	
     3		DATAX,  FILTER,	BNDARR,	GRAL,	GRAP
	real*4	DATARR
	common /A_EXCHAR/
     1		SCALET(NRW), SCALER(NRW),ROUT(NRD,NRW),   OSHIFR(NRW),
     2		GRAL(NRW),   GRAP(NRW),  OSHIFT(NRW),
     3		GDEX(NTARR), GDEY(NTARR),TIMEX(NTARR),   FILTER(NTARR),
     4		TOUTX(NARRX),DATARR(NRDX*NTARR),  XAXES(NRDX,NARRX),
     5		DATAX(NRDX,NARRX), BNDARR(25000)
	character*6
     1		DTNAME,	CFNAME,	PRNAME,	SRNAME,	ARNAME,
     2		NAM7,	EXARNM,	NAMEX
	character*64    LISTSB(NSBMX)
	character*4     NAMET,  NAMER,	TASK
     	common /A_CNTRNM/
     1		DTNAME(24+4*NSBMX),CFNAME(NCONST),
     2		PRNAME(NCONST), SRNAME(NCONST), ARNAME(512), 
     3		NAMER(NRW),     NAMET(NRW),     NAM7(4),
     4		NAMEX(NRW),     EXARNM(NARRX)
	double precision	TIM7, SCM
	integer
     1		NWIND4,	NWIND7,	MODEY,
     2		XSCMAX,	YSCMAX,	DXLET,	DYLET,
     3		IY0,	IYM,	DYM,	DXM,
     4		IDX,	IDT,	KPRI,	NST,
     5		XWX,	XWY,	XWW,	XWH,
     6		IYEAR,	IMONTH,	IDAY,	IHOUR,	IMINUT,	ISEC,
     7		AVERS,  ARLEAS, AEDIT,	NBFLAG
	integer*4 COLTAB(64)	! i*4 needed for compatibility with C
	common /A_MOD4AS/
     1		NWIND4(NRW), NWIND7(NRW),  TIM7(4),	MODEY
	common /A_XWINDOW/
     1		XSCMAX,      YSCMAX,      DXLET,	DYLET,
     2		IY0,         IYM,         DYM,		DXM,
     3		SCM,	     IDX,	  IDT,	  	COLTAB,
     4		XWX,	     XWY,	  XWW,	  	XWH,
     5		KPRI,	     NST,	  NBFLAG
	double precision     CPTOT,	  CPTEQL,	CPTGRA,	
     1			     CPT,	  CPTSBR(NSBMX)
	common /A_STATISTICS/
     1		CPTOT,	     CPTEQL,	  CPTGRA,	CPT,
     2		CPTSBR
	character*80	     VERSION*32,  RUNID,	CODEID,
     1			     AROOT,	  AWD,		WHOME,
     2			     FILEX
	character*40         RDNAME,      EQNAME,       RSNAME,
     1			     TASKID,	  CNFILE,	NBFILE,
     2		     	     ECFILE,	  ICFILE,	MSFILE
	common	/A_RUNID/    AVERS,	  ARLEAS,	AEDIT,
     1		IYEAR,	IMONTH,	IDAY,	IHOUR,	IMINUT,	ISEC,
     2		TASKID,	     VERSION,	  RUNID,	CODEID,
     3		TASK,        RDNAME,      EQNAME,       RSNAME,
     4		CNFILE,	     NBFILE,	  ECFILE,	ICFILE,
     5		MSFILE,	     AROOT,	  AWD,		WHOME,
     6		FILEX
C======================================================================|
