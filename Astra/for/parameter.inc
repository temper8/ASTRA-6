Constants controlling reading data from exp/... and udb/...
C  NCONST   max No. of the simple variables & constants 
C		(both lists are in the file "const.inc")
C	    which can be initialized from a data (exp/...) file
C	    and are recognizable by the Astra Compiler
C		Increase of NCONST beyond 256 requires
C		changes in the Astra Compiler "nambuf.inc"
C  NARRX    No. of arrays readable from a data file (common block A_CARX)
C  NRD	    max No. of radial points in the transport grid
C  NRDX	    max No. of radial points in exp data
C  NTIMES   max No. of stored time slices for output in mode 6
C  NRW	    max No. of radial and time curves for output
	integer NRD,NRW,NCONST,NARRX,NRDX,NTVAR,NTARR,NEQNS,NSBMX
	parameter(NRD=501, NRW=128, NCONST=256, NARRX=67, NSBMX=20)
C  NTVAR total No. of the time slices for all variables in a data file
C  NTARR total No. of the time slices for all arrays in a data file
C  NEQNS max. total No. of eqns in a model, is used in .srv/tmpbuf.inc
	parameter(NRDX=200, NTVAR=25000, NTARR=25000, NEQNS=19)
