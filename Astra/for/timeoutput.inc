C  NTIMES   max No. of stored time slices for output in mode 6
	integer NTIMES
	parameter(NTIMES=1024)
C TOUT   - Time variables output array
C TTOUT  - time-coordinate array for time output [s] TTOUT(1:LTOUT<=NTIMES)
	double precision TTOUT(NTIMES),TOUT(NTIMES,NRW)
	common /A_TIMEOUTPUT/ TTOUT,TOUT
