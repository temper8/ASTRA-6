C IOHM [MA]:	Toroidal ohmic current inside  {O,R}
C			(Pereverzev 17-FEB-00)
	double precision	function	IOHMR(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision	 YR,A(NRD),IINT
	integer  j
	external IINT
	do	j=1,NA1
	   A(j)=CC(j)*ULON(j)/(RTOR*GP2)
	   if (YR .lt. RHO(j)+HRO)	goto	1
	enddo
 1	continue
	IOHMR=IINT(A,YR)
	END
