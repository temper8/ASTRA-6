C IBM [MA]:	Toroidal beam driven current CUBM {O,R} 
C			(Yushmanov 14-DEC-90)
	double precision	FUNCTION	IBMR(YR)
	implicit none
	double precision	IINT,YR
	external IINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	IBMR=IINT(CUBM,YR)
	END
