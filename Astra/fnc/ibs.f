C IBS [MA]:	Toroidal  bootstrep current inside {0,R} 
C			(Yushmanov 14-DEC-90)
	double precision	FUNCTION	IBSR(YR)
	implicit none
	double precision IINT,YR
	external IINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	IBSR=IINT(CUBS,YR)
	END
