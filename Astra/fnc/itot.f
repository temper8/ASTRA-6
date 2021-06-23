C ITOT [MA]:	 Total current inside {O,R}
C			(Yushmanov 14-DEC-90)
	double precision	FUNCTION	ITOTR(YR)
	implicit none
	double precision	IINT,YR
	external IINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	ITOTR=IINT(CU,YR)
	END
