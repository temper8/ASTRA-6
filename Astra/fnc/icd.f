C ICD [MA]:	Toroidal driven current inside  {O,R}
C			(Yushmanov 14-DEC-90)
	double precision	FUNCTION	ICDR(YR)
	implicit none
	double precision 	IINT,YR
	external IINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	ICDR=IINT(CD,YR)
	END
