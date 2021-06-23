C IXR [MA]:	Toroidal driven current inside  {O,R}
C			(Pereverzev 29-NOV-96)
	double precision	FUNCTION	IXR(YR)
	implicit none
	double precision 	IINT,YR
	external IINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	IXR=IINT(CUX,YR)
	END
