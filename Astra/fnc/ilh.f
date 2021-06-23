C ILHR [MA]:	Toroidal driven current inside  {O,R}
C			(Pereverzev 29-NOV-96)
	double precision	FUNCTION	ILHR(YR)
	implicit none
	double precision 	IINT,YR
	external IINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	ILHR=IINT(CULH,YR)
	END
