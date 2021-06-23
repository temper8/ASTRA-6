C QEBM [MW]:	 Integral {0,R} (PEBM) dV
C			(Pereverzev 20-MAY-08)
	double precision FUNCTION QEBMR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QEBMR=VINT(PEBM,YR)
	END
