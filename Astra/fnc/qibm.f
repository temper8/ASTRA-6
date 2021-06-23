C QIBM [MW]:	 Integral {0,R} (PIBM) dV
C			(Pereverzev 20-MAY-08)
	double precision FUNCTION QIBMR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QIBMR=VINT(PIBM,YR)
	END
