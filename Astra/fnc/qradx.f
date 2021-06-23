C QRADX [MW]:	 Integral {0,R} (PRADX) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QRADXR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QRADXR=VINT(PRADX,YR)
	END
