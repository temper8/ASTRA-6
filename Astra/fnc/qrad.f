C QRAD [MW]:	 Integral {0,R} (PRAD) dV
C			(Yushmanov 11-JAN-89)
	double precision function QRADR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QRADR=VINT(PRAD,YR)
	END
