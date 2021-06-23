C QEX [MW]:	 Integral {0,R} (PEX) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QEXR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QEXR=VINT(PEX,YR)
	END
