C QIX [MW]:	 Integral {0,R} (PIX) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QIXR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QIXR=VINT(PIX,YR)
	END
