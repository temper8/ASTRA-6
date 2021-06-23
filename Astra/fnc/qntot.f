C QNTOT [10#19/s]:	 Integral {0,R} (SNTOT) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QNTOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QNTOTR=VINT(SNTOT,YR)
	END
