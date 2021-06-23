C QETOT [MW]:	 Integral {0,R} (PETOT) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QETOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QETOTR=VINT(PETOT,YR)
	END
