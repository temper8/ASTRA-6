C QITOT [MW]:	 Integral {0,R} (PITOT) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QITOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QITOTR=VINT(PITOT,YR)
	END
