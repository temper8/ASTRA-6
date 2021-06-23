C QNX [10#19/s]:	 Integral {0,R} (SNX) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QNXR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	QNXR=VINT(SNX,YR)
	END
