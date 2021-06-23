C QBRAD [MW]:	 Integral {0,R} ( PBRAD ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QBRADR(YR)
	implicit none
	double precision PBRAD
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QBRADR=0.
	DO 1 J=1,JK
	include  'fml/pbrad'
 1	QBRADR=QBRADR+PBRAD*VR(J)
	QBRADR=HRO*(QBRADR-PBRAD*YDR)
	end
