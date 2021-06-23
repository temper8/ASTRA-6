C QEDTF [MW]:	 Integral {0,R} ( PEDTF ) dV
C			(Pereverzev 7-APR-03)
	double precision FUNCTION QEDTFR(YR)
	implicit double precision (y)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QEDTFR=0.
	DO 1 J=1,JK
	include  'fml/pedtf'
 1	QEDTFR=QEDTFR+PEDTF*VR(J)
	QEDTFR=HRO*(QEDTFR-PEDTF*YDR)
	END
