C QIDTF [MW]:	 Integral {0,R} ( PIDTF ) dV
C			(Pereverzev 7-APR-03)
	double precision FUNCTION QIDTFR(YR)
	implicit double precision (a-h,o-z)
	double precision PIDTF
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QIDTFR=0.
	DO 1 J=1,JK
	include  'fml/pidtf'
 1	QIDTFR=QIDTFR+PIDTF*VR(J)
	QIDTFR=HRO*(QIDTFR-PIDTF*YDR)
	END
