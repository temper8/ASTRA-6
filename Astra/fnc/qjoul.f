C QJOUL [MW]:	 Integral {0,R} (PJOUL) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QJOULR(YR)
	implicit none
	double precision PJOUL
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QJOULR=0.
	DO 1 J=1,JK
	include  'fml/pjoul'
 1	QJOULR=QJOULR+PJOUL*VR(J)
	QJOULR=HRO*(QJOULR-PJOUL*YDR)
	END
