C QIDT [MW]:	 Integral {0,R} ( PIDT ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QIDTR(YR)
	implicit double precision (y)
	double precision PIDT,PDT
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QIDTR=0.
	DO 1 J=1,JK
	include  'fml/pidt'
 1	QIDTR=QIDTR+PIDT*VR(J)
	QIDTR=HRO*(QIDTR-PIDT*YDR)
	END
