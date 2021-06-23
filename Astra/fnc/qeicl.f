C QEICL [MW]:	 Integral {0,R} ( PEICL ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QEICLR(YR)
	implicit none
	double precision PEICL,COULG
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	include 'for/yrjkdr.inc'
	QEICLR=0.
	DO 1 J=1,JK
	include  'fml/peicl'
 1	QEICLR=QEICLR+PEICL*VR(J)
	QEICLR=HRO*(QEICLR-PEICL*YDR)
	END
