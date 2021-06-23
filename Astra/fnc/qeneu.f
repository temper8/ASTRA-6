C QENEU [MW]:	 Integral {0,R} ( PENEU ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QENEUR(YR)
	implicit none
	double precision PENEU,SVIE
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QENEUR=0.
	DO 1 J=1,JK
	include  'fml/peneu'
 1	QENEUR=QENEUR+PENEU*VR(J)
	QENEUR=HRO*(QENEUR-PENEU*YDR)
	END
