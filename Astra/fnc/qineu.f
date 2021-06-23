C QINEU [MW]:	 Integral {0,R} ( PINEU ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QINEUR(YR)
	implicit none
	double precision PINEU,SVIE,SVII,SVCX,PICX
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QINEUR=0.
	DO 1 J=1,JK
	include  'fml/pineu'
 1	QINEUR=QINEUR+PINEU*VR(J)
	QINEUR=HRO*(QINEUR-PINEU*YDR)
	END
