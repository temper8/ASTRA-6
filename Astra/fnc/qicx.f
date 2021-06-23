C QICX [MW]:	 Integral {0,R} ( PICX ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QICXR(YR)
	implicit none
	double precision PICX,SVCX
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QICXR=0.
	DO 1 J=1,JK
	include  'fml/picx'
 1	QICXR=QICXR+PICX*VR(J)
	QICXR=HRO*(QICXR-PICX*YDR)
	END
