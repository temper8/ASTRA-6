C WEX [MJ]:	Integral {0:R} ( 3/2*NEX*TEX ) dV
C			(Pereverzev 25-MAR-92)
	double precision FUNCTION WEXR(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	WEXR=0.
	DO 1 J=1,JK
 1	WEXR=WEXR+NEX(J)*TEX(J)*VR(J)
	WEXR=HRO*(WEXR-NEX(JK)*TEX(JK)*YDR)*.0024
	END
