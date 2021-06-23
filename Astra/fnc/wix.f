C WIX [MJ]: Integral {0:R} ( 3/2*NI*TIX ) dV
C			(Pereverzev 25-MAR-92)
	double precision FUNCTION WIXR(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	WIXR=0.
	DO 1 J=1,JK
 1	WIXR=WIXR+NI(J)*TIX(J)*VR(J)
	WIXR=HRO*(WIXR-NI(JK)*TIX(JK)*YDR)*.0024
	END
