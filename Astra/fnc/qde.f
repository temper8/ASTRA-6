C QDE [MW]:	 Integral {0,R} ( PDE ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDER(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDER=0.
	do 1 J=1,JK
 1	QDER=QDER+PDE(J)*VR(J)
	QDER=HRO*(QDER-PDE(J)*YDR)
	end
