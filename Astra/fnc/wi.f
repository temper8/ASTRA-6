C WI [MJ]: Integral {0:R} ( 3/2*NI*TI ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION WIR(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	WIR=0.
	DO 1 J=1,JK
 1	WIR=WIR+NI(J)*TI(J)*VR(J)
	WIR=HRO*(WIR-NI(JK)*TI(JK)*YDR)*.0024
	END
