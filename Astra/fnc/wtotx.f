C WTOTX [MJ]:	Integral {0:R} ( 3/2*NEX*(TEX+TIX) ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION WTOTXR(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	WTOTXR=0.
	DO 1 J=1,JK
 1	WTOTXR=WTOTXR+(NEX(J)*TEX(J)+NI(J)*TIX(J))*VR(J)
	WTOTXR=HRO*(WTOTXR-(NEX(JK)*TEX(JK)+NI(JK)*TIX(JK))*YDR)*.0024
	END
