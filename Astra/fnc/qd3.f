C QD3 [F3/s]:	 Integral {0,R} ( SD3 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD3R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD3R=0.
	do 1 J=1,JK
 1	QD3R=QD3R+SD3(J)*VR(J)
	QD3R=HRO*(QD3R-SD3(J)*YDR)
	end

