C QDV5 [F5/s]:	 Integral {0,R} ( SDV5 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV5R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV5R=0.
	do 1 J=1,JK
 1	QDV5R=QDV5R+SDV5(J)*VR(J)
	QDV5R=HRO*(QDV5R-SDV5(J)*YDR)
	end

