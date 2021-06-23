C QDV3 [F3/s]:	 Integral {0,R} ( SDV3 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV3R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV3R=0.
	do 1 J=1,JK
 1	QDV3R=QDV3R+SDV3(J)*VR(J)
	QDV3R=HRO*(QDV3R-SDV3(J)*YDR)
	end

