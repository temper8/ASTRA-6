C QDV6 [F6/s]:	 Integral {0,R} ( SDV6 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV6R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV6R=0.
	do 1 J=1,JK
 1	QDV6R=QDV6R+SDV6(J)*VR(J)
	QDV6R=HRO*(QDV6R-SDV6(J)*YDR)
	end

