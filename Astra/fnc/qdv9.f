C QDV9 [F9/s]:	 Integral {0,R} ( SDV9 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV9R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV9R=0.
	do 1 J=1,JK
 1	QDV9R=QDV9R+SDV9(J)*VR(J)
	QDV9R=HRO*(QDV9R-SDV9(J)*YDR)
	end

