C QDV4 [F4/s]:	 Integral {0,R} ( SDV4 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV4R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV4R=0.
	do 1 J=1,JK
 1	QDV4R=QDV4R+SDV4(J)*VR(J)
	QDV4R=HRO*(QDV4R-SDV4(J)*YDR)
	end

