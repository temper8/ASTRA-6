C QDV0 [F0/s]:	 Integral {0,R} ( SDV0 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV0R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV0R=0.
	do 1 J=1,JK
 1	QDV0R=QDV0R+SDV0(J)*VR(J)
	QDV0R=HRO*(QDV0R-SDV0(J)*YDR)
	end

