C QDV1 [F1/s]:	 Integral {0,R} ( SDV1 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV1R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV1R=0.
	do 1 J=1,JK
 1	QDV1R=QDV1R+SDV1(J)*VR(J)
	QDV1R=HRO*(QDV1R-SDV1(J)*YDR)
	end

