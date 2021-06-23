C QDV2 [F2/s]:	 Integral {0,R} ( SDV2 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV2R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV2R=0.
	do 1 J=1,JK
 1	QDV2R=QDV2R+SDV2(J)*VR(J)
	QDV2R=HRO*(QDV2R-SDV2(J)*YDR)
	end

