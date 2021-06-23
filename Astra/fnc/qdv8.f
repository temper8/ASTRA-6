C QDV8 [F8/s]:	 Integral {0,R} ( SDV8 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV8R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV8R=0.
	do 1 J=1,JK
 1	QDV8R=QDV8R+SDV8(J)*VR(J)
	QDV8R=HRO*(QDV8R-SDV8(J)*YDR)
	end

