C QDV7 [F7/s]:	 Integral {0,R} ( SDV7 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDV7R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDV7R=0.
	do 1 J=1,JK
 1	QDV7R=QDV7R+SDV7(J)*VR(J)
	QDV7R=HRO*(QDV7R-SDV7(J)*YDR)
	end

