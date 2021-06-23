C QD6 [F6/s]:	 Integral {0,R} ( SD6 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD6R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD6R=0.
	do 1 J=1,JK
 1	QD6R=QD6R+SD6(J)*VR(J)
	QD6R=HRO*(QD6R-SD6(J)*YDR)
	end

