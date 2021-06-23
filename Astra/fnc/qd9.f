C QD9 [F9/s]:	 Integral {0,R} ( SD9 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD9R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD9R=0.
	do 1 J=1,JK
 1	QD9R=QD9R+SD9(J)*VR(J)
	QD9R=HRO*(QD9R-SD9(J)*YDR)
	end

