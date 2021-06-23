C QD4 [F4/s]:	 Integral {0,R} ( SD4 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD4R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD4R=0.
	do 1 J=1,JK
 1	QD4R=QD4R+SD4(J)*VR(J)
	QD4R=HRO*(QD4R-SD4(J)*YDR)
	end

