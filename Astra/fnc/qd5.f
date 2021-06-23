C QD5 [F5/s]:	 Integral {0,R} ( SD5 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD5R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD5R=0.
	do 1 J=1,JK
 1	QD5R=QD5R+SD5(J)*VR(J)
	QD5R=HRO*(QD5R-SD5(J)*YDR)
	end

