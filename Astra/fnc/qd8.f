C QD8 [F8/s]:	 Integral {0,R} ( SD8 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD8R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD8R=0.
	do 1 J=1,JK
 1	QD8R=QD8R+SD8(J)*VR(J)
	QD8R=HRO*(QD8R-SD8(J)*YDR)
	end

