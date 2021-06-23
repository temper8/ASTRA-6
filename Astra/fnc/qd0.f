C QD0 [F0/s]:	 Integral {0,R} ( SD0 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD0R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD0R=0.
	do 1 J=1,JK
 1	QD0R=QD0R+SD0(J)*VR(J)
	QD0R=HRO*(QD0R-SD0(J)*YDR)
	end

