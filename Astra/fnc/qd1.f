C QD1 [F1/s]:	 Integral {0,R} ( SD1 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD1R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD1R=0.
	do 1 J=1,JK
 1	QD1R=QD1R+SD1(J)*VR(J)
	QD1R=HRO*(QD1R-SD1(J)*YDR)
	end

