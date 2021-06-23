C QD2 [F2/s]:	 Integral {0,R} ( SD2 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD2R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD2R=0.
	do 1 J=1,JK
 1	QD2R=QD2R+SD2(J)*VR(J)
	QD2R=HRO*(QD2R-SD2(J)*YDR)
	end

