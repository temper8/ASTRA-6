C QD7 [F7/s]:	 Integral {0,R} ( SD7 ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QD7R(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QD7R=0.
	do 1 J=1,JK
 1	QD7R=QD7R+SD7(J)*VR(J)
	QD7R=HRO*(QD7R-SD7(J)*YDR)
	end

