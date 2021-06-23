C TIXAV [keV]:	Volume average ion temperature (r)
C	Integral {0:R}	(TE/VOL) dV
C			(Perverzev 13-MARCH-91)
	double precision FUNCTION TIXAVR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	TIXAVR=TIX(1)
	DO 1 J=1,JK
	YV=YV+VR(J)
 1	YQ=YQ+TIX(J)*VR(J)
 	YV=YV-YDR
	YQ=YQ-TIX(JK)*YDR
	TIXAVR=YQ/YV
	END
