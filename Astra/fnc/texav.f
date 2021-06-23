C TEXAV [keV]:	Volume average electron temperature (r)
C	Integral {0:R}	(TEX/VOL) dV
C			(Pereverzev 13-MARCH-91)
	double precision FUNCTION TEXAVR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	TEXAVR=TEX(1)
	DO 1 J=1,JK
	YV=YV+VR(J)
 1	YQ=YQ+TEX(J)*VR(J)
 	YV=YV-YDR
	YQ=YQ-TEX(JK)*YDR
	TEXAVR=YQ/YV
	END
