C TIAV [keV]:	Volume average ion temperature (r)
C	Integral {0:R}	(TE/VOL) dV
C			(Yushmanov 11-MAY-87)
	double precision FUNCTION TIAVR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	TIAVR=TI(1)
	DO 1 J=1,JK
	YV=YV+VR(J)
 1	YQ=YQ+TI(J)*VR(J)
 	YV=YV-YDR
	YQ=YQ-TI(JK)*YDR
	TIAVR=YQ/YV
	end
