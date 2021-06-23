C NEXAV [10#19/m#3]:	Volume average density (r)
C	Integral {0:R}	(NEX/VOL) dV
C			(Yushmanov 11-MAY-87)
	double precision FUNCTION NEXAVR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	NEXAVR=NEX(1)
	DO 1 J=1,JK
	YV=YV+VR(J)
 1	YQ=YQ+NEX(J)*VR(J)
 	YV=YV-YDR
	YQ=YQ-NEX(JK)*YDR
	NEXAVR=YQ/YV
	END
