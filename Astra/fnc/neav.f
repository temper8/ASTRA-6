C NEAV [10#19/m#3]:	Volume average density (r)
C	Integral {0:R}	(NE) dV /VOLUME(r)
C			(Yushmanov 11-MAY-87)
	double precision FUNCTION NEAVR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	NEAVR=NE(1)
	do  J=1,JK
	YV=YV+VR(J)
	YQ=YQ+NE(J)*VR(J)
	enddo
 	YV=YV-YDR
	YQ=YQ-NE(JK)*YDR
	NEAVR=YQ/YV
	end
