C QNDNT [10#19/s]: d/dt(Volume integral {0,R}  NE )
C			(Yushmanov 15-FEB-89)
	double precision FUNCTION QNDNTR(YR)
	implicit none
	double precision YQ,YQO
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YQO=0.
	QNDNTR=0.
	DO 1 J=1,JK
	YQO=YQO+NEO(J)*VRO(J)
 1	YQ=YQ+NE(J)*VR(J)
	YQO=YQO-NEO(JK)*YDR
	YQ=YQ-NE(JK)*YDR
	QNDNTR=(YQ-YQO)/TAU*HRO
	END
