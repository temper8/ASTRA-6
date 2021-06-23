C QEDWT [MW]: dWe/dt
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QEDWTR(YR)
	implicit none
	double precision YQ,YQO
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YQO=0.
	QEDWTR=0.
	DO 1 J=1,JK
	YQO=YQO+NEO(J)*TEO(J)*VRO(J)
 1	YQ=YQ+NE(J)*TE(J)*VR(J)
	YQO=YQO-NEO(JK)*TEO(JK)*YDR
	YQ=YQ-NE(JK)*TE(JK)*YDR
	QEDWTR=(YQ-YQO)/TAU*HRO*.0024
	END
