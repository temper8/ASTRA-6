C QIDWT [MW]: dWi/dt
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QIDWTR(YR)
	implicit none
	double precision YQ,YQO
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YQO=0.
	QIDWTR=0.
	DO 1 J=1,JK
	YQO=YQO+NIO(J)*TIO(J)*VRO(J)
 1	YQ=YQ+NI(J)*TI(J)*VR(J)
	YQO=YQO-NIO(JK)*TIO(JK)*YDR
	YQ=YQ-NI(JK)*TI(JK)*YDR
	QIDWTR=(YQ-YQO)/TAU*HRO*.0024
	END
