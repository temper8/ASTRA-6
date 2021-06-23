C TENDN [keV]:	Volume-density average electron temperature (r)
C	<Te*Ne*dV>/<Ne*dV> 
C				(Polevoy 28.09.89)
	double precision function	TENDNR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	TENDNR=TE(1)
	DO 1 J=1,JK
	YV=YV+NE(J)*VR(J)
 1	YQ=YQ+TE(J)*NE(J)*VR(J)
 	YV=YV-YDR*NE(JK)
	YQ=YQ-TE(JK)*YDR*NE(JK)
	TENDNR=YQ/YV
	END
