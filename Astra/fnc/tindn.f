C TINDN [keV]:	Volume-density averaged electron temperature (r)
C	<Ti*Ne*dV>/<Ne*dV>  
C				(Polevoy 28.09.89)
	double precision function	TINDNR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/status.inc'
	include 'for/const.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	TINDNR=TI(1)
	DO 1 J=1,JK
	YV=YV+NE(J)*VR(J)
 1	YQ=YQ+TI(J)*NE(J)*VR(J)
 	YV=YV-YDR*NE(JK)
	YQ=YQ-TI(JK)*YDR*NE(JK)
	TINDNR=YQ/YV
	END
