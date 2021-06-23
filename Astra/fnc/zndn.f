C ZNDN []:	Volume-density averaged Zef (r)
c   <Zef*Ne*dV>/<Ne*dV>
C				 (Polevoy 28.09.89)
	double precision function ZNDNR(YR)
	implicit none
	double precision YQ,YV
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YV=0.
	ZNDNR=ZEF(1)
	DO 1 J=1,JK
	YV=YV+NE(J)*VR(J)
 1	YQ=YQ+ZEF(J)*NE(J)*VR(J)
 	YV=YV-YDR*NE(JK)
	YQ=YQ-ZEF(JK)*YDR*NE(JK)
	ZNDNR=YQ/YV
	END
