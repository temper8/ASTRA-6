C TAUF9 [s]:	Confinement_time@(rho) for the quantity F9
C       TAUF9=Vint(F9)/(Vint(SF9TOT)-dVint(F9)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF9R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF9R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF9TOT(J)*YJ
	   YW = YW+F9(J)*YJ
	   YWO = YWO+F9O(J)*VRO(J)
	enddo
	YQ = YQ-SF9TOT(JK)*YDR
	YW = YW-F9(JK)*YDR
	YWO = YWO-F9O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF9R = 0.
	else
	   TAUF9R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
