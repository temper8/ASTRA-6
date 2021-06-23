C TAUF0 [s]:	Confinement_time@(rho) for the quantity F0
C       TAUF0=Vint(F0)/(Vint(SF0TOT)-dVint(F0)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF0R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF0R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF0TOT(J)*YJ
	   YW = YW+F0(J)*YJ
	   YWO = YWO+F0O(J)*VRO(J)
	enddo
	YQ = YQ-SF0TOT(JK)*YDR
	YW = YW-F0(JK)*YDR
	YWO = YWO-F0O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF0R = 0.
	else
	   TAUF0R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
