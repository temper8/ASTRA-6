C TAUF3 [s]:	Confinement_time@(rho) for the quantity F3
C       TAUF3=Vint(F3)/(Vint(SF3TOT)-dVint(F3)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF3R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF3R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF3TOT(J)*YJ
	   YW = YW+F3(J)*YJ
	   YWO = YWO+F3O(J)*VRO(J)
	enddo
	YQ = YQ-SF3TOT(JK)*YDR
	YW = YW-F3(JK)*YDR
	YWO = YWO-F3O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF3R = 0.
	else
	   TAUF3R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
