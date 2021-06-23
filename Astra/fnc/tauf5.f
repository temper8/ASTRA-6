C TAUF5 [s]:	Confinement_time@(rho) for the quantity F5
C       TAUF5=Vint(F5)/(Vint(SF5TOT)-dVint(F5)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF5R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF5R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF5TOT(J)*YJ
	   YW = YW+F5(J)*YJ
	   YWO = YWO+F5O(J)*VRO(J)
	enddo
	YQ = YQ-SF5TOT(JK)*YDR
	YW = YW-F5(JK)*YDR
	YWO = YWO-F5O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF5R = 0.
	else
	   TAUF5R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
