C TAUF6 [s]:	Confinement_time@(rho) for the quantity F6
C       TAUF6=Vint(F6)/(Vint(SF6TOT)-dVint(F6)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF6R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF6R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF6TOT(J)*YJ
	   YW = YW+F6(J)*YJ
	   YWO = YWO+F6O(J)*VRO(J)
	enddo
	YQ = YQ-SF6TOT(JK)*YDR
	YW = YW-F6(JK)*YDR
	YWO = YWO-F6O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF6R = 0.
	else
	   TAUF6R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
