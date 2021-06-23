C TAUF4 [s]:	Confinement_time@(rho) for the quantity F4
C       TAUF4=Vint(F4)/(Vint(SF4TOT)-dVint(F4)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF4R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF4R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF4TOT(J)*YJ
	   YW = YW+F4(J)*YJ
	   YWO = YWO+F4O(J)*VRO(J)
	enddo
	YQ = YQ-SF4TOT(JK)*YDR
	YW = YW-F4(JK)*YDR
	YWO = YWO-F4O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF4R = 0.
	else
	   TAUF4R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
