C TAUF7 [s]:	Confinement_time@(rho) for the quantity F7
C       TAUF7=Vint(F7)/(Vint(SF7TOT)-dVint(F7)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF7R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF7R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF7TOT(J)*YJ
	   YW = YW+F7(J)*YJ
	   YWO = YWO+F7O(J)*VRO(J)
	enddo
	YQ = YQ-SF7TOT(JK)*YDR
	YW = YW-F7(JK)*YDR
	YWO = YWO-F7O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF7R = 0.
	else
	   TAUF7R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
