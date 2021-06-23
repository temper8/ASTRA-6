C TAUF2 [s]:	Confinement_time@(rho) for the quantity F2
C       TAUF2=Vint(F2)/(Vint(SF2TOT)-dVint(F2)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF2R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF2R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF2TOT(J)*YJ
	   YW = YW+F2(J)*YJ
	   YWO = YWO+F2O(J)*VRO(J)
	enddo
	YQ = YQ-SF2TOT(JK)*YDR
	YW = YW-F2(JK)*YDR
	YWO = YWO-F2O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF2R = 0.
	else
	   TAUF2R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
