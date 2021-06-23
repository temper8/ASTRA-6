C TAUF1 [s]:	Confinement_time@(rho) for the quantity F1
C       TAUF1=Vint(F1)/(Vint(SF1TOT)-dVint(F1)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF1R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF1R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF1TOT(J)*YJ
	   YW = YW+F1(J)*YJ
	   YWO = YWO+F1O(J)*VRO(J)
	enddo
	YQ = YQ-SF1TOT(JK)*YDR
	YW = YW-F1(JK)*YDR
	YWO = YWO-F1O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF1R = 0.
	else
	   TAUF1R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
