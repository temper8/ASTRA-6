C TAUF8 [s]:	Confinement_time@(rho) for the quantity F8
C       TAUF8=Vint(F8)/(Vint(SF8TOT)-dVint(F8)/dt)
C			(Pereverzev 09-OCT-08)
	double precision function TAUF8R(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ = 0.
	YW = 0.
	YWO = 0.
	TAUF8R = 0.
	do	J = 1,JK
	   YJ = VR(J)
	   YQ = YQ+SF8TOT(J)*YJ
	   YW = YW+F8(J)*YJ
	   YWO = YWO+F8O(J)*VRO(J)
	enddo
	YQ = YQ-SF8TOT(JK)*YDR
	YW = YW-F8(JK)*YDR
	YWO = YWO-F8O(JK)*YDR
	if (YQ .eq. 0.d0)	then
	   TAUF8R = 0.
	else
	   TAUF8R = TAU*YW/((YQ*TAU+YWO)-YW)
	endif
	end
