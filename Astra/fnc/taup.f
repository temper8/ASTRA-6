C TAUP [s]:	Particle confinement time at radial position R [m]
C       TAUP=Vint(Ne)/(Vint(SNTOT)-dVint(Ne)/dt)
C			(Yushmanov 11-MAY-87)
	double precision function TAUPR(YR)
	implicit none
	double precision YQ,YW,YWO,YJ
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YW=0.
	YWO=0.
	TAUPR=0.
	do 1 J=1,JK
	YJ=VR(J)
	YQ=YQ+SNTOT(J)*YJ
	YW=YW+NE(J)*YJ
 1	YWO=YWO+NEO(J)*VRO(J)
	YQ=YQ-SNTOT(JK)*YDR
	YW=YW-NE(JK)*YDR
	YWO=YWO-NEO(JK)*YDR
	if(YQ.eq.0.)	then
		TAUPR=0.
			else
		TAUPR=TAU*YW/((YQ*TAU+YWO)-YW)
			endif
	end
