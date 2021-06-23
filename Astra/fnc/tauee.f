C TAUEE [s]: Electron energy confinement time at radial position R [m]
C       TAUEE=We/(Qetot-dWe/dt)
C			(Yushmanov 11-MAY-87)
	double precision function TAUEER(YR)
	implicit none
	double precision YJ,YQ,YW,YWO
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include	'for/yrjkdr.inc'
	YQ=0.
	YW=0.
	YWO=0.
	TAUEER=0.
	do 1 J=1,JK
	YJ=VR(J)
	YQ=YQ+PETOT(J)*YJ
	YW=YW+NE(J)*TE(J)*YJ
 1	YWO=YWO+NEO(J)*TEO(J)*VRO(J)
	YQ=YQ-PETOT(JK)*YDR
	YW=YW-NE(JK)*TE(JK)*YDR
	YWO=YWO-NEO(JK)*TEO(JK)*YDR
	if(YQ.eq.0.)	then
		TAUEER=0.
			else
		TAUEER=TAU*YW/((YQ*TAU*417.+YWO)-YW)
			endif
	end
