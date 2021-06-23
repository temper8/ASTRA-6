C TAUEI [s]: Ion energy confinement time at radial position R [m]
C       TAUEI=Wi/(Qitot-dWi/dt)
C			(Yushmanov 11-MAY-87)
	double precision function TAUEIR(YR)
	implicit none
	double precision YJ,YQ,YW,YWO
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YW=0.
	YWO=0.
	TAUEIR=0.
	do 1 J=1,JK
	YJ=VR(J)
	YQ=YQ+PITOT(J)*YJ
	YW=YW+NI(J)*TI(J)*YJ
 1	YWO=YWO+NIO(J)*TIO(J)*VRO(J)
	YQ=YQ-PITOT(JK)*YDR
	YW=YW-NI(JK)*TI(JK)*YDR
	YWO=YWO-NIO(JK)*TIO(JK)*YDR
	if(YQ.eq.0.)	then
		TAUEIR=0.
			else
		TAUEIR=TAU*YW/((YQ*TAU*417.+YWO)-YW)
			endif
	end
