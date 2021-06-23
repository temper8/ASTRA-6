C TAUIG [s]: Global ion energy confinement time at radial position R [m]
C			(Yushmanov 20-DEC-90)
	double precision function TAUIGR(YR)
	implicit none
	double precision PEICL,YJ,YQ,YW,YWO,COULG
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YW=0.
	YWO=0.
	TAUIGR=0.
	do 1 J=1,JK
	include  'fml/peicl'
	YJ=VR(J)
	YQ=YQ+(PIX(J)+PEICL)*YJ
	YW=YW+NI(J)*TI(J)*YJ
 1	YWO=YWO+NIO(J)*TIO(J)*VRO(J)
	YQ=YQ-(PIX(JK)+PEICL)*YDR
	YW=YW-NI(JK)*TI(JK)*YDR
	YWO=YWO-NIO(JK)*TIO(JK)*YDR
	if(YQ.eq.0.)	then
		TAUIGR=0.
			else
		TAUIGR=TAU*YW/((YQ*TAU*417.+YWO)-YW)
			endif
	end
