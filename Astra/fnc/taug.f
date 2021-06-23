C TAUG [s]: Global energy confinement time at radial position R [m]
C       TAUG=We/(Vint(PEX+PIX+PJOUL)-dWe/dt)
C			(Esipchuk 24-JAN-89)
	double precision function TAUGR(YR)
	implicit none
	double precision YQ,YW,YJ,YWO,PJOUL
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YW=0.
	YWO=0.
	TAUGR=0.
	do 1 J=1,JK
	include  'fml/pjoul'
	YJ=VR(J)
	YQ=YQ+(PEX(J)+PIX(J)+PJOUL)*YJ
	YW=YW+(NE(J)*TE(J)+NI(J)*TI(J))*YJ
 1	YWO=YWO+(NEO(J)*TEO(J)+NIO(J)*TIO(J))*VRO(J)
	YQ=YQ-(PEX(JK)+PIX(JK)+PJOUL)*YDR
	YW=YW-(NE(JK)*TE(JK)+NI(JK)*TI(JK))*YDR
	YWO=YWO-(NEO(JK)*TEO(JK)+NIO(JK)*TIO(JK))*YDR
	if(YQ.eq.0.)	then
		TAUGR=0.
			else
		TAUGR=TAU*YW/((YQ*TAU*417.+YWO)-YW)
			endif
	end
