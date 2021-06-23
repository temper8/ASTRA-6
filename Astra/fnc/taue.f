C TAUE [s]: Full energy confinement time at radial position R [m]
C       TAUE=Wtot/(Qtot-dWtot/dt)
C			(Yushmanov 29-DEC-90)
	double precision function TAUER(YR)
	implicit none
	double precision YJ,YQ,YW,YWO
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	YQ=0.
	YW=0.
	YWO=0.
	TAUER=0.
	do 1 J=1,JK
	YJ=VR(J)
	YQ=YQ+(PETOT(J)+PITOT(J))*YJ
	YW=YW+(NE(J)*TE(J)+NI(J)*TI(J))*YJ
 1	YWO=YWO+(NEO(J)*TEO(J)+NIO(J)*TIO(J))*VRO(J)
	YQ=YQ-(PETOT(JK)+PITOT(JK))*YDR
	YW=YW-(NE(JK)*TE(JK)+NI(J)*TI(JK))*YDR
	YWO=YWO-(NEO(JK)*TEO(JK)+NIO(J)*TIO(JK))*YDR
	if(YQ.eq.0.)	then
		TAUER=0.
			else
		TAUER=TAU*YW/((YQ*TAU*417.+YWO)-YW)
			endif
	end
