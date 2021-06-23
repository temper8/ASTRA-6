C WTOT [MJ]:	Integral {0:R} ( 3/2*NE*(TE+Ti) ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION WTOTR(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	WTOTR=0.
	DO 1 J=1,JK
 1	WTOTR=WTOTR+(NE(J)*TE(J)+NI(J)*TI(J))*VR(J)
	WTOTR=HRO*(WTOTR-(NE(JK)*TE(JK)+NI(JK)*TI(JK))*YDR)*.0024
	END
