C QDVE [MW]:	 Integral {0,R} ( PDVE ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDVER(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDVER=0.
	do 1 J=1,JK
 1	QDVER=QDVER+PDVE(J)*VR(J)
	QDVER=HRO*(QDVER-PDVE(J)*YDR)
	end
