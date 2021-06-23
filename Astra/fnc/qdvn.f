C QDVN [10^19/s]:	 Integral {0,R} ( SDVN ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDVNR(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDVNR=0.
	do 1 J=1,JK
 1	QDVNR=QDVNR+SDVN(J)*VR(J)
	QDVNR=HRO*(QDVNR-SDVN(J)*YDR)
	end

