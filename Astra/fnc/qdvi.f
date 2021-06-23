C QDVI [MW]:	 Integral {0,R} ( PDVI ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDVIR(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDVIR=0.
	do 1 J=1,JK
 1	QDVIR=QDVIR+PDVI(J)*VR(J)
	QDVIR=HRO*(QDVIR-PDVI(J)*YDR)
	end
