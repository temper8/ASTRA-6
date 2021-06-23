C QDN [10^19/s]:	 Integral {0,R} ( SDN ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDNR(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDNR=0.
	do 1 J=1,JK
 1	QDNR=QDNR+SDN(J)*VR(J)
	QDNR=HRO*(QDNR-SDN(J)*YDR)
	end

