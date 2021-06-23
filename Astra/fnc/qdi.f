C QDI [MW]:	 Integral {0,R} ( PDI ) dV
C			(Pereverzev 9-OCT-07)
	double precision FUNCTION QDIR(YR)
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDIR=0.
	do 1 J=1,JK
 1	QDIR=QDIR+PDI(J)*VR(J)
	QDIR=HRO*(QDIR-PDI(J)*YDR)
	end
