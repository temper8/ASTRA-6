C  This subroutine calculates the plasma movement/expansion so that  
C  the minor radius "a"=ABC satisfies the condition q(a)~=QEDGE
C  If OUT>0 then plasma current starts at the outward border
	subroutine LAYER(QEDGE,OUT)
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	double precision QEDGE,OUT,FACT
	if(ELONG.ge.ELONM.and.ABC.ge.AB)	then
		ABC	=AB
		ELONG	=ELONM
		SHIFT	=0.
		return
		endif
	if(1./MU(NA).ge.QEDGE)  return
	if(ABC.ge.AB)	then
		ABC	=AB
		SHIFT	=0.
		ELONG	=min(ELONG*1.01,ELONM)
	else
		FACT	=1.+.1/(NA1-0.5)
		ABC	=min(AB,ABC*FACT)
		SHIFT	=sign(AB-ABC,OUT)
	endif
	TAU	=TAUMIN
	end
