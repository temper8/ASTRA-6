C The subroutine stores the variable Yvar(time=Ytime) in YFIX
	subroutine	STVAR(YVAR,YTIME,YFIX)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	YVAR,YTIME,YFIX
	if (TIME-YTIME) 1,1,2
 1	YFIX	= YVAR
 2	continue
	END
