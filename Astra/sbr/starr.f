C The subroutine stores NAB elements of 
C the array Yvar at the moment time=Ytime in the array YFIX
	subroutine	STARR(YVAR,YTIME,YFIX)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer	j
	double precision	YVAR(*),YTIME,YFIX(*)
	if (TIME-YTIME) 1,1,3
 1	continue
	do	2	j=1,NAB
	YFIX(j)	= YVAR(j)
 2	continue
 3	continue
	end
