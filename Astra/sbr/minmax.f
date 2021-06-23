C The subroutine finds the minimal and maximal values of the array ARREXT
C over the time and returns them in the arrays ARRMIN and ARRMAX so that
C ARRMIN(rho) = min_t{ARREXT(rho,t);  ARRMAX(rho) = max_t{ARREXT(rho,t)}
C
	subroutine	MINMAX(ARREXT,ARRMIN,ARRMAX)
	implicit none
	include	'for/parameter.inc'
	include  'for/const.inc'
	double precision	ARREXT(*),ARRMIN(*),ARRMAX(*)
	integer	JCALL,j
	save	JCALL
	data	JCALL/1/
	goto (10,8),JCALL

 8	do	9	j = 1,NA1
		if (ARREXT(j) .gt. ARRMAX(j))	ARRMAX(j) = ARREXT(j)
		if (ARREXT(j) .lt. ARRMIN(j))	ARRMIN(j) = ARREXT(j)
 9	continue
	return

 10	JCALL = 2
	do	11	j = 1,NA1
		ARRMIN(j)= ARREXT(j)
		ARRMAX(j)= ARREXT(j)
 11	continue
	END
