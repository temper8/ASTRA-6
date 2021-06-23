C 2nd radial derivative of an array defined on the main grid
C YAROUT = d^2(YARIN)/d(rho)^2 		(Pereverzev 20-04-1999)
	subroutine D2ARR(YARIN,YAROUT)
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	integer	j
	double precision	YARIN(*),YAROUT(*),YHRO2
	YHRO2 = 1./HRO**2
	do	j=2,NA-1
	   YAROUT(j) = YHRO2*(YARIN(j-1)+YARIN(j+1)-2.*YARIN(j))
	enddo
	YAROUT(1) = YHRO2*(YARIN(2)-YARIN(1))
	YAROUT(NA) = YAROUT(NA-1)
	YAROUT(NA1) = YAROUT(NA)
	end
