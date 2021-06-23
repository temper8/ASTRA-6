C======================================================================|
C  Subroutine minimizes the value of functional
C  INTEGRAL(alfa*(dU/dx)**2+(U-F)**2)*dx with respect to U(x).
C  FO(1:NA1) is a given array on the grid X(1:NA1)
C  The result is a smoothed array FN(1:NA1) given on the same grid
C	ALFA ~ 0.01*X(NA1)**2  is a regularizator
C	The target function FN obeys the additional conditions:
C	     dFN/dx(x=0)=0 - cylindrical case
C	     FN(XN(NA1))=FO(XO(NA1))
C----------------------------------------------------------------------|
	subroutine	SMEAR1(ALFA,FO,FN,N)
	implicit none
	double precision	ALFA,FO(*),FN(*)
	integer	j,N
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	call	SGLAZH(ALFA,N,FO,RHO,N,FN,RHO)
	if (N .eq. NA1)	return
	do	j=N,NA1
	   FN(j) = FO(j)
	enddo
	end
C======================================================================|
