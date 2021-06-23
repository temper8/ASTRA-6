C======================================================================|
C  	Subroutine minimizes the value of the functional
C		Integral[(j*rho_s*dU/drho)**2+(U-F)**2]drho
C  	with respect to U(rho).
C Input:
C   j [int d/l] sets a characteristic smoothing size (j*rho_s)
C	mapped to the coordinate RHO 
C   F=FO(1:NA1) is a given array on the grid RHO(1:NA1)
C  	The result is a smoothed array
C Output:
C   U=FN(1:NA1) given on the same grid
C  	The target function FN obeys the additional conditions:
C   dFN/dx(x=0)=0 - cylinder-like condition
C   FN(ROC)=FO(ROC)
C----------------------------------------------------------------------|
	subroutine	FEVEN(J,FO,FN)
	implicit none
	integer J
	double precision	FO(*),FN(*)
	call	RHO_SMOOTHING(J,FO,FN)
	end
C======================================================================|
	subroutine	RHO_SMOOTHING(J1,YFO,YFN)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
	integer	J,J1
	double precision YFO(*),YFN(*),YX,YP,YQ,YD,YR,RLS,YS(NRD)
	if (J1 .lt. 1)	then
	   do	j=1,NA1-1
	      YFN(j)= YFO(j)
	   enddo
	   return
	endif
	YP = 0.
	YQ = 0.
	YX = (J1/HRO)**2
	do   1	j=1,NA1-1
	   include 'fml/rls'
	   YR = YX*RLS**2
	   YD = 1.+YP+YR
	   YS(j) = YR/YD
	   YFN(j)= (YFO(j)+YQ)/YD
	   YP	= (1.-YS(j))*YR
	   YQ	= YFN(j)*YR
 1	continue
	YFN(NA1) = YFO(NA1)
	do   2	j=NA1-1,1,-1
	   YFN(j) = YS(j)*YFN(j+1)+YFN(j)
 2	continue
	end
C======================================================================|
