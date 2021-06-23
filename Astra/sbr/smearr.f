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
	subroutine	SMEARR(ALFA,FO,FN)
	implicit none
	double precision	ALFA,FO(*),FN(*)
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
	call	SGLAZH(ALFA,NA1,FO,RHO,NA1,FN,RHO)
	end
C======================================================================|
	subroutine	SGLAZH(ALFA,NO,FO,XO,N,FN,XN) ! same as SMOOTH
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'		! NRD !
	integer	NO,N,J,I
	double precision	ALFA,XO(*),FO(*),XN(*),FN(*),P(NRD)
	double precision	YF,YX,YP,YQ,YD,FJ
	if (N .gt. NRD .or. NO .le. 0)	then
		write(*,*)' >>> SMEARR: array is out of limits'
		stop
	endif
	if (NO .eq. 1)	then
	   do	j=1,N
		FN(j) = FO(1)
	   enddo
	   return
	endif
	if (NO .eq. 2)	then
	  do	j=1,N
	   FN(j)=(FO(2)*(XN(j)-XO(1))-FO(1)*(XN(j)-XO(2)))/(XO(2)-XO(1))
	  enddo
	  return
	endif
	if (N .lt. 2)	then
		write(*,*)' >>> SMEARR: no output grid is provided'
		stop
	endif
	if (abs(XO(NO)-XN(N)) .gt. XN(N)/N)	then
	    write(*,*)'>>> SMEARR: grids are not aligned'
	    write(*,'(1A23,I4,F8.4)')'     Old grid size/edge',NO,XO(NO)
	    write(*,'(1A23,I4,F8.4)')'     New grid size/edge',N,XN(N)
	    stop
	endif
	do	1	j=2,N
	P(j)	=ALFA/(XN(j)-XN(j-1))/XO(NO)**2
C	P(j)	=ALFA*sqrt(XO(j)/XO(NO))/(XN(j)-XN(j-1))/XO(NO)**2
C	P(j)	=ALFA*(XO(j)/XO(NO))/(XN(j)-XN(j-1))/XO(NO)**2
 1	continue
	P(1)	=0.
	FN(1)	=0.
	I	=1
	YF	=(FO(2)-FO(1))/(XO(2)-XO(1))
	YX	=2./(XN(2)+XN(1))
	YP	=0.
	YQ	=0.
	do	5	j=1,N-1
		if(XO(I) .gt. XN(j))	GO TO 4
 3		I	=I+1
		if(I .gt. NO)	I=NO
		if(I .ne. NO .and. XO(I) .lt. XN(j))	GOTO	3
		YF	=(FO(I)-FO(I-1))/(XO(I)-XO(I-1))
 4		FJ	=FO(I)+YF*(XN(j)-XO(I))
		YD=1.+YX*(YP+P(j+1))
		P(j)	=YX*P(j+1)/YD
		FN(j)	=(FJ+YX*YQ)/YD
		YX	=2./(XN(j+2)-XN(j))
		YP	=(1.-P(j))*P(j+1)
		YQ	=FN(j)*P(j+1)
 5	continue
	FN(N)	=FO(NO)
	do	6	j=N-1,1,-1
		FN(j)	=P(j)*FN(j+1)+FN(j)
 6	continue
	end
C======================================================================|
