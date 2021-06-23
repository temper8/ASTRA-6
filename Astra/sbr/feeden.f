C======================================================================|
C The subroutine emulates a Gaussian particle source due to pellet evaporation 
C if an actual average density drops below the required value
!----------------------------------------------------------------------|
! Limitations:	(1) TAU <= YEVTIM
!		(2) The subroutine call should preceed density equation 
!		    call (otherwise, a time delay can cause lost of 
!		    particles). If TAU << YEVTIM the second requirement 
!		    can be abandoned.
!----------------------------------------------------------------------|
C  Input parameters:
C YAVN    [10^19 1/m^3] prescribed average density
C YDEN(*) [10^19 1/m^3] controlled density (eg. NE or F1)
C YNTOT   [10^19]	pellet inventory
C YEVTIM  [s]		evaporation time
C YCENTR  []		(centre) 	parameters for particle
C YWIDTH  []		(width)		deposition function
C  Output parameters:
C YSN(*)  [10^19 1/m^3/s] rhs for use in density eq. (deposition function)
C YFL	  [] 		flag (=1) while pellet is on, (=0) otherwise
C
C Gaussian deposition function P=P0*exp(-((r-r0)/a/width)^2)
C  r    [m]	- coordinate in the equatorial cross-section 
C  r0   [d/l]	- relative position of the centre of distribution 
C width [d/l]	- relative width, normalized to the plasma radius 
C  a	[m]	- minor radius
C
C Usage examples:
C  1)	FEEDEN(ZRD9X,NE,200.,.001,.7,.1,CAR1,CV1):;	NE:EQ;	SN=CAR1;
C  2)	FEEDEN(ZRD9X,NE,200.,.001,.7,.1,SN,CV1):;	NE:EQ;	SN=0.;
C----------------------------------------------------------------------|
	subroutine FEEDEN(YAVN,YDEN,YNTOT,YEVTIM,YCENTR,YWIDTH,YSN,YFL)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	 j,ITSME,IY
	double precision YAVN,YNTOT,YCENTR,YWIDTH,YEVTIM,YDEN(*),YSN(*)
	double precision YPROF(NRD,4),YT(4),Y1,YNUM,YRO,YWH,VINT,YFL
	save	 YT,YPROF,ITSME
	data	 YT/4*-1.E37/ ITSME/0/
	call	getid(YDEN,itsme,IY)
	if (IY .ge. 0 .and. IY .le. 4) goto	1
	if (IY .gt. 4) write(*,*)" >>> FEEDEN >>> too many calls: >",4
	if (IY .eq.-2) write(*,*)"            Calling getID from FEEDEN"
	if (IY .eq.-1) write(*,*)" >>> GETID >>> too many calls: >",2200
	return					! Error !

 1	continue
	Y1 = YT(IY)+YEVTIM-TIME
	if (Y1 .le. 0.)	  goto	2	! -> next pellet
	if (Y1 .ge. TAU)	then	! previous pellet
	   do	j = 1,NA1
	      YSN(j) = YPROF(j,IY)
	   enddo
	else
	   do	j = 1,NA1
	      YSN(j) = YPROF(j,IY)*Y1/TAU
	   enddo
	endif
	return

 2	continue
	Y1 = VINT(YDEN,ROC)/VOLUME
	if (Y1 .gt. YAVN)	then
	   do	j = 1,NA1
	      YSN(j) = 0.
	   enddo
	   YFL = 0.
	   return
	endif
	YRO = YCENTR*ROC
	YWH = YWIDTH*ROC
	do	j = 1,NA1
	   Y1 = (RHO(j)-YRO)/YWH 
	   YPROF(j,IY)= exp(-Y1*Y1)
	enddo
	YNUM = VINT(YPROF(1,IY),ROC)
	Y1 = YNTOT/(YNUM*max(YEVTIM,TAU))
!	if (YEVTIM .lt. TAU)	!warning
	do	j = 1,NA1
	   YSN(j) = Y1*YPROF(j,IY)
	   YPROF(j,IY) = YSN(j)
	enddo
	YT(IY) = TIME
	YFL = 1.
	end
C======================================================================|
