C======================================================================|
	subroutine	BNDRY1(RZPB)
C----------------------------------------------------------------------|
C 1) In case of the plasma boundary defined by 3 moments,
C    this subroutine writes 8 points on the boundary into array BNDARR
C    and into arrays RPB(1:NBND), ZPB(1:NBND)
C 2) If the plasma boundary is defined by a data file then
C    this subroutine uses the array BNDARR as an input and
C    produces output in [time dependent] arrays RPB, ZPB
C----------------------------------------------------------------------|
C NBND     number of points on the plasma vacuum boundary
C NBNT     number of times for the plasma boundary evolution
C  call from ESC:
C		call	BNDRY(RPB,ZPB)
C  call from SPIDER:
C		call	BNDRY(RZPB,RZPB(NBND+1))
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	integer	j,j1,jt
	double precision	ydt,yd1,yd2,yd3
	double precision	RZPB(*)
	if (NBNT .eq. 0)	then	! No boundary points given
	   NBND = 8
	   yd1 = .75d0			! sin^2(pi/3)	fi=pi/3
	   yd2 = .5d0			! cos(pi/3)	
	   yd3 = sqrt(yd1)		! sin(pi/3)
	   yd1 = .5d0			! sin^2(pi/4)	fi=pi/4
	   yd2 = sqrt(yd1)		! cos(pi/4)	
	   yd3 = yd2			! sin(pi/3)
	   RZPB(1) = RTOR+SHIFT+ABC			! R @ 0
	   RZPB(2) = RTOR+SHIFT-ABC*(TRIAN*yd1-yd2)	! R @ fi
	   RZPB(3) = RTOR+SHIFT-ABC*TRIAN		! R @ pi/2
	   RZPB(4) = RTOR+SHIFT-ABC*(TRIAN*yd1+yd2)	! R @ pi-fi
	   RZPB(5) = RTOR+SHIFT-ABC			! R @ pi
	   RZPB(6) = RTOR+SHIFT-ABC*(TRIAN*yd1+yd2)	! R @ pi+fi
	   RZPB(7) = RTOR+SHIFT-ABC*TRIAN		! R @ -pi/2
	   RZPB(8) = RTOR+SHIFT-ABC*(TRIAN*yd1-yd2)	! R @ -fi
	   RZPB(NBND+1) = UPDWN
	   RZPB(NBND+2) = UPDWN+ABC*ELONG*yd3
	   RZPB(NBND+3) = UPDWN+ABC*ELONG
	   RZPB(NBND+4) = UPDWN+ABC*ELONG*yd3
	   RZPB(NBND+5) = UPDWN
	   RZPB(NBND+6) = UPDWN-ABC*ELONG*yd3
	   RZPB(NBND+7) = UPDWN-ABC*ELONG
	   RZPB(NBND+8) = UPDWN-ABC*ELONG*yd3
	   BNDARR(1) = RZPB(1)
	   BNDARR(2) = RZPB(NBND+1)
	   BNDARR(3) = RZPB(2)
	   BNDARR(4) = RZPB(NBND+2)
	   BNDARR(5) = RZPB(3)
	   BNDARR(6) = RZPB(NBND+3)
	   BNDARR(7) = RZPB(4)
	   BNDARR(8) = RZPB(NBND+4)
	   BNDARR(9) = RZPB(5)
	   BNDARR(10) = RZPB(NBND+5)
	   BNDARR(11) = RZPB(6)
	   BNDARR(12) = RZPB(NBND+6)
	   BNDARR(13) = RZPB(7)
	   BNDARR(14) = RZPB(NBND+7)
	   BNDARR(15) = RZPB(8)
	   BNDARR(16) = RZPB(NBND+8)
	   return
	endif
	jt = 1
	if (NBNT .gt. 1)	goto	2
 1	continue
	do	j=1,NBND
	   j1 = jt+(2*j-1)*NBNT
	   RZPB(j) = BNDARR(j1)
	   RZPB(NBND+j) = BNDARR(j1+NBNT)
	enddo
	return
 2	continue
	if (TIME .le. BNDARR(1))	goto	1
	jt = NBNT
	if (TIME .ge. BNDARR(NBNT))	goto	1

	do	j=1,NBNT
	   if (TIME .gt. BNDARR(j)) jt = j
	enddo
	if (jt .eq. NBNT)	
     ,	write(*,*)"last time slice in the plasma boundary in exp/file"
	ydt = BNDARR(jt+1)-BNDARR(jt)
	yd1 = (TIME-BNDARR(jt))/ydt
	yd2 = (TIME-BNDARR(jt+1))/ydt
	RZPB(1:NBND:1) = BNDARR(NBNT+1::NBNT)
	do	j=1,NBND
	   j1 = jt+(2*j-1)*NBNT
	   RZPB(j) = yd1*BNDARR(j1+1)-yd2*BNDARR(j1)
	   j1 = jt+2*j*NBNT
	   RZPB(NBND+j) = yd1*BNDARR(j1+1)-yd2*BNDARR(j1)
	enddo
	end
C======================================================================|
