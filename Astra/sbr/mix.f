	subroutine	MIX(TRCOEF)
	implicit	none
	character*2	TRCOEF
C  Input:	FP(j),MU(j),NA1,HRO,ROC,TRCOEF,CMHD1,  [NE(j),RTOR]
C  Output:	RS1,RS2,RS3,RS4,	Optional: REXT=CMHD4
C			(G.Pereverzev 25-10-90)
C	HE(r):	(1-(r/(2*r0))#2)#2
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	double precision RS0,RS1,RS2,RS3,RS4,RSEXT,AMX
	double precision YY1,YY2,YYY,YRO,YDR,YRS,YCC
	integer	j,IC,JSAWT
	if (TRCOEF .eq. 'HE')	then
	   write(*,*)">>> ERROR >>> formula HMHD1 is obsolete"
	   write(*,*)"              use explicit expression"
	   stop
	endif
	if (TRCOEF .eq. 'XI')	then
	   write(*,*)">>> ERROR >>> formula XMHD1 is obsolete"
	   write(*,*)"              use explicit expression"
	   stop
	endif
	if (TRCOEF .ne. 'CC')	then
	   write(*,*)">>> ERROR calling MIX with wrong argument"
	   stop
	endif
	RS1	=0.
	RS2	=0.
	RS3	=0.
	RS4	=0.
	AMX	=0.
	JSAWT	=0
	do	j=1,NA1-1
	   YY1 = MU(j)-1.
	   YY2 = MU(j+1)-1.
	   AMX = max(YY1,AMX)
	   if (YY1*YY2 .lt. 0.)	then
	      YRS = HRO*(j+YY1/(YY1-YY2)-0.5)
	      if (RS3.ne.0. .and. RS4.eq.0.)	RS4 = YRS
	      if (RS2.ne.0. .and. RS3.eq.0.)	RS3 = YRS
	      if (RS1.ne.0. .and. RS2.eq.0.)	RS2 = YRS
	      if (RS1.eq.0.)	RS1 = YRS
	      JSAWT = max(j,JSAWT)
	   endif
	enddo
	RSEXT = max(RS1,RS2,RS3,RS4)
C	CMHD4 = RSEXT
	RS0   = min(ROC,2.*RSEXT)
	if (RS0 .lt. 0.1*HRO)	return
	AMX = 10.*min(0.1d0,AMX)
C	AMX = 10.*AMX
	IC  = 0
C	write(*,'(3(2F10.5,2X))')RS1,RS2,RSEXT,RS0,AMX,real(JSAWT)
	do	j=NA1-1,1,-1
	   YRO = HRO*J
	   YDR = RSEXT-YRO-HRO
C	   YDR = RSEXT-YRO-CMHD2*HRO
	   if (YDR.gt.0.0 .and. TRCOEF.eq.'CC')	then
	      if (IC.eq.0) YCC = YDR*(CC(J+1)-CC(J))/HRO+CC(J)
	      IC = 1
	      CC(J) = min(CC(J),YCC)
	   endif
	   YYY = max(1.-(YRO/RS0)**2,0.d0)
C	   YHE	=CMHD1*(RSEXT/ROC)**CMHD3
C	   if (TRCOEF.eq.'HE') HE(J) =HE(J)+YHE
	   if (TRCOEF.eq.'HE') HE(J) =HE(J)+CMHD1*AMX*YYY**2
C	   if (TRCOEF.eq.'XI') XI(J) =XI(J)+YHE
	   if (TRCOEF.eq.'XI') XI(J) =XI(J)+CMHD1*AMX*YYY**2
	enddo
	return
	end
