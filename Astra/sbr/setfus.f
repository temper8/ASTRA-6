C SETFUS  Adjust the fusion power at the prescribed level YPSET [MW]
C Input:
C    YPSET  [MW] Required alpha power
C    YTSET  [s]  Stiffness of regulation: large YTSET - soft regulation,
C		  			  small YTSET - stiff regulation
C    YPMAX  [MW] Maximum power in each of two sources
C    YPSTEP [MW] Minimum power step in each source
C Output: YPOUT1,YPOUT2
C						(Pereverzev 04-MAY-01)
C----------------------------------------------------------------------|
	subroutine SETFUS(YPSET,YPOUT1,YPOUT2)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	j
	double precision
     1		YPSET,YP,YSN,YDP,YPOLD,YC,YTOLD,VINT,YPMAX,
     2		PDT,SVDT,YSVDT,YPOUT1,YPOUT2,YTAU,YTSET,YPSTEP
	save	YP,YTOLD
	data	YP /-1./ YTOLD /-1.E10/
	if (TIME .le. YTOLD)	return
	YPOLD = YP
	do   j	=1,NA1
	   include 'fml/pdt'
	   WORK1(j,1) = pdt
	enddo
	YP = VINT(WORK1(1,1),ROC)
	if (YPOLD .lt. 0.)	then
	   YTOLD = TIME
	   return
	endif

	YTSET = 1.				! Stiffnes
	YPMAX = 100.				! 200 MW maximum power
	YPSTEP = 5.				! 5 MW minimum power step

	YDP = (YP-YPOLD)/(TIME-YTOLD)		! dP_fus/dt
	YTAU = (YPSET-YP)/YDP
	YTOLD = TIME
	if (YTAU)	1,1,2
 1	YDP = sign(YPSTEP,YPSET-YP)		! Wrong trend -> unconditional
	goto	3				!	   control needed

 2	continue				! Correct trend but the
	YDP = sign(YPSTEP,YDP)			! approach rate has to be
	if (YTAU .ge. YTSET)	goto	3	! 	-> accelerated 
	YDP =-YDP				! 	-> decelerated 
	if (abs(YPSET-YP)/YPSET.lt..1) YDP = 0.	!	-> retained
 3	continue

	YC = abs(YPSET-YP)/YPSET
	YPOUT1 = YPOUT1+YDP
	YPOUT1 = max(0.d0,min(YPOUT1,YPMAX))
	if (YPOUT1.gt.0. .and. YPOUT1.lt.YPMAX)	return
	YPOUT2 = YPOUT2+YDP
	YPOUT2 = max(0.d0,min(YPOUT2,YPMAX))
	return
	end
