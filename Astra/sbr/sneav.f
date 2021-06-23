C SNEAV  Set NEAV(a) to be equal to the given value GIVNNE [10^19/m^3]
C A magnitude of the neutral density NNCL is adjusted in order to
C increase/decrease the actual volume average density NEAVB 
C GIVNNE[10^19/m^3] pre-set volume average density
C YDT	[s]	characteristic time for density growth  ~ 0.1s
C DELT	[ ]	model stiffness: relative change 	~ 0.001 -> 0.01 
C		in the particle source per time step
C		cannot exceed DELT 
C QNmin	[10^19/s]	minimal neutral source		~ 0.002
C			(Yushmanov 27-2-89)
C----------------------------------------------------------------------|
	subroutine SNEAV(GIVNNE,YDT,DELT,QNmin)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision
     1		GIVNNE,YDT,DELT,QNmin,YAN,YSN,YDN,YANOLD,YC,TOLD,VINT
	save	YAN,TOLD
	data	YAN /-1./ TOLD /-1.E10/
	if (TIME .le. TOLD)	return
	YANOLD	=YAN
	YAN	=VINT(NE,ROC)/VOLUME
	if (YANOLD .lt. 0.)	return
	YSN	=VINT(SNTOT,ROC)/VOLUME
	if (YSN .eq. 0.)	return
	YDN	=1+max(1.d-3,DELT)
	YC = 1.+((GIVNNE-YAN)/YDT-(YAN-YANOLD)/(TIME-TOLD))/YSN
	YC = max(YC,1./YDN)
	YC = min(YC,YDN)
	if (YSN .lt. QNmin)	then
	   NNCL	= NNCL*QNmin/YSN
	else
	   NNCL	= max(NNCL*YC,1.d-5)
	endif
C	if (TIME .gt. 0.5)
C     >write(*,*)NNCL,(GIVNNE-YAN)/YDT,(YAN-YANOLD)/(TIME-TOLD),YSN
	TOLD	=TIME
	end
