C SETQ  Adjust an auxiliary heating power in order to maintain
C	the Q around the given value QREQ
C Input:
C QREQ  required V value,
C YDT   [s] characteristic time of the particle source response or
C	Stiffness of regulation:
C		large YDT - soft regulation, small YDT - stiff regulation
C 		If YDT < 0 then proportional regulation is enabled which
C		works reasonably when <ne> is sufficiently close to QREQ
C Output NNCL
C						(Pereverzev 01-MAR-02)
	subroutine SETQ(QREQ,YDT)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision QREQ,YDT,YAN,YSN,YDN,YANOLD,YC,YTAU,VINT,QDTR
	save	YAN,YTAU
	data	YAN /-1./ YTAU /-1.E10/
	if (TIME .le. YTAU)	return
	YANOLD	=YAN
	YAN = QDTR(ROC)

	YAN	=VINT(NE,ROC)/VOLUME
	if (YANOLD .lt. 0.)	return
	YSN	=VINT(SNTOT,ROC)/VOLUME
	if (YSN .le. 0.)	return
	YDN = (YAN-YANOLD)/(TIME-YTAU)
	YTAU = YAN/(YSN-YDN)
	YC = ((QREQ-YAN)/YTAU-YDN)/YSN
	YSN = (QREQ-YAN)/YAN
	if (YDT .gt. 0)	then	! Stiffness is externally given
	   YC = min(YC,TAU/YDT)
	else			! Proportional stiffness
	   YDN = max(-1.d0,min(1.d0,(QREQ-YAN)/YAN))
	   YC = min(YC,TAU/YTAU*YDN)
	endif
	NNCL = NNCL*(1.+YC)
	NNCL = max(NNCL,1.d-5)
	YTAU	=TIME
	end
