C TUNEN  Tune <ne> to be equal to the given value GIVNNE [10^19/m^3]
C A magnitude of the neutral density NNCL is adjusted in order to
C maintain the volume average density NEAVB about GIVNNE
C GIVNNE [10^19/m^3] pre-set volume average density,
C YDT    [s]	characteristic time of the particle source response
C		Stiffness of regulation:
C		large YDT - stiff regulation, small YDT - soft regulation
C			(Pereverzev 04-MAY-01)
C----------------------------------------------------------------------|
	subroutine TUNEN(GIVNNE,YDT)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision GIVNNE,YDT,YAN,YSN,YDN,YANOLD,YC,YTAU,VINT
	save	YAN,YTAU
	data	YAN /-1./ YTAU /-1.E10/
	if (TIME .le. YTAU)	return
	YANOLD	=YAN
	YAN	=VINT(NE,ROC)/VOLUME
	if (YANOLD .lt. 0.)	return
	YSN	=VINT(SNTOT,ROC)/VOLUME
	if (YSN .le. 0.)	return
	YDN = (YAN-YANOLD)/(TIME-YTAU)
	YTAU = YAN/(YSN-YDN)
	YC = ((GIVNNE-YAN)/YTAU-YDN)/YSN
! Stiffness is externally given
C	YC = min(YC,TAU/YDT)
! Stiffness is proportional to deviation from the given value
	YC = min(YC,TAU/YTAU*(GIVNNE-YAN)/YAN)
	NNCL = NNCL*(1.+YC)
	NNCL = max(NNCL,1.d-5)
	YTAU	=TIME
	end
