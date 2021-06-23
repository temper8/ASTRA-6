C SETNAV  Adjust the density of incoming neutrals NNCL in order to
C	maintain the volume average electron density <ne> around
C	the given value GIVNNE [10^19/m^3]
C Algorithm:      (n0 = GIVNNE)
C       dn         n            n0
C      ---- = S - ---,    S0 = ---- = S(1+y)
C       dt        tau          tau0
C                                              n
C      To find "y" we assume: tau0 = tau = ---------,   then
C                                          S - dn/dt
C                 n0 - n     dn    n0
C          y*S = -------- - ---- = --- - S 
C                   tau      dt    tau
C Input:
C GIVNNE [10^19/m^3] desired volume average density,
C YDT    [s]	characteristic time of the particle source response
C	Stiffness of regulation:
C		large YDT - soft regulation, small YDT - stiff regulation
C 		If YDT < 0 then proportional regulation is enabled which
C		works reasonably when <ne> is sufficiently close to GIVNNE
C Output NNCL
C						(Pereverzev 04-MAY-01)
C----------------------------------------------------------------------|
	subroutine SETNAV(GIVNNE,YDT)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	j
	double precision GIVNNE,YDT,YAN,YSN,YDN,YAO,YC,YTAU,VINT,YM
	save	YAN,YTAU
	data	YAN /-1./ YTAU /-1.E10/
	if (TIME .le. YTAU)	return
	YAO = YAN
	YAN = VINT(NE,ROC)/VOLUME			! <n_e>
	if (YAO .lt. 0.)	return
	YSN    = VINT(SNTOT,ROC)/VOLUME			! <S_e>
	if (YSN .le. 0.)	return
C	if (YM .gt. 2.*GIVNNE)	then
C	endif
	YDN = (YAN-YAO)/(TIME-YTAU)			! d<n_e>/dt
	YAO = 0.d0
	YM = NE(1)
	do	j=2,NA1					! Find max(NE)
	   if (NE(j) .gt. NE(j-1))	then
	      YM = NE(j)
	      if (j .lt. NA)
     >		YAO = (NE(j+1)+NE(j-1)-2.*NE(j))/HRO**2
	   endif
	enddo
	YTAU = YAN/(YSN-YDN)				! tau
	YC = ((GIVNNE-YAN)/YTAU-YDN)/YSN
	if (YDT .gt. 0)	then	! Stiffness is externally given
	   if (GIVNNE.gt.YAN .and. YC.lt.0)	YC = TAU/YDT
	   YC = min(YC,TAU/YDT)
	else			! Proportional stiffness
	   YSN = (GIVNNE-YAN)/YAN
	   YDN = max(-1.d0,min(1.d0,YSN))
	   YC = min(YC,TAU/YTAU*YDN)
	endif
C if |d^2n/dr^2| is too high and <n> is growing then reduce puff
	if (YAO.lt.-1.d1 .and. YDN.gt.0.d0)	then
	   NNCL = 9.d-1*NNCL		! Very hollow profile
	else
	   NNCL = NNCL*(1.+YC)
	endif
	NNCL = max(NNCL,1.d-6)
	if (YM .gt. 1.3*NE(1))	NNCL = 1.d-6	! Very hollow profile
C	write(*,'(6F10.3,1P,E12.3)')
C     >		YAN,YM,1.3*NE(1),1+YC,YAO,YDN,NNCL
	YTAU	=TIME
	end
