C======================================================================|
C Electron flux GNX due to all neutral sources
C GNX	[10^19 particle/m^2/s]
C						(Pereverzev 23-FEB-98)
C			2-NOV-99 GNX redetermined to give the flux density
	subroutine GNEX
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	j
	double precision YSN1,YSN2,YSN3,YH,Y
	double precision SNNEU,SVIE,SVII,SVREC,SNNR,SNNI
	YSN1 = 0.
	YSN2 = 0.
	YSN3 = 0.
	YH = HRO
	do	1	J=1,NA1
	if (j.eq.NA1)	YH = HROA
	include	'fml/snneu'
	YSN1 = YSN1+VR(J)*NE(J)
	YSN2 = YSN2+VRO(J)*NEO(J)
	YSN3 = YSN3+VR(j)*(NE(j)*SNNEU+SNEBM(j))
C Electron sources due to
C  (i)   ionization of the wall neutrals,
C  (ii)  ionization of the beam neutrals
C  (iii) dNE/dt
	GNX(J) = ( YSN3-(YSN1-YSN2)/TAU )*YH/SLAT(J)
 1	continue
	if (NA1.ge.NAB)	return
	do	2	J=NA1+1,NAB
	GNX(J) = GNX(NA1)
 2	continue
	end
C======================================================================|
