C  Electron flux GNX due to the wall (cold + warm) neutral source
C GNX	[10^19 particle/m^2/s]
C			(Pereverzev 25-FEB-88)
C			updated 29-APR-94 to include a real geometry
C			2-NOV-99 GNX redetermined to give the flux density
C					(i.e. G11 -> SLAT)
	subroutine GNXSRC
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	j
	double precision YSE1,YSE2,YSI1,YSI2,YSN3,YE,YI,Y
	double precision SVIE,SVII,SVREC!,SNNEU,SNNR,SNNI
	YSN3	=0.
	YSE1	=0.
	YSE2	=0.
	YSI1	=0.
	YSI2	=0.
	do	1	J=1,NA1
	include	'fml/svie'
	include	'fml/svii'
	include	'fml/svrec'
	YE	=VR(J)*NE(J)
	YI	=VR(J)*NI(J)
	YSN3	=YSN3+NN(J)*(SVIE*YE+SVII*YI)-SVREC*NE(J)*YI
	YSE1	=YSE1+YE
	YSE2	=YSE2+VRO(J)*NEO(J)
C GNX due to  (i) neutral source,  (ii) dNE/dt
	GNX(J)	=( YSN3*(NNCL+NNWM)-(YSE1-YSE2)/TAU )*HRO/SLAT(J)
C
C To get the ion source use the next 3 lines:
C	YSI1	=YSI1+YI
C	YSI2	=YSI2+VRO(J)*NIO(J)
C	GNX(J)	=(YSN3*(NNCL+NNWM)-(YSI1-YSI2)/TAU)*HRO/SLAT(J)
 1	continue
	end
