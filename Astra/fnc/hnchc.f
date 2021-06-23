C HNCHC [m#2/s]:  ion Heat conductivity Neoclassical Chang and Hinton 
C		with a finite orbit width correction according A.Bergmann
C	Phys.Fl.-1493(1982); Phys.Fl.-3314(1986); A.Bergmann EPS-27.
C			(Pereverzev 23-NOV-00)
	double precision function HNCHCR(YR)
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	double precision
     1		YHNCHD,YHNCHC,YHNCHF,YHNCHH,YHNCHM,YHNCHK,YHNCH3,YHNCHR,
     2		YR,RLI,NUI,NUIS,HNCHI,R2POT3,YHNCHE,YHNCHZ,YHNCH2,YHNCHB
	integer	J,JK
C----------------------------------------------------------------------|
	JK = YR/HRO+1.
	if (JK.ge.NA1 .or. JK.le.0)	then
	   HNCHCR = 0.
	   return
	endif
	J = JK
	RLI=0.003235*sqrt(AMAIN(1)*TI(1))/BTOR
	R2POT3 = 135.55*RTOR*(RLI/MU(1))**2 ! (2*Potato_orbit_radius)^3
	if (R2POT3 .lt. AMETR(J)**3)	goto	1
	do J=1,NA1
	   JK = J
	   if (R2POT3 .lt. AMETR(J)**3)	goto	1
	enddo
 1	J = JK
	include	'fml/hnchi'
	HNCHCR = HNCHI
	end
