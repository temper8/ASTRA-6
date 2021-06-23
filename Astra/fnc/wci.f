C WCI [MJ]:	Ion enegry stored in the plasma core (WI - Pedestal)
C		0.0016 * 3/2* Integral {0,R} (Ni*Ti-Nib*Tib) dV
C			(Pereverzev 26-MAR-00)
	double precision function WCIR(YR)
	implicit none
	double precision YR,VINT
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	double precision	YY(NRD)
	integer	j
	do   j = 1,NA1
	   YY(j) = NI(J)*TI(J)-NI(NA1)*TI(NA1)
	enddo
	WCIR=VINT(YY,YR)
	END
