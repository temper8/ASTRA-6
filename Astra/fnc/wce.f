C WCE [MJ]:	Electron enegry stored in the plasma core (WE - Pedestal)
C		0.0016 * 3/2* Integral {0,R} (Ne*Te-Neb*Teb) dV
C			(Pereverzev 26-MAR-00)
	double precision function WCER(YR)
	implicit none
	double precision YR,VINT
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	double precision	YY(NRD)
	integer	j
	do   j = 1,NA1
	   YY(j) = NE(J)*TE(J)-NE(NA1)*TE(NA1)
	enddo
	WCER=VINT(YY,YR)
	END
