C WCIX [MJ]:	Ion thermal enegry stored in the plasma core (WIX - Pedestal)
C		0.0016 * 3/2* Integral {0,R} (Nix*Tix-Nixb*Tixb) dV
C			(Pereverzev 26-MAR-00)
	double precision function WCIXR(YR)
	implicit none
	double precision YR,VINT
	integer  J
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	double precision	YY(NRD)
	do   j = 1,NA1
	   YY(j) = NIX(J)*TIX(J)-NIX(NA1)*TIX(NA1)
	enddo
	WCIXR=VINT(YY,YR)
	end
