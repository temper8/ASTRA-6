c WCX [MJ]:	Exp. thermal enegry stored in the plasma core (WTOTX - Pedestal)
C	0.0016*3/2*Integral {0,R} (Nix*Tix+Nex*Tex-Nixb*Tixb-Nexb*Texb)dV
C			(Pereverzev 26-MAR-00)
	double precision function WCXR(YR)
	implicit none
	double precision YR,VINT
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	double precision	YY(NRD)
	integer	j
	do   j = 1,NA1
	   YY(j) = NIX(J)*TIX(J)-NIX(NA1)*TIX(NA1)
     >		  +NEX(J)*TEX(J)-NEX(NA1)*TEX(NA1)
	enddo
	WCXR=VINT(YY,YR)
	END
