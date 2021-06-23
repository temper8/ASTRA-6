c WC [MJ]:	Exp. thermal enegry stored in the plasma core (WTOTX - Pedestal)
C	0.0016*3/2*Integral {0,R} (Nix*Tix+Nex*Tex-Nixb*Tixb-Nexb*Texb)dV
C			(Pereverzev 26-MAR-00)
	double precision function WCR(YR)
	implicit none
	double precision YR,VINT
	integer  J
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	double precision	YY(NRD)
	do   j = 1,NA1
	YY(j) = NI(J)*TI(J)-NI(NA1)*TI(NA1)+NE(J)*TE(J)-NE(NA1)*TE(NA1)
	enddo
	WCR=VINT(YY,YR)
	END
