C WCEX [MJ]:	Electron enegry stored in the plasma core (WEX - Pedestal)
C		0.0016 * 3/2* Integral {0,R} (Nex*Tex-Nexb*Texb) dV
C			(Pereverzev 26-MAR-00)
	double precision function WCEXR(YR)
	implicit none
	double precision YR,VINT
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	double precision	YY(NRD)
	integer	j
	do   j = 1,NA1
	   YY(j) = NEX(J)*TEX(J)-NEX(NA1)*TEX(NA1)
	enddo
	WCEXR=VINT(YY,YR)
	END
