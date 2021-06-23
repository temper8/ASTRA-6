C======================================================================|
	subroutine CHIMOD
     >		(YTECR,YTICR,YSTE,YSTI,YSM,YSEP,YWID,CH_E,CH_I,YATOP)
C----------------------------------------------------------------------|
C Input:
C	YTECR	[d/l]	grad(TE) - critical gradient T_e
C	YTICR	[d/l]	R*grad(TI)/TI=R/L_Ticr - critical gradient T_i
C	YSTE	[m^2/s]	Stiffness chi_e
C	YSTI	[m^2/s]	Stiffness chi_i
C       YSM	[d/l]	Smoothing factor for chi's
C       YSEP	[d/l]	YSEP*ABC is a separatrix position (SOL: a > YSEP*ABC)
C       YWID	[d/l]	YWID*RLI minimal width of the improved confinement zone
C       1.D-2 - (hidden parameter) time inertia in pedestal top
C       0.1   - (hidden parameter) relative broadening (see P2LIN)
C       0.1   - (hidden parameter) minimal electron chi_e
C Output:
C	CH_E	[m^2/s]	Electron heat conductivity
C	CH_I	[m^2/s]	Ion heat conductivity
C       YATOP	[m]	Pedestal top in the midplane coordinate AMETR
C  						(Pereverzev 22-08-2000)
C----------------------------------------------------------------------|
C              YSTE 
CH_E = chi_e = -----*(grad(TE)-YTECR)
C              YTECR
C
C              TI*YSTI  R*grad_a(TI)
CH_I = chi_i = -------*(------------ - YTICR)
C              R*YTICR       TI
C----------------------------------------------------------------------|
	implicit none
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	integer	j,j1,j2
	double precision YTECR(*),YTICR(*),YSTE(*),YSTI(*),YSM,YSEP,YWID
	double precision CH_E(*),CH_I(*),YATOP,YSE,YSI,RLI
	double precision P2LIN,YRTE,YRTI,YCH_E(NRD),YCH_I(NRD)
	save	j2
	data	j2/1/
	if (j2 .eq. 1)	YATOP = (YSEP-YWID*.001)*ABC
	do  j = NA,1,-1
	   if (AMETR(j) .ge. YSEP*ABC)	then
	      j1 = j
	      j2 = j
	      YCH_E(j) = 0.
	      YCH_I(j) = 0.
	      goto	1
	   endif
	   YRTE = abs(TE(j+1)-TE(j))/HRO
	   YSE = YSTE(j) / YTECR(j)
	   YCH_E(j) = P2LIN(YRTE,YTECR(j),YSE,1.d-1)
	   YRTI = 2.*RTOR/(TI(j+1)+TI(j))
	   YSI = YSTI(j) / YRTI / YTICR(j)
	   YRTI = abs(TI(j+1)-TI(j))/(AMETR(j+1)-AMETR(j)) * YRTI
	   YCH_I(j) = P2LIN(YRTI,YTICR(j),YSI,1.d-1)
	   RLI=0.003235*sqrt(AMAIN(j)*TI(j))/BTOR
	   if (YRTI.le.YTICR(j) .and. AMETR(j2)-AMETR(j).lt.YWID*RLI)
     >		j2 = j
 1	   continue
	enddo

C Forbid too fast (char.time < 10ms) variation of a_top
	YATOP =AMETR(j2)+EXP(-1.D2*TAU)*(YATOP-AMETR(j2))
	do  j = 1,NA1
	   if (AMETR(j) .lt. YATOP)	goto	2
	   YCH_E(j) = 0.
	   YCH_I(j) = 0.
 2	   continue
	enddo
	call	SMEARR(YSM,YCH_E,CH_E)
	call	SMEARR(YSM,YCH_I,CH_I)
	do  j = 1,NA1
	   if (AMETR(j) .lt. YATOP)	goto	3
	   CH_E(j) = 1.d-1
	   CH_I(j) = 0.d0
 3	   continue
	enddo
	end
C======================================================================|
	double precision function P2LIN(X,X_0,K,D)
C piecewise linear function min(0,k*(x-x_0)) is replaced with
C a smoothing function:
C		/ 0			if	x < (1-d)x_0
C		|
C       P2lin = | k*(x-x_0)		if	x > (1+d)x_0
C		|
C		\ 1/4 k/d (x-x_0)^2 /x_0	otherwise
C   K	stiffness
C   D	relative (with respect to x_0) broadening of the transition zone
C----------------------------------------------------------------------|
	implicit none
	double precision	X,X_0,K,D,DX
	if (D .ge. 0.001)	goto	1
	if (X .ge. X_0)		then
	   P2LIN = K*(X-X_0)
	else
	   P2LIN = 0.
	endif
	return
 1	DX = D*X_0
	if     (X .ge. X_0+DX)	then
	   P2LIN = K*(X-X_0)
	elseif (X .le. X_0-DX)	then
	   P2LIN = 0.
	else
	   P2LIN = 0.25*K*(X-X_0+DX)**2/DX
	endif
	end
C======================================================================|
