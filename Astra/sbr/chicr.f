C======================================================================|
	subroutine CHICR(YTECR,YTICR,STE,STI,CHE,CHI)
C----------------------------------------------------------------------|
C Input:
C	YTECR	[d/l]	grad(TE) - critical gradient T_e
C	YTICR	[d/l]	R*grad(TI)/TI=R/L_Ticr - critical gradient T_i
C	STE	[m^2/s]	Stiffness chi_e
C	STI	[m^2/s]	Stiffness chi_i
C Output:
C	CHE	[m^2/s]	Electron heat conductivity
C	CHI	[m^2/s]	Ion heat conductivity
C  						(Pereverzev 24-07-2000)
C----------------------------------------------------------------------|
	implicit none
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	integer	j,j1
	double precision YTECR(*),YTICR(*),STE(*),STI(*),CHE(*),CHI(*)
	double precision P2LIN,YRTE,YRTI,YCHE(NRD),YCHI(NRD)
	double precision RLI,FTAV,YSTE,YSTI,YRTOP
	save	YRTOP
C Find anomalous transport zone: RHO < CV11*ROC <= CF11*ROC
	j1 = NA1
	do  j = 1,NA1
	   j1 = j
	   RLI=0.003235*sqrt(AMAIN(j)*TI(j))/BTOR
	   if (CV12*RLI .ge. (CF11*ROC-RHO(j))/SHEAR(j)**2)  goto	1
	enddo
 1	continue
	CV13 = RHO(j1)/ROC
C Next line is equivalent to 	YRTOP = FTAV(YRTOP,.01)
	YRTOP =RHO(j1)+EXP(-TAU/.01)*(YRTOP-RHO(j1))
	if (YRTOP .gt. CF11*ROC)	YRTOP = CF11*ROC
	CV11 = YRTOP/ROC
	do   j = 1,NA1
	   if (RHO(j) .ge. YRTOP)	goto	2
	   j1 = j
	enddo
 2	continue

	do  j = 1,j1-1
	   YRTE = abs(TE(j+1)-TE(j))/HRO	! /(AMETR(j+1)-AMETR(j))?
	   YSTI = 2.*RTOR/(TI(j+1)+TI(j))
	   YRTI = YSTI*abs(TI(j+1)-TI(j))/(AMETR(j+1)-AMETR(j))
	   YSTI = STI(j)/YSTI
	   YSTE = min(1.d0,((j1-1-j)*CIMP1/j1))	! (100/CIMP1)% inward 
	   YSTI = STI(j)*YSTE
	   YSTE = STE(j)*YSTE
	CAR23(j) = YSTE
	CAR24(j) = YSTI
	   YCHE(j) = P2LIN(YRTE,YTECR(j),YSTE,1.d-1)
	   YCHI(j) = P2LIN(YRTI,YTICR(j),YSTI,1.d-1)
	   CAR32(j) = CV5+P2LIN(RHO(j)/ROC,CF4,CF5,1.d-1)
	enddo
	do   j = j1,NA1
	   YCHE(j) = 0.
	   YCHI(j) = 0.
	CAR23(j) = 0.
	CAR24(j) = 0.
	enddo
	call	SMEARR(CHE3,YCHE,CHE)
	call	SMEARR(CHE3,YCHI,CHI)
	do  j = j1,NA1
	   CHE(j) = 0.
	   CHI(j) = 0.
	enddo
	end
C======================================================================|
	double precision function P2LIN(X,X_0,K,D)
C piecewise linear function min(0,k*(x-x_0)) is replaced with
C a smoothing function:
C	/ 0			if	x < (1-d)x_0
C	|
C	| k*(x-x0)		if	x < (1+d)x_0
C	|
C	\ 1/4 k/d (x-x_0)^2 /x_0	otherwise
C   K	stiffness
C   D	relative (with respect to x_0) broadening of the transition zone
C----------------------------------------------------------------------|
	implicit none
	double precision	X,X_0,K,D,DX
	DX = D*X_0
	if     (X .ge. X_0+DX)	then
	   P2LIN = K*(X-X_0)
	elseif (X .le. X_0-DX)	then
	   P2LIN = 0.
	else
	   P2LIN = 0.25*K*(X-X_0+DX)**2/DX
	endif
	end
C======================================================================|
