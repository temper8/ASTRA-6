C Module call sequence:
C ADCMP -> METRIC -> RHSEQ ->
C                 | -> EQUIL3(E3ASTR,NEWGRD,EDCELL)     |  
C          ->     |                                     |  -> Exit
C                 | -> A2ESC(BNDRY,ESC4A,NEWGRD,EDCELL) |
C
C                 | -> EQGB3
C                 |             | -> EQLVU3
C	E3ASTR -> | -> EQAB3 -> | -> EQK3
C                 |             | -> EQC1
C                 | -> EQPPAB
C----------------------------------------------------------------------|
C  Equilibrium calculation is bypassed when NEQUIL <= 1
C For (NEQUIL = 0 ) metric is prescribed by a simple formula
C For (NEQUIL = -1) metric is taken from a data file
C For (NEQUIL = 1 ) metric is frozen (can be used interactively)
C======================================================================|
	subroutine	METRIC
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	double precision ROC3A
	integer	jn,j
	save	jn
	data	jn/0/
C----------------------------------------------------------------------|
	call	ADDTIME(CPT)
	call	markloc("METRIC"//char(0))
	MEQUIL = 1.d2*(NEQUIL-int(NEQUIL))
C	if (LEQ(5) .ne. -1)	goto	1

C	goto	2
 1	continue
C	write(*,*)"Entry METRIC",NEQUIL,LEQ(5)
	if (LEQ(5).lt.0 .and. NEQUIL.lt.0)	then
C	   Neither NEQUIL nor EquilibriumSolver is defined in the model
	elseif (LEQ(5) .lt. 0)	then
C	   Old type definition: NEQUIL.ge.0 Use NEQUIL
C	   if     (NEQUIL .lt.-1)	then
C	      keep LEQ(5) unchanged and use external equilibrium
	   if (NEQUIL .eq. 0)	then
	      LEQ(5) = 0			! No equilibrium solver
	   elseif (int(NEQUIL) .eq. 42)	then
	      LEQ(5) = 2			! Use ESC solver
	   elseif (NEQUIL.gt.1 .and. MEQUIL.gt.1) then
	      LEQ(5) = 3			! Use SPIDER solver
	   else
	      LEQ(5) = 1			! Use 3m solver
	   endif
	else					! LEQ(5) >= 0
C	   New type definition: EquilibriumSolver is defined in the model
C	   Once LEQ(5) is defined, LEQ(5) it has priority over NEQUIL
	   LEQ(5) = LEQ(5)			! 
	endif
 2	continue
C	write(*,*)LEQ(5),NBND,NBNT,NEQUIL,TIME

C-------------------- Metric replacement ------------------------------|
C Cylindircal case, No equilibrium solver (NEQUIL = -2), No toroidicity
	if (NEQUIL .eq. -2)	then
C	   write(*,*)"Calling EQCYL"
	   call	add2loc("Calling EQCYL"//char(0))
	   call	EQCYL
	   call	RHSEQ
	   goto	9
	endif
C-------------------- Metric initiation -------------------------------|
C No equilibrium solver (NEQUIL<=0) .or. data initiation @ 1st entry
	if (LEQ(5).eq.0 .or. jn.eq.0)	then
C	   write(*,*)"Calling EQGUESS"
	   jn = 1
	   call	add2loc("Calling EQGUESS"//char(0))
	   call	EQGUESS
	   call	add2loc("Calling RHSEQ"//char(0))
	   goto	9
	endif
C-------------------- Pre-calculated metric ---------------------------|
	if (LEQ(5) .eq. -1)	then	! Take metrics from exp/data_file
C	   write(*,*)"Calling SETEXM"
	   call	add2loc("Calling ROC3A"//char(0))
	   ROC = ROC3A(RTOR,SHIFT,ABC,ELONG,TRIAN)
	   call	add2loc("Calling SETGEO"//char(0))
	   call	SETGEO(0)
	   call	add2loc("Calling SETEXM"//char(0))
	   call	SETEXM			! Main grid: (j-00.5d0)*h
	   call	add2loc("Calling RHSEQ"//char(0))
	   call	RHSEQ
	   goto	9
	endif
C-------------------- Equilibrium solver ------------------------------|
	if (TIME-TIMEQL .ge. DTEQL)	then
	   call	RHSEQ			! Define p', FF', j_tor=CUTOR
	   if      (LEQ(5) .eq. 1)	then
C	      write(*,*)"Calling EQUIL3",MU(1)
	      call	EQUIL3
	      call	ADDTIME(CPTEQL)
	   elseif  (LEQ(5) .eq. 2)	then
C	      write(*,*)"Calling A2ESC"
	      call	A2ESC
	   elseif  (LEQ(5) .eq. 3)	then
C	      write(*,*)"Calling a_spid"
C	      call	EQUIL3
C	      call	e_spid
	      call	a_spid
	      call	ADDTIME(CPTEQL)
	   else
	      write(*,*)"Unknown EQLsolver requested"
	      call	a_stop
	   endif
	   TIMEQL = TIME
	   goto	9
	endif
 9	continue
	call	ADDTIME(CPT)
	end
C======================================================================|
	subroutine	EQCYL
C----------------------------------------------------------------------|
C No equilibrium solver (LEQ(5) < 0) .or. data initiation @ 1st entry
C----------------------------------------------------------------------|
C  In:	RTOR,SHIFT,ABC,ELONG,TRIAN,NA1,NB1
C Out:	NA,RHO(j),DRODA,VOLUM,IPOL,G33,GRADRO,G11,G22,SLAT
C
Calling:
C       SETGEO -> SHIF, SHIV, ELONG, TRIA, DRODA, AMETR
C	NEWGRD -> Takes	ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC,  AMETR(NA1)
C		  Returns:
C		  NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, HROA, AMETR(j>NA1)
C	EDCELL -> 
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
	integer	j
C	write(*,*)"Entry EQCYL ",ROC,RHO(NA1),AMETR(NA1),ABC
	VOLUME = GP2*GP*RTOR*ABC*ABC
C       write(*,*)NA1,NB1,VOLUME
	ROC	= ABC
	ROWALL = ABC*1.1d0
	HRO	= ABC/(NA1-0.5d0)
C       HRO	= 1.1d0*ROWALL/(NB1-0.5)
	HROA = HRO
	do	j=1,NB1
	   RHO(J) = (J-0.5d0)*HRO
C          if (RHO(j) .lt. ROC)	NA = j
	   if (RHO(j) .le. ROWALL)	NA = j
	   G22(J) = J*HRO
	   VR(J)  = GP2*GP2*RTOR*RHO(J)
	   VRS(j) = GP2*GP2*RTOR*G22(J)
	   AMETR(J) = RHO(J)
	   SHIF(J) = 0.d0
	   SHIV(J) = 0.d0
	   ELON(J) = 1.d0
	   TRIA(J) = 0.d0
	   IPOL(J) = 1.d0
	   G33(J) = 1.d0
	   G11(J) = VRS(j) 
	   SLAT(J) = VRS(j)
	   BDB02(j) = 1.d0+(RHO(j)*MU(j)/RTOR)**2
	   B0DB2(j) = 1.d0/BDB02(j)
	   BDB0(j) = sqrt(BDB02(j))
	   FOFB(j) = 1.d0
	   BMAXT(j) = BTOR*BDB0(j)
	   BMINT(j) = BMAXT(j)
	   DRODA(j) = 1.d0
	   GRADRO(j) = 1.d0
	enddo
	if (NA .lt. NB1)	NB1 = NA
	NA = NA1-1
	VOLUM(1) = HRO*VR(1)
	do	J=2,NB1
	   VOLUM(J) = VOLUM(J-1)+HRO*VR(J)
	enddo
	VOLUM(NA1) = VOLUM(NA)+(HROA-0.5d0*HRO)*VR(NA)
	VOLUME = VOLUM(NA1)
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	subroutine	EQGUESS
C----------------------------------------------------------------------|
C No equilibrium solver (LEQ(5) < 0) .or. data initiation @ 1st entry
C----------------------------------------------------------------------|
C  In:	RTOR,SHIFT,ABC,ELONG,TRIAN,NA1,NB1
C Out:	NA,RHO(j),DRODA,VOLUM,IPOL,G33,GRADRO,G11,G22,SLAT
C
Calling:
C       SETGEO -> SHIF, SHIV, ELONG, TRIA, DRODA, AMETR
C	NEWGRD -> Takes	ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC,  AMETR(NA1)
C		  Returns:
C		  NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, HROA, AMETR(j>NA1)
C	EDCELL -> 
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
	integer	j,j1,JNA1O
	double precision ROC3A,YA,YAS,YES,YDS,YDV,YIPL,YR1
C----------------------------------------------------------------------|
	ROC = ROC3A(RTOR,SHIFT,ABC,ELONG,TRIAN)
	NA = NA1-1
	do	J=1,NB1
	   RHO(J) = (J-0.5d0)*HRO
	enddo
C Define AMETR,SHIF,ELON,TRIA,SHIV
	call	SETGEO(0)
	YDV = 0.d0
	do J=1,NB1
	   if(j.lt.nb1) then
	      j1=j+1
	   else
	      j1=j
	   endif
	   YA =0.5d0*(AMETR(J)+AMETR(J1))
	   YES=0.5d0*(ELON(J)+ELON(J1))
	   YDS=0.5d0*(TRIA(J)+TRIA(J1))
	   YAS=0.5d0*(SHIF(J)+SHIF(J1))
	   VOLUM(J)=GP*GP2*YA**2*YES*(RTOR+YAS-0.25d0*YA*YDS)
	   VR(J) = (VOLUM(J)-YDV)/HRO
	   YDV = VOLUM(J)
	   YR1 = min(AMETR(j)/2.d0,1.d-2)
	   YAS = ROC3A(RTOR,SHIF(j),AMETR(j),    ELON(j),TRIA(j))
	   YDS = ROC3A(RTOR,SHIF(j),AMETR(j)-YR1,ELON(j),TRIA(j))
	   DRODA(J-1) = (YAS-YDS)/YR1
	enddo
	VOLUM(NA1)=GP*GP2*ABC**2*ELONG*(RTOR+SHIFT-0.25*ABC*TRIAN)
	YIPL = (FP(NA1)-FP(NA)-(FV(NA1)-FV(NA)))/HROA
	JNA1O = NA1
C Input: 	ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC,  AMETR(NA1)
	call	NEWGRD
C Output:	NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, HROA, AMETR(j>NA1)
	if (abs(YIPL-IPL) .lt. .2d0*IPL)	then
	   call	EDCELL(JNA1O,YIPL)
	endif
C If NA1 does not change then EDCELL corrects FP(NA1) only.

!====== Definition --------------------------- Approximation ----------|
!						gradRHO = DRODA
!					<|grad(a)|> = sqrt(G1)/DRODA 
! G1=<(gradRHO)**2>;				G1=DRODA**2
! G11=VR*G1=VR*<(gradRHO)**2>			G11=VR*G1
! G2=<(gradRHO/R)**2>*VR/4/pi**2; 		G22=R*G2/J;
! G22=VR*R*<(gradRHO/R)**2>/(GP2)**2/IPOL;	C22=G11/R/(GP2)**2/IPOL
! SLAT=VR*DRODA*<|gradA|>;			SLAT=VR*sqrt(G1)
	do	J=1,NB1
	   IPOL(J) = 1.d0
	   G33(J) = (RTOR/(RTOR+SHIF(J)))**2
	   GRADRO(J)=DRODA(J)
	   if (j .lt. NB1)	VRS(j)=0.5d0*(VR(J+1)+VR(j))
	   G11(J) = VRS(J)*DRODA(J)**2
	   G22(J) = G11(J)/GP2**2/(RTOR+SHIF(J))
	   SLAT(J) = VRS(J)*DRODA(J)		
	enddo
	VRS(NA)=(VOLUM(NA1)-VOLUM(NA))/(RHO(NA1)-RHO(NA1-1))
	if(NA1.lt.NB1)	
     >		VRS(NA1)=(VOLUM(NA1+1)-VOLUM(NA1))/(RHO(NA1+1)-RHO(NA1))
C   BDB02 - <B**2/B0**2>
C   B0DB2 - <B0**2/B**2>			<(R/R0)^2>
C   BMAXT - BMAXT
C   BMINT - BMINT
C   BDB0  - <B/BTOR>
C   FOFB  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
C   SHEAR -  d[ln(q)]/d[ln(rho)]	(to replace fml/shear)
C Alternative definition for BDB02 (G.Pereverzev 10.02.99)
	do	J=1,NA1
	   if     (j .eq. 1)	then
	      SHEAR(J) = (FP(2)-FP(1))/(2.d0*MU(1)+0.333d0*(MU(1)-MU(2)))
	   elseif (j .lt. NA)	then
	      SHEAR(J) = (FP(j+1)-2.d0*FP(j)+FP(j-1))/(MU(j+1)+MU(j))
	   else
	      SHEAR(J) = HRO*(FP(NA1)-FP(NA))/HROA-FP(NA)+FP(NA-1)
	      SHEAR(J) = 2.d0*SHEAR(J)*(2.d0*HROA-HRO)
	      SHEAR(J) = SHEAR(J)/(2.d0*HROA*MU(NA1)-HRO*MU(NA))
	   endif
	   SHEAR(J) = 1.d0 -SHEAR(J)/(GP*BTOR*HRO**2)
	   YR1      = RHO(j)*G22(j)*(MU(j)/RTOR)**2
	   BDB02(j) = (1.d0+YR1)*G33(J)*IPOL(J)**2
	   B0DB2(j) = ((RTOR+SHIF(J))/RTOR)**2+.75d0*(AMETR(j)/RTOR)**2
	   BMAXT(j) = BTOR*RTOR/(RTOR+SHIF(j)-AMETR(j))
	   BMINT(j) = BTOR*RTOR/(RTOR+SHIF(j)+AMETR(j))
	   BDB0(j)  = (RTOR/(RTOR+SHIF(j)))
C  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
	   YR1     = (RTOR+SHIF(j)-AMETR(j))/(RTOR+SHIF(j))
	   FOFB(j) = 1.d0-sqrt(YR1)*(1.d0+0.5d0*YR1)
	enddo
	do	J=NA1,NAB
	   SHEAR(j) = SHEAR(NA1)
	   BDB02(j) = BDB02(NA1)
	   B0DB2(j) = B0DB2(NA1)
	   BMAXT(j) = BMAXT(NA1)
	   BMINT(j) = BMINT(NA1)
	   BDB0(j) = BDB0(NA1)
	   FOFB(j) = FOFB(NA1)
	enddo
C Volume (on the shifted grid) is calculated using VR:
	VOLUM(1) = HRO*VR(1)
	do	J=2,NA
	   VOLUM(J) = VOLUM(J-1)+HRO*VR(J)
	enddo
	VOLUM(NA1) = VOLUM(NA)+(HROA-0.5d0*HRO)*VR(NA)
	VOLUME = VOLUM(NA1)
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	subroutine SETEXM
C----------------------------------------------------------------------|
C Set external metrics 
C This subroutine is adjusted for external metric calculated for
C TJ-II			(Pereverzev 10.02.2005)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	integer	j
	double precision YNF,YR1
C----------------------------------------------------------------------|
C Piece of BLOCKDATA:
C     43->	'MUX   ','MVX   ','GNX   ','SNX   ','PEX   ',
C     48->	'PIX   ','PRADX ','TEX   ','TIX   ','NEX   ',
C     53->	'CUX   ','ZEFX  ','VRX   ','SHX   ','ELX   ',
C     58->	'TRX   ','G11X  ','G22X  ','G33X  ','DRODAX',
C     63->	'IPOLX ','NIX   ','VPOLX ','VTORX '/
	if (EXARNM(44) .ne. 'MVX')	goto	99
	if (EXARNM(55) .ne. 'VRX')	goto	99
	if (EXARNM(56) .ne. 'SHX')	goto	99
	if (EXARNM(57) .ne. 'ELX')	goto	99
	if (EXARNM(58) .ne. 'TRX')	goto	99
	if (EXARNM(59) .ne. 'G11X')	goto	99
	if (EXARNM(60) .ne. 'G22X')	goto	99
	if (EXARNM(61) .ne. 'G33X')	goto	99
	if (EXARNM(62) .ne. 'DRODAX')	goto	99
	if (EXARNM(63) .ne. 'IPOLX')	goto	99
	if (EXARNM(67) .ne. 'SLATX')	goto	99

C	ROC = ABC
	YNF = RTOR*GP2**2
	do	J=1,NA1
	   if (IFDFAX(56) .gt. 0)	then
	      SHIF(J) = SHX(j)
C	      SHIF(J) = SHIFT*SHX(j)
	   else
	      SHIF(J) = SHIFT
	   endif
	   if (IFDFAX(57) .gt. 0)	then
	      ELON(J) = ELX(j)
C	      ELON(J) = ELONG*ELX(j)
	   else
	      ELON(J) = 1.
	   endif
	   if (IFDFAX(58) .gt. 0)	then
	      TRIA(J) = TRX(j)
C	      TRIA(J) = TRIAN*TRX(j)*RHO(j)**2
	   else
	      TRIA(J) = 0.
	   endif

	   if (IFDFAX(61) .gt. 0)	then
	      G33(J) = G33X(j)
	   else
	      G33(J) = (RTOR/(RTOR+SHIFT))**2
	   endif
	   if (IFDFAX(63) .gt. 0)	then
	      IPOL(J) = IPOLX(j)
	   else
	      IPOL(J) = 1.d0
	   endif
	   if (IFDFAX(55) .gt. 0)	then
	      VR(J) = VRX(j)
	   else
	      VR(J) = YNF*RHO(j)/(IPOL(j)*G33(j))
	   endif
	   AMETR(J) = RHO(J)
C          AMETR(J) = AMX(J)	! AMX not defined
	enddo

C Flux grid: j*h
	do	J=1,NA
	   VRS(j)= 0.5d0*(VR(J+1)+VR(j))
	   if (IFDFAX(67) .gt. 0)	then
	      SLAT(J) = SLATX(j)
	   else
	      SLAT(J) = VRS(j) 
	   endif
	   if (IFDFAX(59) .gt. 0)	then
	      G11(J) = VRS(j)*G11X(j)			! TJ-II convention
C	      G11(J) = G11X(j)
	   else
	      G11(J) = VRS(j) 
	   endif
	   if (IFDFAX(60) .gt. 0)	then
	      G22(J) = RTOR*VRS(j)*G22X(j)/(IPOL(j)*GP2*GP2)
	   else
C?	      G22(J) = RTOR*VRS(j)/(GP2*(RTOR+SHIFT))**2
	      G22(J) = RTOR*VRS(j)/(GP2*(RTOR+SHIFT))**2
	   endif
	   if (IFDFAX(62) .gt. 0)	then
	      DRODA(J) = DRODAX(j)
	   else
	      DRODA(J) = 1.d0
	   endif
	enddo
	SLAT(NA1)= 1.5d0*SLAT(NA)-0.5d0*SLAT(NA-1)
	VRS(NA1)= 1.5d0*VRS(NA)-0.5d0*VRS(NA-1)
	G11(NA1)= 1.5d0*G11(NA)-0.5d0*G11(NA-1)
	G22(NA1)= 1.5d0*G22(NA)-0.5d0*G22(NA-1)
	DRODA(NA1)= 1.5d0*DRODA(NA)-0.5d0*DRODA(NA-1)
	G11(NA1)= 1.5d0*G11(NA)-0.5d0*G11(NA-1)

C   BDB02 - <B**2/B0**2> can be expressed in terms of others, see Eq.(39)
C   B0DB2 - <B0**2/B**2> ~ <(R/R0)^2>
C   BMAXT - BMAXT
C   BMINT - BMINT
C   BDB0  - <B/BTOR>
C   FOFB  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
C Alternative definition for BDB02 (G.Pereverzev 10.02.99)
	do	J=1,NA1
	   if     (j .eq. 1)	then
	      SHEAR(J) = (FP(2)-FP(1))/(2.d0*MU(1)+0.333d0*(MU(1)-MU(2)))
	   elseif (j .lt. NA)	then
	      SHEAR(J) = (FP(j+1)-2.d0*FP(j)+FP(j-1))/(MU(j+1)+MU(j))
	   else
	      SHEAR(J) = HRO*(FP(NA1)-FP(NA))/HROA-FP(NA)+FP(NA-1)
	      SHEAR(J) = 2.d0*SHEAR(J)*(2.d0*HROA-HRO)
	      SHEAR(J) = SHEAR(J)/(2.d0*HROA*MU(NA1)-HRO*MU(NA))
	   endif
	   SHEAR(J) = 1.d0 -SHEAR(J)/(GP*BTOR*HRO**2)
	   GRADRO(j) = DRODA(J)
	   YR1	= RHO(j)*G22(j)*(MU(j)/RTOR)**2
	   BDB02(j) = (1.d0+YR1)*G33(J)*IPOL(J)**2
	   B0DB2(j) = ((RTOR+SHIF(J))/RTOR)**2+0.75d0*(AMETR(j)/RTOR)**2
	   BMAXT(j) = BTOR*RTOR/(RTOR+SHIF(j)-AMETR(j))
	   BMINT(j) = BTOR*RTOR/(RTOR+SHIF(j)+AMETR(j))
C  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
	   BDB0(j) = (RTOR/(RTOR+SHIF(j)))
	   YR1     = (RTOR+SHIF(j)-AMETR(j))/(RTOR+SHIF(j))
	   FOFB(j) = 1.d0-sqrt(YR1)*(1.d0+0.5d0*YR1)
	enddo
	if (NA1 .ge. NAB)	return
	do	J=NA1+1,NAB
	   SHIF(J) = SHIFT
	   ELON(J) = 1.d0
	   TRIA(J) = 0.d0
	   G33(J) = (RTOR/(RTOR+SHIFT))**2
	   IPOL(J) = 1.d0
	   VR(J) = YNF*RHO(j)/(IPOL(j)*G33(j))
	   AMETR(J) = RHO(J)
	   VRS(j)= 0.5d0*(VR(J+1)+VR(j))
	   G11(J) = VRS(j)
	   G22(J) = RTOR*VRS(j)/(RTOR+SHIFT)**2
	   DRODA(J) = 1.d0
	   SLAT(J) = VRS(J)*DRODA(J)
	   SHEAR(j) = SHEAR(NA1)
	   BDB02(j) = BDB02(NA1)
	   B0DB2(j) = B0DB2(NA1)
	   BMAXT(j) = BMAXT(NA1)
	   BMINT(j) = BMINT(NA1)
	   BDB0(j) = BDB0(NA1)
	   FOFB(j) = FOFB(NA1)
	enddo
C Volume (on the shifted grid) is calculated using VR:
	VOLUM(1) = HRO*VR(1)
	do	J=2,NA
	   VOLUM(J) = VOLUM(J-1)+HRO*VR(J)
	enddo
	VOLUM(NA1) = VOLUM(NA)+(HROA-0.5d0*HRO)*VR(NA)
	VOLUME = VOLUM(NA1)
	return
 99	continue
	write(*,*)">>> Illegal changes in the Astra BLOCKDATA"
	call	a_stop
	end
C======================================================================|
C ROC3A [m]:	Analytical formula for the "rho_tor" in vacuum
C				 (Pereverzev 24.03.00)
	double precision function ROC3A(YR,YD,YA,YE,YT)
C----------------------------------------------------------------------|
C Input: YR - Major radius	[m]
C	 YD - Shafranov shift	[m]
C	 YA - Minor radius	[m]
C	 YE - Elongation	[d/l]
C	 YT - Triangularity	[d/l]
C Output:
C	 ROC3A - Dimensional toroidal "rho" [m]
C----------------------------------------------------------------------|
	implicit none
	double precision YR,YD,YA,YE,YT,YD1,YT2,YGE,Y1,Y2,Y3,Y4,Y5,Y6
	YT2 = 2.d0*YT
	YGE = YA/(YR+YD)
	if (abs(YGE) .gt. 1.d0)	goto	8
	if (YGE) 9,5,1
 1	continue
	YD1 = 1.d0+YT2*(YT2-2./YGE)
	if (YD1) 3,2,2
 2	continue
	YD1 = sqrt(YD1)
	Y1 = (2.d0-YT2*YGE)*YD1/(1.d0+YD1)
	Y2 = YGE-YT2
	Y3 = sqrt(1.d0-YT2*YGE+YGE*YD1)
	Y4 = sqrt(1.d0-YT2*YGE-YGE*YD1)
	Y5 = sqrt(1.d0+YGE)
	Y6 = sqrt(1.d0-YGE)
	Y1 = (Y1+Y2)/(Y3+(1.d0-YT2)*Y5)-(Y1-Y2)/(Y4+(1.d0+YT2)*Y6)
	Y1 = Y1/sqrt(2.d0*(1.d0-YT2*YGE+Y5*Y6))
	goto	4
 3	continue
	Y5 = sqrt(1.d0+YGE)
	Y6 = sqrt(1.d0-YGE)
	Y1 = (1.d0-YT2)*Y5+(1.d0+YT2)*Y6
	Y1 = (1.d0-Y1/sqrt(2.d0*(1.d0-YT2*YGE+Y5*Y6)))/YT2
 4	continue
	ROC3A = 2.d0*sqrt(YR*YA*YE*Y1)
	return
 5	continue
	ROC3A = 0.d0
	return
 8	write(*,*)" >>> Error >>> ROC3A >>> Illegal input: R+Delta < a"
	write(*,'(1P,16X,2(A,E10.3))')"R+Delta =",YR+YD,",    a =",YA
	call	a_stop
	return
 9	write(*,*)" >>> Error >>> ROC3A: Illegal input"
	write(*,*)YR,YD,YA,YE,YT
	call	a_stop
	end
C======================================================================|
	subroutine	A2ESC
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	integer Nbndmx,Nescmx,Nesc,JNA1O,Fail,J,Fjob,Mpol,kAst,ICALL
	parameter(Nbndmx=121, Nescmx=257)
	double precision
     1		Zrb(Nbndmx),Zzb(Nbndmx),ZRX,ZRBtor,Zrho,Zgrid(Nescmx),
     2		ZPres(Nescmx),ZJpar(Nescmx+1),ZgY(Nescmx),XTR(NRD),
     3		XEQ(NRD),ZP(NRD),ZA(NRD),ZB(NRD),ZC(NRD),ZJ(NRD),CJ,CP,
     4		ALFA,mu0,Y1,Y2,Y3,Y4,YF,YIPL,Zzm,Zrm,ZgFtor,ROC3A
	double precision coulg,ccsp
	save	Zrm,ZRBtor,mu0,ICALL
	data	mu0/1.25663706d0/ ICALL/0/

	if (ICALL.eq.0)	then
	   call	markloc("Initializing ESC"//char(0))
C          Machine='ASDEX'//char(0)
	   j = 1			! Launch ESC in interactive mode
	   j = 0			! Launch ESC in background
Control via "j" is disabled! Environment variable XESC is used
	   call	esc0(j)
	endif
	ICALL = ICALL+1

	call	markloc("Preparing ESC"//char(0))
C Equilibrium calculation start:
C Note: j_{||Z}=j_{||A}*mu0*RTOR/(RTOR+SHIF(1))/IPOL(j)/G33(j)
	CJ = mu0*RTOR/(RTOR+SHIF(1))
	CP = 1.6d-3*mu0
	do	10	J=1,NA1
	   ZP(J) = CP*
     & (PFAST(J)+NE(J)*TE(J)+NI(J)*TI(J)+
     & 0.5d0*NB2EQL*(PBLON(J)+PBPER(J)))
	   if (j .eq. NA1)	goto	10
	   ZJ(J) = CJ*CU(j)/(IPOL(j)*G33(j))
 10	continue
	ZJ(NA1) = ZJ(NA) ! added 23-11-2006, otherwise j(Nesc) is undefined

	if (NA1 .gt. Nescmx-2)	then
	   Nesc = Nescmx-2
	   do	J=1,Nesc
	      XEQ(J) = (J-1.d0)/(Nesc-1.d0)
	   enddo
	   print *, 'A2ESC 576'
	   print *, XEQ(1),XEQ(2),XEQ(3)
	   do	J = 1,NA1
	      XTR(J) = RHO(J)/ROC
	   enddo
	   ALFA = 0.001d0
	   call	SMOOTH(ALFA,NA1,ZP,XTR,Nesc,ZA,XEQ)
	   call	SMOOTH(ALFA,NA1,ZJ,XTR,Nesc,ZB,XEQ)
	   call	SMOOTH(ALFA,NA1,CC,XTR,Nesc,ZC,XEQ)
	   do	J=1,Nesc
	      Zgrid(J) = XEQ(J)*ROC
	      ZPres(J) = ZA(J)
	      ZJpar(J) = ZB(J)
	   enddo
	   call	SMOOTH(ALFA,NA1,FP,XTR,Nesc,ZB,XEQ)
	   do	J=1,Nesc
	      ZgY(J) = ZB(J)
	   enddo
	else
	   Nesc = NA1
	   do	J=1,Nesc
	      Zgrid(J) = RHO(J)
	      ZPres(J) = ZP(J)
	      ZJpar(J) = ZJ(J)
	      ZgY(J) = FP(J)
	      ZC(J) = CC(J)
	   enddo
	endif

	ZgY(Nesc+1)=Ipl*mu0*RTOR/(IPOL(NA1)*G22(NA))
	ZgY(Nesc+2)=CU(NA1)*mu0*RTOR/(RTOR+SHIF(1))/IPOL(NA1)/G33(NA1)
	ZRX = RTOR
	ZRBtor = RTOR*BTOR
	ZgFtor = GP*BTOR*ROC**2

	if (  Nbnd .gt. Nbndmx)	then
	    write(*,*)">>> ESC: Boundary grid out of limits"
	    write(*,*)"         Call ignored."
	    return
	endif
	call BNDRY(Zrb,Zzb)

	if (ICALL .eq. 1)	then
	   kAst = 221
	   do j=1,Nesc
C	      ZgY(j) = ZJpar(j)
	      ZgY(j) = 1.	! 1st call with homogeneous current
	   enddo
	   ZgY(Nesc+1) = IPL
	else
	   kAst = 1228
	endif
	if (ICALL .le. 4)	then
C	   Y1 = exp(-3.d0+ICALL)
	   do j=1,Nesc
	      Zpres(j) = 0.1d0*Zpres(j)
	   enddo
	endif
	if (ICALL .le. 4)	then
C	   Y1 = exp(-3.d0+ICALL)
	   do j=1,Nesc
	      Zpres(j) = 0.1d0*Zpres(j)
	   enddo
	endif
	Y1 = TIME
C	if (ICALL .ge. 21)	then
C	   kAst = 1229			!z'1229'
C	   Y1 = TIME+TAU
C	   write(*,*)ICALL,TIME,kAst
C	   do j=1,Nesc
C	      ZgY(j) = ZC(j)
C	   enddo
C	endif
C----------------------------------------------------------------------|
C ESC Input:
C   double ZPres(1:Nesc) - pressure on the grid Zgrid,
C   double ZgY  (1:Nesc) - poloidal flux on the grid Zgrid,
C   double Zgrid(1:Nesc) - radial grid
C   int    Nesc - number of radial grid points
C   int    kAst - index of input profiles/a now it is:
c       12-$\sqrt(\bgF/\gp B_0)$, 2-$\bsp$, 1-$\bsj_\parallel$
c       i.e.   12 - $\sqrt(\Phi/\pi B_0)$,
c               2 - $\bsp$,
c               1 - $\bsj_\parallel$
C   double Zrb - (Zrb,Zzb) array of boundary points
C   double Zzb
C   int    Nbnd   - number of boundary points
C   double ZRBtor - RTOR*BTOR
C   double ZRX
C	double	*pProf	/* p[] - pressure related profile */
C	double	*jProf	/* j[] - current related profile */
C	double	*aProf	/* a[] - normalized sqrt(Phi) square root 
C				 of toroidal flux.
C				 If a=NULL (dropped in FORTRAN) - uniform 
C				 grid from 0 till 1,
C				 otherwise it is considered as an array of 
C				 grid points */
C	int *nProf /* number <= 129 of profile points including \
C		  	   magnetic axis and plasma edge */
C	int *kAst  /* First digit
C	   	 	n - 0 - normalized radial coordinate (0-1)
C	   	 	    1 - absolute radial coordinate;
C	   	 	r - 0 - b (vertical semiaxis)
C	   	 		    as the radial coordinate,
C	   	 	    1 - V (volume),
C	   	 	    2 - gF (toroidal flux),
C	   	 	    3 - gY (poloidal flux, DO NOT use),
C	   	 	h specifies the meaning of pProf:
C	   	 	    0 - j_p;
C	   	 	    1 - P =dp/dpsi;
C	   	 	    2 -p [MPa].
C	   	 	Second digit H specifies the meaning
C	   	 		of jProf:
C	   	 	    0 - j_s [MA/m^2];
C	   	 	    1- j_|| [MA/m^2]=(j dot B)/(R_0 B grad phi),
C	   	 		                R_0 = magnetic axis;
C	   	 	    2- j_||R_0 [MA/m] = (j dot B)/(B grad phi);
C	   	 	    3 - T=FF';
C	   	 	    6 - q;
C	   	 	    7 - 1/q;
C	   	 	    8 - \gY;
C			    9 - \sigma, conductivity [MS]
C	   	 	For Example:
C	   	 		26 - p[] and q[] profiles are supplied;
C	   	 		21 - p[] and j||[] profiles are supplied;
C	   	 		0 -  jp[] and js[] profiles are supplied.
C
C	   	 		Possible combinations are limited to:
C	   	 		0,1,2,6,7
C	   	 		10,11,12,16,17
C	       	 	21,22,26,27   */
C	double	*Rpv	/* R[m]- plasma boundary */
C	double	*Zpv	/* Z[m]- plasma boundary */
C	int	*Npv	/* number <= 257 of the plasma-vacuum points.
C			   If *Npv >12, the first and last points coincide */
C	double	*RBtor	/* RBtor [m Tesla] outside the plasma */
C	double	*ZRX	/* Reference major radius [m] */
C----------------------------------------------------------------------|
	call	markloc("Filling ShMem"//char(0))
C	write(*,*)NESC,kAst,Nbnd,ZRBtor,ZRX
C	write(*,*)(Zpres(j),j=1,Nesc)
C	write(*,*)(ZgY(j),j=1,Nesc+2)
C,ZgY,Zgrid,Nesc,kAst,Zrb,Zzb,Nbnd,ZRBtor,ZRX)
	call	escin(Zpres,ZgY,Zgrid,Nesc,kAst,Zrb,Zzb,Nbnd,ZRBtor,ZRX)
	Fail = 0
	Fjob = 0
	Mpol = 4
	call	markloc("Calling ESC"//char(0))
	call	esc(Fail,Fjob,Mpol,Y1)
	CPTEQL = Y1

C	call	markloc("Post-ESC processing 1"//char(0))
C	call	esigetmagfluxes(Y2,Y1,1.d0)
C	ROC = sqrt(Y2/GP/BTOR)
C	write(*,'(2(A,F10.6))')'Y2 =',Y2,',    Ftor =',ZgFtor
	call	markloc("Calling ESCVOUT"//char(0))
	call	escvout(Zrm,Zzm,ZgFtor)	! (ZRm,ZZm) - the magnetic axis
	ROC = sqrt(ZgFtor/GP/BTOR)	! ZgFtor - \Phi at the plasma edge
	call	add2loc("Calling get3ma"//char(0))
	call	get3ma (1.d0,ABC,Y1,Y2,Y3,Y4)
C	write(*,'(3(A,F10.6))')'ROC =',ROC,',    t =',TIME
C	write(*,'(3(A,F10.6))')'AB =',AB,'   ABC =',ABC,'   h_a =',HROA
	if (ABC .gt. AB)	then
	   write(*,*)
	   write(*,*)">>> Inconsistency in boundary settings:"
	   write(*,'(2(A,F6.3))')
     >		"     ESC determines ABC = ",ABC,"  >  AB =",AB
	   write(*,*)
     >		"    Please increase parameter AB and restart ASTRA"
	   write(*,*)
	   call	a_stop
	endif

	call	add2loc("Calling NEWGRD"//char(0))
	JNA1O = NA1		! save current NA1 for later use in EDCELL 
	YIPL = HROA
C input: 	ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC
C Output:	NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, HROA, AMETR(j>=NA1)
	call	NEWGRD
C	write(*,'(3(A,I3))')'JNA1O -> NA1:',JNA1O,' ->',NA1,',  NAB=',NAB
C	write(*,'(3(A,F8.6))')'h_a -> h_a:  ',HROA/HRO,' -> ',YIPL/HRO

	if (ICALL .eq. 1)	then
	   call	add2loc("Calling esiPlCurrent"//char(0))
	   call	esiplcurrent(YIPL)
	   call	add2loc("Calling esiGetMagFluxes"//char(0))
	   do	j=1,NA1
	      Y1 = RHO(j)/ROC
	      call	esigetmagfluxes(Y2,YF,Y1)	! Y2 = \Phi(j)
	      FP(j)=-IPL*YF/YIPL
	      FP(j)=-YF
	   enddo
	   call	add2loc("Calling CUOFP"//char(0))
	   call	CUOFP
	   Y1 = GP2*RTOR
	   do	J = 1,NA1
	      include 'fml/ccsp'
	      CC(j) = CCSP
	      if (CC(J).gt.0.0) then
		 ULON(J)=Y1*CU(J)/CC(J)
	      else
		 ULON(J)=0.d0
	      endif
	      UPL(J)   =ULON(J)/(IPOL(J)*G33(J))
	   enddo
	   ULON(NA1) = ULON(NA)
	   UPL(NA1)  = ULON(NA1)/(IPOL(NA1)*G33(NA1))
	else		! save current dPsi/drhoIPL for use in EDCELL 
	   YIPL = (FP(NA1)-FP(NA)-(FV(NA1)-FV(NA)))/YIPL
	endif
C	write(*,'(3(A,F8.3))')'I_pl -> I_pl:  ',IPL,' ->',
C     ,		YIPL*IPOL(NA1)*G22(NA)/mu0/RTOR

	call	markloc("Post-ESC processing 2"//char(0))
	call	add2loc("Calling EDCELL"//char(0))
	call	EDCELL(JNA1O,YIPL)
	call	markloc("Post-ESC processing 3"//char(0))

C Define main transport (half-integer) grid:
	call	add2loc("Calling ESC3Mout"//char(0))
	do  j = 1,NA1
	   Zrho = RHO(j)/ROC
	   call	get3ma (Zrho,AMETR(j),SHIF(j),ELON(j),TRIA(j),SHIV(j))
	   call escmout(Zrho,VR(j),G33(j),IPOL(j),SHEAR(j))
	   VR(j) = VR(j)/ROC
	   G33(j) = G33(j)*RTOR*RTOR
	   IPOL(j) = IPOL(j)/RTOR/BTOR
	   SHIF(j) = SHIF(j)-RTOR
	enddo
	if (NA1 .eq. NAB)	goto	12
	do  j = NA1+1,NAB
	   SHIF(j) = SHIF(NA1)
	   ELON(j) = ELON(NA1)
	   TRIA(j) = TRIA(NA1)
	   AMETR(j) = (AB*(j-NA1)+ABC*(NAB-j))/(NAB-NA1)
	enddo
 12	continue

C Define auxiliary shifted (integer) grid:
	Y4 = HRO/ROC
	Zrho = 0.d0
	call	add2loc("Calling ESCaout"//char(0))
	do  j = 1,NA1
	   Zrho = Zrho+Y4
	   if (Zrho .gt. 1.d0)	then
	      Zrho = 1.d0
	      if (j .le. NA)	write(*,*)"Out of grid"
	   endif
	   if (j .eq. NA1)	then
	      Zrho = 1.d0
	      DRODA(j) = HROA/(AMETR(NA1)-AMETR(NA))	! d\rho/da
	   elseif  (j .eq. 1)	then
	      DRODA(j) = HRO/(AMETR(2)-AMETR(1))	! d\rho/da
	   else
	      DRODA(j) = HRO/(AMETR(j)-AMETR(j-1))	! d\rho/da
	   endif
	   call escaout(Zrho,VRS(j),SLAT(j),Y1,Y2,Y3,	! Y3 = 1/q
     &		BMINT(j),BMAXT(j),BDB0(j),BDB02(j),B0DB2(j))
	   VRS(j) = VRS(j)/ROC
	   G11(j) = Y1*ROC
	   G22(j) = Y2*RTOR*RTOR*BTOR*ROC
	   GRADRO(j) = SLAT(j)/VRS(j)			! <|\nabla\rho|>
	   BDB0(j)  = BDB0(j)/BTOR			! <B/B0>
	   BDB02(j) = BDB02(j)/(BTOR*BTOR)		! <(B/B0)^2>
	   B0DB2(j) = B0DB2(j)*BTOR*BTOR		! <(B0/B)^2>
	enddo
C Volume (on the shifted grid) is calculated using VR:
	VOLUM(1) = HRO*VR(1)
	do	J=2,NA
	   VOLUM(J) = VOLUM(J-1)+HRO*VR(J)
	enddo
	VOLUM(NA1) = VOLUM(NA)+(HROA-0.5d0*HRO)*VR(NA)
	VOLUME = VOLUM(NA1)

	do	J=NA1,NAB
	   SHEAR(j) = SHEAR(NA1)
	   BDB02(j) = BDB02(NA1)
	   B0DB2(j) = B0DB2(NA1)
	   BMAXT(j) = BMAXT(NA1)
	   BMINT(j) = BMINT(NA1)
	   BDB0(j) = BDB0(NA1)
	   FOFB(j) = FOFB(NA1)
	enddo
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	subroutine	BNDRY(RPB,ZPB)
C----------------------------------------------------------------------|
C If 3M solver is used the subroutine is not called.
C Otherwise, if a general equilibrium solver, ESC or SPIDER, is called
C then 
C 1) In case of the plasma boundary defined by 3 moments,
C    this subroutine writes 8 points on the boundary into array BNDARR
C    and into arrays RPB(1:NBND), ZPB(1:NBND)
C 2) If the plasma boundary is defined by a data file then
C    this subroutine uses the array BNDARR as an input and
C    produces output in [time dependent] arrays RPB, ZPB
C----------------------------------------------------------------------|
C NBND     number of points on the plasma vacuum boundary
C NBNT     number of times for the plasma boundary evolution
C  call from ESC:
C		call	BNDRY(RPB,ZPB)
C  call from SPIDER:
C		call	BNDRY(RZPB,RZPB(NBND+1))
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	integer	j,j1,jt
	double precision	ydt,yd1,yd2
	double precision	RPB(*),ZPB(*)
	if (NBNT .gt. 1)	goto	2
	if (NBNT .eq. 0)	then	! No boundary points given
	   NBND = 8
	   yd1 = .75d0			! sin^2(pi/3)
	   yd2 = .5d0		! cos(pi/3)	
	   ydt = sqrt(yd1)		! sin(pi/3)
	   BNDARR(1) = RTOR+SHIFT-ABC*TRIAN
	   BNDARR(2) = UPDWN+ABC*ELONG
	   BNDARR(3) = RTOR+SHIFT-ABC*TRIAN
	   BNDARR(4) = UPDWN-ABC*ELONG
	   BNDARR(5) = RTOR+SHIFT-ABC
	   BNDARR(6) = UPDWN
	   BNDARR(7) = RTOR+SHIFT+ABC
	   BNDARR(8) = UPDWN
	   BNDARR(9) = RTOR+SHIFT-ABC*(TRIAN*yd1+yd2)
	   BNDARR(10) = UPDWN+ABC*ELONG*ydt
	   BNDARR(11) = RTOR+SHIFT-ABC*(TRIAN*yd1-yd2)
	   BNDARR(12) = UPDWN+ABC*ELONG*ydt
	   BNDARR(13) = RTOR+SHIFT-ABC*(TRIAN*yd1-yd2)
	   BNDARR(14) = UPDWN-ABC*ELONG*ydt
	   BNDARR(15) = RTOR+SHIFT-ABC*(TRIAN*yd1+yd2)
	   BNDARR(16) = UPDWN-ABC*ELONG*ydt
	endif
 1	jt = max(1,NBNT)
	do	j=1,NBND
	   j1 = NBNT+(2*j-1)*jt
	   RPB(j) = BNDARR(j1)
	   ZPB(j) = BNDARR(j1+jt)
	enddo
	return
 2	continue
C Input order:
C	t_1		t_2		t_3
C	r_1(t_1)	r_1(t_2)	r_1(t_3)	
C	z_1(t_1)	z_1(t_2)	z_1(t_3)	
C	r_2(t_1)	r_2(t_2)	r_2(t_3)	
C	z_2(t_1)	z_2(t_2)	z_2(t_3)	
	if (TIME .le. BNDARR(1))	goto	1
	if (TIME .ge. BNDARR(NBNT))	goto	1

	do	j=1,NBNT
	   if (TIME .gt. BNDARR(j)) jt = j
	enddo
	if (jt .eq. NBNT)	write(*,*)"OGOGO"
	ydt = BNDARR(jt+1)-BNDARR(jt)
	yd1 = (TIME-BNDARR(jt))/ydt
	yd2 = (TIME-BNDARR(jt+1))/ydt
	do	j=1,NBND
	   j1 = jt+(2*j-1)*NBNT
	   RPB(j) = yd1*BNDARR(j1+1)-yd2*BNDARR(j1)
	   j1 = jt+2*j*NBNT
	   ZPB(j) = yd1*BNDARR(j1+1)-yd2*BNDARR(j1)
	enddo
	if (NBND .gt. 12)	return
C The order is essential
C  Top(1), bottom(2), inward(3), outward(4), ...
	YDT = -1.d0
	do	j=1,NBND
	   if (ZPB(j) .gt. YDT)	then
	      YDT = ZPB(j)
	      j1 = j
	   endif
	enddo
	YD1 = RPB(1)
	YD2 = ZPB(1)
	RPB(1) = RPB(j1)
	ZPB(1) = ZPB(j1)
	RPB(j1) = YD1
	ZPB(j1) = YD2
C Bottom(2)
	YDT = 1.d0
	do	j=2,NBND
	   if (ZPB(j) .lt. YDT)	then
	      YDT = ZPB(j)
	      j1 = j
	   endif
	enddo
	YD1 = RPB(2)
	YD2 = ZPB(2)
	RPB(2) = RPB(j1)
	ZPB(2) = ZPB(j1)
	RPB(j1) = YD1
	ZPB(j1) = YD2
C Inward(3)
	YDT = 1000.d0
	do	j=3,NBND
	   if (RPB(j) .lt. YDT)	then
	      YDT = RPB(j)
	      j1 = j
	   endif
	enddo
	YD1 = RPB(3)
	YD2 = ZPB(3)
	RPB(3) = RPB(j1)
	ZPB(3) = ZPB(j1)
	RPB(j1) = YD1
	ZPB(j1) = YD2
C Outward(4)
	YDT = -1.d0
	do	j=4,NBND
	   if (RPB(j) .gt. YDT)	then
	      YDT = RPB(j)
	      j1 = j
	   endif
	enddo
	YD1 = RPB(4)
	YD2 = ZPB(4)
	RPB(4) = RPB(j1)
	ZPB(4) = ZPB(j1)
	RPB(j1) = YD1
	ZPB(j1) = YD2

C Alternative boundary setting:
C	RPB(1) = ZRD17X
C	RPB(2) = ZRD19X
C	RPB(3) = ZRD21X
C	RPB(4) = ZRD23X
C	RPB(5) = ZRD25X
C	RPB(6) = ZRD27X
C	RPB(7) = ZRD29X
C	RPB(8) = ZRD31X
C	ZPB(9) = ZRD18X
C	ZPB(10) = ZRD20X
C	ZPB(11) = ZRD22X
C	ZPB(12) = ZRD24X
C	ZPB(13) = ZRD26X
C	ZPB(14) = ZRD28X
C	ZPB(15) = ZRD30X
C	ZPB(16) = ZRD32X
	end
C======================================================================|
	subroutine	RHSEQ
C-----------------------------------------------------------------------
C Input: RTOR,BTOR,NA,NA1,HRO,HROA,   NB2EQL,
C	 NE,NI,TE,TI,MU,CU,AMETR,RHO, PBLON,PBPER,G22,G33,IPOL
C Output:
C       EQPF
C       EQFF
C       CUTOR
C Both quantities EQPF (~p') and EQFF (~II') are given in [MA/m^2]
C	EQPF = -1.6E-3*(2*\pi*R_0)\prti{n_13*T_keV}{\psi[Vs=T*m^2]}
C	     = -1.E-6*/(2*\pi*R_0)\prti{p[J/m^3=Pascal]}{\psi[Vs]}
C	EQFF = -1.E-6*2*\pi/(R_0*\mu_0)*I*\prti{I}{\psi}
C	     = -5./R_0*I*\prti{I}{\psi}
C Local toroidal current density j[MA/m^2] is EQPF*r/R_0+EQFF*R_0/r, i.e.
C       j(r,z) = r*(\vec j\cdot\nabla\zeta) = EQPF*r/R_0+EQFF*R_0/r ,
C ASTRA average toroidal current density is
C	  R_0*<\vec j\cdot\nabla\zeta> = EQPF+EQFF*<R_0^2/r^2>
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	j
	double precision	YCB,YG,YTH2
C Preparing input for the 3M equilibrium solver:
	YCB = 1.6E-3*RTOR/(BTOR*HRO*HRO)
	do	J=2,NA
	   if (j.eq.NA) YCB = YCB*HRO/HROA
C	   EQFF(J) = -YCB/(MU(J))*(
C     +		.5d0*(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))/(J+1)
C     +	        +( (NE(J+1)*TE(J+1)-NE(J)*TE(J))+
C     +		   (NI(J+1)*TI(J+1)-NI(J)*TI(J)) )/J )
	   EQFF(J) = ( (NE(J+1)*TE(J+1)-NE(J)*TE(J))+
     +		       (NI(J+1)*TI(J+1)-NI(J)*TI(J)) )/J
	   EQFF(J) = EQFF(J)+0.5d0*NB2EQL*
     +			(PBLON(J+1)-PBLON(J)+PBPER(J+1)-PBPER(J))/J
	   EQFF(J) = EQFF(J)+(PFAST(J+1)-PFAST(J))/J
	   EQFF(J) = -YCB*EQFF(J)/(MU(J))
	enddo
	EQFF(1)	=EQFF(2)
	EQFF(NA1)=EQFF(NA)+(EQFF(NA)-EQFF(NA-1))
     .		*(AMETR(NA)-AMETR(NA-1))/(AMETR(NA1)-AMETR(NA))
	EQFF(NA1)=EQFF(NA)+(EQFF(NA)-EQFF(NA-1))/HRO*HROA
	do	J=1,NA1
	   EQPF(j) = EQFF(j)
	   YTH2	= RHO(j)*G22(J)*(MU(J)/RTOR)**2
	   YG	= (1.d0+YTH2)*G33(J)
	   EQFF(J) = (CU(J)/IPOL(J)-EQPF(J))/YG
C	   if (J.eq.NA1) EQFF(J) = -EQPF(J) /YG
	   CUTOR(J) = (CU(J)/IPOL(J)+YTH2*EQPF(J))/(1.d0+YTH2)
	enddo
	return
C	if(TIME.gt.0.3)	then
C	write(*,*)
C	write(*,*)"From RHS"
C	write(*,*)"PBLON"
C	write(*,100)(PBPER(j),j=1,6)
C	write(*,*)"PBPER"
C	write(*,100)(NE(j),j=NA1-5,NA1)
C	write(*,*)"TE"
C	write(*,100)(TE(j),j=NA1-5,NA1)
C	write(*,*)"TI"
C	write(*,100)(TI(j),j=NA1-5,NA1)
C	write(*,*)"dPe/dPsi"
C	write(*,100)(YCB/(MU(J)*J)*(NE(j+1)*TE(j+1)-NE(j)*TE(j)),j=NA1-5,NA1)
C	write(*,*)"dPi/dPsi"
C	write(*,100)(YCB/(MU(J)*J)*(NI(j+1)*TI(j+1)-NI(j)*TI(j)),j=NA1-5,NA1)
C	write(*,*)"dPlong/dPsi"
C	write(*,100)(0.5d0*YCB/(MU(J)*(J+0.5))*(PBLON(j+1)-PBLON(j)),j=1,6)
C	write(*,*)"dPperp/dPsi"
C	write(*,100)(0.5d0*YCB/(MU(J)*(J+0.5))*(PBPER(j+1)-PBPER(j)),j=1,6)
C	write(*,*)"CU"
C	write(*,100)(CU(j),j=NA1-5,NA1)
C	write(*,*)"MU"
C	write(*,100)(MU(j),j=NA1-5,NA1)
C	write(*,*)"IPOL"
C	write(*,100)(IPOL(j),j=NA1-5,NA1)
C	write(*,*)"RHO"
C	write(*,100)(RHO(j),j=NA1-5,NA1)
C	write(*,*)"G22"
C	write(*,100)(G22(j),j=NA1-5,NA1)
C	write(*,*)"G33"
C	write(*,100)(G33(j),j=NA1-5,NA1)
C	write(*,*)"CUBS"
C	write(*,100)(CUBS(j),j=NA1-5,NA1)
C	write(*,*)"CD"
C	write(*,100)(CD(j),j=NA1-5,NA1)
C	write(*,*)"EQPF"
C	write(*,100)(EQPF(j),j=NA1-5,NA1)
C	write(*,*)"EQFF"
C	write(*,100)(EQFF(j),j=NA1-5,NA1)
C	EQFF(NA1)=EQFF(NA)
C	endif
 101	format(1P,5E13.3)
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	subroutine	CUOFP
C----------------------------------------------------------------------|
C Compute CU(rho) and MU(rho) from FP(rho)
C----------------------------------------------------------------------|
C	Input:	HRO	- radial step (m)
C		HROA	- edge radial step (m)
C		NA1	- number of grid points
C		G22(1:NA)	- <g22/g>*............
C		G33(1:NA)	- 
C		IPOL(1:NA)	- 
C		CD(1:NA)	- external (+bootstrap) current
C		FP(1:NA1)	- poloidal flux
C	Output:	CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
C		MU(1:NA1) - (1/rho)dF/d(rho)	    rotational transform
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	integer	j
	double precision	HH,YAJ,YCJ,ARRNA1
	HH = HRO*HRO
	YAJ = 0.d0
C	write(*,*)"CU"
C	write(*,100)(CU(j),j=NA1-5,NA1)
	do	1	J=1,NA
	   YCJ = YAJ
	   if (j .lt. NA)	then
	      YAJ = (FP(j+1)-FP(j))/HH
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
	      CU(j) = (YAJ-YCJ)/HRO
	   else
	      YAJ = (FP(j+1)-FP(j))/HRO/HROA
	      MU(j) = YAJ/j
	      YAJ = G22(j)*YAJ
C30-11	      CU(j) = 2.d0*(YAJ-YCJ)/(HRO+HROA)
	      CU(j) = (YAJ-YCJ)/HRO
	   endif
	   CU(j) = CU(j)/(j-0.5d0)
 1	continue
C	CU(NA1) = ARRNA1(CU(NA),HROA/HRO)
	MU(NA1) = ARRNA1(MU(NA),HROA/HRO)		! See DEFARR
C	MU(NA1) = MU(NA)*NA/(NA-0.5+HROA/HRO)
	YCJ = 1.25d0/(GP*GP*RTOR)
	YAJ = 0.5d0/(GP*BTOR)
	do	2	J=1,NA1
	   CU(j) = YCJ*CU(j)*G33(J)*IPOL(J)**3
           MU(J) = YAJ*MU(j)
 2	continue
C	CU(NA)=CUBS(NA)+CD(NA)+CC(NA)*ULON(NA)/(RTOR*GP2)
C	CU(NA1)=CUBS(NA1)+CD(NA1)+CC(NA1)*ULON(NA1)/(RTOR*GP2)
C	write(*,100)Time
C	write(*,100)(CU(j),j=NA1-5,NA1)
C	write(*,100)(MU(j),j=NA1-5,NA1)
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	subroutine	CUOFMU
C----------------------------------------------------------------------|
C Compute CU(rho) and FP(rho) from MU(rho)
C----------------------------------------------------------------------|
C	Input:	RHO - radial grid (m)
C		NA1 - number of grid points
C		GP2 - 2\pi
C		RTOR - R_0 [m]
C		BTOR - B_0 [T]
C		G22(1:NA1-1)	- <g22/g>*............
C		G33(1:NA1-1)	- 
C		IPOL(1:NA1-1)	- 
C		CD(1:NA1-1)	- external (+bootstrap) current
C		MU(1:NA1) - (1/rho)dF/d(rho)	    rotational transform
C	Output:	CU(1:NA1) - (1/rho)d{K*dF/d(rho)}/d(rho) current density
C		FP(1:NA1) - poloidal flux [Vs]
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	integer	j
	double precision YH,YM,YM1,YM2,YC,YF
	YC = .2d0*GP2*RTOR/BTOR
	YH = RHO(2)-RHO(1)
	YF = GP2*YH*YH*BTOR
	YM = 0.d0
	do	J=1,NA1
	   YM1 = MU(j)*j
	   YM2 = YM1*G22(j)
	   if (j .lt. NA)	then
	      CU(j) = (YM2-YM)*G33(j)*IPOL(j)**3/YC/RHO(j)
	      FP(j+1) = FP(j)+YF*YM1
	   elseif (j .eq. NA)	then
	      CU(j) = (YM2-YM)*G33(j)*IPOL(j)**3/YC/RHO(j)
	      FP(j+1) = FP(j)+GP2*BTOR*YM1*YH*(RHO(NA1)-RHO(NA))
	   else
! Be careful because of instability in the loop: j -> FF' -> G22 -> j ->
	      CU(NA1) = CU(NA)
	   endif
	   YM=YM2
	enddo
C	CU(NA)=CUBS(NA)+CD(NA)+CC(NA)*ULON(NA)/(RTOR*GP2)
C	CU(NA1)=CUBS(NA1)+CD(NA1)+CC(NA1)*ULON(NA1)/(RTOR*GP2)
C	write(*,100)Time
C	write(*,100)(CU(j),j=NA1-5,NA1)
C	write(*,100)(MU(j),j=NA1-5,NA1)
 100	format(3(2F12.6,2X))
	end
C======================================================================|
	subroutine	EQUIL3
C The parameter NEQL =< NP  (see also the file "for/emeq.inc")
C and		NEQL =< NRD (see the file "for/parameter.inc")
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	NEQL,j,jp,jt,jc,jna1o,jcall
	double precision
     1		ACEQLB,ALFA,Y1,YM,YY,YIPL,YDA,GPP4,YRO,YCB,TRIABC,IINT,
     2		BA(NRD),BB(NRD),GR(NRD),GBD(NRD),GL(NRD),GSD(NRD),
     3		A(NRD),B(NRD),C(NRD),D(NRD),BC(NRD),BD(NRD),XTR(NRD),
     4		XEQ(NRD),B2B0EQ(NRD),B0B2EQ(NRD),BMAXEQ(NRD),BMINEQ(NRD)
     5		,BMODEQ(NRD),FOFBEQ(NRD),GRDAEQ(NRD),betrr,lintr,ARRNA1
	character*80	STRI
	external IINT
	save	jcall,ACEQLB,NEQL
	data	jcall/0/ ACEQLB/1.d-3/
C----------------------------------------------------------------------|
C Prepare input data for the 3M equilibrium solver:
	call	markloc("Calling EQUIL3"//char(0))
C	write(*,'(2(A,I3,5F8.3))')"Call No.",jcall,NEQUIL
	if (jcall .eq. 0)	then
	   NEQL = NEQUIL+1.d-3
	   open(2,file='for/emeq.inc')
 1	   continue
	   read(2,fmt='(A79)',end=2)STRI
	   j = index(STRI,"NP=")
	   if (j .eq. 0)	goto	1
	   jt = index(STRI(j+3:),")")-2
	   read(STRI(j+3:j+3+jt),*)j
	   close(2)
	   if (NEQL .gt. j)	write(*,'(A,I4/A,I4)')
     .	      " >>> Warning >>> Maximum size of the equilibrium grid is"
     .	      ,j,"                 The grid will be reduced to NP =",j
	   if (NEQL .gt. NA1)	write(*,'(1X,A,I4)')
     .">>> Warning >>> 3M Equilibrium grid will be reduced to NA1 =",NA1
	   NEQL = min(NA1,j)
	endif
 2	continue			! Now NEQL is set to min(NA1,NP)
	jt = 5
	TRIABC = TRIAN*ABC
	if (jcall .le. jt)	then
	   TRIABC = jcall*TRIAN*ABC/jt
C	   write(*,'(7F10.3)')TRIABC,TRIAN*ABC
	endif

	jp = 10
	if (jcall .le. jp)	then
	   do	J=1,NA1
	      A(J) = TE(J)
	      B(J) = TI(J)
	      C(J) = CU(J)
	      TE(J) = jcall*TE(J)/jp
	      TI(J) = jcall*TI(J)/jp
	      CU(J) = CU(J)+(jp-jcall)*(1-(RHO(J)/ROC)**2)/jp
	   enddo
	   Y1 = BETRR(ROC)
	   YY = LINTR(ROC)
	   YCB = IINT(CU,ROC)
	   do	J=1,NA1
	      CU(J) = CU(J)*IPL/YCB
	   enddo
	   YM = CU(1)
	   call	RHSEQ
	   do	J=1,NA1
	      TE(J) = A(J)
	      TI(J) = B(J)
	      CU(J) = C(J)
	   enddo
C	   write(*,'(I4,7F10.3)')jcall,		!(TE(j),j=1,1),
C     >		Y1,BETRR(ROC),YY,LINTR(ROC),SHIF(1)/ABC,YM,CU(1)
	   jcall = jcall+1
	endif

	YCB	=RTOR/(RTOR+SHIFT)
	do	J=1,NA1
	   XTR(J) = AMETR(J)/AMETR(NA1)
	   B(J)	= EQPF(J)/YCB
	   A(J)	= B(J)+YCB*EQFF(J)
C	   CAR5(J) = A(J)
C	   CAR6(J) = B(J)
	enddo
	do	J=1,NEQL
		XEQ(J) = (j-1.d0)/(NEQL-1.d0)
	enddo
	print *, 'EQUIL3 1310'
	print *, XEQ(1),XEQ(2),XEQ(3)
C From transport grid in "a" to equidistant grid in "a"
	ALFA	=.001d0
	call	SMOOTH(ALFA,NA1,A,XTR,NEQL,BA,XEQ)
	call	SMOOTH(ALFA,NA1,B,XTR,NEQL,BB,XEQ)

	if (TIME .gt. 1.d10)	then
C	if (jcall .gt. 1)	then
C	   open(1,file="m.out")
C	   write(1,*)NA1,TIME,NEQL,"    EQPF  EQFF  CU  CUTOR  A  B"
C	   write(1,'(3(2F10.5,2X))')RTOR+SHIFT,ABC,ELONG,TRIABC
C	   write(1,'(3(2F10.5,2X))')BTOR*RTOR/(RTOR+SHIFT),IPL
C	   do j=1,na1
C	      write(1,'(I4,7F10.5))')j,A(j),B(j)
C	   enddo
C	   close(1)
C	do j=1,neql
C	write(*,'(I4,7F10.5))')j,BA(j),BB(j)
C	enddo
C	write(*,'(3(2F10.5,2X))')(EQPF(j),j=NA1-5,NA1)
C	write(*,'(3(2F10.5,2X))')(EQFF(j),j=NA1-5,NA1)
C	write(*,'(3(2F10.5,2X))')(CU(j)/IPOL(j),j=NA1-5,NA1)
C	write(*,'(3(2F10.5,2X))')(CUTOR(j),j=NA1-5,NA1)
C	write(*,'(3(2F10.5,2X))')(2.d0*(CU(j)/IPOL(j)-CUTOR(j))
C     >	/(CU(j)/IPOL(j)+CUTOR(j)),j=1,NA1)
C	write(*,'(3(2F10.5,2X))')((EQPF(j)-CUTOR(j))
C     >			*RHO(j)*G22(J)*(MU(J)/RTOR)**2,j=NA1-5,NA1)
C	write(*,'(3(2F10.5,2X))')(BA(j),j=1,NEQL)
C	write(*,'(3(2F10.5,2X))')(BB(j),j=1,NEQL)
	endif
C	write(*,*)
C	write(*,*)IPL,
C     +		(FP(NA1)-FP(NA)-(FV(NA1)-FV(NA)))/(1.25663706d0*HROA)*
C     +		G22(NA)*IPOL(NA1)/RTOR,ROC
	call	markloc(" 3-moment solver"//char(0))
C	call	add2loc("User drawing mode"//char(0))
        call E3ASTR(BA,BB,RTOR+SHIFT,ABC,ELONG,TRIABC,NEQL,ACEQLB,
     .  GR,GBD,GL,GSD,A,BD,B,BA,BB,BC,C,D,BTOR*RTOR/(RTOR+SHIFT),IPL,
     .  B2B0EQ,B0B2EQ,BMAXEQ,BMINEQ,BMODEQ,FOFBEQ,GRDAEQ,TIME)
C-------------------------------------------------------------------------
C	Output:	GR  - rho
C		GBD - shift
C		GL  - elongation
C		GSD - triangularity	= \delta^{Zakh}	= a*\delta^{Astra}
C		A   - <g^{11}>		= <[nabla(a)]^2>
C		BD  - <sqrt[g^{11}]>	= <|nabla(a)|> lateral surface
C		B   - <g^{11}g^{33}>	= <[nabla(a)/r]^2>
C		BA  - <g^{33}>		= <1/r^2>=G33/RTOR**2
C		BB  - IPOL*RTOR*BTOR	= I
C		BC  - d\rho/da
C		C   - (dV/da)/(4\pi^2)
C		D   - V(a)
C               B2B0EQ - <B**2/B0**2>
C               B0B2EQ - <B0**2/B**2>
C               BMAXEQ - BMAXT
C               BMINEQ - BMINT
C               BMODEQ - <B/BTOR>
C               FOFBEQ - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
C               GRDAEQ - <grad a>
C-------------------------------------------------------------------------
	call	markloc("3-moment solver"//char(0))
	YRO = sqrt(RTOR/(RTOR+SHIFT))
C	if(TIME.gt.0.38)	then
C	write(*,*)ROC,YRO*GR(NEQL)
C	write(*,'(3(2F10.5,2X))')BTOR*RTOR/(RTOR+SHIFT),IPL,GR(NEQL)
C	write(*,'(3(2F10.5,2X))')(BB(j),j=1,NEQL)
C	write(*,'(3(2F10.5,2X))')(GBD(j),j=1,NEQL,NEQL-1)
C	write(*,'(3(2F10.5,2X))')(GR(j),j=1,NEQL)
C	write(*,'(3(2F10.5,2X))')(BC(j),j=1,NEQL)
C	write(*,'(3(2F10.5,2X))')(GL(j),j=1,NEQL)
C	write(*,'(3(2F10.5,2X))')(GSD(j),j=1,NEQL)
C	call	a_stop
C	endif
C--------------  NEQL <= 1 can be returned by EQAB3 via E3ASTR
	if (NEQL .le. 10)  then
	   NEQUIL = 1
	   write(*,*)"Time =",TIME,",  Continue without MES"
	   return
	endif

C Define a new "rho" grid:
C	write(*,*)"ROC =",ROC,GR(NEQL),NEQL
	ROC = YRO*GR(NEQL)

	YIPL = (FP(NA1)-FP(NA)-(FV(NA1)-FV(NA)))/HROA
C	write(*,*)IPL,(FP(NA1)-FP(NA))/HROA*G22(NA)*IPOL(NA1)/.4/GP/RTOR
C	write(*,*)NA1,HRO,HROA
	JNA1O = NA1

	call	NEWGRD
C	write(*,*)IPL,
C     +		(FP(NA1)-FP(NA)-(FV(NA1)-FV(NA)))/(1.25663706d0*HROA)*
C     +		G22(NA)*IPOL(NA1)/RTOR,ROC
	call	EDCELL(JNA1O,YIPL)
C	write(*,*)NA1,HRO,HROA
C	write(*,*)IPL,(FP(NA1)-FP(NA))/HROA*G22(NA)*IPOL(NA1)/.4/GP/RTOR

C Define quantities attributed to the auxiliary grid:
	GPP4 = GP2*GP2
C A -> <[nabla(rho)]^2>;     B -> <[nabla(rho)/r]^2>;    C -> dV/d(rho)
C 			    BD -> <|nabla(rho)|>
	do	J=1,NEQL
	   BC(J) = YRO*BC(J)
	   A(J) = A(J)*BC(J)**2
	   B(J) = B(J)*BC(J)**2
	   C(J) = GPP4*C(J)/BC(J)
	   BA(J) = BA(J)*RTOR*RTOR
	   BB(J) = BB(J)/RTOR/BTOR
	   BD(J) = BD(J)*BC(J)
	   XEQ(J) = GR(J)/GR(NEQL)
	enddo
	print *, 'EQUIL3 1422 XEQ(J) = GR(J)/GR(NEQL)', NEQL
	print *, GR(1), GR(2), GR(3)	
	print *, XEQ(1), XEQ(2), XEQ(3)	
C Define auxiliary (shifted) grid:
	do	J=1,NA
	   XTR(J) = J*HRO/ROC
	enddo
	XTR(NA1) = 1.d0

C -> dV/d(rho);  BA -> G33;  BB -> IPOL;  BC -> DRODA  (to auxiliary grid)
C Lateral surface: GRADRO*VRS=V'*<|grad(rho)|>:
C <[nabla(ro)]^2> = <g^{11}> = <g_{22}*r^2/g> :
C <[nabla(ro)/r]^2> = <g^{11}/r^2> = <g_{22}/g> :
	call	TRANSF(     NEQL,C, XEQ,NA1,VRS,XTR)
	call	SMOOTH(ALFA,NEQL,A, XEQ,NA1,G11,XTR)
	call	SMOOTH(ALFA,NEQL,B, XEQ,NA1,G22,XTR)
	call	SMOOTH(ALFA,NEQL,BA,XEQ,NA1,G33,XTR)
	call	SMOOTH(ALFA,NEQL,BB,XEQ,NA1,IPOL,XTR)
	call	SMOOTH(ALFA,NEQL,BC,XEQ,NA1,DRODA,XTR)
	call	SMOOTH(ALFA,NEQL,BD,XEQ,NA1,GRADRO,XTR)

C Arrays are artificially defined at NA1 to get smooth GRAD
	YRO = HROA/ROC/(XEQ(NEQL)-XEQ(NEQL-1))
C	VRS(NA1) = VRS(NA)+(C(NEQL)-C(NEQL-1))*YRO
C	VRS(NA1) = ARRNA1(VRS(NA),(ROC-NA*HRO)/HRO)
C	write(*,*)(NA-1)*HRO,NA*HRO,ROC
C	write(*,*)VRS(NA-1),VRS(NA),VRS(NA1)
C      write(*,*)(VRS(NA)-VRS(NA-1))/HRO,(VRS(NA1)-VRS(NA))/(ROC-NA*HRO)
C	write(*,*)2.d0*VRS(NA-1)-VRS(NA-2)-VRS(NA)
C     >		,VRS(NA-1)-VRS(NA)+HRO/(ROC-j*HRO)*(VRS(NA1)-VRS(NA))
C	G11(NA1) = G11(NA)+(A(NEQL)-A(NEQL-1))*YRO
C	G22(NA1) = G22(NA)+(B(NEQL)-B(NEQL-1))*YRO
C	G33(NA1) = G33(NA)+(BA(NEQL)-BA(NEQL-1))*YRO
C	IPOL(NA1) = IPOL(NA)+(BB(NEQL)-BB(NEQL-1))*YRO
C	DRODA(NA1) = DRODA(NA)+(BC(NEQL)-BC(NEQL-1))*YRO
C	GRADRO(NA1) = GRADRO(NA)+(BD(NEQL)-BD(NEQL-1))*YRO
	do	J=1,NA1
	   SLAT(J) = GRADRO(J)*VRS(J)
	   G22(J) = G22(J)/G33(j)*(RTOR/IPOL(j))**2
	   if (j.lt.NA1)	then
	      G22(J) = G22(J)*j*HRO
	   else
	      G22(J) = G22(J)*(NA*HRO+HROA)
	   endif
	   G11(J) = G11(J)*VRS(j)
	enddo
C-------------------------------------------------------------------------
C Determine IPOL, G33, AMETR, SHIF, ELON, TRIA, VR 
C 	    on the main transport grid:
	do	J = 1,NA1
	   XTR(J) = RHO(J)/ROC	
	enddo
C	XTR(NA) = (HRO*(NA-1)+0.5d0*HROA)/ROC
C -> VR = dV/d(rho), \Delta, \lambda, G33, IPOL, <|\na\rho|>
	call	TRANSF(NEQL,C,XEQ,NA1,VR,XTR)
	call	SMOOTH(ALFA,NEQL,GBD,XEQ,NA1,SHIF, XTR)
	call	SMOOTH(ALFA,NEQL,GL, XEQ,NA1,ELON, XTR)
	call	SMOOTH(ALFA,NEQL,BA, XEQ,NA1,G33,XTR)
	call	SMOOTH(ALFA,NEQL,BB, XEQ,NA1,IPOL, XTR)
C (option 1) Volume as taken from the equilibrium solver:
C	call	SMOOTH(ALFA,NEQL,D,XEQ,NA1,VOLUM, XTR)

C GL -> a
C GSD -> \delta^{ASTRA} (dimensionless)
	YDA = ABC/(NEQL-1.d0)
	do	j =1,NEQL
	   A(j) = YDA*(J-1.d0)
	   if (j.ne.1) B(J) = GSD(J)/A(J)
	enddo
	B(1)=0.d0
	print *, '----'
	print *, A(1),A(2),A(3)
	print *, XEQ(1),XEQ(2),XEQ(3)
	call	TRANSF(NEQL,B,XEQ,NA1,TRIA,XTR)
	call	TRANSF(NEQL,A,XEQ,NA1,AMETR,XTR)
	print *, AMETR(1),AMETR(2),AMETR(3)
	print *, XTR(1),XTR(2),XTR(3)
	pause
	do	J=1,NA1
C (option 2a, default) Volume is calculated using VR (shifted grid):
C (option 2b) Volume is calculated using VINT(1) (main grid):
C	   B(j) = 1.
C	   VOLUM(j) = VINT(B,RHO(j))
C (option 3) Volume is calculated using analytic formula
C	   VOLUM(J) = GP*GP2*ELON(J)*AMETR(J)**2
C     *		    * (RTOR+SHIF(J)-0.25*TRIA(J)*AMETR(J))
	   SHIF(J) = SHIFT+SHIF(J)
	   if     (j .eq. 1)	then
	      SHEAR(J) = (FP(2)-FP(1))/(2.d0*MU(1)+0.333d0*(MU(1)-MU(2)))
	   elseif (j .lt. NA)	then
	      SHEAR(J) = (FP(j+1)-2.d0*FP(j)+FP(j-1))/(MU(j+1)+MU(j))
	   elseif (j .eq. NA)	then
	      SHEAR(J) = HRO*(FP(j+1)-FP(j))/HROA-FP(j)+FP(j-1)
	      SHEAR(J) = SHEAR(J)/(MU(j+1)+MU(j))
C	   else
C	      SHEAR(J) = HRO*(FP(NA1)-FP(NA))/HROA-FP(NA)+FP(NA-1)
C	      SHEAR(J) = 2.d0*SHEAR(J)*(2.d0*HROA-HRO)
C	      SHEAR(J) = SHEAR(J)/(2.d0*HROA*MU(NA1)-HRO*MU(NA))
	   endif
	   SHEAR(J) = 1.d0 -SHEAR(J)/(GP*BTOR*HRO**2)
	enddo
	SHEAR(NA1) = SHEAR(NA)
C	write(*,'(4f9.3,3X,4f9.3)')
C     +		(SHEAR(j),j=NA1-3,NA1),(MU(j),j=NA1-3,NA1)
	do	J = 1,NAB
	   SHIV(J) = UPDWN	
	enddo

C Alternative definition for BDB02 (G.Pereverzev 10.02.99)
C	do	J=1,NA1
C		YTH2	=RHO(j)*G22(j)*(MU(j)/RTOR)**2
C		BDB02(J)=(1.+YTH2)*G33(J)*IPOL(J)**2
C	enddo
C
CMR additional quantities
C Note!	   BDB02(j) = G33(j)*IPOL(j)**2+(RHO(j)*MU(j)/RTOR)**2
C   BDB02 - <B**2/B0**2>
C   B0DB2 - <B0**2/B**2>
C   BMAXT - BMAXT
C   BMINT - BMINT
C   BDB0  - <B/BTOR>
C   FOFB  - <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>

        call TRANSF(NEQL,B2B0EQ,XEQ,NA1,BDB02,XTR)
        call TRANSF(NEQL,B0B2EQ,XEQ,NA1,B0DB2,XTR)
        call TRANSF(NEQL,BMAXEQ,XEQ,NA1,BMAXT,XTR)
        call TRANSF(NEQL,BMINEQ,XEQ,NA1,BMINT,XTR)
        call TRANSF(NEQL,BMODEQ,XEQ,NA1,BDB0,XTR)
        call TRANSF(NEQL,FOFBEQ,XEQ,NA1,FOFB,XTR)
	do	J=NA1,NAB
	   SHEAR(j) = SHEAR(NA1)
	   BDB02(j) = BDB02(NA1)
	   B0DB2(j) = B0DB2(NA1)
	   BMAXT(j) = BMAXT(NA1)
	   BMINT(j) = BMINT(NA1)
	   BDB0(j) = BDB0(NA1)
	   FOFB(j) = FOFB(NA1)
	enddo
C Volume (on the shifted grid) is calculated using VR:
	VOLUM(1) = HRO*VR(1)
	do	J=2,NA
	   VOLUM(J) = VOLUM(J-1)+HRO*VR(J)
	enddo
	VOLUM(NA1) = VOLUM(NA)+(HROA-0.5d0*HRO)*VR(NA)
	VOLUME = VOLUM(NA1)
	end
C======================================================================|
	subroutine	NEWGRD
C----------------------------------------------------------------------|
C input: 	ROC, HRO, NB1, NA1, NA=NA1-1, AB, ABC,  AMETR(NA1)
C Output:	NA1, NA=NA1-1, NAB, RHO(NA1)=ROC, HROA, AMETR(j>NA1)
C----------------------------------------------------------------------|
C Define size of the edge cell to be:   .6 <= HROA/HRO < 1.8
C	with a hysteresis of 0.2d0*HRO, so that
C	   if (HROA<=0.6d0*HRO)  jumps to  HROA<=1.6d0*HRO ! remove one cell
C	   if (HROA>=1.8d0*HRO)  jumps to  HROA>=0.8d0*HRO !   add  one cell
C	write(*,*)NA1,TIME,RHO(NA1)," -> ",ROC
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	j
	double precision YDA
C	write(*,*)NA1,ROC,HRO*(NA+0.1),HRO*(NA1+0.3),HRO*(NB1-0.49)
C	write(*,'(1P,6E13.5)')(RHO(j),j=NA-2,NA1+1)
C	write(*,'(1P,6E13.5)')(AMETR(j),j=NA-2,NA1+1)
	if (ROC .gt. HRO*(NB1-0.49d0))	then
C rho_edge gets out of the grid
	   write(*,*)'>>> WARNING: the allocated grid is too small'
	   write(*,*)'Time =',TIME,'   Rho_b =',ROC,
     >		'   Rho_max = ',HRO*(NB1-0.5d0)
	   write(*,*)'Try to increase AWALL'
	   write(*,*)'If this does not help inspect equilibrium input'
	   NA = NB1-1
	elseif (ROC .le. HRO*(NA+0.1d0))	then	! 0.1 <=> 0.6-0.5
!	    if (HROA < 0.6d0*HRO) then   reduce NA
C	   write(*,*)" <- ",ROC,HRO*NA,RHO(NA1),HRO*NA1,NA1
	   do	j=NA-1,1,-1
		NA = j
C	   	write(*,*)j,NA*HRO,ROC,(NA+1)*HRO,ROC/HRO-(NA-0.5)
		if (ROC .gt. HRO*(j+0.1d0))	goto	10
	   enddo
	elseif (ROC .gt. HRO*(NA1+0.3))	then	! 0.3 <=> 1.8-1.5
!	    if (HROA > 1.8d0*HRO) then   increase NA
C	   write(*,*)" -> ",ROC,HRO*NA,RHO(NA1),HRO*NA1,NA1
	   do	j=NA1,NB1
		NA = j
C	   	write(*,*)j,NA*HRO,ROC,(NA+1)*HRO
		if (ROC .le. HRO*(j+1.3d0))	goto	10
	   enddo
	endif
 10	continue

C HRO*(NA+0.6) < ROC <= HRO*(NA+1.8)
C	write(*,*)AB,ABC,TIME,TSTART
C	write(*,*)NA+1,NA+1-NA1,HRO*(NA+0.6),ROC,HRO*(NA+1.8)
	if (abs(NA+1-NA1).gt.1 .and. TIME.gt.TSTART)
     >	   write(*,'(A,F6.3,/,2A)')
     >     ">>> Warning: Too fast volume variation at t =",TIME,
     >	   "             Reduce TAUMAX or NB1.",char(7)
	if (NA .lt. 10) write(*,*)
     >	  ">>> Warning: Spatial grid < 10 points. Increase NB1.",char(7)
	NA1 = NA+1
	do	J=1,NB1
	   RHO(J)=(J-0.5d0)*HRO
	enddo
	HROA = ROC-RHO(NA)
	RHO(NA1) = ROC
	AMETR(NA1) = ABC

	if (HROA/HRO .lt. 0.59d0)	write(*,*)	! Should never happen?
     >	'Warning: edge cell is too small  ',HROA/HRO,ROC/HRO-(NA-0.5)
 	NAB = NA1
	ROB = RHO(NAB)
	if ( 2.d0*abs(AB-ABC) .lt. AMETR(NA1)-AMETR(NA) )	return
	NAB = NA1/ABC*AB
	ROB = RHO(NAB)
C	write(*,'(A,2I5,4f7.4)')"NAB->",NA1,NAB,ROC,ROB,ABC,AB
C	if (NA1 .gt. NAB)	write(*,*)"AMETR grid error"

C "a" is calculated assuming "rho" given
C	call	SETGEO(NA1)
	if (NA1+1.gt.NB1 .or. NA1.eq.NAB)	return

C	write(*,'(4F7.3,3I5)')ABC,AB,ROC,ROB,NA1,NAB,NB1
	YDA = (AB-ABC)/(NAB-NA1)
	do	j=NA1+1,NB1
	   AMETR(j) = ABC+YDA*(j-NA1)
C	   DRODA(j) = YDA/(YR2-YR1)
	enddo

C	write(*,'(1P,6E13.5)')(  RHO(j),j=NA,NAB)
C	write(*,'(1P,6E13.5)')(AMETR(j),j=NA,NAB)
C	write(*,'(1P,6E13.5)')(  RHO(j)-RHO(j-1),j=NA,NAB)
C	write(*,'(1P,6E13.5)')(AMETR(j)-AMETR(j-1),j=NA,NAB)
	if (NA1 .lt. NB1)	then
	   do	j=NA1+1,NB1
		if (AMETR(j) .lt. AB)   NAB = j
	   enddo
	   if (NAB .lt. NB1)	NAB=NAB+1
	endif
	ROB = RHO(NAB)
	AMETR(NAB) = AB
C	write(*,*)"Done"
C	write(*,*)NA1
C	write(*,'(1P,6E13.5)')(AMETR(j),j=NA-2,NA1+1)
C	write(*,*)
C	write(*,*)NA1,NAB,NB1,HA,ABC-AMETR(NA)
C	write(*,*)ROC,RHO(NA1),ROB,RHO(NAB)
C	write(*,'(10f7.4)')(RHO(j),j=NA-1,NB1)
C	write(*,*)AMETR(NA),ABC,AB,AMETR(NAB)
C	write(*,'(10f7.4)')(AMETR(j),j=NA-1,NB1)
C "rho" is given, "a" calculated
C	NAB = NA1
C	Y1 = sqrt((RTOR+SHIFT)**2-AB**2)
C	Y2 = sqrt((RTOR+SHIFT)**2-ABC**2)
C	ROB = sqrt(ROC**2+2.d0*ELONG*RTOR*(AB**2-ABC**2)/(Y1+Y2))
C	if (NA1 .lt. NB1)	then
C	   Y1 = 1./(ELONG*RTOR)
C	   do	j=NA1+1,NB1
C		Y3 = Y1*(RHO(j)**2-ROC**2)
C		AMETR(j) = sqrt(ABC**2+Y3*(Y2-0.25*Y3))
C		if (AB .ge. AMETR(j))	NAB = j
C	   enddo
C	endif
C	if (NAB .lt. NB1)	then
C	   if (AB-AMETR(NAB) .gt. AMETR(NAB+1)-AB)	NAB = NAB+1
C	endif
C	AMETR(NAB) = AB
C	RHO(NAB) = ROB
	end
C======================================================================|
	subroutine	SETGEO(jst)
C----------------------------------------------------------------------|
C Presently the subroutine is called with jst=0 only
C----------------------------------------------------------------------|
C input: 	jst,NB1,HRO,ROC,AB,RHO(j)
C Output:	SHIF(jst:NB1),  ELON(jst:NB1),  TRIA(jst:NB1),
C		AMETR(jst:NB1), DRODA(jst:NB1)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	jst,j
	double precision	YDA,YA,YR1,YR2,ROC3A
C----------------------------------------------------------------------|
	if (jst+1 .gt. NB1)	return
	do	j=jst+1,NB1
	    YR2		= min(1.d0,(RHO(J)/ROC)**2)
	    SHIF(J)	= SHIFT
	    SHIV(J)	= UPDWN
	    ELON(J)	= 0.5d0*(1.d0+ELONG+(ELONG-1.d0)*YR2)
C	    ELON(J)	= (0.75+0.25*YR2)*ELONG
	    TRIA(J)	= TRIAN*YR2
	enddo
	YDA = .1d0*AB/NB1

	YR1 = .0d0
	if (jst .eq. 0)	then
	   YA = 0.d0
	   YR2 = 0.d0
	else
	   YA = AMETR(jst+1)
	   YR2 = RHO(jst+1)
	endif
	do	j=jst+1,NB1
 2	   if (YR2 .gt. RHO(j))	goto	3
	   YR1 = YR2
	   YA = YA+YDA
	   YR2 = ROC3A(RTOR,SHIF(j),YA,ELON(j),TRIA(j))
	   if (YR2 .le. RHO(j))	goto	2
	   DRODA(j) = YDA/(YR2-YR1)
 3	   AMETR(j) = YA-YDA+(RHO(j)-YR1)*DRODA(j)
	enddo
	return
	end
C======================================================================|
	subroutine	EDCELL(JNA1O,YIPL)
C----------------------------------------------------------------------|
C Input:
C	JNA1O - old grid size
C	YIPL = (FP(NA1)-FP(NA))/HROA
C Output:
C	if (NA1.eq.JNA1O) then FP(NA1) 
C       	all main arrays otherwise
C Note! FP defined in this subroutine will be assigned to FPO in OLDNEW
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	JNA1O,JNAO,j
	double precision	YIPL,YL,YR,YD,ARRNA1
	if (NA1 .gt. JNA1O)	then	! NA1 increases: JNA1O >= NA
	   JNAO = JNA1O-1
	   TE(NA1) = TE(JNA1O)
	   TI(NA1) = TI(JNA1O)
	   NE(NA1) = NE(JNA1O)
	   NI(NA1) = NI(JNA1O)
	   NIZ1(NA1) = NIZ1(JNA1O)
	   NIZ2(NA1) = NIZ2(JNA1O)
	   NIZ3(NA1) = NIZ3(JNA1O)
	   NALF(NA1) = NALF(JNA1O)
	   NHE3(NA1) = NHE3(JNA1O)
	   NHYDR(NA1) = NHYDR(JNA1O)
	   NDEUT(NA1) = NDEUT(JNA1O)
	   NTRIT(NA1) = NTRIT(JNA1O)
	   F0(NA1) = F0(JNA1O)
	   F1(NA1) = F1(JNA1O)
	   F2(NA1) = F2(JNA1O)
	   F3(NA1) = F3(JNA1O)
	   F4(NA1) = F4(JNA1O)
	   F5(NA1) = F5(JNA1O)
	   F6(NA1) = F6(JNA1O)
	   F7(NA1) = F7(JNA1O)
	   F8(NA1) = F8(JNA1O)
	   F9(NA1) = F9(JNA1O)
	   YD = RHO(NA1)-RHO(JNAO)
	   FP(NA1) = FP(JNAO)+(FV(NA1)-FV(NA))+YD*YIPL
C	YIPL = (FP(NA1)-FP(NA)-(FV(NA1)-FV(NA)))/HROA
	   do	j=JNA1O,NA ! Linear interpolation into the second last node
	      YR = (RHO(j)-RHO(JNAO))/YD
	      YL = (RHO(NA1)-RHO(j))/YD
	      TE(j) = TE(NA1)*YR+TE(JNAO)*YL
	      TI(j) = TI(NA1)*YR+TI(JNAO)*YL
	      NE(j) = NE(NA1)*YR+NE(JNAO)*YL
	      NI(j) = NI(NA1)*YR+NI(JNAO)*YL
	      F1(j) = F1(NA1)*YR+F1(JNAO)*YL
	      F2(j) = F2(NA1)*YR+F2(JNAO)*YL
	      F3(j) = F3(NA1)*YR+F3(JNAO)*YL
	      FP(j) = FP(JNAO)+(FV(j)-FV(JNAO))+HRO*(j-JNAO)*YIPL
	      CU(j) = CU(JNAO)
	   enddo
C Define FP(NA) from sigma*E=j_OH at j=NA
C	   j = NA-1
C	   YD = VR(j)*CD(j)+GP2*RHO(j)*CC(j)*(FP(j)-FPO(j))/TAU
C	   YD = 0.2d0*HRO*HRO*YD/IPOL(j)**2
C	   FP(NA) = FP(j)+(G22(j)*(FP(j)-FP(j-1))+YD)/G22(NA)
C	   FP(NA1) = FP(NA)+.4d0*GP*RTOR*IPLOLD*HROA/(G22(NA1)*IPOL(NA1))
	   CU(NA1) = 0.d0
	   YD = RTOR*IPL/IPOL(JNAO)/G22(JNAO)
	   MU(NA1) = YD/(5.d0*BTOR*HRO*(JNAO))
	elseif (NA1 .lt. JNA1O)	then		! NA1 decreases
	   TE(NA1) = TE(JNA1O)
	   TI(NA1) = TI(JNA1O)
	   NE(NA1) = NE(JNA1O)
	   NI(NA1) = NI(JNA1O)
	   NIZ1(NA1) = NIZ1(JNA1O)
	   NIZ2(NA1) = NIZ2(JNA1O)
	   NIZ3(NA1) = NIZ3(JNA1O)
	   NALF(NA1) = NALF(JNA1O)
	   NHE3(NA1) = NHE3(JNA1O)
	   NHYDR(NA1) = NHYDR(JNA1O)
	   NDEUT(NA1) = NDEUT(JNA1O)
	   NTRIT(NA1) = NTRIT(JNA1O)
	   F0(NA1) = F0(JNA1O)
	   F1(NA1) = F1(JNA1O)
	   F2(NA1) = F2(JNA1O)
	   F3(NA1) = F3(JNA1O)
	   F4(NA1) = F4(JNA1O)
	   F5(NA1) = F5(JNA1O)
	   F6(NA1) = F6(JNA1O)
	   F7(NA1) = F7(JNA1O)
	   F8(NA1) = F8(JNA1O)
	   F9(NA1) = F9(JNA1O)
C	   MU(NA1) = RTOR*IPL/(5.d0*BTOR*HRO*NA*G22(NA)*IPOL(NA))
C	   MU(NA1) = ARRNA1(MU(NA),HROA/HRO)
	   MU(NA1) = RTOR*IPL/(5.d0*BTOR*ROC*G22(NA)*IPOL(NA1))
	   FP(NA1) = FP(NA)+(FV(NA1)-FV(NA))+HROA*YIPL
	else					! NA1 doesn't change
	   FP(NA1) = FP(NA)+(FV(NA1)-FV(NA))+HROA*YIPL
	endif
	end
C======================================================================|
	subroutine	ADCMP
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision	FACTOR,XST(NRD),XFN(NRD+1),A(NRD+1)
	double precision	ALFA,YN,YT,YM,YM1,YCM,YCU,CMU,BTOLD,YA
	integer	NFN,NF11,J
	save	NFN,BTOLD
	data	BTOLD/-1.d0/
	call	markloc("ADCMP"//char(0))
	if (BTOLD .lt. 0.d0)	BTOLD = BTOR
	YA = (BTOR-BTOLD)/(TAU*(BTOR+BTOLD))		! 1/2 d[ln(B)]/dt
	do	j=1,NA1
	   VPFP(j) = YA*RHO(j)
	enddo
	if (abs(BTOR-BTOLD)/BTOR .lt. 1.d-6)	return

	FACTOR = SQRT(BTOR/BTOLD)
	BTOLD = BTOR
	ALFA	=1.d-4
	do	1	J=1,NA1
	   XST(J) = (J-0.5d0)/(NA1-0.5d0)
	   XFN(J) = XST(J)/FACTOR
1	continue
	do	2	J=1,NA1
	if (XFN(J) .lt. 1.d0)	goto	2
C			Decompression (FACTOR < 1) :
	if (J.ne.NA1)	then
	   NE(J) = NE(NA1)
	   TE(J) = TE(NA1)
	   TI(J) = TI(NA1)
	endif
	NFN	=J
	goto	3
2	continue
C			Compression (FACTOR > 1) :
	NFN	=NA1+1
3	continue
	XFN(NFN)=1.d0
	YN	=FACTOR*FACTOR
	YT	=FACTOR**1.333d0
	do	4	J=1,NA1-1
	   if(J.ge.NFN)	goto	41
	   NE(J) = NE(J)*YN
	   TE(J) = TE(J)*YT
	   TI(J) = TI(J)*YT
4	continue
41	continue
	do	5	J=1,NA1
	   A(J)	= NE(J)
5	continue
	A(NA1+1)=NE(NA1)
	call	SMOOTH(ALFA,NFN,A,XFN,NA1,NE,XST)
	do	6	J=1,NA1
		A(J)	=TE(J)
6	continue
	A(NA1+1)=TE(NA1)
	call	SMOOTH(ALFA,NFN,A,XFN,NA1,TE,XST)
	do	7	J=1,NA1
		A(J)	=TI(J)
7	continue
	A(NA1+1)=TI(NA1)
	call	SMOOTH(ALFA,NFN,A,XFN,NA1,TI,XST)
C	call	TRANSF(NFN,A,XFN,NA1,TI,XST)
	do	71	J=1,NA
		A(J)	=MU(J)
		XST(J)	=J/(NA1-0.5d0)
		XFN(J)	=XST(J)/FACTOR
71	continue
	if(FACTOR.lt.1.d0)	then
	   NFN	=NA
	   NF11	=NA-1
	   call	TRANSF(NA,A,XFN,NF11,MU,XST)
	   MU(NA)	=MU(NA)/FACTOR**2
	else
	   do	72	J=1,NA1
	      if(XST(J).le.XFN(NA))	NFN	=J
 72	   continue
	   call	TRANSF(NA,A,XFN,NFN,MU,XST)
	endif
C	CMU	=MU(NFN)*(NFN)*G22(NFN)
	CMU	=MU(NFN)*NFN*NFN
	do	9	J=NFN+1,NA1
C		MU(J)	=CMU/(J*G22(J))
		MU(J)	=CMU/(J*J)
9	continue
	YM=0.d0
	YCM	=GP2*BTOR*HRO**2
	YCU=2.5d0*BTOR/(GP*RTOR*HRO)
	do	91	J=1,NA
		FP(J+1)	=FP(J)+YCM*J*MU(J)
		YM1=YM
		YM=J*G22(J)*MU(J)
		CU(J)=YCU*G33(J)*IPOL(J)**3*(YM-YM1)/(J-0.5d0)
91	continue
	FP(NA1)	= FP(NA)+GP2*BTOR*HRO*NA*MU(NA)*HROA
	CU(NA1) = CU(NA)+(CU(NA)-CU(NA-1))/HRO*HROA
C	CU(NA1) = CU(NA)
99	end
C======================================================================|
C Input:
C BA(1:NEQL)[MA/m^2] - flux function in front of R/r on the rhs of G-Sh. eq.
C BB(1:NEQL)[MA/m^2] - flux function in front of (r/R-R/r) - " - " - " - " -
C BR00 [m] - R the geom. center of the edge magnetic surface
C SA0  [m] - a_{edge} - "a" for the edge magnetic surface
C GL0  [d/l] - elongation of the edge magnetic surface
C GD30 [m]   - triangularity of the edge magnetic surface (\delta(a_{edge}))
C NEQL [d/l]- # of grid points
C B0T  [T]   - vacuum toroidal magnetic field at the point r=R
C PLCUR[MA]  - plasma current
C
C Output:
C GR(1:NEQL) [m]	- \rho(a)
C GBD(1:NEQL) [m] - \Delta(a)
C GL(1:NEQL) [d/l] - \lambda(a)
C GSD(1:NEQL) [d/l] - \delta(a)/[a_{edge}(a/a_{edge})^2]
C GRA(1:NEQL) [d/l] - <(\nabla a)^2>
C GRAR(1:NEQL) [1/m^2] - <[(\nabla a)/r]^2>
C AVR2(1:NEQL) [1/m^2] - <[1/r]^2>
C AI0(1:NEQL) [mT] - I = r B_{tor}
C Everywhere above 
C (i)  the same equidistant grid with respect to the "radial" 
C      variable "a" is assumed to be used both for input and output.
C (ii) average <f> is an [0.5/\pi] integral over the poloidal angle \tau 
C      of the quantity [f\sqrt{g}], where g_{ij} is the metric tensor 
C      of {a,\tau,\zeta}.
C (iii) r is the current polar radius (the distance to the major torus axis)
C

C - BR00,SA0 - MAJOR & MINOR RADII /METER/
C - GLO,GD3O - ELONGATION AND TRIANGULARITY /BOUNDARY VALUE/
C - GSD = TRIANGULARITY

C=======================================================================
C	Input:	BA [MA/m^2] Zakharov's input function on the rhs current  
C		BB [MA/m^2] --"--"--"--"--"--"--"--"--"--- presure 
C		BR00  [m]   major radius of the edge magnetic surface
C		SA0   [m]   minor radius in the equatorial plane
C		GL0  [d/l]  edge elongation
C		GD30  [m]   edge triangularity
C		NA1         number of radial points
C		B0T   [T]   vacuum magnetic field at the position BR00 
C		PLCUR [MA]  total plasma current
C	Output:	GR  - rho
C		GBD - shift
C		GL  - elongation
C		GSD - triangularity
C		GRA - <G11>
C		SQGRA - <sqrt(G11)>
C		GRAR - <G11/R#2>
C		AVR2 - <1/R#2>
C		AI0  - IPOL*RTOR*BTOR
C======================================================================|
      subroutine E3ASTR
C Input:
     &	(BA,BB	! j_zeta = BA*(R_0/r) + BB*(r/R_0-R_0/r)
     &	,BR00	! R_0+\Delta_edge		! RTOR+SHIFT
     &	,SA0	! a_edge			! ABC
     &	,GL0	! \lambda_edge			! ELONG
     &	,GD30	! \delta_edge			! TRIAN*ABC
     &	,NA1	! radial grid point No.		! NEQL
     &	,ACC	! radial grid point No.		! ACEQLB
C Output:
     &	,GR,GBD,GL,GSD	! \Rho(a),\Delta(a),\lambda(a),\delta(a)
     &	,GRA	! <g^{11}>		= <[nabla(a)]^2>
     &	,SQGRA	! <\sqrt{|g^{11}|}>	= <|nabla(a)|>
     &	,GRAR	! <g^{11}g^{33}>	= <[nabla(a)/r]^2>
     &	,AVR2	! <g^{33}>		= <1/r^2>=G33/RTOR**2
     &	,AI0	! I				-> IPOL*RTOR*BTOR
     &	,dgrda	! \prti{\rho}{a}
     &	,avsqg	! \prti{V}{a}/(4\pi^2)
     &	,Vol	! V(a)
C Input:
     &	,B0T	! B_0 at the magnetic axis	! BTOR*RTOR/(RTOR+SHIFT)
     &	,PLCUR	! Total plasma current		! IPL
C Output
     &  ,B2B0EQ ! <B**2/B0**2>
     &  ,B0B2EQ ! <B0**2/B**2>
     &  ,BMAXEQ ! BMAXT
     &  ,BMINEQ ! BMINT
     &  ,BMODEQ ! <B/BTOR>
     &  ,FOFBEQ ! <(BTOR/B)**2*(1.-SQRT(1-B/Bmax)*(1+.5B/Bmax))>
     &  ,GRDAEQ,TIME)! <grad a>
C - BR00,SA0 - MAJOR & MINOR RADII /METER/
C - GLO,GD3O - ELONGATION AND TRIANGULARITY /BOUNDARY VALUE/
C - GSD = TRIANGULARITY
C	\Delta' =WDSD1(I)*WSA(I)
C	\lambda'=WDGL(I)*WGL(I)*WSAA(I)
C	\delta' =WDSD3(I)*WSA(I)

	implicit none
	include 'for/emeq.inc'
	integer	NA1,NAOLD,NA,NT,NT1,I,I1,J,K,NITER
	double precision	BR00,SA0,GL0,GD30,ACC,B0T,PLCUR,TIME
	double precision	BA(NP),BB(NP),GR(NP),GBD(NP),GL(NP),
     &		GSD(NP),GRA(NP),SQGRA(NP),GRAR(NP),AVR2(NP),AI0(NP),
     &		dgrda(NP),avsqg(NP),Vol(NP),GR2AUX(NP)
	double precision	B2B0EQ(*),B0B2EQ(*),BMAXEQ(*),BMINEQ(*),
     1		BMODEQ(*),FOFBEQ(*),GRDAEQ(*),
     2		A,AA,C,CC,S,SS,SR,SX,SX1,T,Y,Y1,SDT,SDT0,
     3		DRDA,DZDA,DRDT,DZDT,DMETR,DA2,DGR2,FI,FJ,D0,CGP,
     4		GP,GP2,GR2,GLOLD,G3DOLD,AOLD,G22A2,SKGGG,
     5		SKDR,SKGA,SQG,SQG22R,YLIN,YVOL,YMIN,YMAX
      save AOLD,GLOLD,G3DOLD,NAOLD,cgp,NITER
      common /EMEQMR/SKDR(NP),SKGA(NP),SQG22R(NP)
      data AOLD/0.d0/GLOLD/0.d0/G3DOLD/0.d0/NAOLD/1/
      data cgp/3.14159265359d0/	NITER/60/

*************************************************
C	if (TIME .gt. .1)	then
C	write(*,*)BR00	! R_0+\Delta_edge		! RTOR+SHIFT
C	write(*,*)SA0	! a_edge			! ABC	   
C	write(*,*)GL0	! \lambda_edge			! ELONG	   
C	write(*,*)GD30	! \delta_edge			! TRIAN*ABC0
C	write(*,*)NA1
C	write(*,*)B0T	! B_0 at the magnetic axis ! B_0*R_0/(R_0+SHIFT)
C	write(*,*)PLCUR	! Total plasma current		! IPL
C	write(*,'(2(E14.6))')(BA(j),BB(j),j=1,NA1)
C	stop
C	endif
**************************************************

      WBBS0=B0T
      WBJ0=0.2d0*PLCUR
      NA=NA1-1
      NT=12*MAX(1,NA1/8)
C 		Set initial conditions / zero iteration
      if (               NAOLD .eq. NA
     &	.and. abs(AOLD-SA0)    .lt. 1.E-6
     &	.and. abs(GLOLD-GL0)   .lt. 1.E-6
     &	.and. abs(G3DOLD-GD30) .lt. 1.E-6
     *   )	goto	4
      call add2loc("Calling EQGB3"//char(0))
      call EQGB3(BR00,SA0,GL0,GD30,NA)
      AOLD=SA0
      GLOLD=GL0
      G3DOLD=GD30
      NAOLD=NA
      WBR00=BR00
	  print *, 'WBR00=BR00'
 4    continue
      do  I=1,NA1
         WSJP(I)=BA(I)
         WSP(I)=BB(I)
      enddo
      NITER=30				! Use 60 for the 1st entry only
      call add2loc("Calling EQAB3"//char(0))
      call EQAB3(NA,NT,NITER,ACC)		! Call MEM equil solver
      if (NA .eq. 0)	then
	 j = 7
	 write(*,'(2A,F10.6)')char(j),
     >	 ">>> No convergence in the 3M equilibrium solver at t =",TIME
	 call	a_stop
      endif
      if (NA .le. 1)	then
	 NA1 = NA
	 return
      endif
      call add2loc("Calling EQPPAB"//char(0))
	  call EQPPAB(NA)
	  print *, 'WBR00', WBR00	  
	  print *, 'wbr0', wbr0	  
      D0=WSD1(NA1)*WSAA(NA1)
      GR(1)=0.d0
      Vol(1)=0.d0
      GR2=0.d0
      fi=0.d0
      s=4.d0*cgp*cgp
      do I=1,NA1
         Vol(i)=s*WSAA(i)*(WBR0*wsl0(i)+WSAA(i)*wsl1(i))
         j=i-1
         fj=fi
         AI0(I)=WBF(I)*WBR0
         sqg=wbg0(i)*wgl(i)*wbr0
         avsqg(i)=WSA(i)*sqg
         GRA(I)	=SKGA(I)/sqg
         SQGRA(I)=SQG22R(I)/sqg
         GRAR(I)=wbg22(i)/(wgl(i)*wbr0*sqg)
         AVR2(I)=wbg33(i)*wgl(i)/(wbr0*sqg)
         GL(I)=WGL(I)
         GBD(I)=D0-WSD1(I)*WSAA(I)
         GSD(I)=WSD3(I)*WSAA(i)
         fi=wbg33(i)*wgl(i)*AI0(I)/(wbr0*WBBS0)
         if (I.gt.1) then
            DGR2=fj*WSCJ1(j)+fi*WSCI1(j)
            GR2=GR2+DGR2
	    if (DGR2.le.0.)	then
		write(*,*)"NA1 = 0",GR2,DGR2
	 	NA1 = 0
	 	return
	    endif
            GR(I)=SQRT(GR2*2.d0)
            dgrda(i)=fi*WSA(i)/gr(i)
         else
            dgrda(i)=sqrt(fi)
         endif
      enddo
	  print *, 'GR', GR(1),GR(2),GR(3)
	  print *, 'wbr0', wbr0
CMR extra quantities
      GR2AUX(1)=0.d0
      GP=3.1415926d0
      GP2=2.d0*GP
      NT1=NT+1
      SDT0=1.d0/NT
	YLIN=0.d0
	YVOL=0.d0
      do 10 I=1,NA1
         A=WSA(I)
         AA=A*A
         B2B0EQ(I)	=0.d0
         B0B2EQ(I) 	=0.d0
         BMODEQ(I) 	=0.d0
         FOFBEQ(I) 	=0.d0
         GRDAEQ(I)	=0.d0
         SKGGG		=0.d0
         YMIN		=99999.d0
         YMAX		=0.
         AI0(I)=WBF(I)*WBR0
         IF(I.EQ.1) GOTO 3
         K=I-1
         GR2AUX(I)=GR2AUX(K)+(SKDR(K)*AI0(K)*WSCJ1(K)+
     +                 SKDR(I)*AI0(I)*WSCI1(K))*2.d0/WBBS0
 3       CONTINUE
C Need to calculate Ymax and Ymin, before doing flux surface averages
C                                        ( Y=B^2 )
         DO 50 K=1,NT1
            T=SDT0*GP*(K-1)
            SDT=SDT0
            IF(K.EQ.1.OR.K.EQ.NT1) SDT=SDT0*0.5d0
            C=COS(T)
            S=SIN(T)
            SS=S*S
            CC=1.d0-SS
            SX1=-WSD1(I)-WSD3(I)*SS
            SX=C+A*SX1
            SR=WBR0+A*SX
C R,Z derivatives
            DRDA=-WDSD1(I)*A+C-WDSD3(I)*A*SS
            DZDA=S*WGL(I)*(AA*WDGL(I)+1)
            DRDT=-A*S-2.d0*AA*WSD3(I)*C*S	
            DZDT=WGL(I)*A*C
C metric tensor components
            DMETR=DRDA*DZDT-DRDT*DZDA
            SKGGG=SKGGG+DMETR*SR*SDT
            if(I.NE.1)then
               I1=I-1
               Y=(AI0(I)/SR/WBBS0)**2 + (DRDT**2+DZDT**2)*
     *	        ((GR2AUX(I)-GR2AUX(I1))/(WSA(I)-WSA(I1))/
     /	        (WSQ(I)+WSQ(I1))/SR/DMETR)**2
            else
               Y=(AI0(I)/SR/WBBS0)**2
            endif
            YMAX=MAX(Y,YMAX)
            YMIN=MIN(Y,YMIN)
50       continue
         DO 40 K=1,NT1
            T=SDT0*GP*(K-1)
            SDT=SDT0
            IF(K.EQ.1.OR.K.EQ.NT1) SDT=SDT0*0.5d0
            C=COS(T)
            S=SIN(T)
            SS=S*S
            CC=1.d0-SS
            SX1=-WSD1(I)-WSD3(I)*SS
            SX=C+A*SX1
            SR=WBR0+A*SX
C R,Z derivatives
            DRDA=-WDSD1(I)*A+C-WDSD3(I)*A*SS
            DZDA=S*WGL(I)*(AA*WDGL(I)+1)
            DRDT=-A*S-2.d0*AA*WSD3(I)*C*S
            DZDT=WGL(I)*A*C
c G22/A**2
            G22A2=SS+4.d0*A*WSD3(I)*SS*C+(2.d0*A*WSD3(I)*S*C)**2+
     +            (WGL(I)*C)**2
C D**2/A**2
            DA2=(C-A*(WDSD1(I)+WDSD3(I)*SS))*WGL(I)*C+
     +             WGL(I)*SS*(WDGL(I)*AA+1.d0)*(1.d0+2.d0*A*C*WSD3(I))
            DA2=DA2*DA2
C metric tensor components
            DMETR =DRDA*DZDT-DRDT*DZDA
            if(I.NE.1)then
               I1=I-1
               Y=(AI0(I)/SR/WBBS0)**2 + (DRDT**2+DZDT**2)*
     *	          ((GR2AUX(I)-GR2AUX(I1))/(WSA(I)-WSA(I1))/
     /            (WSQ(I)+WSQ(I1))/SR/DMETR)**2
               B2B0EQ(I)=B2B0EQ(I)+ Y*DMETR*SR*SDT
               B0B2EQ(I)=B0B2EQ(I)+ DMETR*SR*SDT/Y
               BMODEQ(I)=BMODEQ(I)+ SQRT(Y)*DMETR*SR*SDT
               IF (Y.GT.YMAX) YMAX=Y
	       Y1 = SQRT(Y/YMAX)
               FOFBEQ(I)=FOFBEQ(I)+1.d0/Y*(1.d0-SQRT(abs(1.d0-Y1))
     +          *(1.d0+.5d0*Y1))*DMETR*SR*SDT
               GRDAEQ(I)=GRDAEQ(I)+ SQRT(G22A2/DA2)*DMETR*SR*SDT
	       YLIN	=YLIN+(DRDT**2+DZDT**2)*
     *	          ((GR2AUX(I)-GR2AUX(I1))/(WSA(I)-WSA(I1))/
     /            (WSQ(I)+WSQ(I1))/SR/DMETR)**2*DMETR*SR*SDT
	yvol = yvol + DMETR*SR*SDT
            else
                Y=(AI0(I)/SR/WBBS0)**2
                B2B0EQ(I)=B2B0EQ(I)+Y*SDT
                B0B2EQ(I)=B0B2EQ(I)+ SDT/Y
                BMODEQ(I)=BMODEQ(I)+ SQRT(Y)*SDT
                IF (Y.GT.YMAX) YMAX=Y
		Y1 = SQRT(Y/YMAX)
        Y1 = 1.d0/Y*(1.d0-SQRT(abs(1.d0-Y1))*(1.d0+.5d0*Y1))*SDT
                FOFBEQ(I)=FOFBEQ(I)+Y1
                GRDAEQ(I)=GRDAEQ(I)+SQRT(G22A2/DA2)*SDT
             endif
40        continue
          BMAXEQ(I)=SQRT(YMAX)*WBBS0
          BMINEQ(I)=SQRT(YMIN)*WBBS0
          if(I.NE.1)then
             B0B2EQ(I)=B0B2EQ(I)/SKGGG
             B2B0EQ(I)=B2B0EQ(I)/SKGGG
             BMODEQ(I)=BMODEQ(I)/SKGGG
             FOFBEQ(I)=FOFBEQ(I)/SKGGG
             GRDAEQ(I)=GRDAEQ(I)/SKGGG
          endif
 10   continue
CMR
	YLIN=YLIN*WBBS0*WBBS0*A*2./(.4d0*cgp*PLCUR)**2/BR00
C	open(33,file='dat/lin3')
C	write(33,*) YLIN
!,yvol*A
C	close(33)
 30   format(5(f10.4))
      END

C - 3 MOMENT EQUILIBRIUM SOLVER			AUGUST 17,1988 
C NITER - max number of iterations
C ACC   - relative tolerance parameter
      subroutine EQAB3(NA,NT,NITER,ACC)
	implicit none
	include	'for/emeq.inc'
	integer	NA,NT,NITER,NA1,I,J,I1,JJ
	double precision
     1		ACC,AA,AAI,A4AI,A6AI,S,SC2,W0,W1,W2,W3,BI0,BP1,BP3,
     2		FU0,FU0I,FU0J,FU1I,FU1J,FU2I,FU2J,FU3I,FU3J,
     3		FV0I,FV0J,FV1I,FV1J,FV2I,FV2J,FV3I,FV3J
	NA1=	NA+1
	ITER=	0
  200	ITER=	ITER+1
	WSAC1=	WDSD1(NA1)
	WSAC2=	WDGL(NA1)
	WSAC3=	WDSD3(NA1)
	print *, 'EQAB3 NA1=', NA1
	print *, WBR00, WSD1(NA1), WSAA(NA1)
	WBR0=	WBR00+WSD1(NA1)*WSAA(NA1)
      call add2loc("Calling EQK3"//char(0))
      CALL EQK3(NA,NT)
      call add2loc("Calling EQALVU3"//char(0))
      CALL EQLVU3(NA)
      call add2loc("Calling EQC1"//char(0))
      CALL EQC1(NA)
C  GMC,SW0,SDD1,GDL
	WGMC(1)=WBA(1)*WSL0(1)/WBK0(1)
      BI0=	0.d0
	W0=	0.d0
	FU0I=	WDBA(1)*WSL0(1)
	FV0I=	WDBB(1)*WSV0(1)
	FU0=	WGMC(1)*WSU0(1)
	WSW0(1)=0.25d0*(FU0-FU0I)
      W2=	0.d0
	FU2I=	WGMC(1)*WSU2(1)-WDBA(1)*WSL2(1)
	FV2I=	WDBB(1)*WSV2(1)
	SC2=	0.25d0*(WGL(1)**2-1.d0)
	WBS2(1)=((WBA(1)*WSL22(1)+WBB(1)*(WSV2(1)+SC2*WSV0(1))+
     ,		FU2I/6.d0+SC2*WSW0(1))/WGMC(1)-WBK20(1))/WBK22(1)
	WDGL(1)=0.d0
      W1=	0.d0
	FU1I=	WGMC(1)*WSU1(1)
	FV1I=	WDBA(1)*WSL1(1)+WDBB(1)*WSV1(1)
	BP1=	(WBA(1)*WSL1(1)+WBB(1)*WSV1(1)+
     *	 	0.25d0*FU1I)/WGMC(1)-WBK10(1)
      W3=	0.d0
	FU3I=	WGMC(1)*WSU3(1)
	FV3I=	WDBA(1)*WSL3(1)+WDBB(1)*WSV3(1)
	BP3=	(WBA(1)*WSL3(1)+WBB(1)*WSV3(1)+
     *	 	FU3I/6.d0)/WGMC(1)-WBK30(1)
	S=	1.d0/(WBK11(1)*WBK33(1)-WBK13(1)*WBK31(1))
	WBS1(1)=(BP1*WBK33(1)-BP3*WBK13(1))*S
	WBS3(1)=(BP3*WBK11(1)-BP1*WBK31(1))*S
	WDSD1(1)=WBS1(1)
      DO 1 I=2,NA1
	J=	I-1
	AA=	WSAA(I)
	AAI=	1.d0/AA
	A4AI=	AAI*AAI
	A6AI=	A4AI*AAI
      FU0J=	FU0I
	FU0I=	WDBA(I)*WSL0(I)
	FV0J=	FV0I
	FV0I=	WDBB(I)*WSV0(I)
	BI0=	BI0+WSCJ3(J)*FU0J+WSCI3(J)*FU0I+
     * 	    	WSCJ5(J)*FV0J+WSCI5(J)*FV0I
	W0=	W0+WSCJ3(J)*FU0
	WGMC(I)=(WBA(I)*WSL0(I)+WBB(I)*AA*WSV0(I)+(W0-BI0)*AAI)
     ,		/(WBK0(I)-WSCI3(J)*WSU0(I)*AAI)
	FU0=	WGMC(I)*WSU0(I)
	W0=	W0+WSCI3(J)*FU0
	WSW0(I)=(W0-BI0)*A4AI
      FU2J=	FU2I
	FU2I=	WGMC(I)*WSU2(I)-WDBA(I)*WSL2(I)
	FV2J=	FV2I
	FV2I=	WDBB(I)*WSV2(I)
	W2=	W2+FU2J*WSCJ5(J)+FU2I*WSCI5(J)-
     *		FV2J*WSCJ7(J)-FV2I*WSCI7(J)
	SC2=	0.25d0*(WGL(I)**2-1.d0)
	WBS2(I)=((WBA(I)*WSL22(I)+WBB(I)*(WSV2(I)+SC2*WSV0(I))+
     ,		W2*A6AI+SC2*WSW0(I))/WGMC(I)-WBK20(I))/WBK22(I)
        WDGL(I)=WDGL(J)+(WBS2(J)+WBS2(I))*WSCJ1(J)
      FV1J=	FV1I
	FV1I=	WDBA(I)*WSL1(I)+WDBB(I)*WSV1(I)
	FU1J=	FU1I
	FU1I=	WGMC(I)*WSU1(I)
	W1=	W1+FU1J*WSCJ3(J)+FU1I*WSCI3(J)-
     *		FV1J*WSCJ5(J)-FV1I*WSCI5(J)
	BP1=	(WBA(I)*WSL1(I)+WBB(I)*WSV1(I)+
     *	 	W1*A4AI)/WGMC(I)-WBK10(I)
      FU3J=	FU3I
	FU3I=	WGMC(I)*WSU3(I)
	FV3J=	FV3I
	FV3I=	WDBA(I)*WSL3(I)+WDBB(I)*WSV3(I)
	W3=	W3+FU3J*WSCJ5(J)+FU3I*WSCI5(J)-
     *		FV3J*WSCJ7(J)-FV3I*WSCI7(J)
	BP3=	(WBA(I)*WSL3(I)+WBB(I)*WSV3(I)+
     *	 	W3*A6AI)/WGMC(I)-WBK30(I)
	S=	1.d0/(WBK11(I)*WBK33(I)-WBK13(I)*WBK31(I))
	WBS1(I)=(BP1*WBK33(I)-BP3*WBK13(I))*S
	WBS3(I)=(BP3*WBK11(I)-BP1*WBK31(I))*S
	WDSD1(I)=WBS1(I)
 1	CONTINUE
C  SD1,GL
      print *, 'WDSD1 1 NA1', WDSD1(1), WDSD1(NA1)
	W1=	0.d0
	WSD1(1)=0.5d0*WDSD1(1)
	S=	WGL(NA1)*EXP(-WDGL(NA1))
C	write(*,*)S,WGL(1),WGL(NA1),WDGL(NA1),EXP(-WDGL(NA1)),WDSD1(1)
	WGL(1)=	S
	WDGL(1)=WBS2(1)
      DO 2 I=2,NA1
	J=	I-1
        WGL(I)=	S*EXP(WDGL(I))
	WDGL(I)=WBS2(I)
	W1=	W1+(WDSD1(J)+WDSD1(I))*WSCJ1(J)
	WSD1(I)=W1/WSAA(I)
 2	CONTINUE
C  SD3,SDD3
     	FU3J=	1.d0+WDGL(NA1)*WSAA(NA1)
	WDSD3(NA1)=WBS3(NA1)+2.d0*WSD3(NA1)*FU3J
      DO 3 I1=1,NA
	J=	NA1-I1
	I=	J+1
	AA=	WSAA(J)
	FU3I=	FU3J
     	FU3J=	1.d0+WDGL(J)*AA
	WSD3(J)=(WSD3(I)*(WSAA(I)-2.d0*WSCI1(J)*FU3I)-(WBS3(I)+
     ,		WBS3(J))*WSCJ1(J))/(AA+2.d0*WSCJ1(J)*FU3J)
	WDSD3(J)=WBS3(J)+2.d0*WSD3(J)*FU3J
 3	CONTINUE
      print *, 'WDSD1 1 NA1',WDSD1(1),WDSD1(NA1)
	WSAC1=	WSA(NA1)*ABS(WDSD1(NA1)-WSAC1)
	WSAC2=	WSAA(NA1)*ABS(WDGL(NA1)-WSAC2)
	WSAC3=	WSA(NA1)*ABS(WDSD3(NA1)-WSAC3)
	WSACC=	WSAC1+WSAC2+WSAC3
      print *, 'acc', WSACC, ACC, ITER
      IF(ITER.LT.NITER.AND.WSACC.GT.ACC) GO TO 200

C	write(*,'(2I5,5F10.7)')ITER,NITER,WSACC,WSAC1,WSAC2,WSAC3
        if (ITER.EQ.NITER)	then
	   NA = 0	! Suppress optional output
	   return
	   write(*,'(1A)')"        You can: (0) Continue anyway"
	   write(*,'(1A)')"                 (1) Try more iterations"
	   write(*,'(2A)')"                 (2) Freeze the current ",
     >			"equilibrium and continue"
	   write(*,'(1A)')"                 (3) Stop the run"
	   write(*,'(A,$)')"     Your selection > "
	   read(*,*)jj
	   if (jj .ge. 3) call	a_stop
	   if (jj .eq. 0)	goto	4
	   ITER = 0
	   if (jj .eq. 1)	goto	200
	   if (jj .eq. 2)	NA = 0
	   return
	endif

C  SD2D1,GD2L,SD2D3
 4	WD2SD1(1)=WDSD1(1)
	FU1I=	0.d0
	WD2GL(1)=WDGL(1)*WGL(1)
	FU2I=	-WGL(1)
	WD2SD3(1)=WDSD3(1)
	FU3I=	0.d0
      DO 5 I=2,NA1
	J=	I-1
	AA=	WSAA(I)
	FU1J=	FU1I
	FU1I=	AA*(WDSD1(I)-WSD1(I))
	WD2SD1(I)=(FU1I-FU1J)/WSCI1(J)-WD2SD1(J)
	FU2J=	FU2I
	FU2I=	WGL(I)*(AA*WDGL(I)-1.d0)
	WD2GL(I)=(FU2I-FU2J)/WSCI1(J)-WD2GL(J)
	FU3J=	FU3I
	FU3I=	AA*(WDSD3(I)-WSD3(I))
	WD2SD3(I)=(FU3I-FU3J)/WSCI1(J)-WD2SD3(J)
 5	CONTINUE
	RETURN
	END
C======================================================================|

C EQGB3.FOR - INITIAL CONDITIONS FOR THREE MOMENT EQUILIBRIA	AUGUST 17,1988

	BLOCK DATA GB3
	implicit none
	include	'for/emeq.inc'
	data	WSD1/NP*0.d0/,WDSD1/NP*0.d0/,WDGL/NP*0.d0/
	data	WBS1/NP*0.d0/,WBS2/NP*0.d0/,WBS3/NP*0.d0/
	end

      subroutine EQGB3(BR00,SA0,GL0,GD30,NA)
C----------------------------------------------------------------------|
C Parameters used:
C Input:
C	Variables:	NA,BR00,SA0,GL0,GD30
C	Arrays:		
C Output
C	Arrays defined here and used elsewhere
C   WSA  = j*SA0/NA	current minor radius, a
C   WSAA = WSA**2	a^2
C   WSCI1 = 
C   WSCI3 = 
C   WSCI5 = 
C   WSCI7 = 
C   WSCJ1 = 
C   WSCJ3 = 
C   WSCJ5 = 
C   WSCJ7 = 
C	Arrays first defined here and modified in EQAB3
C   WGL  = \lambda
C   WSD3 = \delta/a^2
C   WDSD3
C----------------------------------------------------------------------|
	implicit none
	include	'for/emeq.inc'
	integer	NA,NA1,I,J
	double precision
     1	BR00,SA0,GL0,GD30,DA,D3,D32,AAI,AAJ,H,C0,CI,CJ,CII,CJJ
C	double precision	WBR0,WBR00
C	double precision	WSA(1),WSAA(1),WGL(1),WSD3(1),WDSD3(1)
C	double precision	WSCI1(1),WSCI3(1),WSCI5(1),WSCI7(1)
C	double precision	WSCJ1(1),WSCJ3(1),WSCJ5(1),WSCJ7(1)
	WBR00=BR00
	WBR0=BR00
        DA=SA0/NA
        D3=GD30/SA0**2
	D32=2.d0*D3
        NA1=NA+1
      DO 1 I=1,NA1
	J=I-1
	WSA(I)=DA*J
	WSAA(I)=WSA(I)**2
	WGL(I)=GL0
	WSD3(I)=D3
	WDSD3(I)=D32
C---C(J,N),C(I,N)
      IF(J.EQ.0) GO TO 1
	AAI=WSAA(I)
	AAJ=WSAA(J)
	H=(AAI-AAJ)*0.5d0
	WSCJ1(J)=H*0.5d0
	WSCI1(J)=WSCJ1(J)
	CJ=AAI+2.d0*AAJ
	CI=AAJ+2.d0*AAI
	C0=H/6.
	WSCJ3(J)=C0*CJ
	WSCI3(J)=C0*CI
C---CJJ=AJ**2N
	CJJ=AAJ*AAJ
C---CII=AI**2N
	CII=AAI*AAI
	CJ=AAI*CJ+3.d0*CJJ
	CI=AAJ*CI+3.d0*CII
	C0=C0*0.5d0
	WSCJ5(J)=C0*CJ
	WSCI5(J)=C0*CI
	C0=H*0.05d0
	WSCJ7(J)=C0*(AAI*CJ+4.d0*CJJ*AAJ)
	WSCI7(J)=C0*(AAJ*CI+4.d0*CII*AAI)
 1	CONTINUE
	RETURN
	END


C - EQK3.FOR - LEFT HAND SIDE AVERAGER		AUGUST 17,1988

      subroutine EQK3(NA,NT)
C----------------------------------------------------------------------|
C Parameters used:
C Input:
C	Variables:	NA,WBR0
C	Arrays:		
C Output (arrays):
C			
C----------------------------------------------------------------------|
	implicit none
	include	'for/emeq.inc'
	integer	NA,NT,NA1,NT1,I,K
	double precision
     1		SKDR,SKGA,SQG22R,CGP,A,AA,ASPA,EE,SC1,SC2,T,DEN,
     2		DA1,DA2,RA,RI,XR,G22A2,R0RD,SF2,SF20,SF21,SF3,SF30,SF31,
     3		XRXR,SK10,SK11,SK13,SK0,SK01,SK2,SR,SX,SXX,SX1,SXX1,SZZ,
     4		S,SS,SY1,C,CC,BM1,BM2,BMI,BN,BN0,BN1,BN2,BN12,SDT,SDT0
	common /EMEQMR/SKDR(NP),SKGA(NP),SQG22R(NP)
	save	CGP
      data CGP/3.14159265359d0/
	NA1=	NA+1
	NT1=	NT+1
	SDT0=	1.d0/NT
	print *, 'SDT0=', SDT0
	RI=	1.d0/WBR0
      DO 1 I=1,NA1
	A=	WSA(I)
	AA=	A*A
	ASPA=	AA*RI
	EE=	WGL(I)**2
	SC1=	0.25d0*(EE-1.d0)
	SC2=	0.5d0*(EE+1.d0)
	WBK02(I)=0.d0
	WBG332(I)=0.d0
	WBG222(I)=0.d0
	WSU1(I)=0.d0
	WSU2(I)=0.d0
	WBK20(I)=0.d0
	WBK22(I)=0.d0
	WBK10(I)=0.d0
	WBK11(I)=0.d0
	WBK13(I)=0.d0
	WBK30(I)=0.d0
	WBK31(I)=0.d0
	WBK33(I)=0.d0
C------------------------------------------
C KONOVALOV LIKES TO DO SOMETHING BY HIMSELF
C------------------------------------------
	SKGA(I)=0.d0
        SKDR(I)=0.d0
	SQG22R(I)=0.d0
C------------------------------------------
      DO 2 K=1,NT1
	T=	SDT0*CGP*(K-1)
	SDT=	SDT0
	IF(K.EQ.1.OR.K.EQ.NT1) SDT=SDT0*0.5d0
C--------------------------------------
	C=	COS(T)
	S=      SIN(T)
	print *, 'cos sin', C, S, T
	SS=	S*S
	CC=	C*C
	SX1=	-WSD1(I)-WSD3(I)*SS
	SXX1=	SX1*SX1
	SZZ=	EE*SS
	SX=	C+A*SX1
	SXX=	SX*SX
	SR=	WBR0+A*SX
	SY1=	2.d0*C*WSD3(I)
C---BN0=N0*DT
	BN0=	(EE*CC+SS)*SDT
C---BN1=N1*DT
	BN1=	2.d0*SY1*SS*SDT
C---BN2=N2*DT
	BN2=	SS*SY1*SY1*SDT
	BN12=	BN1+A*BN2
	BN=	BN0+A*BN12
	BM1=	-(WBS1(I)+WBS3(I)*SS)*C
	BM2=	WBS2(I)*SS
	BMI=	1.d0/(1.d0+A*BM1+AA*BM2)
C  EVEN EQUATIONS
	DEN=	1.d0/(1.d0+A*BM1)
	SK01=	BN1-BN0*BM1
	SK0=	(BN2-BN1*(BM1+A*BM2))*BMI+BN0*BM1*BM1*DEN
	SK2=	-BN0*SS*DEN*BMI
	SF20=	CC-SZZ+SC1
	SF21=	2.d0*C*SX1
	SF2=	SF20+A*SF21+AA*SXX1
	WBK02(I)=WBK02(I)+SK0+WBS2(I)*SK2
	WBK20(I)=WBK20(I)+BN0*SXX1+SK01*SF21+SK0*SF2
	WBK22(I)=WBK22(I)+SK2*SF2
C  ODD EQUATIONS
	DEN=	1.d0/(1.d0+AA*BM2)
	SK0=	BN0*DEN
	SK10=	BN12*BMI
	SK11=	SK0*C*BMI
	SK13=	SK11*SS
	SF3=	CC-3.d0*SZZ
	SF30=	C*SF3
	SF31=	SX1*(SXX+C*SX+SF3)
	SF3=	SF30+A*SF31
	WBK10(I)=WBK10(I)+SK0*SX1+SK10*SX
	WBK30(I)=WBK30(I)+SK0*SF31+SK10*SF3
	WBK11(I)=WBK11(I)+SK11*SX
	WBK31(I)=WBK31(I)+SK11*SF3
	WBK13(I)=WBK13(I)+SK13*SX
	WBK33(I)=WBK33(I)+SK13*SF3
        XR=	SX*RI
	XRXR=	XR*XR
	R0RD=	WBR0/SR
	WBG222(I)=WBG222(I)+BN*BMI*XRXR*R0RD
	R0RD=R0RD*SDT
	WBG332(I)=WBG332(I)+(XRXR+BM2-BM1*XR)*R0RD
	R0RD=	R0RD*SX*C
	WSU1(I)=WSU1(I)+R0RD*SZZ
	WSU2(I)=WSU2(I)+R0RD*SXX
C HERE ZAKHAROV'S AUTHOR RIGHTS ARE CANCELED WITHOUT ANY CEREMONY
c	DRDA=-WDSD1(I)*A+C-WDSD3(I)*A*SS
c	DZDA=S*WGL(I)*(AA*WDGL(I)+1.d0)
c	DRDT=-A*S-2.d0*AA*WSD3(I)*C*S	
c	DZDT=WGL(I)*A*C
C METRIC (SUBSRIPT) TENSOR COMPONENTS
c	G11=DRDA**2+DZDA**2
c	G22=DRDT**2+DZDT**2
c	G12=DRDA*DRDT+DZDA*DZDT
c G22/A**2
        G22A2=SS+4.d0*A*WSD3(I)*SS*C+
     &	(2.d0*A*WSD3(I)*S*C)**2+(WGL(I)*C)**2
C D/A
        DA1=WGL(I)*(CC-A*(WDSD1(I)+WDSD3(I)*SS)*C+
     +       SS*(WDGL(I)*AA+1.d0)*(1.d0+2.d0*A*C*WSD3(I)))
C---- ASTRA METRIC COMBINATIONS -------
        SKGA(I)=SKGA(I)+G22A2*SR*SDT/DA1
        SQG22R(I)=SQG22R(I)+sqrt(G22A2)*SR*SDT
C D**2/A**2
        DA2=(C-A*(WDSD1(I)+WDSD3(I)*SS))*WGL(I)*C+
     +  WGL(I)*SS*(WDGL(I)*AA+1.d0)*(1.d0+2.d0*A*C*WSD3(I))
        SKDR(I)=SKDR(I)+SDT*DA2/SR      
 2    CONTINUE
      WBG332(I)=WBG332(I)+(WSD1(I)+0.5d0*WSD3(I))*RI
	WBG222(I)=WBG222(I)+WBK02(I)-(WBK10(I)+
     *	WBS1(I)*WBK11(I)+WBS3(I)*WBK13(I))*RI
	WBD02(I)= 0.5d0*WBS2(I)
	RA=	  1.d0-WSD1(I)*AA*RI
	WBG02(I)=0.5d0*WBS2(I)*(RA-0.75d0*WSD3(I)*ASPA)-
     *	(WSD1(I)+0.5d0*(WSD3(I)+WBS1(I)+0.25d0*WBS3(I)))*RI
	WBD12(I)= WBG02(I)-WBG332(I)
	WBK0(I)=  SC2+AA*WBK02(I)
	WDBK00(I)=EE*WDGL(I)
	WBD0(I)=  1.d0+AA*WBD02(I)
	WBG0(I)=  1.d0+AA*WBG02(I)
	WBG33(I)= 1.d0+AA*WBG332(I)
	WBG22(I)= SC2+AA*WBG222(I)
 1	CONTINUE
	END
C - EQLVU3.FOR - RIGHT HAND SIDE AVERAGER	AUGUST 17,1988

 	subroutine EQLVU3(NA)
C----------------------------------------------------------------------|
C Parameters used:
C Input:
C	Variables:	NA,WBR0
C	Arrays:	
C   WSAA(*)
C   WGL(*)
C   WSD1(*)
C   WSD3(*)
C Output
C	Arrays defined here and used elsewhere
C   WSL0
C   WSL1
C   WSL2
C   WSL3
C   WSL22
C   WSV0
C   WSV1
C   WSV2
C   WSV3
C	Arrays modified here, first defined in EQK3, used in EQAB3
C   WSU0
C   WSU1
C   WSU2
C   WSU3
C----------------------------------------------------------------------|
	implicit none
	include	'for/emeq.inc'
C	double precision	WBR0,WSAA(1),WGL(1),WSD1(1),WSD3(1)
C	double precision	WSL0(1),WSL1(1),WSL2(1),WSL3(1),WSL22(1)
C	double precision	WSV0(1),WSV1(1),WSV2(1),WSV3(1)
C	double precision	WSU0(1),WSU1(1),WSU2(1),WSU3(1)
	integer	NA,NA1,I
	double precision
     1		CR2,CR3,CR4,CR6,CR8,AA,E,EE,T1,T3,T5,RBR0,
     2		S,S1,S2,X20,X22,X30,X32,X40,X42,X50,X60,X5X1,X6X1,
     3		Y0,Y2,Y4,UX0,UX2,UX12,UX30,UT0,UT1,UT2,UT4,UT6,
     4		UY0,UY1,UY2,UY3,UY4,UZ1,UZ3
	data CR2 /5.d-1/,CR3/3.33333333d-1/,CR6/1.66666666d-1/,
     ,		CR4/2.5d-1/,CR8/1.25d-1/
	NA1=	NA+1
	RBR0=1.d0/WBR0
      do 1 I=1,NA1
        AA=	WSAA(I)
        E=	WGL(I)
        EE=	E*E
	UX2=	CR2*WSD3(I)
	UX0=	-WSD1(I)-UX2
C-------------------------------
	UY4=	CR2*UX2*UX2
	UY0=	UX0*UX0+UY4
	Y0=	CR2+UY0*AA
	UY1=	2.d0*UX0+UX2
	UY2=	2.d0*UX0*UX2
	Y2=	CR2+UY2*AA
	UY3=	UX2
	Y4=	UY4*AA
	UT0=	UX0*Y0+CR2*(UY1+UX2*Y2)
	UT1=	UX0*UY1+UY0+CR2*(UY2+UX2*(UY3+UY1))
	T1=	0.75d0+UT1*AA
	UT2=	UX0*Y2+Y0*UX2+CR2*(UY3+UY1+UX2*Y4)
	T3=	AA*UX0*UY3+CR2*(Y2+Y4+AA*UX2*UY1)
	UT4=	UX0*Y4+CR2*(UY3+UX2*Y2)
	T5=	CR2*(Y4+AA*UX2*UY3)
	UT6=	CR2*UX2*Y4
	UZ1=	UY1*(2.d0*Y0+Y2)+UY3*(Y2+Y4)
	UZ3=	UY3*2.d0*Y0+UY1*(Y4+Y2)
	X5X1=	T1*Y0+AA*UT0*UY1+CR2*(Y2*(T1+T3)+
     ,		Y4*(T3+T5)+AA*(UT2*UY1+UY3*(UT2+UT4)))
	X6X1=	T1*(2.d0*UT0+UT2)+T3*(UT2+UT4)+T5*(UT4+UT6)
      S=	EE*CR8
	X20=	CR4*UY1
	X22=	S*(UY1-UY3)*CR2
	X30=	CR6*T1
	X32=	S*(T1-T3)*CR3
	UX30=	CR6*UT1
	X40=	CR8*UZ1
	X42=	S*(UZ1-UZ3)*0.25d0
	X50=	X5X1*0.1d0
	X60=	CR2*X6X1*CR6
C-------------------------------
	S1=(1.d0-EE)*CR4
        WSL0(I)=	E*CR2
	WSL22(I)=	E*UX30
	WSL2(I)=	WSL0(I)*S1+AA*WSL22(I)
	WSL1(I)=	E*X20
	WSL3(I)=	E*(X40-3.d0*X22)
	S=	2.d0*E*RBR0
	S1=	CR2*RBR0
	S2=	AA*S1
        WSV0(I)=	S*(X20+S1*X30)
        WSV1(I)=	S*(X30+S2*X40)
	WSV2(I)=	S*(X40-X22+S1*(X50-X32))
	WSV3(I)=	S*(X50-3.d0*X32+S2*(X60-3.d0*X42))

	UX30=	WSU2(I)
        UX12=	WSU1(I)
	S=	EE*RBR0
	S2=	AA*RBR0
        WSU1(I)=	S*(CR2+S2*(UX30*RBR0-2.d0*X20))
        WSU0(I)=	-RBR0*WSU1(I)
        WSU2(I)=	S*(2.d0*X20+RBR0*(UX12-UX30))
        WSU3(I)=	S*(UX30-3.d0*UX12)
 1	CONTINUE
	RETURN
	END

	subroutine EQC1(NA)
C----------------------------------------------------------------------|
C Parameters used:
C Input:
C	Variables:	NA,WBR0,WBR00
C	Arrays:		WSA,WSP,WSJP
C Output:
C      Arrays:
C   WBA   Zakharov's function A
C   WBB   Zakharov's function B
C   WDBA  h/a*dA/da
C   WDBB  
C----------------------------------------------------------------------|
	implicit none
	include 'for/emeq.inc'
C	double precision WBR0,WBR00
C	double precision WSA(1),WSP(1),WSJP(1),WBA(1),WBB(1),WDBA(1),WDBB(1)

	integer	NA,NA1,I
	double precision	H, S, RS, SS
	NA1=NA+1
	WDBA(1)=2.d0*(WSJP(2)-WSJP(1))/WSA(2)**2
	WBA(1)=WSJP(1)
	WBB(1)=WSP(1)
	WDBB(1)=2.d0*(WSP(2)-WSP(1))/WSA(2)**2
	  H=0.5d0/WSA(2)
	  do  I=2,NA
	    WBA(I)=WSJP(I)
	    WBB(I)=WSP(I)
	    WDBB(I)=(WSP(I+1)-WSP(I-1))*H/WSA(I)
    	WDBA(I)=(WSJP(I+1)-WSJP(I-1))*H/WSA(I)
      end do
	WBA(NA1)=WSJP(NA1)
	WDBA(NA1)=(WSJP(NA1)-WSJP(NA))/WSA(2)/WSA(NA1)
	WBB(NA1)=WSP(NA1)
	WDBB(NA1)=(WSP(NA1)-WSP(NA))/WSA(2)/WSA(NA1)
	S=WBR0/WBR00
	RS=1.d0/S
	SS=S-RS
      do I=1,NA1
	    WBA(I)=WBA(I)*RS+WBB(I)*SS
	    WBB(I)=WBB(I)*S
	    WDBA(I)=WDBA(I)*RS+WDBB(I)*SS
    	WDBB(I)=WDBB(I)*S
      end do
	end

	subroutine EQPPAB(NA)
C----------------------------------------------------------------------|
C	gm0 = 0.4d0*cgp
C	WSA(I) = a[m], WSAA(I) = a**2
C	WBR00 = (RTOR+SHIFT) = R[m]
C	WBR0 = R0[m]
C	WBJ0 = 0.2d0*I[MA]
C	WBBS0 = Bs[T],
C	WBF(I)=0.2d0*F/R0[MA/m], WBFF(I) = WBF(I)**2
C	WSJ(I)=gm0*<j>[MA/m**2], WSJP(I) = <jB>/B
C	WSJSL(I) = j(r,gt=cgp), WSJSR(I) = j(r,gt=0)
C	WSP(I) = gm0*p[MJ/m**3]
C	WGP(I) = gP[VS]/(2*cgp*R0)
C	WDSQRQ(I) = aq'/q
C   WGPINT, WGPRES, WSLI are defined but not used (not present in commons
C----------------------------------------------------------------------|
C Parameters used:
C Input:
C	Variables:	NA,WBR0,WBR00,WBJ0,WBBS0
C	Arrays:	
C   WSA(*)
C   WSAA(*)
C   WGL(*)
C   WSD1(*)
C   WSD3(*)
C Output
C     Variables:
C   WGB
C   WGB0
C   WGBD
C   WGBJ
C   WGBST
C   WGMJ
C   WGMJEX
C   WSQC
C   WGPINT
C   WGPRES
C   WSLI
C	Arrays:
C   WSL0
C   WSL1
C   WSL2
C   WSL3
C   WSL22
C   WSV0
C   WSV1
C   WSV2
C   WSV3
C----------------------------------------------------------------------|
	implicit none
	include	'for/emeq.inc'
C	double precision WBR0,WBR00,WBJ0,WBBS0
C	double precision WGB,WGB0,WGBD,WGBJ,WGBST,WGMJ,WGMJEX,WSQC
C	double precision WSP(1),WSA(1),WSAA(1),WGL(1),WGP(1),WSJ(1),WSJP(1)
C	double precision WBA(1),WBB(1),WBF(1),WBFF(1),WSD1(1)
C	double precision WDBA(1),WDBB(1),WDGL(1),WBK0(1),WBK02(1),WDBK00(1)
C	double precision WDSP(1),WDSJ(1),WSQ(1),WDSQ(1),WDSQRQ(1),WDSJP(1)
C	double precision WBD0(1),WBD12(1),WBD02(1),WBG22(1),WBG33(1),WBG332(1)
C	double precision WSCI3(1),WSCJ3(1),WSCJ1(1)
C	double precision WSL0(1),WSU0(1),WSW0(1),WSV0(1),WGMC(1),WDGMC(1)
C	double precision WSJSL(1),WSJSR(1)
	integer	NA,NA1,I,J
	double precision
     1		S,SLI,SLJ,G33I,G33J,CGP,DG33,BJ,BJI,BJJ,RR,RL,
     2		DL0,DL0I,DL0J,D12I,D12J,DD12,DV0,FF0,FFI,FFJ,V0I,V0J,
     3		GB,GBI,GBJ,GBS,GBSI,GBSJ,GBD,GBDI,GBDJ,GPI,GPJ,GMJI,GMJJ,
     4		WGPINT,WGPRES,WSLI
        save CGP
	data CGP/3.14159265359d0/
	NA1=NA+1
C --- GP,SJ,DSJ,BJ,SJL,SJR

	G33I=0.d0
	DG33=2.d0*(WBA(1)*(WBG332(1)-WBD02(1))+WBB(1)*WBD12(1))
	WSJ(1)=WBA(1)
	WDSJ(1)=WDBA(1)+DG33
	BJI=WSJ(1)*WGL(1)
	BJ=0.d0

	DO 1 I=2,NA1
	J=I-1
	G33J=G33I
	G33I=(WBA(I)*(WBG332(I)-WBD02(I))+WBB(I)*WBD12(I))*
     ,	WSAA(I)/WBD0(I)
	DG33=(G33I-G33J)/WSCJ1(J)-DG33
	WSJ(I)=WBA(I)+G33I
	WDSJ(I)=WDBA(I)+DG33
	BJJ=BJI
	BJI=WSJ(I)*WBD0(I)*WGL(I)
	BJ=BJ+(BJJ+BJI)*WSCJ1(J)
1	CONTINUE

	S=WBJ0/BJ

	DO 2 I=1,NA1
	WBA(I)	=WBA(I)*S
	WDBA(I)	=WDBA(I)*S
	WBB(I)	=WBB(I)*S
	WDBB(I)	=WDBB(I)*S
	WSJ(I)	=WSJ(I)*S
	WDSJ(I)	=WDSJ(I)*S
	WGMC(I)	=WGMC(I)*S
	WSW0(I)	=WSW0(I)*S

	RR=WBR0-WSD1(I)*WSAA(I)
	RL=WBR0/(RR-WSA(I))
	RR=WBR0/(RR+WSA(I))
	WSJSL(I)=WBA(I)*RL+WBB(I)*(1.d0/RL-RL)
	WSJSR(I)=WBA(I)*RR+WBB(I)*(1.d0/RR-RR)

	WDSP(I)=-WBB(I)*WGMC(I)*WGL(I)
2	CONTINUE

	GPI=WGMC(1)*WGL(1)
	WGP(1)=0.d0
	SLI=WGP(1)*WSJ(1)*WGL(1)*WBD0(1)
	WSLI=0.d0

	FFI=(WBA(1)-WBB(1))*GPI
	WBFF(1)=0.d0
	WSP(1)=0.d0

	GBI=WBB(1)*GPI*WSL0(1)
	GB=0.d0
	GBSI=0.d0
	GBS=0.d0
	GBDI=(WBA(1)-WBB(1))*WGMC(1)*WSU0(1)
	GBD=0.d0

	DO 3 I=2,NA1
	J=I-1

	GPJ=GPI
	GPI=WGMC(I)*WGL(I)
	WGP(I)=WGP(J)-(GPJ+GPI)*WSCJ1(J)
	SLJ=SLI
	SLI=WGP(I)*WSJ(I)*WGL(I)*WBD0(I)
	WSLI=WSLI+(SLJ+SLI)*WSCJ1(J)

	FFJ=FFI
	FFI=(WBA(I)-WBB(I))*GPI
	WBFF(I)=WBFF(J)+(FFJ+FFI)*WSCJ1(J)
	WSP(I)=WSP(J)+(WDSP(J)+WDSP(I))*WSCJ1(J)

	GBJ=GBI
	GBI=WBB(I)*GPI*WSL0(I)
	GB=GB+GBJ*WSCJ3(J)+GBI*WSCI3(J)

	GBSJ=GBSI
	GBSI=WSP(I)*GBI
	GBS=GBS+GBSJ*WSCJ3(J)+GBSI*WSCI3(J)

	GBDJ=GBDI
	GBDI=(WBA(I)-WBB(I))*WGMC(I)*WSU0(I)
	GBD=GBD+GBDJ*WSCJ3(J)+GBDI*WSCI3(J)
3	CONTINUE
	WGPINT=-2.d0*CGP*WBR0*WGP(NA1)
c
c inernal iductance
c
	WSLI=2.d0*WBR0/(WBR00*WBJ0**2)*(WSLI-WGP(NA1)*WBJ0)
	WGPRES=CGP*WBR00*WSLI*WBJ0
c
c betta j
c
	WGBJ=4.d0*GB/(WBJ0**2)
	S=WSL0(NA1)*WSAA(NA1)
c
c betta
c
	WGB=2.d0*GB/(WBBS0**2*S)
c
c betta*
c
	WGBST=2.d0*SQRT(2.d0*(GBS-WSP(NA1)*GB)/S)/WBBS0**2
	WGBD=2.d0*GB/(WBR0**2*(2.d0*GBD-FF0*WSU0(NA1)*WSAA(NA1)/WGL(NA1)))
	WGMJ=4.d0*GBD*(WBR0/WBJ0)**2
	WSP(1)=-WSP(NA1)
	FF0=(WBBS0*WBR00/WBR0)**2
	WBFF(1)=FF0+2.d0*WBFF(NA1)
	WBF(1)=SQRT(WBFF(1))
	WGB0=2.d0*WSP(1)/WBBS0**2
	WSQC=WBBS0*WSAA(NA1)*(WGL(NA1)**2+1.d0)*0.5d0/(WBJ0*WBR00)

C --- SJP,SDJP,SP,BFF,DGMC,SQ,SDQ,GMJ,GMJEX
	DV0=2.d0*(WBB(1)*WSV0(1)+WSW0(1)-WGMC(1)*WBK02(1))
	V0I=0.d0
	DL0I=(WDGL(1)-2.d0*WBD02(1))*WGL(1)
	DL0=0.d0
	WDGMC(1)=(WDBA(1)*WSL0(1)+WBA(1)*(WGL(1)*WBD02(1)+.25*DL0I)-
     ,	WGMC(1)*WDBK00(1)+DV0)/WBK0(1)

	FFI=(WBA(1)-WBB(1))*WGMC(1)*WGL(1)
	WSQ(1)=WBF(1)/(WBR0*WGMC(1))
	DG33=2.d0*WBG332(1)*WSQ(1)
	G33I=0.d0
	WDSQ(1)=WSQ(1)*(-FFI/WBFF(1)-WDGMC(1)/WGMC(1))+DG33
	WDSQRQ(1)=0.d0

	DD12=2.d0*((WBA(1)-WBB(1))*WBG22(1)/(WBR0*WSQ(1))**2+
     ,	WBB(1)*WBD12(1))
	D12I=0.d0
	WSJP(1)=WBA(1)
	WDSJP(1)=WDBA(1)+DD12

	BJI=WGL(1)
	BJ=0.d0
	GMJI=FFI/WBF(1)*BJI*0.5d0
	WGMJEX=0.d0

	DO 4 I=2,NA1
	J=I-1
	WSP(I)=WSP(I)-WSP(NA1)
C --- DERIVATIVES
	S=1.d0/WSCJ1(J)
	GPI=WGMC(I)*WGL(I)

	WBFF(I)=WBFF(1)-2.d0*WBFF(I)
	WBF(I)=SQRT(WBFF(I))

	DL0J=DL0I
	DL0I=(WDGL(I)-2.d0*WBD02(I))*WGL(I)
	DL0=DL0+DL0J*WSCJ3(J)+DL0I*WSCI3(J)
	V0J=V0I
	V0I=(WBB(I)*WSV0(I)+WSW0(I)-WGMC(I)*WBK02(I))*WSAA(I)
	DV0=(V0I-V0J)*S-DV0
	WDGMC(I)=(WBA(I)*(WGL(I)*WBD02(I)+DL0/WSAA(I)**2)+
     ,	WDBA(I)*WSL0(I)-WGMC(I)*WDBK00(I)+DV0)/(WGL(I)**2+1.d0)*2.

	FFI=(WBA(I)-WBB(I))*WGL(I)*WGMC(I)
	G33J=G33I
	WSQ(I)=WBF(I)/(WBR0*WGMC(I))
	G33I=WSQ(I)*WBG332(I)*WSAA(I)
	DG33=(G33I-G33J)*S-DG33
	WDSQ(I)=WSQ(I)*(-FFI/WBFF(I)-WDGMC(I)/WGMC(I))+DG33
	WSQ(I)=WSQ(I)+G33I
	WDSQRQ(I)=WSAA(I)*WDSQ(I)/WSQ(I)

	D12J=D12I
	D12I=((WBA(I)-WBB(I))*WBG22(I)*WBG33(I)/(WBR0*WSQ(I))**2+
     ,	WBB(I)*WBD12(I)/WBG33(I))*WSAA(I)
	DD12=(D12I-D12J)*S-DD12
	WSJP(I)=WBA(I)+D12I
	WDSJP(I)=WDBA(I)+DD12

	BJJ=BJI
	BJI=WBG33(I)*WGL(I)
	BJ=BJ+(BJJ+BJI)*WSCJ1(J)
	GMJJ=GMJI
	GMJI=FFI*BJ/(WBF(I)*WSAA(I))
	WGMJEX=WGMJEX+GMJJ*WSCJ3(J)+GMJI*WSCI3(J)

4	CONTINUE

	WGMJEX=4.d0*WGMJEX*WBBS0/WBJ0**2

	END
C======================================================================|
