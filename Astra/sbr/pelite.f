	Subroutine pelite(YAM,YVP,YRP)
c (not for distribution) contact: polevoa@itergps.naka.jaeri.go.jp
c Simplified Mass Ablation and Relocation Treatment	Polevoy 22-09-2000
corrected 17-01-2001, 11-MAY-2001(for G.Pacher), 18-FEB-2002 (for GVP)
c
c Polevoi,Shimada, PPCF 43 (2001) pp.1525-1533 
cRestriction: 	H,D,T (or two of them(mixed)) pellet only
c		HFS(LFS) injection in the mid plain only
c	How to use:
C		call 
C	pelite(YAM,YVP,YRP):TIME1:TIME2:
C		
C	I used YAM = 2.5 (but see quotation from PDD below), YRP = 0.5 cm, YVP = 0.5-0.3
c	it is better to choose 
C		YVP < 0.5 km/s (0.3-0.5 km/s is better since pellet is not been 
c				defragmented yet in experiments)
C		YRP < 0.3 cm (if it is not a killer pellet, then more than doubling 
c						of number of particles is not good)
c in Plant Description Document Chapt. 2.7 p. 4 it is written:
C The pellet injection system will provide a fuelling rate of 50 Pam3/s with
c 90%T/10%D pellets and 100 Pam3/s for other Hydrogenic species in the form of 3-6 mm
c diameter, sized to limit density and fusion power excursions to < 10%. The maximum 
c repetition rate varies between 7 Hz for 6 mm pellets to 50 Hz for 3 mm pellets,
c and is available for pulse length up to 3,000 s. Pellet speed is up to 0.5 km/s,
c from the HFS, are considered necessary to achieve a penetration beyond the ELM 
c affected area (~15% of minor radius)
c
C		call after!!! TE:eq,	TI:EQ,	NE:eq, NI=... (not between)
c=========================================================================
C Pellet ablation model by B.Kuteev, NF, 35 (1995) 431
c Cloud size by P.B.Parks, Phys. of Plasmas, 7 (2000) 1968
c Mass relocation model by H.R.Strauss, 5 (1998) 2676
c=========================================================================
c	Sperical/Cylindrical/Cubic pellet
c	Input:
c 		YAM 1-3 pellet atomic mass H,D,T or HD,DT pellets only
c		YVP > 0 [km/s] pellet velocity
c		YRP > 0 [cm] Rp initial radius of a pellet
c		YCOS = 1 (-1) for High Field Side (LFS) injection
c		YEFF > 0 fuelling efficiency
c		YSHAPE = 1(Spherical d=2Rp)
c			 2(Cylinrical d=2Rp h=2Rp)
c			 3(Cubic h=2Rp) 
c	Ablation rate is represented in the form:
c	dN/dt[atoms/s] =
c	coeff * Ne**an * Te**at *rp**ap * Mi**am,
c	where 	rp is the pellet radius in cm,
c		Mi is the mass of the pellet material in atomic units
c		Te is plasma electron temperature in eV
c		Ne is the plasma density in cm-3
C 	Output:
C		Te,Ti,Nhydr,(Ndeut,Ntrit) after injection
C
c-----------------------------------------------------------------
C Warnings 	NHYDR,NDEUT,NTRIT densities must be prescribed in ASTRA
c 		model explicily if proper pellet is used
c	Subroutine SMART should be called after NE,TE,TI transport
c	equations
c=======================================================================
c	The values of the coefficients:
c	Parks' model:
c	an=0.333 at=1.64 ap=1.333 am=-0.333
c	Hydrogen: coeff = 1.12e16
c	Krypton:  coeff = 2.1e15
c
c	Kuteev model:
c	Deuterium: an=0.453 at=1.72 ap=1.443 am=-0.283
C	coeff=3.46e14 (2e14 New)
c------------------------------------------------------------------
c=================================================================
c	Input parameters:
c		Np  is the relative mass of the pellet in g/cm+3
c		Vb  is the pellet velosity in km/s
c		rp0 is the pellet radius in cm
c		Mi  is the mass of the pellet material in atomic units
c		cf = coeff/10e15
c		Zp  is the charge number of the pellet material
c		Eion is energy required for full ionization 
c		of a pellet atom 13ev +molecule dissosiation 2.2 eV
c=====================================================================
	implicit none
	include 'for/parameter.inc'
	include 'for/status.inc'
	include 'for/const.inc'
	integer	JS1,J0,JJ,JABS,JDEL,JBEG,JEND,JS,J
	double precision DNI(NRD),DNE(NRD),YKCR(NRD),YSFT(NRD),JJFP(NRD)
     ,	,YTE(NRD),YTI(NRD),YNTE(NRD),YNTI(NRD),YXJ(NRD),YX12(NRD),ALFA
	double precision YAM,YCN,YCE,YCI,YSTNE,YSTNI,YSTNE1,YSTNI1,YTST
	double precision YSDNE,YSDNI,YX,YDX,YF,YDF,YNE,YNI,YDNE,YDNI,YDL
	double precision YR,YDA,YDV,YLC,YFI,YPSI,YBETB,YR1,YR2,YF1,YF2
	double precision YRP1,YRP2,YA1,YA2,YA3,Y16,YPP,YP23,YCOS,YCOSA
	double precision YRP,YVP,YKC1,YEFF,YCOS0,YNp,YSIG,YHRO
	double precision Zp,an,Eion,at,ap,am,coeff,ycoef,yshape
C----------------------------------------------------------------------|
	if(YAM.lt.1..or.YAM.gt.3.) then
	write(*,*) 'Warnig from SMART: wrong input data'
	write(*,*) '	pellet mass must be H/D/T/HD/DT only'
	return
	endif
	if((YAM*YVP*YRP*(yshape-.999)*YEFF).le.0..or.abs(YCOS0).gt.1.) 
     .	then
	write(*,*) 'Warnig from SMART: wrong input parameters'
	return
	endif
Constants for ITER FEAT:
	YCOS0	 =0.7
	YEFF	 =1.
	YSHAPE	 =2.
Constants for Kuteev's model
	Eion=(13.+2.2)*1.e-3
	Zp=1.
	an=0.453
	at=1.72
	ap=1.443
	am=-0.283
	YNp=2.98	!e22 molecules/cm3 of H/D/T
	coeff=2.	!3.46	!e14 coeff for dN/dt
cShape 
	if(yshape.le.1.) then
	write(*,*)	'spherical pellet'
	ycoef=1
	endif
	if(yshape.le.2..and.yshape.gt.1.) then
	write(*,*)	'cylindrical pellet'
	ycoef=1.5
	endif
	if(yshape.gt.2.) then
	write(*,*)	'cubic pellet'
	ycoef=6./GP
	endif

c=====================================================================
Constants
C effective opacity for inclined pellet
C*NEW	vvvvvvv
	YDL=1.
	YCOSA=abs(YCOS0)
	if(YCOSA.gt.1.e-6) then
		YSIG=YCOS0/YCOSA
	else
		YSIG=1
	endif
		YCOS=YSIG
c YCOS - normal to the magn. surf. is parallel to the midplain
C*NEW	^^^^^^^^

	YA1=4.*GP*2*YNP
C 10**11=10**14*100[m->cm]/10**5[km/s->cm/s]/10**22
C /2 Simpson integration
	YA2=coeff*(3-ap)*1.e-11*(1.e13)**an*1000.**at*YAM**am
	YA2=YA2/YVP/YA1/2*YDL
C 1000=10**22/10**19 * shape correction (1 for sphere)
	YA3=1000/3*YA1*YEFF	*ycoef
	YP23=2./3./(3.-ap)
	YPP=3./(3.-ap)
C for Kc = Rperp/Rp (Parks) cm->m
	Y16	=(1/6.)
C*NEW	vvvvvvv
c TSTAR =2 eV
	YTST	=2.
	YKC1	=1.54*(10./YAM)**Y16*(.01)**.667
     .	*sqrt((4*5/3*YTST+13.6+2.2)/(1.-0.5))		/1.4
C*NEW	/2.
C*NEW	^^^^^^^^

	do J=1,NA1
		DNI(J)=0.
		DNE(J)=0.
	YTE(J)=0.
	YTI(J)=0.
	YNTE(J)=0.
	YNTI(J)=0.

		YKCR(J)=0.
		YSFT(J)=0.
	YXJ(J)=(j-1)*HRO/ROC
	YX12(J)=(j-.5)*HRO/ROC

cc		CAR9(J)=0.
cc		CAR8(J)=0.
cc		CAR12(J)=0.
cc		CAR13(J)=0.
cc		CAR14(J)=0.
cc		CAR15(J)=0.
cc		CAR16(J)=0.
	enddo
	YXJ(NA1)=1.
	YX12(NA1)=1.
	
	JS=1
	J=NA1
	YR1=SHIFT-YSIG*JS*ABC
	YF1=NE(NA1)**an*TE(NA1)**at
	YRP1=YRP**(3.-ap)

c	write(*,*) YA1,YA2,YA3
 1	J=J-JS
	YR2=YR1
	YF2=YF1
	YRP2=YRP1
	YR1=SHIF(J)-YSIG*JS*AMETR(J)
	YF1=NE(J)**an*TE(J)**at
	YRP1=YRP2 - YA2*abs(YR2-YR1)*(YF1+YF2)
	
	YKCR(J)=YKC1*TE(J)**Y16/(NE(J)*LOG(2000.*TE(J)/7.5))**.3333
     .	*YRP2**YP23

c	write(*,*) YF1,YF2,YR1,YR2,YRP1,YRP2,J
	if(YRP1.gt.0.)	then
	DNI(J)= YA3*(YRP2**YPP-YRP1**YPP)+DNI(J)

	else
	DNI(J)= YA3*YRP2**YPP+DNI(J)

	endif
C*NEW
	JABS=J

	if(YRP1.le.0.) goto 2

	if(J.eq.1) then
		J=0
		JS=-1
		goto 1
	endif
	if(J.lt.NA1) goto 1	

 2	continue
c*NEW-1 vvvvvvvvvvvvv
	YKCR(NA1)=YKCR(NA)
	do	J=NA1,JABS,-1
		JDEL=2*YKCR(J)/HRO+1
		JBEG=JDEL/2
		JBEG=J+JBEG
		YHRO=2*YKCR(J)/JDEL
	if(JBEG.gt.NA1) JBEG=NA1
		JEND=JBEG-JDEL
	if(JEND.lt.1) JEND=1
		YDV=0.
		YDNE=0.
	do	JJ=JBEG,JEND,-1
	if(JJ.eq.1)	then
		YDV=VOLUM(1)
	else
		YDV=VR(JJ)*YHRO+YDV
CVOLUM(J)-VOLUM(J-1)
	endif
		YDNE=YDNE+DNI(JJ)
	enddo
	YDNI = YDNE/YDV
	YNE	=NE(J)+YDNI*Zp
	YNI	=NI(J)+YDNI
C*NEW	vvvvvvv
C
	YR=RTOR+SHIF(J)
C	YLC=SQRT(2.*YKCR(J)*YR)
	YDA=YKCR(J)
	YLC=SQRT(YKCR(J)*YR)
	YFI=YLC/(YR/MU(J))
c	YFI=YLC/YR
	YBETB=4.e-3*(TE(J)*NE(J)+TI(J)*NI(J))
	YPSI=-YR/MU(J)*YBETB/(BTOR)*YDNE/YLC/
     .  AMETR(J)/YNE/(1+YDA/AMETR(J)/YFI)*YCOS
cc		CAR14X(J)=YDNE
cc		CAR15X(J)=YNE
cc		CAR16X(J)=YDNE/YNE

c*YSIG
C*NEW	^^^^^^^^^^^^^
	YSFT(J)=YPSI

	enddo

c*NEW-1 ^^^^^^^^^^^^^^
	do J=1,NA1
	if(J.eq.1)	then
		YDV=VOLUM(1)
	else
		YDV=VR(J)*HRO
	endif
	YDNI = DNI(J)/YDV
	YNE	=NE(J)+YDNI*Zp
	YNI	=NI(J)+YDNI
Ctemporary output
	YTE(J)=(TE(J)-0.667*YDNI/NE(J)*Eion)*NE(J)/YNE
	YTI(J)=TI(J)*NI(J)/YNI
cc	CAR12(J)=YKCR(J)
cc	CAR13(J)=YSFT(J)/(FP(NA1)-FP(1))
cc	CAR16(J)	=YDNI

	enddo
c			write(*,*) 'TYTA'
c index for density shift

	YDX=1./NRD
	J=2
	YDF=FP(NA1)-FP(1)
	YF=(FP(J)-FP(1))/YDF

	do JJ=1,NRD
	YX=YDX*JJ
	if(YX.gt.YF) then
		J=J+1
		YF=(FP(J)-FP(1))/YDF
	endif
		JJFP(JJ)=J-1
	enddo		
		JJFP(NRD)=NA1
C density, energy shift	
	do j=NA1,1,-1

	if(DNI(J).ne.0.) then
		JJ=NRD*(FP(J)+YSFT(J)-FP(1))/YDF + 1	
	if(JJ.lt.0) JJ=-JJ
	if(JJ.lt.NRD) then

		J0=JJFP(JJ)
	JDEL=2*YKCR(J)/HRO+1
	JBEG=JDEL/2
	JBEG=J0-JBEG
	do JS1=1,JDEL
		JS=JBEG+JS1-1
		if(JS.lt.0)JS=-JS
		if(JS.eq.0.) JS=1
	 	if(JS.le.NA1) then
	DNE(JS)=DNI(J)/JDEL+DNE(JS)
	YNTE(JS)=YTE(J)*DNI(J)/JDEL+YNTE(JS)
	YNTI(JS)=YTI(J)*DNI(J)/JDEL+YNTI(JS)
	endif
	enddo
	endif
	endif
	enddo
C smoothing with energy/particle conservation
	ALFA = 0.001
	call	SMOOTH(ALFA,NA1,DNE,YXJ,NA1,DNI,YX12)
	call	SMOOTH(ALFA,NA1,YNTE,YXJ,NA1,YKCR,YX12)
	call	SMOOTH(ALFA,NA1,YNTI,YXJ,NA1,YSFT,YX12)
		YSDNE=0.
		YSDNI=0.
		YSTNE=0.
		YSTNI=0.
		YSTNE1=0.
		YSTNI1=0.
	do J=1,NA1
		YSDNE=YSDNE+DNE(J)
		YSDNI=YSDNE+DNI(J)
		YSTNE=YSTNE+YNTE(J)
		YSTNI=YSTNI+YNTI(J)
		YSTNE1=YSTNE1+YKCR(J)
		YSTNI1=YSTNI1+YSFT(J)
	YNTE(J)=YKCR(J)
	YNTI(J)=YSFT(J)
	enddo

	YCN=YSDNE/YSDNI
	YCE=YSTNE/YSTNE1
	YCI=YSTNI/YSTNI1		
	DO J=1,NA1
cc	if(DNE(J).gt.0.) then
	if(J.eq.1)	then
		YDV=VOLUM(1)
	else
		YDV=VR(J)*HRO

	endif
	YDNI = DNI(J)/YDV*YCN
	YNE	=NE(J)+YDNI*Zp
	YNI	=NI(J)+YDNI
Ctemporary output
c Te
	TE(J)	=(YTE(J)*NE(J)+YNTE(J)/YDV*YCE)/YNE
cc	CAR14(J)=(YTE(J)*NE(J)+YNTE(J)/YDV*YCE)/YNE
c Ti
	TI(J)	=(YTI(J)*NI(J)+YNTI(J)/YDV*YCI)/YNI
cc	CAR15(J)=(YTI(J)*NI(J)+YNTI(J)/YDV*YCI)/YNI
c Ne
	NE(J)	=NE(J)+YDNI*Zp
cc	CAR8(J)	=CAR8(J)+YDNI*Zp
c Ni
c
cc	goto 777
	if(YAM.eq.1.)	then
	NI(J)=NI(J)-NHYDR(J) 
	NHYDR(J)=NHYDR(J)+YDNI
	NI(J)=NI(J)+NHYDR(J)
	endif
	if(YAM.gt.1.and.YAM.lt.2.)	then
	NI(J)=NI(J)-NHYDR(J)-NDEUT(J) 
	NHYDR(J)=NHYDR(J)+(2.-YAM)*YDNI
	NDEUT(J)=NDEUT(J)+(YAM-1.)*YDNI
	NI(J)=NI(J)+NHYDR(J)+NDEUT(J)
	endif	
	if(YAM.eq.2.)	then
	NI(J)=NI(J)-NDEUT(J) 
	NDEUT(J)=NDEUT(J)+YDNI
	NI(J)=NI(J)+NDEUT(J)
	endif
	if(YAM.eq.3.)	then
	NI(J)=NI(J)-NTRIT(J) 
	NTRIT(J)=NTRIT(J)+YDNI
	NI(J)=NI(J)+NTRIT(J)
	endif
	if(YAM.gt.2.and.YAM.lt.3.)	then
	NI(J)=NI(J)-NTRIT(J)-NDEUT(J) 
	NTRIT(J)=NTRIT(J)+(YAM-2.)*YDNI
	NDEUT(J)=NDEUT(J)+(3.-YAM)*YDNI
	NI(J)=NI(J)+NTRIT(J)+NDEUT(J)
	endif
	F3(j)=NDEUT(j)
	F2(j)=NTRIT(j)
cc 777	continue
c
cc	CAR9(J)	=CAR9(J)+YDNI

	enddo
	return
	end
