C=======================================================================
	SUBROUTINE	MIXQ(OPTION,RECOND,QRES)
C-----------------------------------------------------------------------
C					(G.Pereverzev 17-NOV-94)
C					(     updated 14-JUL-99)
C-----------------------------------------------------------------------
C ! Note !   The Astra compiler places this subroutine
C		AFTER all the transport equations solved
C-----------------------------------------------------------------------
C B.B.Kadomtsev, Disruptive instability in tokamaks,
C        Sov. J. Plasma Phys. Vol.1, No.5, Sept.-Oct. (1975) pp.389-391.
C V.V.Parail, G.V.Pereverzev, Internal disruption in a tokamak,
C        Sov. J. Plasma Phys. Vol.6, No.1, Jan-Feb. (1980) pp.14-17.
C-----------------------------------------------------------------------
C   OPTION = 0,10,20,30
C		Kadomtsev-type reconnection for a single q=1 resonance 
C		surface. Reconnection condition: dt > RECOND
C   RECOND[sec] determines a minimal time interval
C		between the successive disruptions.
C		For this type of reconnection it is additionally 
C		required that q = 1 surface exists in the plasma 
C
C   OPTION = 1,11,21,31	
C		Kadomtsev-type reconnection for a single q=1 resonance 
C		surface. Reconnection condition: rho_s > RECOND*ROC
C   RECOND[d/l] minimal resonance radius
C
C		In both cases above, when multiple resonance encountered
C		the outermost resonance surface is taken into account
C		only. This can cause a contradiction when
C		(a) a resonance surface first appears long after
C		    the prescribed time interval (OPTION = 0)
C		(b) a resonance surface appears far outside
C		    rho = RECOND*ROC (OPTION = 1)
C		In those cases, message(*) can be printed (see below).
C
C   OPTION = 2,12,22,32
C		q = 1 double tearing mode central reconnection.
C   RECOND[s]	is ignored when two resonance surfaces occur.
C		If only one resonance surface is found then
C		the OPTION = 0 type reconnection is done with
C		the reconnection condition saw_period >= RECOND.
C		This option is useful for the first call.
C		In this case, message(*) can be printed (see below).
C
C	Note:	The subroutine affects control times TAUMIN and TAUMAX 
C		so that	TAUMAX = min(TAUMAX,PERIOD/5.)
C			TAUMIN = min(TAUMIN,TAUMAX/2.)
C		where	PERIOD is the recent sawtooth duration
C
C	To print energy conservation check enable lines marked with "Co"
C
C     (*) To enable message prints use 
C   OPTION = 10, 11, 12, respectively.
C
C	Additional information (involving CFs and CARs ) is available for 
C   OPTION = 20, 21, 22.
C		See the examples below:
C		CAR1 gives psi* quantity before reconnection
C		CAR2 gives psi* quantity after reconnection
C		CMHD1 - relative radius of the 1st resonance surface
C		CMHD2 - relative radius of the 2nd resonance surface
C		CF3 - relative radius of the 3rd resonance surface
C		CMHD3 - relative radius of the reconnection region
C		CMHD4 - psi* maximal value
C
C	Options 30, 31, 32 include (10+20), (11+21), (12+22).
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	double precision RHOS(3),PSIM(3),FPSTAR(NRD),WBPOLR,OPTION,QRES
	integer	IS(3),IOPT,KOPT,JNRES,j,jj,IMIX
	double precision TMIX,RECOND,NEAVR,YMURES,YMUC,YXS,YRS,YWP,YX,
     1		HFP,YSIGN,YSIGNO,YFPC,YCM,YCM1,D1PS,D2PS,RHOMIX,YV0,DPS,
     2		YV2,YV3,YV4,YV5,YV6,YNE,YNI,YTE,YTI,YX2,YCU,YM,YM1,YT,
     3		DEMAX,DIMAX,TEMAX,TIMAX,DEMIN,DIMIN,TEMIN,TIMIN,DDE,DDI,
     4		DECONT,DICONT,WECONT,WICONT,YTE0,YTI0,YNE0,YDTE,YDTI
	character*132	STRI
	save	TMIX
	data	TMIX/-99999./
C Store the ASTRA start time:
	if (TMIX .lt. -99998.) TMIX = TSTART
C or the time of the 1st call:
C	if (TMIX .lt. -99998.) TMIX = TIME
	IOPT = OPTION+0.001
C KOPT = 0,1,2,3 - output option
	KOPT = IOPT/10
C IOPT = 1,2,3 - reconnection option
	IOPT = IOPT-10*KOPT
C Other value "n/m" can be used instead of 1 for the double tearing mode
	YMURES = 1./QRES
	YMUC = (4.*MU(1)-MU(2))/3.
	YFPC = (9.*FP(1)-FP(2))/8.
	YCM = GP2*BTOR*HRO*HRO*YMURES
	YCM1 = 0.5*YCM
	JNRES = 0
	FPSTAR(1) = FP(1)-YFPC-0.125*YCM
	if (KOPT .ge. 2) CAR1(1) = FPSTAR(1)
	if (KOPT .ge. 2) CAR2(1) = FPSTAR(1)
	YSIGNO = sign(1.d0,YMUC-YMURES)
	do  10  j=2,NB1
C "Psi" for a homogeneous current can be calculated as
C	   HFP(1) = 0.125*YMURES*YCM
C	   HFP(j) = HFP(j-1)+YMURES*YCM*(j-1)
C or as
	   HFP = YCM1*(j-0.5)**2
	   FPSTAR(j) = FP(j)-YFPC-HFP
	   if (KOPT .ge. 2) CAR1(j) = FPSTAR(j)
	   if (KOPT .ge. 2) CAR2(j) = FPSTAR(j)
	   YSIGN  = sign(1.d0,MU(j)-YMURES)
C           if (JNRES.eq.0 .and. j.gt.NA1)	return
           if (YSIGN.eq.YSIGNO)	goto	10
	   YSIGNO = YSIGN
	   if (j .eq. 2)	goto	10
	   JNRES = JNRES+1
	   if (JNRES .gt. 3)	goto	10
	   YRS = (j-1+(YMURES-MU(j-1))/(MU(j)-MU(j-1)))*HRO
	   RHOS(JNRES) = YRS
	   IS(JNRES) = j
 10	continue
	if (FPSTAR(NB1).gt.0. .and. (KOPT.eq.1 .or. KOPT.eq.3))   then
	   write(*,*)MU
	   pause'Inverse q-profile'
	   endif
	if ( KOPT .ge. 2)   then
	   CMHD1 = 0.
	   CMHD2 = 0.
C	   CF3 = 0.
	   endif
	if (JNRES .eq. 0)   return		! No resonance found
	if ( KOPT .ge. 2)   CMHD1 = RHOS(1)/ROC

C One resonance surface
	if (JNRES .gt. 1)   goto   20
C The two lines allow to switch the reconnection condition between
C	the OPTION = 1 type with the condition rho_s > RECOND*ROC and
C	the OPTION = 0 type with the condition saw_period > RECOND.
C	if ( IOPT .gt. 1 .and. RHOS(JNRES) .gt. RECOND*ROC)   then
	if ( IOPT .gt. 1 .and. TIME-TMIX+.5*TAU .gt. RECOND)   then
	   if (KOPT .eq. 1 .or. KOPT .eq. 3)  write(*,*)
     +	   'One resonance',RHOS(1)/ROC,' found instead of ',
     +	   'two expected. ',"Kadomtsev's reconnection done."
	   goto  12
	endif
 11	if ( IOPT .eq. 1 .and. RHOS(JNRES) .gt. RECOND*ROC .or.
     +	     IOPT .eq. 0 .and. TIME-TMIX+.5*TAU .gt. RECOND)   goto   12
	return

 12	if (JNRES .eq. 2 .and. (KOPT.eq.1 .or. KOPT.eq.3))
     +		write(*,*)'>>> MIXQ >>> Double tearing mode ignored'
	j = IS(JNRES)
	YXS = RHOS(JNRES)/HRO-j+0.5
	D1PS = FPSTAR(j+1)-FPSTAR(j-1)
	D2PS = FPSTAR(j+1)-2.*FPSTAR(j)+FPSTAR(j-1)
	PSIM(JNRES) = 0.5*YXS*(D1PS+YXS*D2PS)+FPSTAR(j)
	do   14   j=IS(JNRES)-1,NA1
	   if (FPSTAR(j) .gt. 0.)   goto   14
	   RHOMIX=(j-1.5-FPSTAR(j-1)/(FPSTAR(j)-FPSTAR(j-1)))*HRO
	   IMIX = j-1
	   goto  15
 14	continue
	if (FPSTAR(NA1) .gt. 0.)   return
 15	continue

	if (KOPT .ge. 2)  CMHD4 = PSIM(JNRES)
	if (KOPT .ge. 2)  CMHD3 = RHOMIX/ROC

Conserving quantities:
Co	YN = NEAVR(ROC)
Co	YWE = WER(ROC)
Co	YWI = WIR(ROC)

	YWP = WBPOLR(ROC)
	YTE0 = TE(1)
	YTI0 = TI(1)
	YNE0 = NE(1)
	YV0=0.
	DECONT=0.
	DICONT=0.
	WECONT=0.
	WICONT=0.
	do   16   j=1,IMIX
	   FPSTAR(j) = PSIM(JNRES)*(1.-((j-0.5)*HRO/RHOMIX)**4)
	   FP(j) = FPSTAR(j)+YFPC+YCM1*(j-0.5)**2
	   if (KOPT .ge. 2) CAR2(j) = FPSTAR(j)
	   YV0  = YV0+VR(J)
	   DECONT = DECONT+NE(J)*VR(J)
	   DICONT = DICONT+NI(J)*VR(J)
	   WECONT = WECONT+NE(J)*TE(J)*VR(J)
	   WICONT = WICONT+NI(J)*TI(J)*VR(J)
 16	continue
	do   17   j=1,IMIX
	   if(LEQ(1).gt.0)	NE(j) = DECONT/YV0
	   if(LEQ(1).gt.0)	NI(j) = DICONT/YV0
	   if(LEQ(2).gt.0)	TE(j) = WECONT/DECONT
	   if(LEQ(3).gt.0)	TI(j) = WICONT/DICONT
 17	continue
	goto   50
	     
C Two resonance surfaces
 20	if ( KOPT .ge. 2)   CMHD2 = RHOS(2)/ROC
	if (JNRES .gt. 2)   goto   30
	do  21  jj =1,JNRES
	    j = IS(jj)
	    YXS = RHOS(jj)/HRO-j+0.5
	    D1PS = FPSTAR(j+1)-FPSTAR(j-1)
	    D2PS = FPSTAR(j+1)-2.*FPSTAR(j)+FPSTAR(j-1)
	    PSIM(jj) = 0.5*YXS*(D1PS+YXS*D2PS)+FPSTAR(j)
 21	continue

	if (PSIM(2) .lt. 0.)   return

C Enforced switch to the OPTION = 0 or 1 type reconnection.
C The outermost resonance is taking into account only.
	if (IOPT .lt. 2)   goto   11
		
C This "goto" allows OPTION = 0 or 1 type reconnection.
C The reconnection takes into account the outermost resonance only.
C	if (IOPT .eq. 2 .and. abs(PSIM(1)) .lt. .1*PSIM(2))  goto  12

	do  23  j=IS(JNRES)-1,NA1
	    if (FPSTAR(j) .gt. PSIM(1))   goto   23
	    DPS = (PSIM(1)-FPSTAR(j-1))/(FPSTAR(j)-FPSTAR(j-1))
	    RHOMIX=(j-1.5-DPS)*HRO
	    IMIX = j-1
	    goto  24
 23	continue
 24	continue

	if (KOPT .ge. 2)  CMHD4 = PSIM(JNRES)
	if (KOPT .ge. 2)  CMHD3 = RHOMIX/ROC

	do  25  jj =1,JNRES
	    j = IS(jj)
	    YXS = RHOS(jj)/HRO-j+0.5
	    if (jj .eq. 1)   then
	       DEMAX = NE(j-1)+YXS*(NE(j)-NE(j-1))
	       DIMAX = NI(j-1)+YXS*(NI(j)-NI(j-1))
	       TEMAX = TE(j-1)+YXS*(TE(j)-TE(j-1))
	       TIMAX = TI(j-1)+YXS*(TI(j)-TI(j-1))
	    endif
	    if (jj .eq. 2)   then
	       DEMIN = NE(j-1)+YXS*(NE(j)-NE(j-1))
	       DIMIN = NI(j-1)+YXS*(NI(j)-NI(j-1))
	       TEMIN = TE(j-1)+YXS*(TE(j)-TE(j-1))
	       TIMIN = TI(j-1)+YXS*(TI(j)-TI(j-1))
	    endif
 25	continue

Conserving quantities:
Co	YN = NEAVR(ROC)
Co	YWE = WER(ROC)
Co	YWI = WIR(ROC)

	YWP = WBPOLR(ROC)
	YTE0 = TE(1)
	YTI0 = TI(1)
	YNE0 = NE(1)
	YV0=0.
	YV2=0.
	YV3=0.
	YV4=0.
	YV5=0.
	YV6=0.
	DECONT=0.
	DICONT=0.
	WECONT=0.
	WICONT=0.
	do   26   j=1,IMIX
	   YX = (j-0.5)*HRO/RHOMIX
C	   FPSTAR(j) = PSIM(1)*(1.-(1.-(YX)**4)**2)
	   FPSTAR(j) = (PSIM(1)-PSIM(2))*(1.-(1.-(YX)**4)**2)+PSIM(2)
	   FP(j) = FPSTAR(j)+YFPC+YCM1*(j-0.5)**2
	   if (KOPT .ge. 2) CAR2(j) = FPSTAR(j)
	   YV0 = YV0+VR(J)
	   YV2 = YV2+VR(J)*YX**2
	   YV3 = YV3+VR(J)*YX**3
	   YV4 = YV4+VR(J)*YX**4
	   YV5 = YV5+VR(J)*YX**5
	   YV6 = YV6+VR(J)*YX**6
	   DECONT = DECONT+NE(J)*VR(J)
	   DICONT = DICONT+NI(J)*VR(J)
	   WECONT = WECONT+NE(J)*TE(J)*VR(J)
	   WICONT = WICONT+NI(J)*TI(J)*VR(J)
 26	continue
	DDE = DEMAX-DEMIN
	DDI = DIMAX-DIMIN
	YDTE = TEMAX-TEMIN
	YDTI = TIMAX-TIMIN
	YNE = (DECONT-DEMIN*YV0-DDE*YV2)/(YV2-YV3)
	YNI = (DICONT-DIMIN*YV0-DDI*YV2)/(YV2-YV3)
	YTE = (WECONT-DEMIN*TEMIN*YV0-(DEMIN*YDTE+TEMIN*DDE)*YV2
     +	      -YNE*(TEMIN*(YV2-YV3)+YDTE*(YV4-YV5))-DDE*YDTE*YV4)/
     +	      (DEMIN*(YV2-YV3)+(DDE+YNE)*(YV4-YV5)-YNE*(YV5-YV6))
	YTI = (WICONT-DIMIN*TIMIN*YV0-(DIMIN*YDTI+TIMIN*DDI)*YV2
     +	      -YNI*(TIMIN*(YV2-YV3)+YDTI*(YV4-YV5))-DDI*YDTI*YV4)/
     +	      (DIMIN*(YV2-YV3)+(DDI+YNI)*(YV4-YV5)-YNI*(YV5-YV6))
	do   27   j=1,IMIX
	   YX = (j-0.5)*HRO/RHOMIX
	   YX2 = YX*YX
	   if(LEQ(1).gt.0)   NE(j) = DEMIN+(DDE+YNE*(1.-YX))*YX2
	   if(LEQ(1).gt.0)   NI(j) = DIMIN+(DDI+YNI*(1.-YX))*YX2
	   if(LEQ(2).gt.0)   TE(j) = TEMIN+(YDTE+YTE*(1.-YX))*YX2
	   if(LEQ(3).gt.0)   TI(j) = TIMIN+(YDTI+YTI*(1.-YX))*YX2
C Homogeneous profiles:
C	   if(LEQ(1).gt.0)   NE(j) = DECONT/YV0
C	   if(LEQ(1).gt.0)   NI(j) = DICONT/YV0
C	   if(LEQ(2).gt.0)   TE(j) = WECONT/DECONT
C	   if(LEQ(3).gt.0)   TI(j) = WICONT/DICONT
 27	continue
	goto   50

C Three resonance surfaces
 30	continue
C	if (KOPT .ge. 2)   CF3 = RHOS(3)/ROC
	if (JNRES .gt. 3)   goto   40
C	if (IOPT  .lt. 2 .and. (KOPT.eq.1 .or. KOPT.eq.3))
C     +			write(*,*)'Triple tearing mode encountered'
	if (IOPT .lt. 2)   then
C This "goto" switches to the OPTION = 0 or 1 type reconnection so that
C the reconnection takes into account the outermost resonance only.
		goto   11
	endif
	return

 40	if (JNRES .ge. 4 .and. (KOPT.eq.2 .or. KOPT.eq.3))
     +		 write(*,*)'>>> MIXER: too many resonant surfaces'
	return

 50	YCU=2.5*BTOR/(GP*RTOR*HRO)
	YM=0.
C	write(*,*)(MU(j),j=1,IMIX)
	DO 51 J=1,IMIX
	MU(J)=min(YMURES,(FP(J+1)-FP(J))/(J*YCM))
	YM1=YM
	YM=J*G22(J)*MU(J)
	CU(J)=YCU*G33(J)*IPOL(J)**3*(YM-YM1)/(J-0.5)
C	UPL(J) = (FP(J)-FPO(J))/TAU
C	ULON(J) = IPOL(J)*G33(J)*UPL(J)
 51	continue
C	write(*,*)(MU(j),j=1,IMIX)

Check of energy conservation
Co	write(*,*)(YN-NEAVR(ROC))/(YV0*HRO)
Co     +	,1000*(YWE-WER(ROC))/(0.0024*DECONT*HRO)
Co     +	,1000*(YWI-WIR(ROC))/(0.0024*DECONT*HRO),1000.*YDTE
	YDTE = max(0.d0,(YWP-WBPOLR(ROC))/(0.0024*DECONT*HRO))
	do   53   j=1,IMIX
	   if(LEQ(2).gt.0)   TE(j) = TE(j)+YDTE
 53	continue
C	write(STRI,
C     +	 '(1A9,1F7.2,1A3,1A17,1F6.1,1A3,1A9,1F5.2,1A20,1F4.2,1A16)')
	YT = TIME-TMIX
	if (YT .lt. 0.75)	then
	write(STRI,'(1A9,1F7.2,1A4)')' Period =',1000.*YT,' ms;'
	else
	write(STRI,'(1A9,1F7.3,1A4)')' Period =',YT,' s; '
	endif
C-----------------------------------------------------------------------
	YT = YTE0-TE(1)
	if (YT .lt. 0.75)	then
	write(STRI(21:),'(1A,1F6.1,1A,1F5.2)')
     +	'  Te amplitude: ',1000.*YT,' eV;   r0/a =',RHOMIX/ROC
	else
	write(STRI(21:),'(1A,1F6.3,1A,1F5.2)')
     +	'  Te amplitude: ',YT,' keV;  r0/a =',RHOMIX/ROC
	endif
C     +	,'  Magnetic energy: ',1000.*YDTE,' eV per particle'
C	59 = 9+7+3+17+6+3+9+5
	if (KOPT .ge. 1)write(*,'(A)')STRI
	YT = max(1.d-4,(TIME-TMIX)/5.)
C	TAUMAX = min(TAUMAX,YT)
C	TAUMIN = min(TAUMIN,TAUMAX/2.)
	TAU = TAUMIN
	TMIX = TIME
	end
C=======================================================================
