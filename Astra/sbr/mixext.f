C=======================================================================
	SUBROUTINE	MIXEXT(OPTION,RECOND)
C-----------------------------------------------------------------------
C					(G.Pereverzev 15-NOV-94)
C-----------------------------------------------------------------------
C   OPTION = 0	Kadomtsev-type reconnection for a single q=1 resonance 
C		surface. Reconnection condition: dt > RECOND
C   RECOND[sec] should determine minimal time interval
C		between successive disruptions.
C		For this type of reconnection the additional condition
C		that q = 1 exists in the plasma is required
C	Note:	This option affects TAUMIN := min(TAUMIN,RECOND/5.)
C   OPTION = 1	Kadomtsev-type reconnection for a single q=1 resonance 
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
C		In those cases, message is printed.
C
C   OPTION = 2	q = 1 double tearing mode central reconnection.
C   RECOND[d/l]	is ignored when two resonance surfaces occur.
C		If only one resonance surface is found then
C		the OPTION = 1 type reconnection is done with
C		the reconnection condition rho_s > RECOND*ROC
C		In this case, message is printed.
C
C   OPTION = 3	q = n/m double tearing mode edge reconnection.
C		Not implemented.
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	integer	IS(3),IOPT,j,jj,IMIX,JNRES
	double precision RHOS(3),PSIM(3),FPSTAR(NRD),HFP(NRD)
	double precision TMIX,POLNUM,TORNUM,OPTION,RECOND,NEAVR
	double precision YNE,YNI,YDE,YDI,YTE,YTI,YX,YX2,YM,YM1,YCU
	double precision YV0,YV2,YV3,YV4,YV5,YV6,YDTE,YDTI,YTE0,YTI0
	double precision TEMIN,TIMIN,DEMIN,DIMIN,TEMAX,TIMAX,DEMAX,DIMAX
	double precision DECONT,DICONT,WECONT,WICONT,DPS,D1PS,D2PS
	double precision RHOMIX,WER,WIR,WBPOLR,YN,YWE,YWI,YWP,YNE0
	double precision YSIGN,YSIGNO,YFPC,YMUC,YMURES,YCM,YC1,YRS,YXS
	save	TMIX
	data	TMIX/-99999./	POLNUM/1./	TORNUM/1./
C-----------------------------------------------------------------------
C Store the ASTRA start time:
	if (TMIX .lt. -99998.) TMIX = TINIT
C or the time of the 1st call:
C	if (TMIX .lt. -99998.) TMIX = TIME
	IOPT = OPTION+0.001
	YMURES = TORNUM/POLNUM
	if (IOPT .le. 2)   YMURES = 1.
	YMUC = (4.*MU(1)-MU(2))/3.
	YFPC = (9.*FP(1)-FP(2))/8.
	YCM = GP2*BTOR*HRO*HRO
	YC1 = YMURES*YCM
	JNRES = 0
	HFP(1) = 0.125*YC1
	FPSTAR(1) = FP(1)-YFPC-HFP(1)
	CAR1(1) = FPSTAR(1)
	CAR3(1) = FPSTAR(1)
	YSIGNO = sign(1.d0,YMUC-YMURES)
	do   10   j=2,NB1
	     HFP(j) = HFP(j-1)+YC1*(j-1)
C	     HFP(j) = 0.5*YCM*(j-0.5)**2
	     FPSTAR(j) = FP(j)-YFPC-HFP(j)
	CAR1(j) = FPSTAR(j)
	CAR3(j) = FPSTAR(j)
	     YSIGN  = sign(1.d0,MU(j)-YMURES)
             if(YSIGN.eq.YSIGNO)   goto   10
	     YSIGNO = YSIGN
	     if (j .eq. 2)         goto   10
	     JNRES = JNRES+1
	     YRS = (j-1+(YMURES-MU(j-1))/(MU(j)-MU(j-1)))*HRO
	     RHOS(JNRES) = YRS
	     IS(JNRES) = j
 10	continue
	if (FPSTAR(NB1) .gt. 0.)   then
	   write(*,*)MU
	   pause'Inverse q-profile'
	   endif
	if (JNRES .eq. 0)   return
	CF1 = RHOS(1)/ROC

C One resonance surface
	if (JNRES .gt. 1)   goto   20
	if ( IOPT .gt. 1 .and. RHOS(JNRES) .gt. RECOND*ROC)   then
	   write(*,*)'One resonance',RHOS(1)/ROC,' found instead of ',
     +		'two expected. ',"Kadomtsev's reconnection done."
	   goto  12
	endif
 11	if ( IOPT .eq. 1 .and. RHOS(JNRES) .gt. RECOND*ROC .or.
     +	     IOPT .eq. 0 .and.  TIME-TMIX  .gt. RECOND)   goto   12
	return

 12	if (JNRES .eq. 2) write(*,*)'Double tearing mode ignored'
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
 15	continue

	CF4 = PSIM(JNRES)
	CF5 = RHOMIX/ROC

Conserving quantities:
C	write(*,*)NEAVR(ROC),WER(ROC),WIR(ROC)
	YN = NEAVR(ROC)
	YWE = WER(ROC)
	YWI = WIR(ROC)

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
	   FP(j) = FPSTAR(j)+YFPC+HFP(j)
	CAR3(j) = FPSTAR(j)
	   YV0  = YV0+VR(J)
	   DECONT = DECONT+NE(J)*VR(J)
	   DICONT = DICONT+NI(J)*VR(J)
	   WECONT = WECONT+NE(J)*TE(J)*VR(J)
	   WICONT = WICONT+NI(J)*TI(J)*VR(J)
 16	continue
	do   17   j=1,IMIX
	   NE(j) = DECONT/YV0
	   NI(j) = DICONT/YV0
	   TE(j) = WECONT/DECONT
	   TI(j) = WICONT/DICONT
 17	continue
	   if ( IOPT .eq. 0)  TAUMIN = min(TAUMIN,RECOND/5.)
	goto   50
	     
C Two resonance surfaces
 20	if (JNRES .gt. 2)   goto   30
	CF2 = RHOS(2)/ROC
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

	CF4 = PSIM(JNRES)
	CF5 = RHOMIX/ROC

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
	YN = NEAVR(ROC)
	YWE = WER(ROC)
	YWI = WIR(ROC)

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
	   FPSTAR(j) = PSIM(1)*(1.-(1.-(YX)**4)**2)
	   FP(j) = FPSTAR(j)+YFPC+HFP(j)
	CAR3(j) = FPSTAR(j)
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
	YDE = DEMAX-DEMIN
	YDI = DIMAX-DIMIN
	YDTE = TEMAX-TEMIN
	YDTI = TIMAX-TIMIN
	YNE = (DECONT-DEMIN*YV0-YDE*YV2)/(YV2-YV3)
	YNI = (DICONT-DIMIN*YV0-YDI*YV2)/(YV2-YV3)
	YTE = (WECONT-DEMIN*TEMIN*YV0-(DEMIN*YDTE+TEMIN*YDE)*YV2
     +	      -YNE*(TEMIN*(YV2-YV3)+YDTE*(YV4-YV5))-YDE*YDTE*YV4)/
     +	      (DEMIN*(YV2-YV3)+(YDE+YNE)*(YV4-YV5)-YNE*(YV5-YV6))
	YTI = (WICONT-DIMIN*TIMIN*YV0-(DIMIN*YDTI+TIMIN*YDI)*YV2
     +	      -YNI*(TIMIN*(YV2-YV3)+YDTI*(YV4-YV5))-YDI*YDTI*YV4)/
     +	      (DIMIN*(YV2-YV3)+(YDI+YNI)*(YV4-YV5)-YNI*(YV5-YV6))
	do   27   j=1,IMIX
	   YX = (j-0.5)*HRO/RHOMIX
	   YX2 = YX*YX
	   NE(j) = DEMIN+(YDE+YNE*(1.-YX))*YX2
	   NI(j) = DIMIN+(YDI+YNI*(1.-YX))*YX2
	   TE(j) = TEMIN+(YDTE+YTE*(1.-YX))*YX2
	   TI(j) = TIMIN+(YDTI+YTI*(1.-YX))*YX2
C Homogeneous profiles:
C	   NE(j) = DECONT/YV0
C	   NI(j) = DICONT/YV0
C	   TE(j) = WECONT/DECONT
C	   TI(j) = WICONT/DICONT
 27	continue
	TAUMAX = min(TAUMAX,(TIME-TMIX)/5.)
	TAUMIN = min(TAUMIN,TAUMAX/2.)
	goto   50

C Three resonance surfaces
 30	if (JNRES .gt. 3)   goto   40
	if (IOPT  .lt. 2) write(*,*)'Triple tearing mode encountered'
	CF3 = RHOS(3)/ROC
	if (IOPT .lt. 2)   then
C This "goto" enables switching to the OPTION = 0 or 1 type reconnection
C then the reconnection takes into account the outermost resonance only.
C		write(*,*)'Double tearing mode ignored'
		goto   11
	endif
C Mixing suppressed:
C	goto   50
	return

 40	if(JNRES.ge.4) write(*,*)'>>> MIXER: too many resonant surfaces'
	return
C-----------------------------------------------------------------------
 50	YCU=2.5*BTOR/(GP*RTOR*HRO)
	YM=0.
	DO 51 J=1,NB1
	if(J.lt.NB1)	then
	MU(J)=min(.99999d0,(FP(J+1)-FP(J))/(J*YCM))
			else
	MU(j)=min(.99999d0,MU(NB1-1)*(NB1-1)/NB1*G22(NB1-1)/G22(NB1))
			endif
	YM1=YM
	YM=J*G22(J)*MU(J)
	CU(J)=YCU*G33(J)*IPOL(J)**3*(YM-YM1)/(J-0.5)
	UPL(J) = (FP(J)-FPO(J))/TAU
	ULON(J) = IPOL(J)*G33(J)*UPL(J)
 51	continue
	TAU = TAUMIN
	YDTE = max(0.d0,(YWP-WBPOLR(ROC))/(0.0024*DECONT*HRO))
	do   52   j=1,IMIX
	   TE(j) = TE(j)+YDTE
 52	continue
Check of energy conservation
C	write(*,*)(YN-NEAVR(ROC))/(YV0*HRO)
C     +	,1000*(YWE-WER(ROC))/(0.0024*DECONT*HRO)
C     +	,1000*(YWI-WIR(ROC))/(0.0024*DECONT*HRO),1000.*YDTE
	write(*,'(1A17,1F5.2,1A4,1A22,1F5.1,1A4,/1A32,1F4.2,1A16)')
     +	   ' Time interval = ',1000.*(TIME-TMIX),' ms;'
     +	,'   Te saw amplitude:  ',1000.*(YTE0-TE(1)),' eV;'
C     +	,'   Magnetic energy dissipated: ',1000.*YDTE,' eV per particle'
	TMIX = TIME
	end
C=======================================================================
