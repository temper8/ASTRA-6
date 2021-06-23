C=======================================================================
	SUBROUTINE	SMASH(QTHR,QOUT)
C-----------------------------------------------------------------------
C					(G.Pereverzev 26-JUN-03)
C-----------------------------------------------------------------------
C ! Note !   The subroutine should be called
C		AFTER all the transport equations solved
C Input parameters: 
C   QTHR = q_tr -> if q_min drops below q_tr re-distribution turns on
C   QOUT = q after redistribution
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	double precision	RHOS(3),WBPOLR,NEAVR,FRMIN,RFMIN,RADIAL
	integer	IS(3),JNRES,j,IMIX
	double precision	QTHR,QOUT,YMUOUT,YMUC,YWP,FPMIX,FPHMIX,
     1		YSIGN,YSIGNO,YCM,RHOMIX,YV0,YCU,YM,YM1,YT,
     4		DECONT,DICONT,WECONT,WICONT,YTE0,YDTE
	character*132	STRI
	if ( FRMIN(MU) .gt. 1./QTHR )	return
	if ( RFMIN(MU) .gt. 0.7*ROC )	return

	YMUOUT = 1./QOUT
	YMUC = (4.*MU(1)-MU(2))/3.		! \mu(0)
	JNRES = 0
	YSIGNO = sign(1.d0,YMUC-YMUOUT)
	do  10  j=2,NA1/2
	   YSIGN  = sign(1.d0,MU(j)-YMUOUT)
           if (YSIGN .eq. YSIGNO)	goto	10
	   YSIGNO = YSIGN
C	   if (j .eq. 2)	goto	10
	   JNRES = JNRES+1
	   if (JNRES .gt. 3)	goto	10
	   IS(JNRES) = j
	   RHOS(JNRES) = (j-1+(YMUOUT-MU(j-1))/(MU(j)-MU(j-1)))*HRO
 10	continue
	if (JNRES .eq. 0)	then
	   if (MU(NA1) .lt. YMUOUT)	return
	   write(*,*)"ERROR: mu=q_max/2 is not found inside ROC/2"
	endif
	IMIX = IS(JNRES)
	if (IMIX .gt. NA1/2)	return
	RHOMIX = RHOS(JNRES)

Conserving quantities:
Co	YN = NEAVR(ROC)
Co	YWE = WER(ROC)
Co	YWI = WIR(ROC)
	YWP = WBPOLR(ROC)

	YCM    = GP*BTOR*YMUOUT
	FPHMIX = YCM*RHOMIX**2
	FPMIX  = RADIAL(FP,RHOMIX)
C	write(*,*)IMIX
C	write(*,*)(FP(j),j=1,IMIX)
	do  15  j=1,IMIX			! 1,IMIX-1
	   FP(j) = FPMIX+YCM*RHO(j)**2-FPHMIX
 15	continue
C	write(*,*)(FP(j),j=1,IMIX)
	YCU=2.5*BTOR/(GP*RTOR*HRO)
	YM=0.
C	write(*,*)(MU(j),j=1,IMIX)
	YCM = GP2*BTOR*HRO*HRO
	do	J=1,IMIX
	   MU(J)=(FP(J+1)-FP(J))/(J*YCM)
	   YM1=YM
	   YM=J*G22(J)*MU(J)
	   CU(J)=YCU*G33(J)*IPOL(J)**3*(YM-YM1)/(J-0.5)
C	   UPL(J) = (FP(J)-FPO(J))/TAU
C	   ULON(J) = IPOL(J)*G33(J)*UPL(J)
	enddo
C	write(*,*)(MU(j),j=1,IMIX)

	YTE0  = TE(1)
	YV0   =0.
	DECONT=0.
	DICONT=0.
	WECONT=0.
	WICONT=0.
	do   16   j=1,IMIX
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
	TAU = TAUMIN
C 100	format(1P,5E15.6)
	end
C=======================================================================
