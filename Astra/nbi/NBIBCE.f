C All subroutines for time dependent Fokker-Planck calculations
C	AERF, CERF, NBCOEF, NBMESH, NBPOMU, NBPOVE, NBIONR
C========================================== for time-dependent FP
 	double precision	function CERF(X)
C----------------------------------------------------- 17.11.89
C                     ____
C CERF =( exp(-x2)/x/V pi + erf(x)*(1-.5/x2))/2x
C	
C-------------------------------------------------------------
	implicit none
	double precision	X,X2,EX2,T,T2,T3,ERF
	IF(X.LT.0.06)	THEN
	CERF	=0.3761264
	RETURN
	ENDIF
		X2	=X*X
		EX2	=EXP(-X2)
		T	=1./(1.+0.3275911*X)
		T2	=T*T
		T3	=T2*T
	ERF	=(1.-(0.254829592*T-0.284496736*T2+1.421413741*T3-
     -	1.453152027*T2*T2+1.061405429*T2*T3)*EX2)
	CERF=(.56418959*EX2+ERF*(X-.5/X))/(2.*X2)
	return
	END
	double precision	function	AERF(X)
C----------------------------------------------------- 17.11.89
C			 		     ___
C		AERF=erf(x) - 2.*x*exp(-x2)/V pi
C       X2
C ERF =[S  exp(-s2)ds] 2/SQ pi (1.5E-7) Abramowitz,Stegun p.122 
C	0
C-------------------------------------------------------------
	implicit none
	double precision	X,EX2,T,T2,T3,ERF
	IF(X.LT.0.055)	THEN
	AERF	=0.7522528*X**3
	RETURN
	ENDIF
		EX2	=exp(-X*X)
		T	=1./(1.+0.3275911*X)
		T2	=T*T
		T3	=T2*T
	ERF	=(1.-(0.254829592*T-0.284496736*T2+1.421413741*T3-
     -	1.453152027*T2*T2+1.061405429*T2*T3)*EX2)
	AERF	=ERF	-1.1283792*X*EX2
	return
	end
	
C============================================= F-P Coefficients
	SUBROUTINE NBCOEF
	implicit none
	INCLUDE 'nbi/nbfp.inc'
	integer	ISPEND,ISPE,IV1,I,J,I1
	common /CNBSP/ ISPEND,ISPE(9)
	double precision	PI,PI2,SQPI,VB2,YE,YA,YB,YD,YI,YDEL
	double precision	YEXPI,AERF,CERF,YHDVB,YIP,YAB
	PI	=3.1415926536
	PI2	=PI*PI
	SQPI	=SQRT(PI)
	IV1	=IV+1
	DO 1	I=1,IV1
		A1(I)	=0.
		B(I)	=0.
		A(I)	=((I+1)*I/(HV*(I+0.5))**3)
		B1(I)	=0.
 1		D(I)	=0.
C....Coefficients for plasma species Aelectron,Ap,Ad
	DO 2	J=1,ISPEND
		VB2	=VB(J)*VB(J)
		YE	=HV*HV/EB(J)
		YDEL	=ZB(J)**2*RNB(J)/RNB(1)
		YB	=.5*VB2*YDEL*EXP(0.75*YE)
		YA	=EXP(YE)
		YD	=YDEL/(HM*HM)/VB(J)
		YHDVB	=HV/VB(J)
		DO 21	I=1,IV1
			I1	=I+1
			YI	=YHDVB*I
			YIP	=YHDVB*(I+.5)
			YEXPI	=EXP(I*YE)
			YAB	=YB*A(I)*AERF(YIP)
			B(I)	=B(I)+YAB*YEXPI
			A1(I)	=A1(I)+YAB/(YEXPI*YA)
 21			D(I)	=D(I)+YD*CERF(YI)
	IF(J.EQ.1)	THEN
	do	22	I=1,IV1
			AE(I)	=A1(I)
 22			BE(I)	=B(I)
	ENDIF
 2	CONTINUE
		A(1)	=0.
	DO 30	I=1,IV1
			AI(I)	=A1(I)-AE(I)
 30			BI(I)	=B(I)-BE(I)
	DO 3	I=1,IV
		I1	=I+1
		D(I)	=D(I)*DV2(I)
		A(I1)	=A1(I)*DV2(I1)
		A1(I)	=A1(I)*DV2(I)
		B1(I1)	=B(I)*DV2(I1)
 3		B(I)	=B(I)*DV2(I)
		B(IV1)	=B(IV1)*DV2(IV1)
		D(IV1)	=D(IV1)*DV2(IV1)
		A1(IV1)	=A1(IV1)*DV2(IV1)
	RETURN
	END
C======================================================== V,MU mesh
	SUBROUTINE NBMESH
	implicit none
	INCLUDE 'nbi/nbfp.inc'
	double precision	PI,PI2,SQPI
	integer	IV1,IT1,JV,JT
	PI	=3.1415926536
	PI2	=PI*PI
	SQPI	=SQRT(PI)
C...Vi=i*HV i=1,2,...,IV,IV1;	MUj=-1+(j-1/2)*Hmu j=1,2,..,IT
	IV1	=IV+1
	IT1	=IT+1
	HM	=2./IT
	HV	=4./3./IV
C...YM2j = 1 - (Mj+1/2)**2
	DO 1	JT=1,IT
		YM1(JT)	=-1.+HM*(JT-.5)
 1		YM2(JT)	=JT*HM*(2.-JT*HM)
C...DV2i = 1/Vi**2
	DO 2	JV	=1,IV1
 2		DV2(JV)	=1./(HV*JV)**2
	RETURN
	END
C======================================================== MU sweep
	SUBROUTINE NBPOMU
	implicit none
	INCLUDE 'nbi/nbfp.inc'
	double precision	PI,PI2,SQPI,AJ,BJ,CJ,DJ,Y
	integer	IV1,I,ITM1,J,J1,JM
	PI	=3.1415926536
	PI2	=PI*PI
	SQPI	=SQRT(PI)
	IV1	=IV+1
	ITM1	=IT-1
C....Calculations of ALFA,BETA
	DO 1	I=1,IV
		BJ	=D(I)*YM2(1)
cr		CJ	=BJ+DT		+RMN(I,1)
cr		DJ	=FSRS(I,1)+DT*FVM(I,1)
		CJ	=BJ+DT
crr		DJ	=FSRS(I,1)+(DT		-RMN(I,1))*FVM(I,1)
		DJ	=FSRS(I,1)+DT*FVM(I,1)
		AL(2)	=BJ/CJ
		BT(2)	=DJ/CJ
	DO 11	J=2,ITM1
		J1	=J+1
cr		DJ	=FSRS(I,J)+DT*FVM(I,J)
crr		DJ	=FSRS(I,J)+(DT-RMN(I,J))*FVM(I,J)
		DJ	=FSRS(I,J)+DT*FVM(I,J)
		AJ	=BJ
		BJ	=D(I)*YM2(J)
cr		CJ	=AJ+BJ+DT	+RMN(I,J)
		CJ	=AJ+BJ+DT
		Y	=(CJ-AL(J)*AJ)
		AL(J1)	=BJ/Y
 11		BT(J1)	=(DJ+AJ*BT(J))/Y
C....Sweeping
		DJ	=FSRS(I,IT)+DT*FVM(I,IT)
		AJ	=BJ
crr		CJ	=AJ+DT		+RMN(I,IT)
		CJ	=AJ+DT
		FVM(I,IT)=(DJ+AJ*BT(IT))/(CJ-AL(IT)*AJ)
	DO	1	JM=1,ITM1
		J	=IT-JM
		J1	=J+1
 1		FVM(I,J)	=AL(J1)*FVM(I,J1)+BT(J1)
	RETURN
	END
C========================================================= V sweep
	SUBROUTINE NBPOVE
	implicit none
	INCLUDE 'nbi/nbfp.inc'
	double precision	PI,PI2,SQPI,Y
	integer	IV1,J,JP,I,I1,IM
	PI	=3.1415926536
	PI2	=PI*PI
	SQPI	=SQRT(PI)
	IV1	=IV+1
	DO 1 J=1,IT
C.....Calculations of alfa, beta
		JP	=J+1
		Y	=DT+A1(1)
     .	+RMN(1,J)
crr HERBHO ucTO4HuK HA BEPXHEM C/\OE
cr	+RMN(1,J)
		AL(2)	=B(1)/Y
     		BT(2)	=(FSRS(1,J)+FVM(1,J)*DT)/Y
C....Sweeping
	DO 11	I	=2,IV
		I1	=I+1
		Y	=DT+A1(I)+B1(I)-AL(I)*A(I)
     .	+RMN(I,J)
crr HERBHO ucTO4HuK HA BEPXHEM C/\OE

cr	+RMN(I,J)
		AL(I1)	=B(I)/Y
		BT(I1)	=(A(I)*BT(I)+FSRS(I,J)+FVM(I,J)*DT)/Y
 11	CONTINUE
		FVM(IV1,J)=BT(IV1)*A1(IV)/(B(IV)-A1(IV)*AL(IV1))
	DO 1	IM	=1,IV
		I	=IV-IM+1
		I1	=I+1
 1		FVM(I,J)	=AL(I1)*FVM(I1,J)+BT(I1)
	RETURN
	end
C==================================================================
        SUBROUTINE NBIONR(EBEAM,ABEAM,AMJ,RTOR,NA1,TAU,NNCL,NNWM,
     >		CBM1,CBM2,CBM3,CBM4,CBMI1,CBMI2,CBMI3,CBMI4)
C====================================================== 22-MAR-99+30-MAY-06
C	2D (v, cosTET) Fokker - Plank solver for different
c	magnetic surfaces, ripple losses
C======================================================== Polevoy
	implicit none
	INCLUDE  'for/parameter.inc'
        double precision	EBEAM,ABEAM,AMJ,RTOR,TAU,NNCL,NNWM,
     >		CBM1,CBM2,CBM3,CBM4,CBMI1,CBMI2,CBMI3,CBMI4
C Arrays used
C	double precision PEBM(1),PIBM(1),CUBM(1),CUFI(1),PBPER(1),NIBM(1)
C	double precision PBLON(1),SNEBM(1),SNNBM(1)
C	double precision AMAIN(1),NE(1),TE(1),TI(1),ZEF(1),SHIF(1),AMETR(1)
C	double precision EXTARR(1,1),NN(1),TN(1),NNBM1(1),NNBM2(1),NNBM3(1)
	INCLUDE  'for/status.inc'
	INCLUDE	'nbi/nbfp.inc'
	integer	JFPBEG,JN1OLD,ISPEND,ISPE
	common /CNBSP/ ISPEND,ISPE(9)
	double precision YASBA,PBCX,YFCUR,YFI,CNSFI,YLNI,YLNE,YLNZ
	common /CSRS/	YASBA(3,NRD,51)
	common/CION2/	PBCX(NRD),YFCUR(NRD),YFI(3),CNSFI(3)
     .	,YLNI(NRD),YLNE(NRD),YLNZ(NRD)
	CHARACTER	YSRSNM*12
	double precision FNBF,SPEX,PBEIE,PBICX,GP,YET,CBMI33,YM2F
	double precision SQPI,DTION,YEV21,YEV22,YEV23,YJ2,YEPS,T0,TSNBI
	double precision CNSTN,CNSTQ,CNSTP,CNSTC,CNSTE,CNSTT,CNSTQT
	double precision DTAU,EBDTI,EBDTI0,CNSTE0,CNSNN,CNSNN0,YNN0,YV2
	double precision YRMN,YSRSE,YEBEAM,YPB,YIP,YV4,FVMMIN,F0J,YEXARG
	double precision YE,YI,YPEBM,YPIBM,YPBPER,YPBLON,YNB,YCUFI,YMF
	double precision PTHBM,PTHERM,YDELPE,YDELPI,YDELMI,YDELME
	integer	N,N1,ITC,IT1,IV1,NTET1,IEB,J1BEG,J1END,I,J,JT,JV,JN
	integer	JE,JBMS4,JNA,JNAC,J2,JSP,ISP,ITRAP,JTDTS,ITIME
	integer	JN22,IE,IVE,J1,JTIME,I1,I2,NA1,JLREC,JSRREC,JSRNUM,JDBL
	SAVE JFPBEG,JN1OLD
	DATA JFPBEG/0/
	DATA JN1OLD/1/
c	for double precision JDBL=2 (single precision 1) 30-MAY-06
	JDBL=2
	GP	=3.1415926536
	SQPI	=SQRT(GP)
c...Time step
	if(TAU.le.0.) goto 998
	DTION	=CBMI2*TAU
C	ITIME	=CBMI33
C	TAUDT0	=1./CBMI33
c...Flux surf.
	N1	=(NA1-1)/CBMI3+1
	N	=N1-1
c...Pitch angle
	ITC	=50/CBMI1
	NTET1	=ITC+1
	IT	=2*ITC
	IT1	=IT+1
C...Velocity	NBFP.OLD
c	IV	=160/CBMI4
	IV	=160
C...Velocity	NBFP.INC
C	IV	=80/CBMI4
	IV1	=IV+1
C...Beam energy components
	YFI(3)	=(IV*3)/4
	YFI(2)	=(YFI(3)*100)/142
	YFI(1)	=(YFI(3)*100)/173
	IEB	=3
C	JEB	=4-CBMS3
	JSRREC	=JDBL*4*(3*N*NTET1+1)
C...Clean previous distributions
	DO 1001	JN	=1,NA1
		PEBM(JN)=0.
		PIBM(JN)=0.
		CUBM(JN)=0.
		CUFI(JN)=0.
		PBPER(JN)=0.
		NIBM(JN)=0.
		PBLON(JN)=0.
 1001	CONTINUE
		JLREC	=JDBL*4*IV1*IT
		OPEN(31,FILE='dat/fij.dat',form='unformatted',
     ,		ACCESS='DIRECT',RECL=JLREC)
C....Zero values for initial conditions for distribution function F
	IF(N1.gt.JN1OLD) then
	   J1BEG   =JN1OLD
	   J1END   =N1
	else
	   J1END   =JN1OLD
	   J1BEG   =N1
	endif
	IF(J1BEG.ne.J1END) then

	DO 1	JN=J1BEG,J1END
		DO 10	I=1,IV1
		DO 10	J=1,IT
 10		FVM(I,J)	=0.
 1	WRITE(31,REC=JN,ERR=999) ((FVM(JV,JT),JV=1,IV1),JT=1,IT)
	endif

	JFPBEG	=1
	JN1OLD  =N1

C... Input sourse for CO & CONTR injection
	DO 11	JT=1,NTET1
	DO 11	JN=1,N1
	DO 11	JE=1,3
 11	YASBA(JE,JN,JT)=0.
	JBMS4	=CBM1
	CALL NBMESH
		EZ	=0.
		CNSTN	=HM*HV
		CNSTQ	=3.2E-3*EBEAM*HV**2*CNSTN
		CNSTP	=EBEAM*CNSTN*HV*HV
		CNSTC	=0.7*SQRT(EBEAM/ABEAM)*HV*CNSTN
		CNSTE	=1.6E-3*CNSTP/DTION
		CNSTT	=3.5E-2*EBEAM*SQRT(EBEAM*ABEAM)
		YEV21	=EBEAM/ABEAM
		YEV22	=YEV21/2.
		YEV23	=YEV21/3.
C...Change to program units
	YJ2	=0.
C...Radial distribution cicle
	OPEN(33,FILE='dat/scon.dat',form='unformatted',
     .		ACCESS='DIRECT',RECL=JSRREC)

	OPEN(34,FILE='dat/sctr.dat',form='unformatted',
     .		ACCESS='DIRECT',RECL=JSRREC)
	DO 2	JN	=1,N
		JNA	=1+CBMI3*(JN-1)
		JNAC	=JNA-1+CBMI3
		J2	=JNA-CBMI3*YJ2
	if(JN.EQ.N) JNAC=NA1-1
	IF(CBMI3.GT.1.d0)	YJ2	=0.5d0
		YLNE(JN)=15.85+LOG(TE(J2)/SQRT(NE(J2)))
	IF(EBEAM.GT.100.*ABEAM)	THEN
		YLNI(JN)=23.7+LOG(AMAIN(J2)/(AMAIN(J2)+ABEAM)*
     *             SQRT(1.E-3*ABEAM*EBEAM*TE(J2)/NE(J2)))
				ELSE
		YLNI(JN)=25.4+LOG(1.E-3*EBEAM*AMAIN(J2)/
     /		  (AMAIN(J2)+ABEAM)*SQRT(TE(J2)/NE(J2)))
				ENDIF
		YLNZ(JN)=YLNI(JN)
		YEPS	=AMETR(J2)/(RTOR+SHIF(J2))
		ITRAP	=(1.-SQRT(2.*YEPS/(1.+YEPS)))/HM
C=====Ripple losses
c				yyyyy=RIPCOS(ZUPDWN,JNA)
c				IRIP=RIPCOS(ZUPDWN,JNA)/HM
c				car16(jna)=yyyyy
C=============================
 		YFCUR(JN)=(1.-FNBF(ZEF(J2),YEPS)/ZEF(J2))
		PBCX(JN)=0.
		RNB(1)	=NE(J2)*YLNE(JN)
		EB(1)	=TE(J2)/EBEAM
		VB(1)	=SQRT(EB(1)*ABEAM/RMB(1))
	DO 20	JSP	=2,ISPEND
		ISP	=ISPE(JSP)
		RNB(JSP)=EXTARR(J2,ISP)*YLNI(JN)
		EB(JSP)	=TI(J2)/EBEAM
 20		VB(JSP)	=SQRT(EB(JSP)*ABEAM/RMB(JSP))
C...Time step DTAU[s]
		T0	=CNSTT/(NE(J2)*YLNE(JN))
C...Beam prtcls. slowing down time
		TSNBI	=2.*ABEAM/AMJ*SQRT(TE(J2))*TE(J2)/
     .			(14.78+LOG(TE(J2)))/NE(J2)
CBMI33 new is now const (21-FEB-97):
		CBMI33	=2.
		JTDTS	=CBMI33*DTION/TSNBI
	IF(JTDTS.GT.1)	THEN
		ITIME	=JTDTS
				ELSE
		ITIME	=1
				ENDIF
		DTAU	=DTION/ITIME
		DT	=T0/DTAU
C...Constants
		CNSTQT	=CNSTQ/T0
		EBDTI	=1./EB(2)
		EBDTI0	=EBDTI/DV2(1)
		CNSTE0	=EB(2)*SQRT(EB(2))*SQPI/2.
C		CNSTF0	=exp(EBDTI/DV2(1))
C...Sources and losses distributions...............................
cr						*2
cr		YNN0	=(NNCL+NNWM)*NN(J2)	*4.373E7*T0/2.
		CNSNN	=4.373E7*T0*CBM3        *.5d0
		CNSNN0	=(NNCL+NNWM)*4.373d7*T0*CBM4  *.5d0
	DO 	JV	=1,IV1
	DO 	JT	=1,IT
		RMN(JV,JT)=0.d0
 		FSRS(JV,JT)=0.d0
	ENDDO
	ENDDO
C...Fast ions CX due to cold neutrals
C...Fast ions CX due to NB neutrals
	IF(CBM4.GT.0.d0.or.CBM3.GT.0.d0) 		then
		YNN0	=CNSNN0*NN(J2)	+ 
     +			CNSNN*(NNBM1(J2)+NNBM2(J2)+NNBM3(J2))
	DO 21	JV	=1,IV1
		YV2	=YEV21/DV2(JV)
		YET	=MAX(TN(J2),YV2)
		YRMN	=YNN0*SPEX(YET)*SQRT(YV2)
	DO 21	JT	=1,IT
 21		RMN(JV,JT)=YRMN
 				endif
C...End of CX losses..................................................

		YSRSE	=0.
	DO 22 JSRNUM=1,JBMS4
	READ(33,REC=JSRNUM,ERR=211) 
     .	YEBEAM,(((YASBA(JE,JN22,JT),JE=1,3),JN22=1,N),JT=NTET1,IT1)
 211	READ(34,REC=JSRNUM,ERR=212)
     .	YEBEAM,(((YASBA(JE,JN22,JT),JE=1,3),JN22=1,N),JT=1,NTET1)
 212	continue
	DO 22	JT	=1,IT
	DO 22	IE	=1,3
	IF(YASBA(IE,JN,JT).GT.0.)	THEN
		IVE	=YFI(IE)*sqrt(YEBEAM/EBEAM)
		CNSFI(IE)	=DV2(IVE)/CNSTN*0.5
     .			*(YFI(IE)/IVE)**2*YEBEAM/EBEAM
C******************* CORRECTION OF POWER BALANCE
		FSRS(IVE,JT)	=YASBA(IE,JN,JT)*CNSFI(IE)*T0
     .				+ FSRS(IVE,JT)
		YSRSE	=YSRSE+YASBA(IE,JN,JT)*T0/CNSTN*0.5*IVE**2
					ENDIF
 22	continue
	DO 221	IE	=1,3
	IF(YASBA(IE,JN,JT).GT.0.)	THEN
		IVE	=YFI(IE)*sqrt(YEBEAM/EBEAM)
		CNSFI(IE)	=DV2(IVE)/CNSTN*0.5
     .			*(YFI(IE)/IVE)**2*YEBEAM/EBEAM
C******************* CORRECTION OF POWER BALANCE
		FSRS(IVE,IT)	=YASBA(IE,JN,IT1)*CNSFI(IE)*T0
     .				+ FSRS(IVE,IT)
	YSRSE	=YSRSE+YASBA(IE,JN,IT1)*T0/CNSTN*0.5*IVE**2
					ENDIF
 221	continue
		YSRSE	=YSRSE*CNSTQT
C...Coefficients
	CALL	NBCOEF
	READ(31,REC=JN) ((FVM(JV,JT),JV=1,IV1),JT=1,IT)
C...For dPb/dt
		YPB	=0.
	DO 224	I=1,IV
		YIP	=0.
C		YV4	=(I-.5)**2/DV2(I)
		YV4	=I**2/DV2(I)
	DO 2240	J=1,ITC
		J1	=IT-J+1
 2240		YIP	=YIP+FVM(I,J)+FVM(I,J1)
 224		YPB	=YPB+YV4*YIP
C...Fij	SWEEPING
	DO 99	JTIME=1,ITIME
		CALL	NBPOMU
 		CALL	NBPOVE
 99	CONTINUE
C...Fij linearization Fij b =Fij - Foj exp(-Ei/T)
		FVMMIN	=FVM(1,1)
	DO	990 J=2,IT
 		FVMMIN	=MIN(FVMMIN,FVM(1,J))
 990	CONTINUE
CCCC	IF(FVMMIN.LT.0.) WRITE(*,*) 'FVM<0,JN', JN,FVMMIN
		F0J	=FVMMIN*0.999999
		PTHBM	=F0J*TI(J2)*CNSTE0	*1.5
		PTHERM	=1.6E-3*PTHBM/DTION
		YE	=0.
C		IF(F0J.LT.1.E-15)	F0J=0.
	DO 991	JV	=1,IV1
	YE=0.
		YEXARG	=EBDTI0-EBDTI/DV2(JV)
	IF(YEXARG.GT.-30.)YE	=EXP(YEXARG)*F0J
	DO 991	JT	=1,IT
 		FVM(JV,JT)=FVM(JV,JT)-YE
CCCC	IF(FVM(JV,JT).LT.0.)WRITE(*,*) 'FVM<0,JN',FVM(JV,JT),YE,JN,JT
 991	CONTINUE
c===================================================== Ripple cone
c	if(IRIP.gt.0)	then
c	DO	JT	=ITC-IRIP,ITC+IRIP
c	DO	JV	=1,IV1
c 		FVM(JV,JT)=0.
c	ENDDO
c	ENDDO
c			endif
c===================================================== Ripple cone end
	WRITE(31,REC=JN) ((FVM(JV,JT),JV=1,IV1),JT=1,IT)
C...Power to plasma, beam pressure, density, current
		YIP	=0.
		YI	=0.
		YPEBM	=0.
		YPIBM	=0.
		YPBPER	=0.
		YPBLON	=0.
		YNB	=0.
		YCUFI	=0.
		YMF	=0.
		YM2F	=0.
	DO 23	J=1,ITC
		J1	=IT-J+1
		YM2F	=YM2F+YM2(J)*(FVM(1,J)+FVM(1,J1))
 23		YIP	=YIP+FVM(1,J)+FVM(1,J1)
	DO 24	J=1,ITRAP
		J1	=IT-J+1
 24		YMF	=YMF+YM1(J)*(FVM(1,J)-FVM(1,J1))
		YNB	=YNB+YIP/DV2(1)
C		YV4	=0.25/DV2(1)
		YV4	=1./DV2(1)
		YPBPER	=YPBPER+YV4*YM2F
		YPBLON	=YPBLON+YV4*(YIP-YM2F)
		YCUFI	=YCUFI+YMF/DV2(1)
			
		YDELPE	=0.
		YDELPI	=0.

	DO 25	I	=1,IV
		I1	=I+1
		I2	=I*I
		YI	=YIP
		YIP	=0.
		YMF	=0.
		YM2F	=0.
C		YV4	=(I1-.5)**2/DV2(I1)
		YV4	=(I1**2/DV2(I1))
	DO 250	J	=1,ITC
		J1	=IT-J+1
		YM2F	=YM2F+YM2(J)*(FVM(I1,J)+FVM(I1,J1))
 250 		YIP	=YIP+FVM(I1,J)+FVM(I1,J1)
	DO 251	J	=1,ITRAP
		J1	=IT-J+1
 251		YMF	=YMF+YM1(J)*(FVM(I1,J)-FVM(I1,J1))
		YDELMI	=YDELPI
		YDELPI	=BI(I)*YIP-AI(I)*YI
		YDELME	=YDELPE
		YDELPE	=BE(I)*YIP-AE(I)*YI
		YNB	=YNB+YIP/DV2(I1)
		YCUFI	=YCUFI+YMF*I1/DV2(I1)
		YPBPER	=YPBPER+YV4*YM2F
		YPBLON	=YPBLON+YV4*(YIP-YM2F)
cc		YPEBM	=YPEBM+I*(BE(I)*YIP-AE(I)*YI)
cc 25		YPIBM	=YPIBM+I*(BI(I)*YIP-AI(I)*YI)
		YPEBM	=YPEBM-I*I*(YDELPE-YDELME)
 25		YPIBM	=YPIBM-I*I*(YDELPI-YDELMI)
	do 26	J	=JNA,JNAC
C...Beam power distribution
	include 'fml/pbeie'
		PEBM(J)		=YPEBM*CNSTQT	/2.	-PBEIE
		PIBM(J)		=YPIBM*CNSTQT	/2.	+PTHERM
	if(CBM2.gt.0.)		then
	include 'fml/pbicx'
		PIBM(J)		=PIBM(J)-PBICX
				endif
cc		PEBM(J)		=YPEBM*CNSTQT
cc		PIBM(J)		=YPIBM*CNSTQT	+PTHERM
c...Beam energy Wbeam = Pbper + Pblon/2
c...Beam perpendicular pressure <Mb Vort2/2> [10*19 keV/m3]
		PBPER(J)	=YPBPER*CNSTP
c...Beam parallel pressure 	  <Mb Vpar2> [10*19 keV/m3]
		PBLON(J)	=2.*YPBLON*CNSTP
c...Beam density		 	      [10*19 /m3]
		NIBM(J)		=YNB*CNSTN
c...Fast ion current without trapping correction for Eb [MA/m2]
C		CAR10(J)	=T0
C		CAR11(J)	=YSRSE-(YPBLON+YPBPER-YPB)*CNSTE
cc		CAR12(j)	=ysrse
cc		car15(j)	=yasb
C		CAR13(J)	=PTHERM
 		CUFI(J)		=YCUFI*CNSTC
26	 	CUBM(J)		=CUFI(J)*YFCUR(JN)
2	CONTINUE
	CLOSE(31)
	CLOSE(33)
	CLOSE(34)

ctest
c	OPEN(35,FILE='DAT\YAS')
c	WRITE(35,*) (YASBA(1,1,JT),JT=1,IT1)
c		CLOSE(35)
c9997	open (47,file='a:\iters.dat')
c	write(47,*) ((YASBA(3,JN,JT),JN=1,N),JT=1,IT)
c	close (47)
 998	return
 999	write(*,*) 'error NBION2'
	return
	end
