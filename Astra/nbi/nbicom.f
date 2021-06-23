C====================================================================10-DEC-07
	subroutine	NB0(ADQB,JSRC,YHM,CH1,CH2,CR1,CR2,CS3)
C==== input beam footprint distribution  ================== 03-MAR-99 
C	Upward/Downward shift ZUPDWN
C	Perpendicular injection
C============================================================ 10-DEC-07
C       Extension to the rare mesh N1 
C============================================================ Polevoy
	implicit none
        include 'for/parameter.inc'
        include 'nbi/nbicom.inc'
        double precision	YHM,CH1,CH2,CR1,CR2,CS3,YS,Z,ZY,Y1,ZY12
        double precision	ADQB(*),nbfhz,nbfry,YR1,YR2,Y,YC,YD,YE
        double precision	YDR,YR0,YZ1,YDS,YSB,YHBD2,YHOCU,YNOCU
        double precision	YHTOP,YHTOPA,YHEDGE,YHEDGA,YD05H,YDZOUT
        double precision	YDHZ,YDZJH
	integer JSRC,J,JH,JR,JHH,JE,N,JH1,JHB,JHB05,JEDGE
        	N	=N1-1
		
		YR1	=RBMIN1/A
		YR2	=RBMAX1/A
		if(YR2.le.YR1) then
		   write(*,*) 
     & 'for NBI source N ',JSRC,' RMAX = ',RBMAX1,' < RMIN = ',RBMIN1
		   write(*,*) 'continue with (RMAX - RMIN)/A = 1.d-5'
		   YR2	=YR1+1.d-5	   
		   endif
		YC	=R/A
		Y	=YC+1.d0

		YHBD2	=abs(HB)/2.d0
		if(YHBD2.eq.0.d0) then
		   write(*,*) 'NBI source N ',JSRC,' vertical size = 0'
		   write(*,*) 'continue with vertical size = 1.e-5'
		   YHBD2	=1.d-5	   
		   endif
	
CCC	IF(YR1.LT.-Y) YR1	=-Y
CCC	IF(YR2.GT.Y)  YR2	=Y
		YDR	=(YR2-YR1)/IRB
		YD	=2.d0/(abs(YR2-YR1))
		YR0	=0.5d0*(YR1+YR2)
C		Hocu	=(HBEAM-ZUPDWN)*100.
		Hocu	=YHM
	do      JH      =1,N1
	        YDZ(JH) =0.d0
	enddo
cccc	if(JTANG.EQ.0)		then
C...Perpendicular NBI: calculation of the NBI height
		Y	=Hocu/A
c17JUN02vvvvvvvvvvvvvvvvvv
		YHocu=Y
	if(YR0.lt.(YC+DX(1))) then
		J	=1
 1		YHocu	=Y + CS3*sqrt((YC+DX(J)+YR0)*(YC+DX(J)-YR0))
	if(abs(YHocu).lt.X(J)*ELON1(J).and.J.lt.N1)	then
		J	=J+1
		goto	1
							endif
	endif
c17JUN02^^^^^^^^^^^^^^^^^^^^
		Hocu	=A*YHocu
cccc				endif			

		YS	=0.d0
	do 10	JR	=1,IRB
		AZ(JR)	=(JR-0.5d0)*YDR+YR1
 		Y	=(ABS(AZ(JR)-YR0)+1.d-7)*YD
     		YDRY(JR)=NBFRY(Y,JSRC,CR1,CR2)
 10		YS	=YS+YDRY(JR)
		YS	=(IRB/YS)
	do 11	JR	=1,IRB
 11		YDRY(JR)=YDRY(JR)*YS

		YHocu	=abs(Hocu)
		YHTOP	=YHocu+YHBD2
		YHTOPA  =YHTOP/A
		YNocu	=YHocu/YHBD2
C*FEB*99		YNTOP	=YNocu+1.

		YHEDGA   =X(N1)*ELON1(N1)
		YHEDGE  =A*YHEDGA
		JHB05   =0
	if(YHEDGE.lt.(YHocu-YHBD2)) then
C NBI out of plasma
	        IHB     =0
	   RETURN
	                            endif
		YD05H	=A/YHBD2
		YS	=0.d0


	if(YHTOP.gt.YHEDGE) then
C NBI is partly out of plasma
	        JEDGE   =(YHTOP/YHEDGE-1.d0)*(N1-1)

CW 
CW		write(*,*) 'N1,IRB,JEDGE', N1,IRB,JEDGE

	if(JEDGE.lt.1) JEDGE =1
	        YDZOUT  =(YHTOP-YHEDGE)/(JEDGE*A)
		ZY      =YHTOPA
	do 	JHH	=1,JEDGE
		YZ1	=ZY
ccc		JH1	=JH-1
		ZY	=YZ1-YDZOUT
		ZY12    =0.5d0*(ZY+YZ1)
		Z	=abs(ZY12*YD05H-YNocu)
		YDHZ	=0.d0
	IF(Z.LE.1.d0)		THEN
		YDHZ	=NBFHZ(Z,JSRC,CH1,CH2)
		Z	=abs(ZY12*YD05H+YNocu)
	IF(Z.LE.1.d0)	YDHZ	=YDHZ+NBFHZ(Z,JSRC,CH1,CH2)

				ENDIF
	YDZJH	=(YZ1-ZY)*YDHZ
		YS	=YS+YDZJH
	enddo
	        JHB     =N1
	        JHB05   =N1-1
	endif
	
			      
C NBI completely inside the plasma
C looking for the NBI edge
		YZ1     =YHEDGA
		ZY      =YZ1
		JHB      =N1

 3		YZ1	=X(JHB)*ELON1(JHB)
	if(YZ1.gt.YHTOPA.and.JHB.gt.1) then
	        JHB     =JHB-1
		ZY      =YZ1
		goto 3
	endif
	        JHB05   =JHB
C*DEC*07		JHB     =JHB05+1

	IF(JHB.GT.1)THEN

	do 	JHH	=2,JHB05+1
	        JH1     =JHB05-JHH+1
		YZ1	=ZY
		JH	=JH1+1
		ZY	=X(JH)*ELON1(JH)
ccTEST	if(JH.le.0) write(*,*) JHH,JHB,JH
		ZY12    =0.5d0*(ZY+YZ1)
		Z	=abs(ZY12*YD05H-YNocu)
		YDHZ	=0.d0
	IF(Z.LE.1.d0)		THEN
		YDHZ	=NBFHZ(Z,JSRC,CH1,CH2)
		Z	=abs(ZY12*YD05H+YNocu)
	IF(Z.LE.1.d0)	YDHZ	=YDHZ+NBFHZ(Z,JSRC,CH1,CH2)

				ENDIF
		YDZ(JH)	=(YZ1-ZY)*YDHZ
		YS	=YS+YDZ(JH)
	enddo
	ENDIF
C*DEC*07
	if(YS.le.0.d0) then
	   Z	=0.d0
	   JH	=JHB05
	   YDHZ	=YDHZ+NBFHZ(Z,JSRC,CH1,CH2)
		YDZ(JH)	=YDHZ*2.d0*YHBD2
		YS	=YS+YDZ(JH)
		endif
C*DEC*07
		YDS	=YHBD2/(A*YS)
	do 	JH	=1,JHB05
 		YDZ(JH)	=YDZ(JH)*YDS
		enddo

		IHB	=2*JHB05
C*DEC*07		YSB	=(YR2-YR1)*A*HB
		YSB	=(YR2-YR1)*A*2.d0*YHBD2 
		Y	=2.d0*451.9d0*QB/YSB
		Y1	=EB/PB
	do 	JE	=JEB,IEB
		YE	=Y1/(IEB-JE+1)
 		YAQB(JE)=Y*ADQB(JE)/(YE*PB*dSQRT(YE))
	enddo
	return
	END
C====================================================================
        subroutine NB1TR(YAQBP,CONTR,AR,JSRNUM)
C====== Neutral beam ionization ================ Version by 08-JUL-97
C coinjection  with ion's trapping and orbital losses
C  AQBA[MW],ACBA[MA m/s],ANBA[10*13 prtcls],
C  ASBA[10#19 prtcl/s]*Dcos(JT)
C  YQSHth [MW] shine through power in the file dat\shth.dat
C============================================================ Polevoy
	implicit none
	double precision YAQBP(*),AR(*),DRL(3),YS(3),YCU(3),YCT2(3)
        double precision	PLEJ2(51),CONTR,YRIPLR,YDEDJ,YVE(3)
	include	 'for/parameter.inc'
	include 'nbi/nbicom.inc'
c	COMMON /CNBFOA/YPSI(3,201),RLM(3,201),RJM(201),
c     ,	RCR,RJ,RJT,JIM0R,JIM0L,RM0R,RM0L,YFM0R,YFM0L,YFM0T,YF0,
c     ,	YDYH,JNR,JNL,JE,II(201),jfig,ixyu,SQRCR,JCENTR
	include 'nbi/nbfoa.inc'
c	include 'nblej2.inc'
	character*12 SHNAME
c=====Ripple	normalized radius of the ripple boundary
c==== banana with RTRAP > YRIPLR is lost 
	common	/CRIPL/	YRIPLR(NRD)
        double precision Y0,Y1,Y2,Y,YE,YH,DS,DV,YQSHTH,YLOSS,YVEDE,YDYDT
        double precision YE3,RJ0,YR,YDH,XJH,ZJH,YZJH,Y12,YR1,YR2,YC1,YC2
        double precision YDT,YDY,YCDV,YDYS,YDEX,YQBP,YFA,YDDD,YD,YDCOS,
     ,	YTCOS
	integer JSRNUM,N,NTET,JI,JN,JT,JR,JJN,JJ,JHB05,JN1,IN,IN1
	integer JJH,JH,JH1,JJR,JBB,JII,JNN,JHP1,JETH,JT1,JT2
	DATA SHNAME/'dat/shth.dat'/
	DATA PLEJ2/
     .-5.000000d-01,-4.994000d-01,-4.976000d-01,-4.946000d-01,
     .-4.904000d-01,-4.850000d-01,-4.784000d-01,-4.706000d-01,
     .-4.616000d-01,-4.514000d-01,-4.400000d-01,-4.274000d-01,
     .-4.136000d-01,-3.986000d-01,-3.824000d-01,-3.650000d-01,
     .-3.464000d-01,-3.266000d-01,-3.056000d-01,-2.834000d-01,
     .-2.600000d-01,-2.354000d-01,-2.096000d-01,-1.826000d-01,
     .-1.544000d-01,-1.250000d-01,-9.440002d-02,-6.260002d-02,
     .-2.960002d-02, 4.599977d-03, 3.999998d-02, 7.659998d-02,
     . 1.144000d-01, 1.534000d-01, 1.936000d-01, 2.350000d-01,
     . 2.776000d-01, 3.214000d-01, 3.664000d-01, 4.126000d-01,
     . 4.599999d-01, 5.085999d-01, 5.584000d-01, 6.094000d-01,
     . 6.615999d-01, 7.150000d-01, 7.695999d-01, 8.253999d-01,
     . 8.823999d-01, 9.405999d-01, 9.999999d-01/

c=====Ripple
	if(JSRNUM.gt.99) then
	 write(*,*) '>>> Too many NBI sources > 99 (see NBI1TR)'
		return
			endif			
	if(JSRNUM.lt.10) write(SHNAME(12:12),98) JSRNUM
	if(JSRNUM.gt.9.and.JSRNUM.lt.100) write(SHNAME(11:12),99) JSRNUM
	if(JSRNUM.gt.99) then
	 write(*,*) '>>> Too many NBI sources > 99 (see NBI1TR)'
		return
			endif			
 98	FORMAT(1I1)
 99	FORMAT(1I2)
 
			open(38,FILE=SHNAME,STATUS='UNKNOWN')
			YQSHTH=0.d0
cc			write(1,*) 'R[m]      Z[m]      Q/S[W/cmZ]'
    		NTET	=NTET1-1
		N	=N1-1
		Y0	=R/A
		DS	=ABS(RBMAX1-RBMIN1)*A/(N*IRB)
		DV	=DS*A
		YLOSS	=0.d0
Clean the sources
	do 1	JN	=1,N1
	do 10	JE	=JEB,IEB
		YANBA(JE,JN)=0.d0
		YACBA(JE,JN)=0.d0
		YAQBA(JE,JN)=0.d0
		YATBA(JE,JN)=0.d0
	do 10	JT	=1,NTET1
		YASBA(JE,JN,JT)=0.d0
 10	continue
		RC(JN)	=Y0+DX(JN)
 1	continue
c*TRACE		write(*,'(20I4)')IHB,YH
c*TRACE		read(*,*)
	if(IHB.eq.0) return
c NBI out of plasma
		JHB05	=IHB/2
		YE3	=EB/PB
cc	include 'nbtrap.inc'
c for calculation of surface index II=JN(X(Rcrit)) in trapping analysis
c YR=(RJ-Rii/AB) normalised distance from ext. boundary RJ
 		RJ	=Y0+1.d0
		RJT	=Y0-1.d0
		RJ0	=Y0+DX(1)
		YDYH	=100.0d0
		YH	=0.01d0
		JN1	=N1
	DO 8	JI	=1,200
		YR	=YH*(JI-1)
 81	IF(JN1.GT.1)	THEN
		JN	=JN1-1
	 IF(YR.GE.(X(N1)+DX(N1)-DX(JN1)-X(JN1))
     .		.AND.YR.LT.(X(N1)+DX(N1)-DX(JN)-X(JN))) THEN
		II(JI)	=JN
	 ELSE
	 	JN1	=JN
		go to	81
	 ENDIF
	ELSE
	 IF(JN1.EQ.1)	JN1=-1
 82		IN1	=-JN1
		IN	=IN1+1
	 IF(YR.GE.(X(N1)+DX(N1)+X(IN1)-DX(IN1))
     .		.AND.YR.LT.(X(N1)+DX(N1)+X(IN)-DX(IN))) THEN
		II(JI)	=IN1
	 ELSE
	 	JN1	=JN1-1
		go to 82
	 ENDIF
	ENDIF
 8	CONTINUE
		II(201) =N
c*TRACE		write(*,'(20I4)')(j,II(j),j=1,201)
c*TRACE		read(*,*)
	do 80	JE	=JEB,IEB
		DRL(JE)=1./(0.0144d0*SQRT(EB*PB/(IEB-JE+1))/(BZ*A))
	do 80	JN	=1,N1
		ARD(JE,JN)=AR(JN)*DRL(JE)
 80	continue
C===End trapping arrays
cc	include 'nbtrag.inc'

	do 9	JE	=JEB,IEB
		YE	=YE3/(IEB-JE+1)
  		YVEDE	=DS*SQRT(YE)/451.9d0
		YVE(JE)	=YVEDE*PB*YE*1.d-3	
 		YCU(JE)	=-CONTR*VNB(JE)*YVEDE*5.d-6
 9	CONTINUE

CHH 2 - loop beam`s height H = X(JN)+dX/2  HHHHHHHHHHHHHHHHH
	do 2	JJH	=1,JHB05
		JH	=JHB05-JJH+1
		JH1	=JH+1
		Y	=ELON1(JH1)*X(JH1)-ELON1(JH)*X(JH)
	if(JH.EQ.JHB05)	THEN
		YDH	=N*YDZ(JH)
		XJH	=X(JH)+0.5d0*YDH/(N*ELON1(JH))
		RE(JH)	=RC(JH)-(RC(JH1)-RC(JH))*(XJH-X(JH))/
     /			 (X(JH1)-X(JH))-TRIA1(JH)
			ELSE
		XJH	=X(JH)+0.5d0*Y/ELON1(JH)
		YDH	=N*YDZ(JH)
		RE(JH)	=0.5d0*(RC(JH1)+RC(JH)-TRIA1(JH)-TRIA1(JH1)) 
			ENDIF
		ZJH	=XJH*ELON1(JH)
		RI(JH)	=RE(JH)
	do 21	JJR	=JH,N
		JN	=N1-JJR+JH
		YZJH	=ZJH/ELON1(JN)
		Y	=SQRT((X(JN)-YZJH)*(X(JN)+YZJH))
		RE(JN)	=RC(JN)+Y-TRIA1(JN)*(YZJH/X(JN))**2
 21		RI(JN)	=RC(JN)-Y-TRIA1(JN)*(YZJH/X(JN))**2
C*NEW BCTABKA
	do JN=JH+1,N1    
	if(jn.gt.1.and.RI(JN-1).lt.RI(JN)) write(*,*) JN,JH,JHB05,N1
	enddo 
c 21	continue
	do 20	JN	=JH,N1
		DRE(JN)	=1.d0/RE(JN)
		DRI(JN)	=1.d0/RI(JN)
 20	continue
CZZ 22 - loop beam`s radius R = Z(JR)    ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	do 22	JR	=1,IRB
	if(YAQB(3).le.0.)	goto 23
C---- no power in this pencil ------------------ go to 23 -------------
		JBB	=JEB
	do 221	JE	=JBB,IEB
		YS(JE)	=0.d0
 		YAQBP(JE)=YAQB(JE)*YDH*YDRY(JR)
 221	continue
CYY 220  integral along beam  Y		 YYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
C<<<<<<<<<  R  decreases <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
		Y12	=(RE(N1)-AZ(JR))*(RE(N1)+AZ(JR))
	if(Y12.LE.0.)	GO TO 23
C---- beam beyond tokamak ---------------------- go to 23 -------------
		Y2	=-SQRT(Y12)
		YC2	=AZ(JR)*DRE(N1)
		JJN	=N1-JH
		RCR	=Y12*DRE(N1)
		YR2	=CONTR*SQRT(RJ*(RJ-RCR))
	if(RCR.LT.RJT)	THEN
		YR1	=CONTR*SQRT(RJT*(RJT-RCR))
	ELSE
		YR1	=0.d0
	ENDIF
	do 2200	JE	=JEB,IEB
		YCT2(JE)=0.d0
 2200	continue
	IF(CONTR.LT.0.)	THEN
c	include 'nbtcon.inc'
	include 'nbi/nbco_0.rip'
	ELSE
	include 'nbi/nbctr_0.rip'
	ENDIF
ctest	write(*,*) 'passed'
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 220	continue
CYYYYY end of integral along beam YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
cc	Shinethrough power distribution R[m],Z[m],Q/S[W/cm2]	
	do	JETH=JEB,IEB
		YQSHTH	=YQSHTH+YAQBP(JETH)*YVE(JETH)
	enddo
cc			YYYD	=(YAQBP(1)*YVE(1)+YAQBP(2)*YVE(2)
cc     +	+YAQBP(3)*YVE(3))*0.5E4/(DS*YDH)
cc			write(1,*) AZ(JR)*A*1.d-2,ZJH*A*1.d-2,YYYD
 22	continue
 23	continue
CZZZZZ end of loop on beam`s radius R = Z(JR) ZZZZZZZZZZZZZZZZZZZZZZZZZ
 2	continue
CHHHHH end of loop on beam`s height HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
	DO 3	JE=JEB,IEB
	DO 3	JN=1,N1
		YAQBA(JE,JN)	=YAQBA(JE,JN)*YVE(JE)
		YANBA(JE,JN)	=YANBA(JE,JN)*YVE(JE)
 		YACBA(JE,JN)	=YACBA(JE,JN)*YCU(JE)
 		YATBA(JE,JN)	=YATBA(JE,JN)*YCU(JE)
 3	continue
		YDEDJ	=1.d6/(1.6d0*EB)
	DO 4	JE=JEB,IEB
		Y	=(IEB-JE+1)*YDEDJ*YVE(JE)
	DO 4	JN=1,N1
		YSLEJ0(JE,JN)=0.0d0
		YSLEJ2(JE,JN)=0.0d0
	DO 4	JT=1,NTET1
	IF(YASBA(JE,JN,JT).GT.0.)	THEN
 		YASBA(JE,JN,JT)	=YASBA(JE,JN,JT)*Y
     	YSLEJ2(JE,JN)=YSLEJ2(JE,JN)+YASBA(JE,JN,JT)*2.5d0*PLEJ2(JT)
     	YSLEJ0(JE,JN)=YSLEJ0(JE,JN)+0.5d0*YASBA(JE,JN,JT)
					ENDIF
 4	CONTINUE
cc>>>>>Shinethrough
			
	write(38,*) YQSHTH,' Q shine through [MW] for source ',JSRNUM
			CLOSE(38)
cc>>>>>Shinethrough end
	return
	END
C=======================================================================
	subroutine NBSRSR(JSRNUM,YCONTR,
     ,		NA1,	RTOR,	SHIFT,	AB,	BTOR,
     ,		NNCL,	NNWM,	HRO,	YHM,
     ,		CBMH1, 	CBMH2,	CBMS1, 	CBMS2,	CBMS3, 	CBMS4,
     ,		CBMR1, 	CBMR2, 	CBMI3, 	CBMI1,
     ,		EBEAM,	DBM1,	DBM2,	DBM3,
     ,		ABEAM,	CONTR,	QBEAM,	RBMAX,	RBMIN)
C=========================================================== 06-APR-98
C fast ion's sourses (for multi sources) + ripple	    
C double precision 					     30-MAY-06
C--------------------------------------------------------------Polevoy
C	entry:	AMETR,SHIF,NA1,RTOR,AB,BTOR,NI,HBEAM,RBMIN,RBMAX,
C     		CBMS2,CBMI3,EBEAM,QBEAM,ABEAM,DBM1,DBM2,DBM3,
C	     	MU,TE,TI,NE,NN,NNB,ZEF,AMAIN,PBEAM,SCUBM,CONTR,ELON,
C	CBMS2	number of 'pencils'
C	CBMI3	number of internal mesh points
C			 41 (CBMI3=1),21 (CBMI3=2)
C	DBM3,2,1 power fraction of energy comps.
C		    	 3(EB,EB/2,EB/3),2(EB,EB/2),1(EB)
C	CONTR	Qcontr/Qbeam
c	JSRNUM 	Number of the current hot ion source
c	JSRREC 	Length of the hot ion source record
C	exit:   COMMONS for HOTION
C	exit:	PBEAM,SCUBM,SNEBM,SNNBM		for MAIN
C	SCUBM	Toroidal pulse [kg*m/s2/m3]	05-AUG-96
c======================================================================
	implicit none
	double precision
     ,		RTOR,	SHIFT,	AB,	BTOR,
     ,		NNCL,	NNWM,	HRO,	YHM,
     ,		CBMH1, 	CBMH2,	CBMS1, 	CBMS2,	CBMS3, 	CBMS4,
     ,		CBMR1, 	CBMR2, 	CBMI3, 	CBMI1,
     ,		EBEAM,	DBM1,	DBM2,	DBM3,
     ,		ABEAM,	CONTR,	QBEAM,	RBMAX,	RBMIN
	include 'for/parameter.inc'
	include 'for/status.inc'
	include 'nbi/nbicom.inc'
        double precision	yVB,RNB,ZB,RMB,yEB,DT,EZ,YCONTR,YSIMPI
	COMMON /CNBCH/yVB(9),RNB(9),ZB(9),RMB(9),yEB(9),DT,EZ
	INTEGER JSRREC,JSRNUM,ISPEND,ISPE,JDBL
	COMMON /CNBSP/ ISPEND,ISPE(9)
	COMMON /CNBSPI/ YSIMPI(3,9)
	double precision	NNB,AQBP(3),ADQB(3),AR(NRD),GP2,YE2,YCX
	double precision	YZERO(3,NRD,51),SEIV,YSCU1,YFVDA,YNHDT
        double precision	YDBM,YD,YJN,YSCU,YPOW,YDV,YDDV,YJE,YE1
        integer	N,J,NTET,JN,JN1,JNA,JNAX,JE,JT,JNA1,JNAC,JS,JSP,NA1
c	for double precision JDBL=2 (single precision 1) 30-MAY-06
	JDBL=2

	GP2 = 6.283185d0
	N1	=(NA1-1)/CBMI3+1
	N	=N1-1
  	if(RBMAX.LE.RBMIN)	then
	write(*,*) 'ILLEGAL: NBI source N',
     ,					JSRNUM,' RBMAX <= RBMIN !!!'
			return
			endif
  	if( EBEAM*ABEAM.eq.0.)	then
	write(*,*) 'ILLEGAL: NBI source N',
     ,					JSRNUM,' EBEAM*ABEAM = 0 !!!'
			return
			endif
		IEB	=3
        	IRB	=CBMS2
C...IF DBM1,2,3 - SOURCE CURRENT FRACTIONS
C		Z	=1./(DBM1+0.5*DBM2+DBM3/3.)
C		ADQB(1)	=Z*DBM3/3.
C		ADQB(2)	=Z*DBM2*0.5
C		ADQB(3)	=Z*DBM1
		YDBM	=DBM1+DBM2+DBM3
C*23MAY00 Filling by default 0.
	do j=1,3
		ADQB(J) =0.d0
	enddo
  	if(DBM1.LE.0.d0)	then
	write(*,*) 'ILLEGAL: NBI source N',JSRNUM,' DBM1 <= 0 !!!'
			return
			endif
		ADQB(3)	=DBM1/YDBM
		JEB	=IEB
  	if(DBM2.GT.0.d0)	then
		ADQB(2)	=DBM2/YDBM
		JEB	=IEB-1
			endif
  	if(DBM3.GT.0.d0)	then
  		ADQB(1)	=DBM3/YDBM
		JEB	=IEB-2
			endif
		NTET1	=50/CBMI1+1
		NTET	=NTET1-1
CW	write(*,*) 'nbsrs in',YDBM,NTET

	if(CBMI1.GT.1.d0)				THEN
 		JSRREC	=JDBL*4*(3*N*NTET1+1)
	OPEN(33,FILE='dat/scon.dat',FORM='UNFORMATTED',STATUS='UNKNOWN',
     .		ACCESS='DIRECT',RECL=JSRREC)
	OPEN(34,FILE='dat/sctr.dat',FORM='UNFORMATTED',STATUS='UNKNOWN',
     .		ACCESS='DIRECT',RECL=JSRREC)
						ENDIF
		X(1)	=0.d0
		XJ(1)	=0.d0
		XJ(NA1)	=1.d0
		AR(1)	=0.d0
        	EB	=EBEAM*1000.d0
	       	QB	=QBEAM*1000.d0
        	PB	=ABEAM
	      	R	=(RTOR+SHIFT)*100.d0
        	A	=50.d0*(AMETR(NA1-1)+AMETR(NA1))
	if(((RBMAX+RBMIN)/2.d0).gt.(RTOR-AB))	then
C== Tangentional NBI
		JTANG	=1
						else
C== Perpendicular NBI
		JTANG	=0
						endif						
        	BZ	=BTOR/(1+SHIFT/RTOR)
C*FEB*99
        	HB	=CBMS4*(RBMAX-RBMIN)*100.*sqrt(1.+CBMS3*CBMS3)
	if(HB.le.0.)then 
	write(*,*) 'NBI input for the source ',JSRNUM,' is not correct'
        write(*,*) 'please use: RBMAX > RBMIN and CBMS4 > 0'
	endif
		NNB	=NNCL+NNWM
cc	IF(HB.GT.2.*A*ELONG)WRITE(*,*)' Vertical NBI aperture>2*A*Elon'
        	RBMIN1	=RBMIN*100.
        	RBMAX1	=RBMAX*100.
cc		YD	=1./(BZ*GP2*ABC*ABC)
		YD	=10000.d0/(BZ*GP2*A*A)
	do 10	JN	=2,NA1
		JN1	=JN-1
		AR(JN)	=(FP(JN)-FP(1))*YD
 10 		XJ(JN)	=50.d0*(AMETR(JN)+AMETR(JN1))/A
			YJN		=0.d0
        do 1		JN		=1,N1
			JNAX		=1+CBMI3*(JN-1)
			JNA		=JNAX-CBMI3*YJN
	IF(CBMI3.GT.1.d0)	YJN		=.5d0
        		PM(JN)		=max(AMAIN(JNA),1.d0)
        		PLI(JN)		=max(NI(JNA),1.d-4)
			ELON1(JN)	=max(ELON(JNAX),1.d0)
			TRIA1(JN)	=TRIA(JNAX)*XJ(JNAX)
			X(JN)		=XJ(JNAX)
			AR(JN)		=AR(JNA)
cc        		DX(JN)		=SHIF(JNAX)/AB
          		DX(JN)		=100d0*(SHIF(JNAX)-SHIFT)/A
        		AMU(JN)		=max(MU(JNA),1.d-4)
        		PL(JN)		=max(NE(JNA),1.d-4)
        		PN0(JN)		=max(NN(JNA)*(NNWM+NNCL),0.d0)
	do 1		JE		=1,3
CJEB,IEB
 			YAQBA(JE,JN)	=0.d0
 			YACBA(JE,JN)	=0.d0
 			YATBA(JE,JN)	=0.d0
 			YANBA(JE,JN)	=0.d0
			YSLEJ0(JE,JN)	=0.0d0
			YSLEJ2(JE,JN)	=0.0d0
	do 1		JT		=1,NTET1
			YASBA(JE,JN,JT)	=0.d0
			YZERO(JE,JN,JT)	=0.d0
 1	CONTINUE

		X(N1)	=1.0d0
		DX(N1)	=0.d0
C...IF DBM1,2,3 - SOURCE CURRENT FRACTIONS
C		Z	=1.d0/(DBM1+0.5d0*DBM2+DBM3/3.d0)
C		ADQB(1)	=Z*DBM3/3.d0
C		ADQB(2)	=Z*DBM2*0.5d0
C		ADQB(3)	=Z*DBM1
	if(QB.LE.0.d0)	goto 999

	call	NBSISN(NA1,EBEAM/ABEAM,CBMI3)
C	YHM = (HBEAM-UPDWN*(1.-CBMH4)-CBMH3)*100.
	call	NB0(ADQB,JSRNUM,YHM,CBMH1,CBMH2,CBMR1,CBMR2,CBMS3)

	if (CBMS1.lt.1.d0)	then
	call	NB1TR(AQBP,YCONTR,AR,JSRNUM)
				else
	call	NBTRAG(AQBP,YCONTR,AR,JSRNUM,CBMS1)
				endif				
		YSCU	=2.d3/(9.79d0*DSQRT(2000.d0*EBEAM/ABEAM))
	YPOW=0.d0
	
 	do 3	JN	=2,N1
		JN1	=JN-1
		JNA1	=1+CBMI3*(JN1-1)
		JNA	=JNA1+1
		JNAC	=JNA1-1+CBMI3
C*06APR00vvvvvvvvvvvvvv
	if(JN.eq.N1)then
		JNAC=NA1-1
		YDV=VR(JNAC)*(RHO(NA1)-JNAC*HRO)
	else
		YDV=0.d0
	endif
		JNAX	=JNAC
C*06APR00^^^^^^^^^^^^^^^
	do 30	J	=JNA1,JNAC
  30	 	YDV	=VR(J)*HRO+YDV
		YDDV	=1./YDV
	do 31	JE	=JEB,IEB
		YSLEJ0(JE,JN1)	=YSLEJ0(JE,JN)*YDDV
		YSLEJ2(JE,JN1)	=YSLEJ2(JE,JN)*YDDV
 		YAQBA(JE,JN1)	=YAQBA(JE,JN)*YDDV
 		YACBA(JE,JN1)	=YACBA(JE,JN)*YDDV
 		YATBA(JE,JN1)	=YATBA(JE,JN)*YDDV
		YJE		=JE
c		YSCU1		=YSCU*SQRT(YJE)
c*26FEB2003vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
		YSCU1		=ABEAM*.0209d0*.5d0
cYSCU*SQRT(YJE)^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
 	DO 311	JT=1,NTET1
CC	IF(YASBA(JE,JN,JT).GT.0.)	THEN
		YASBA(JE,JN1,JT)=YASBA(JE,JN,JT)*YDDV
CC					ENDIF
 311	CONTINUE
	do 310	J	=JNA1,JNAC
		PBEAM(J)=PBEAM(J)+YAQBA(JE,JN1)
  310		SCUBM(J)=SCUBM(J)+YATBA(JE,JN1)*YSCU1
	YANBA(JE,JN1)=YANBA(JE,JN)*YDDV/EBEAM*625.d0*(IEB-JE+1)
  31	CONTINUE
	do 32	J	=JNA1,JNAC
		YFVDA	=F(1,JN1)*VNB(1)/A
ccc	write(*,*) YFDVA
  	if(YFVDA.NE.0.d0)
     . 		NNBM1(J)= NNBM1(J)+YANBA(3,JN1)/YFVDA
		YFVDA	=F(2,JN1)*VNB(2)/A
  	if(YFVDA.NE.0.d0)
     .		NNBM2(J)= NNBM3(J)+YANBA(2,JN1)/YFVDA
		YFVDA	=F(3,JN1)*VNB(3)/A
  	if(YFVDA.NE.0.d0)
     .  	NNBM3(J)= NNBM3(J)+YANBA(1,JN1)/YFVDA
C...Calc. of total proton content:	YNHDT
		YNHDT	=0.d0
	DO 320	JS	=2,ISPEND
		JSP	=ISPE(JS)
	IF(ZB(JS).EQ.1.)
     .		YNHDT	=YNHDT	+EXTARR(JNA,JSP)
 320	CONTINUE
C...Calc. total particle sourse:	SNEBM
C 		YE1	=NE(JNAX)*SEIV(TE(JNAX))
 		YE1	=NE(JNA)*SEIV(TE(JNA))
 	DO 321	JE	=JEB,IEB
		YE2	=0.d0
	DO 3210	JS	=2,ISPEND
		JSP	=ISPE(JS)
 		YE2	=YE2+
     +	 EXTARR(JNA,JSP)*YSIMPI(JE,JS)*VNB(JE)
 3210	CONTINUE
C...Calc. thermal neutral  sourCe:	SNNBM
		YCX	=YNHDT*SVEX(4,JE)
C*15-DEC
	SNNBM(J)=SNNBM(J)+YCX/(YCX+YE1+YE2)*YANBA(JE,JN1)
c*3-SEP-03		SNNBM(J)=SNNBM(J)+YCX/(YE1+YE2)*YANBA(JE,JN1)
C...Calc. electron particle sourCe:	SNEBM
c*3-SEP-03		SNEBM(J)=SNEBM(J)+(1-YCX/(YE1+YE2))*YANBA(JE,JN1)
C*15-DEC
	SNEBM(J)=SNEBM(J)+(1d0-YCX/(YCX+YE1+YE2))*YANBA(JE,JN1)
 321	CONTINUE
 32		CONTINUE
  3	continue

 999	CONTINUE
	if(CBMI1.GT.1.d0)				THEN
		if(YCONTR.GT.0.d0)			then
		WRITE(34,REC=JSRNUM) EBEAM,(((YASBA(JE,JN,JT),
     ,				JE=1,3),JN=1,N),JT=NTET1,1,-1)
	if(CONTR.eq.1.d0)WRITE(33,REC=JSRNUM) EBEAM,(((YZERO(JE,JN,JT),
     ,				JE=1,3),JN=1,N),JT=1,NTET1)
    						else
		WRITE(33,REC=JSRNUM) EBEAM,(((YASBA(JE,JN,JT),
     ,				JE=1,3),JN=1,N),JT=1,NTET1)
	if(CONTR.eq.0.d0)WRITE(34,REC=JSRNUM) EBEAM,(((YZERO(JE,JN,JT),
     ,				JE=1,3),JN=1,N),JT=NTET1,1,-1)
						endif
		close(33)
		close(34)
						ENDIF
	return
 	END
C======================================================================
	subroutine	NBION0(NA1,NNCL,NNWM,ABEAM,EBEAM,RTOR,CBMI3)
c================================================ Version by: 21-FEB-97
C Steady State (1+2D:(x+MU,V)) Fokker-Plank Solver 
C	PEBM,PIBM(X)-power to electrons,ions [MW/m3]	
C	CUFI,CUBM(X)- fast ion's and NBI driven currents [MA/m2]
C	NIBM(x) - fast ion's density [10^19/m3]
C============================================================== Polevoi
C	common/CION0/
	implicit none
	include  'for/parameter.inc'
	include  'for/status.inc'
	include  'nbi/nbicom.inc'
	double precision PBCX(NRD),YZ2D3(NRD),YTSE(NRD),YFCUR(NRD)
	double precision YLNI(NRD),YLNE(NRD),YLNZ(NRD),FNBF,FNBI1,FNB2
	double precision YEBDEC(NRD),RTOR,CBMI3,YC1,YC2,PBEIE,FNBP
        double precision yVB,RNB,ZB,RMB,yEB,DT,EZ,NNCL,NNWM,ABEAM,EBEAM
	COMMON /CNBCH/yVB(9),RNB(9),ZB(9),RMB(9),yEB(9),DT,EZ
        integer	ISPEND,ISPE,JN,JS,JSP,JN1,JNA,JNA1,JNAC,J,J2,JE,NA1
	common /CNBSP/ ISPEND,ISPE(9)
        double precision YS,YSTE,YEPS,YN0,YDN,Y,X1,X2,X3,YA,YB,YC,YC0
        double precision Y1,Y12,YA1,YD1,YI0,YI2,YP11,YCRNT,YPIDPB,STSD3
	do 1	JN	=1,NA1
		YSTE	=SQRT(TE(JN))
	IF(EBEAM.GT.100.*ABEAM)	THEN
		YLNI(JN)=23.7d0+DLOG(AMAIN(JN)/(AMAIN(JN)+ABEAM)*
     *             SQRT(1.d-3*ABEAM*EBEAM*TE(JN)/NE(JN)))
				ELSE
		YLNI(JN)=25.4d0+DLOG(1.d-3*EBEAM*AMAIN(JN)/
     /		  (AMAIN(JN)+ABEAM)*SQRT(TE(JN)/NE(JN)))
				ENDIF
		YLNE(JN)=15.85d0+LOG(TE(JN)/SQRT(NE(JN)))
		YLNZ(JN)=YLNI(JN)
		YS	=0.
	DO	10	JS	=2,ISPEND
		JSP	=ISPE(JS)
 10		YS	=YS+EXTARR(JN,JSP)*ZB(JS)**2*ABEAM/RMB(JS)
		YS	=YLNI(JN)*YS
		YEBDEC(JN)=
     &           EBEAM/(TE(JN)*ABEAM*18.3d0*(0.726d0*YS/ABEAM/
     /				(NE(JN)*YLNE(JN)))**0.6667)
		YZ2D3(JN)=ZEF(JN)*NE(JN)*YLNI(JN)/(YS*3.0d0)
		YTSE(JN)=2.d0*ABEAM*YSTE*TE(JN)/(NE(JN)*YLNE(JN))
		YEPS	=AMETR(JN)/(RTOR+SHIF(JN))
		YFCUR(JN)=(1.d0-FNBF(ZEF(JN),YEPS)/ZEF(JN))
		PBCX(JN)=0.d0
		YN0	= (NNWM+NNCL)*NN(JN)
		YDN	=0.d0
 1	continue
 	do 3	JN	=2,N1
		JN1	=JN-1
		JNA1	=1+CBMI3*(JN1-1)
		JNA	=JNA1+1
		JNAC	=JNA1-1+CBMI3
	if(JN.eq.N1)JNAC=NA1-1
		J2	=1+CBMI3*(JN1-1)-CBMI3*YDN
	IF(CBMI3.GT.1.d0)	YDN	=.5d0
 	do 30	JE	=JEB,IEB
		Y	=1.d0/(IEB-JE+1)
	do 300	J	=JNA1,JNAC
		X2	=YEBDEC(J2)*Y
		X1	=SQRT(X2)
		X3	=X2*X1
		YPIDPB	=2.0d0*FNB2(X2)/X2
		STSD3	=208.33d0*YTSE(J2)*YAQBA(JE,JN1)
C...Tail correction for current
 		YB	=1.d0+1.d0/X3
		YA=0.5*TE(J2)/(EBEAM*Y)*(1.d0+TI(J2)/(TE(J2)*X3))
		YC0	=YTSE(J2)*DTCX(JE,J2)-3.d0
		YC	=YZ2D3(J2)*3.d0/X3
		YC1	=YC0+YC
		YC2	=YC1+2.d0*YC
		Y1	=4.d0*YA*YC1/(YB*YB)
		Y12	=1.d0+SQRT(1.d0+Y1)
		YA1	=0.5d0*YB*Y12/YA
		YD1	=2.d0/(Y12+0.5d0*Y1)
		YI0	=(0.5d0-FNB2(X2)/X2)
		YI2	=FNBP(X1,YZ2D3(J2))
c...Beam energy Wbeam[10#19keV/m3] = Pbper + Pblon/2 
	YP11=0.6666667d0*(YI0*YSLEJ0(JE,JN1)+0.4d0*YI2*YSLEJ2(JE,JN1))
c...Beam perpendicular pressure <Mb Vper2>/2 [10*19 keV/m3]
	PBPER(J)=PBPER(J)+EBEAM*Y*YTSE(J2)*(2.d0*YI0*YSLEJ0(JE,JN1)-YP11)
c...Beam parallel pressure 	  <Mb Vpar2> [10*19 keV/m3]
		PBLON(J)=PBLON(J)+EBEAM*Y*YTSE(J2)*2.d0*YP11
c...Beam density
		NIBM(J)	=NIBM(J)+STSD3*DLOG(1.d0+X3)/(EBEAM*Y)
c...Fast ion current without trapping correction for Eb
		YCRNT	=YACBA(JE,JN1)*YTSE(J2)*YD1*FNBI1(X2,YZ2D3(J2))
C +X3/((1+X3)*YA1**4)*(3.*YA1*(YA1+2.)+7.))
		CUFI(J)	=CUFI(J)+YCRNT
		PIBM(J)	=PIBM(J)+YAQBA(JE,JN1)*YPIDPB
	include 'fml/pbeie'
		PEBM(J)	=PEBM(J)+YAQBA(JE,JN1)*(1.0d0-YPIDPB)-PBEIE
 300		CUBM(J)	=CUBM(J)+YCRNT*YFCUR(J2)
 30	continue
 3	continue
	return
	END
C=======================================================================
	subroutine NBSISN(NA1,YEOA,CBMI3)
C=========================================================== 16-DEC-96
C	Fij = A/L(ri,Ej) - local normalized inverse Neutral Beam
C	mean free path for Ej energy beam component, Ej beam power 
C	absorption along the beam line :
C	Ij(x)	=Ij0*exp(-S Fj dx),	dx=dL/A
C========================================================= Polevoy A.R.
	implicit none
	include  'for/parameter.inc'
	include  'for/status.inc'
	include  'nbi/nbicom.inc'
        double precision	yVB,RNB,ZB,RMB,yEB,DT,EZ,YSIMPI,YEOA
	COMMON /CNBCH/yVB(9),RNB(9),ZB(9),RMB(9),yEB(9),DT,EZ
        integer	ISPEND,ISPE,J,JN,JSP,JE,JS,NA1
	COMMON /CNBSP/ ISPEND,ISPE(9)
	COMMON /CNBSPI/ YSIMPI(3,9)
        double precision	CBMI3,SIMPI,STOTQ,SEIV
        double precision	YY,Y1,Y2,Y3,Y12,Y13,Y23,YJE,SPII,SPEX
C		Y3	=EBEAM/ABEAM
		Y3	=YEOA
		Y1	=Y3/3.d0
		Y2	=Y3/2.d0
	DO 1	JE	=1,IEB
		YJE	=Y3/(IEB-JE+1)
		VNB(JE)	=1.383d6*DSQRT(1.d3*YJE)
		SVII(JE)=SPII(YJE)
		SVEX(4,JE)=SPEX(YJE)
	DO 10	JS	=2,ISPEND
	IF(ZB(JS).GT.1.)	THEN
		YSIMPI(JE,JS)	=SIMPI(YJE,ZB(JS))
	ELSE
		YSIMPI(JE,JS)	=SVEX(4,JE)+SVII(JE)
	ENDIF
 10	CONTINUE
		SVEX(4,JE)=SPEX(YJE)*VNB(JE)
 1	CONTINUE
		SVEX(4,4)=SPEX(1.d0)
		Y12	=Y1+Y2-2.*SQRT(Y1*Y2)
		Y13	=Y1+Y3-2.*SQRT(Y1*Y3)
		Y23	=Y2+Y3-2.*SQRT(Y2*Y3)
		SVEX(1,2)=SPEX(Y12)
		SVEX(1,3)=SPEX(Y13)
		SVEX(2,3)=SPEX(Y23)
	DO 2	JE	=1,IEB
		YY	=Y3/(IEB-JE+1)
	DO 20	JN	=1,N1
		J	=(JN-1)*CBMI3+1
C==Janev,Boley (for Eb>0.1 MeV)
	IF(YY.GE.100.)	THEN
		F(JE,JN)=0.d0
	DO 200	JS	=2,ISPEND
		JSP	=ISPE(JS)
 		F(JE,JN)=F(JE,JN)+
     +	EXTARR(J,JSP)*STOTQ(YY,NE(J),TE(J),ZB(JS),RMB(JS))
 200	CONTINUE
	ELSE
		F(JE,JN)=NE(J)*SEIV(TE(J))/VNB(JE)
	DO 201	JS	=2,ISPEND
		JSP	=ISPE(JS)
 		F(JE,JN)=F(JE,JN)+ EXTARR(J,JSP)*YSIMPI(JE,JS)
 201	CONTINUE
	ENDIF
		F(JE,JN)=F(JE,JN)*(AMETR(NA1-1)+AMETR(NA1))*50.d0
C==Charge-exchange transparancy
C	if(CCD2.EQ.2.)F(JE,JN)=F(JE,JN)*(1.d-7+DTI(JE,J)+DTE(JE,J))/
C     /			(1.d-7+DTI(JE,J)+DTE(JE,J)+DTCX(JE,J))
 20	continue
 2	continue

c	do 3	j=1,n1
c	jn	=2*(j-1)
c	jn1	=jn+1
c 	car10(jn)	=f(3,j)
c3	car10(jn1)	=f(3,j)

	return
	END
C====================================================================
        SUBROUTINE NBSPEC(AIM1,AIM2,AIM3,AMJ,ZMJ,NA1,NB1)
C.....................................RENEWED C*27APR95 (Polevoy)
C Input (arrays only)
C   NHYDR, NDEUT, NTRIT, NHE3, NALF, NIZ1, NIZ2, NIZ3, NI
C   EXTARR(,)
C Output:
C   ISPEND    = 1 + total number of ion species
C	Arrays: 1(e), p, d, t, He3, He4, Imp1, Imp2, Imp3 
C   RMB(9)    m/m_p
C   ZB(9)     
C   ISPE(9)   number of ion species in the EXTARR
C   EXTARR	! One of arrays for (p,d,t,He3) is spoiled !
	implicit none
	include  'for/parameter.inc'
	include  'for/status.inc'
        double precision	AMJ,ZMJ
        double precision	yVB,RNB,ZB,RMB,yEB,DT,EZ,AIM1,AIM2,AIM3
	COMMON /CNBCH/yVB(9),RNB(9),ZB(9),RMB(9),yEB(9),DT,EZ
        integer	ISPEND,ISPE,JIHYDR,JIDEUT,JITRIT,JIHE3,JIALF,JN
	integer	JIZ1,JIZ2,JIZ3,J,JSP,NA1,NB1
	common /CNBSP/ ISPEND,ISPE(9)
C...Identification of plasma species for NBI
C...Polevoy A.R.	18.03.91
	JIHYDR	=0
	JIDEUT	=0
	JITRIT	=0
	JIHE3	=0
	JIALF	=0
	JIZ1	=0
	JIZ2	=0
	JIZ3	=0
	DO	1	JN=1,NA1
	IF(NHYDR(JN).GT.0.d0)	JIHYDR=1
	IF(NDEUT(JN).GT.0.d0)	JIDEUT=1
	IF(NTRIT(JN).GT.0.d0)	JITRIT=1
	IF(NHE3(JN).GT.0.d0)	JIHE3=1
	IF(NALF(JN).GT.0.d0)	JIALF=1
	IF(NIZ1(JN).GT.0.d0)	JIZ1=1
	IF(NIZ2(JN).GT.0.d0)	JIZ2=1
	IF(NIZ3(JN).GT.0.d0)	JIZ3=1
 1	CONTINUE
 	ISPEND	=1
	RMB(1)	=1./1836.d0 
	ZB(1)	=1.
	IF(JIHYDR.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=9
	RMB(ISPEND)	=1.d0 
	ZB(ISPEND)	=1.d0
	ENDIF
	IF(JIDEUT.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=10
	RMB(ISPEND)	=2.d0 
	ZB(ISPEND)	=1.d0
	ENDIF
	IF(JITRIT.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=11
	RMB(ISPEND)	=3.d0 
	ZB(ISPEND)	=1.d0
	ENDIF
	IF(JIALF.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=6
	RMB(ISPEND)	=4.d0 
	ZB(ISPEND)	=2.d0
	ENDIF
	IF(JIHE3.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=12
	RMB(ISPEND)	=3.d0 
	ZB(ISPEND)	=2.d0
	ENDIF
C...Main ion specia identification 
	IF(ISPEND.EQ.1)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=(AMJ+ZMJ+7.00001)
	IF(AMJ.EQ.4.)		ISPE(ISPEND)	=6
C*27APR95 below
	JSP=ISPE(ISPEND)
	DO 99 J=1,NB1
 99	EXTARR(J,JSP)	=NI(J)
C*27APR95 above
	RMB(ISPEND)	=AMJ 
	ZB(ISPEND)	=ZMJ
	ENDIF
C...Impurity identification 
	IF(JIZ1.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=3
	RMB(ISPEND)	=AIM1 
	ZB(ISPEND)	=EXTARR(1,13)
	ENDIF
	IF(JIZ2.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=4
	RMB(ISPEND)	=AIM2 
	ZB(ISPEND)	=EXTARR(1,14)
	ENDIF
	IF(JIZ3.NE.0)	THEN
	ISPEND=ISPEND+1
	ISPE(ISPEND)	=5
	RMB(ISPEND)	=AIM3 
	ZB(ISPEND)	=EXTARR(1,15)
	ENDIF
	return
	end
C==============================================================16-MAR-98
C 	with Ar stopping crossection by Leonov
C======================================================================|
	double precision	FUNCTION	FNBI1(X2,Y)
C------------------------------------------  former NBI
C                                Y X          Y+1
C aproximation for 1/X((1+X3)/X3)  S(U3/(1+U3))    dU
C    0<Y<4                         0
C Mikkelsen D.R.,Singer C.E.//Nucl.Tecnol./Fus.,V.4,Sept.1983,PP.237-252
	implicit none
	double precision	X,X2,Y
		X	=DSQRT(X2)
		FNBI1	=X2*X/(4.d0+3.d0*Y+X2*(X+1.39d0+0.61d0*Y**0.7))
	RETURN
	END	
C======================================================================|
	double precision	FUNCTION	FNB2(X2)
C-----------------------------------------22.07.89
C           X     
C FNB2=      S udu/(1+u3)
C           0
	implicit none
	double precision	X,X2
	X	=DSQRT(X2)
	FNB2	=(0.166666667d0*DLOG((1d0-X+X2)/(1.d0+2.d0*X+X2))
     .	+0.57735026d0*(DATAN(0.57735026d0*(2.d0*X-1.d0))+0.52359874d0))
	RETURN
	END
C======================================================================|
	double precision	FUNCTION	FNBF(Z,E)
C toroidal correction for NB driven current
C Kim Y.B.,Callen J.D.,Hamnen H.//Nucl.Fus.,(1988)
C Neocl.Cur.And Transp.In Aux.Heat.Tokamaks
	implicit none
	double precision	Z,E,FT,G,Z2
		FT	=DSQRT(E)*(1.46d0-0.46d0*E)
		G	=FT/(1.d0-FT)
		Z2	=Z*Z
		FNBF=((Z2+1.41d0*Z)+(Z2+0.45d0*Z)*G)/((Z2+1.41d0*Z)+
     +		(2.d0*Z2+2.66d0*Z+0.75d0)*G+(Z2+1.24*Z+0.35)*G*G)
	RETURN
	END
C======================================================================|
	double precision FUNCTION	FNBP(X,B)
C----------------------------------------------------
C	       3b    X             3b+1       
C FNBP  =(1+1/x3)/x2* S z*(z3/(1+z3)) dz 
C	             0			step 0.05
C----------------------------Polevoy 18.05.89 -------
	implicit none
	double precision	X,B,YA,A,YS,Z,Z1,Z3,Y,YS1,YS2
	integer	JEND,J
	YA	=3.d0*B
	A	=YA+1.d0
	JEND	=20.000001*X
	IF(JEND.LT.7)	THEN
	IF(X.LE.0.0d0)	THEN
	FNBP	=0.0d0
	RETURN
	ENDIF
	FNBP	=X**(3*A+2.)/(3*A+2)*(1.d0+1./(X**3))**YA/(X*X)
	RETURN
	ENDIF
	YS	=0.0d0
	Z	=0.0d0
	IF(JEND.GE.100)	THEN
		DO	1 J=1,100
		Z	=Z+0.050d0
		Y	=Z*(1.d0-1.d0/(1.d0+Z**3))**A
 1		YS	=YS+Y
	FNBP	=((YS-0.5d0*Y)*0.050d0+0.5d0*(X-5.0d0)*(X+5.0d0))
     .		*(1.d0+1.d0/(X**3))**YA/(X*X)
	RETURN
	ELSE
	DO 2	J=1,JEND
		Z	=Z+0.050d0
		Z3	=Z**3
		Y	=Z*(1.d0-1.d0/(1.d0+Z3))**A
 2		YS	=YS+Y
		YS1	=(YS-0.5d0*Y)*(1.+1./Z3)**YA/(Z*Z)
		Z1	=Z
		Z	=Z+0.05d0
		YS2	=(YS+0.5d0*Z*(1.d0-1.d0/(1.+Z**3))**A)
     .		*(1.d0+1./(Z**3))**YA/(Z*Z)
	FNBP	=0.05d0*YS1+(YS2-YS1)*(X-Z1)
	ENDIF
	RETURN
	END
C======================================================================|
	double precision	function	SIMPI(EKEV,Zq)
C----------------------------------------------- Simpi [10#-13 cm2]
C	Neutral beam impurity impact ionization cross section for
C	H(D,T) beam	(3He,4He,C,O,Fe)	impurity
C	Janev R.K.,Boley C.D.,Post D.E.,"Penetration of Neutral
C	Beams into Fusion Plasma",Nucl.Fusion,Vol.29.,No.12,(1989),
C 	pp. 2125-2137.
C	Ebeam [KeV]=0-10000,
C	Simpi	=Zq*C1*[1/(1+C2*E)+C3*ln(1+C5*E)/(C4+E)]
C	E [keV]	=Eb/Mb/Zq
C--------------------------------------------Polevoy A.R. 21.03.91
	implicit none
	double precision	EKEV,Zq,E
	E	=EKEV/ZQ
	SIMPI	=7.457d-3*ZQ*(1./(1.+.08095*E)+
     +			2.754*LOG(1.+1.27*E)/(64.58+E))
	RETURN
	end
C======================================================================|
	double precision	function	SEIV(TE)
C------------------------------------------ <S*Ve> [10#-13 cm3/s]
C	The electron velocity weighted Maxwellian average for
C	neutral beam ionization cross section by electrons
C	TE [keV] -  electron temperature 
C------------------------------------------- Polevoy A.R. 21.03.91
	implicit none
	double precision	TE,Y
		Y	=13.6d-3/TE
	if(TE.LT..02)	THEN
		SEIV	=0.61d6*EXP(-Y)/(SQRT((1.+Y)/Y)*(Y+0.73d0))
			ELSE
		Y	=LOG10(TE)+3.d0
		SEIV	=10.**(7.769d0-0.5151d0*Y-2.563d0/Y)
			ENDIF
	return
	END
C======================================================================|
	double precision	function	SPII(E)
C---------------------------------------------- Spii [10#-13 cm2]
C	Neutral beam ionization cross section by proton impact
C	Reviere A.C //Nucl.Fusion.v.11(1971).p.363
C	E [kev] = Ebeam*Mp/Mi
C------------------------------------------- Polevoy A.R. 21.03.91
	implicit none
	double precision	E,Y1,Y2
		Y1	=LOG10(E)+3.d0
		Y2	=0.d0
	IF(E.GT.3.d0.AND.E.LT.150.d0) 
     &    Y2=10**((-0.8712d0*Y1+8.156d0)*Y1-21.833d0)
	IF(E.GE.150.d0)	Y2=36.d-3*(LOG10(0.1666d0*E)+3.d0)/E
			SPII	=Y2
	return
	END
C======================================================================|
	double precision	function	SPEX(E)
C---------------------------------------------- Spex [10#-13 cm2]
C	Neutral beam ionization cross section by proton impact
C	Reviere A.C //Nucl.Fusion.v.11(1971).p.363
C	E [kev] = Ebeam*Mp/Mi
C------------------------------------------- Polevoy A.R. 21.03.91
	implicit none
	double precision	E,Y
	Y	=1.d3*E
	SPEX	=0.06937d0*(1.d0-0.155d0*LOG10(Y))**2/
     &   (1.d0+0.1112d-14*Y**3.3)
	return
	END
C======================================================================|
	double precision	FUNCTION	YEXP(X)
C----------------------------------------	Artificial exponent
	implicit none
	double precision	X
	YEXP	=0.d0
	IF(X.GT.-30.d0)	YEXP=EXP(X)
	RETURN
	END
C======================================================================|
	double precision function	STOTQ(EKEV,NE19,TEKEV,Zq,Aq)
C----------------------------------------------- Stotq [10#-13 cm2]
C	Partial neutral beam stopping cross section for
C	H(D,T) beam	H(D,T,3He,4He,C,O,Fe)	plasma
C	Janev R.K.,Boley C.D.,Post D.E.,"Penetration of Neutral
C	Beams into Fusion Plasma",Nucl.Fusion,Vol.29.,No.12,(1989),
C 	pp. 2125-2137.
C	Ne [10#19 m-3]=0.1-100, Ebeam [KeV]=100-10000, Te [KeV]=1-50
C	E [keV]=Ebeam*Mi/Mb
C	Stotq	=Zq*(1+Sq*(Zq-1))*exp[S1(E,ne,Te)]/E,
C	n*Stot	=nq*Stotq+... , q=1,2,...,N for all ion species
C--------------------------------------------- Polevoy A.R. 21.03.91
C---------------------------Ne,Ar added by Leonov V.M. in 1999
	implicit none
	double precision EKEV,TEKEV,Zq,Aq,SQ,S1,T,ALT,ALN,ALE,AT,AN,AE
	integer	I,J,K
	double precision	NE19,A(2,3,2),B(3,2,2)
	IF(EKEV.LT.100.)
     .	WRITE(*,*) 'Illegal use of STOT: Eb<100 keV'
 	A(1,1,1)=4.4d0
	A(2,1,1)=2.3d-1
	A(1,1,2)=-2.49d-2
	A(2,1,2)=-1.15d-2
	A(1,2,1)=7.46d-2
	A(1,2,2)=2.27d-3
	A(1,3,1)=3.16d-3
	A(1,3,2)=-2.78d-5
	A(2,2,1)=-2.55d-3
	A(2,2,2)=-6.2d-4
	A(2,3,1)=1.32d-3
	A(2,3,2)=3.38d-5
	IF(ZQ.GT.1.) GO TO 91
	DO	100	I=1,3
	DO	100	J=1,2
	DO	100	K=1,2
 100	B(I,J,K)	=0.d0
	GOTO 99
 91	IF(AQ.LT.11.d0) GOTO 4
	IF(AQ.GT.11.d0.AND.AQ.LT.13.d0) GOTO 12
	IF(AQ.GT.15.d0.AND.AQ.LT.17.d0) GOTO 16
	IF(AQ.GT.19.d0.AND.AQ.LT.21.d0) GOTO 20
	IF(AQ.GT.39.d0.AND.AQ.LT.41.d0) GOTO 40
	IF(AQ.GT.51.d0.AND.AQ.LT.54.d0) GOTO 52
	WRITE(*,*)	'NBI: No stop cross section DATA for Aimp=',AQ
	WRITE(*,*)	'    cross section for Aimp=52 is substituted insted'
	GOTO	52
 4	B(1,1,1)=-2.36d0
	B(1,1,2)=1.85d-1
	B(1,2,1)=-2.5d-1
	B(1,2,2)=-3.81d-2
	B(2,1,1)=8.49d-1
	B(2,1,2)=-4.78d-2
	B(2,2,1)=6.77d-2
	B(2,2,2)=1.05d-2
	B(3,1,1)=-5.88d-2
	B(3,1,2)=4.34d-3
	B(3,2,1)=-4.48d-3
	B(3,2,2)=-6.76d-4
	GOTO 99
 12	B(1,1,1)=-1.49d0
	B(1,1,2)=-1.54d-2
	B(1,2,1)=-1.19d-1
	B(1,2,2)=-1.50d-2
	B(2,1,1)=5.18d-1
	B(2,1,2)=7.18d-3
	B(2,2,1)=2.92d-2
	B(2,2,2)=3.66d-3
	B(3,1,1)=-3.36d-2
	B(3,1,2)=3.41d-4
	B(3,2,1)=-1.79d-3
	B(3,2,2)=-2.04d-4
	GOTO 99
 16	B(1,1,1)=-1.41d0
	B(1,1,2)=-4.08d-4
	B(1,2,1)=-1.08d-1
	B(1,2,2)=-1.38d-2
	B(2,1,1)=4.77d-1
	B(2,1,2)=1.57d-3
	B(2,2,1)=2.59d-2
	B(2,2,2)=3.33d-3
	B(3,1,1)=-3.05d-2
	B(3,1,2)=7.35d-4
	B(3,2,1)=-1.57d-3
	B(3,2,2)=-1.86d-4
	GOTO 99
 20	B(1,1,1)=-1.32d0
	B(1,1,2)=1.8d-2
	B(1,2,1)=-9.7d-2
	B(1,2,2)=-1.13d-2
	B(2,1,1)=4.45d-1
	B(2,1,2)=-4.5d-3
	B(2,2,1)=2.3d-2
	B(2,2,2)=2.85d-3
	B(3,1,1)=-2.78d-2
	B(3,1,2)=2.05d-3
	B(3,2,1)=-1.4d-3
	B(3,2,2)=-1.7d-4
	GOTO 99
 40	B(1,1,1)=-1.12d0
	B(1,1,2)=7.3d-2
	B(1,2,1)=-7.0d-3
	B(1,2,2)=-6.7d-3
	B(2,1,1)=3.62d-1
	B(2,1,2)=-2.4d-2
	B(2,2,1)=1.6d-2
	B(2,2,2)=1.68d-3
	B(3,1,1)=-2.18d-2
	B(3,1,2)=2.92d-3
	B(3,2,1)=-9.5d-4
	B(3,2,2)=-9.5d-5
	GOTO 99
 52	B(1,1,1)=-1.03d0
	B(1,1,2)=1.06d-1
	B(1,2,1)=-5.58d-2
	B(1,2,2)=-3.72d-3
	B(2,1,1)=3.22d-1
	B(2,1,2)=-3.75d-2
	B(2,2,1)=1.24d-2
	B(2,2,2)=8.61d-4
	B(3,1,1)=-1.87d-2
	B(3,1,2)=3.53d-3
	B(3,2,1)=-7.43d-4
	B(3,2,2)=-5.12d-5
	GOTO 99
 99	CONTINUE
	SQ=0.
	S1=0.0
	T=MAX(1.d0,TEKEV)
	ALT=LOG(T)
	ALN=LOG(NE19)
	ALE=LOG(EKEV)
		AT	=1.0d0
	DO	1	K=1,2
		AN	=AT
	DO	2	J=1,2
		AE	=AN
	DO	3	I=1,2
			S1	=A(I,J,K)*AE+S1
			SQ	=B(I,J,K)*AE+SQ
 3		AE	=AE*ALE
			SQ	=B(3,J,K)*AE+SQ
 2		AN	=AN*ALN
 1		AT	=AT*ALT
		AT	=ALN*ALN
	DO	10	K=1,2
		AE	=AT
	DO	30	I=1,2
			S1	=A(I,3,K)*AE+S1
 30		AE	=AE*ALE
 10		AT	=AT*ALT
	STOTQ	=1.d-3*EXP(S1)*ZQ*(1.+(ZQ-1.0d0)*SQ)/EKEV
	RETURN
	end
C======================================================================|
