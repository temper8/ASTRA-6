C====================================================================
        subroutine NBTRAG(YAQBP,CONTR,AR,JSRNUM,CBMS1)
C====== Neutral beam ionization ================ Version by 08-JUL-97+ 24-MAY-06
C coinjection  with ion's trapping and orbital losses
C  AQBA[MW],ACBA[MA m/s],ANBA[10*13 prtcls],
C  ASBA[10#19 prtcl/s]*Dcos(JT)
C  YQSHth [MW] shine through power in the file dat\shth.dat
C============================================================ Polevoy
	implicit none
	double precision YAQBP(*),AR(*),DRL(3),YS(3),YCU(3),YCT2(3)
     2 ,     YVE(3),PLEJ2(51),CONTR,CBMS1
	include	 'for/parameter.inc' 
	include 'nbi/nbicom.inc' 
c	COMMON /CNBFOA/YPSI(3,201),RLM(3,201),RJM(201),II(201),
c     ,	RCR,RJ,RJT,JIM0R,JIM0L,RM0R,RM0L,YFM0R,YFM0L,YFM0T,YF0,
c     ,	YDYH,JNR,JNL,JE,jfig,ixyu,SQRCR,JCENTR
	include 'nbi/nbfoa.inc'
	integer JSRNUM,JLOSS
c=====Ripple	normalized radius of the ripple boundary
c==== banana with RTRAP > YRIPLR is lost 
        double precision	YRIPLR
	common	/CRIPL/	YRIPLR(NRD) 
        double precision	YQSHTH,Y0,DS,DV,YLOSS,YE3,RJ0,YH,YR,YRN
        double precision	YE,YVEDE,Y,YDH,XJH,ZJH,YZJH,Y12,Y2,YC2
        double precision	YRN1,YR2,YR1,YQBP,YDEDJ
        double precision	Y1,YC1,YDT,YDY,YCDV,YSQR,YSQL,YDYS,YDEX
        double precision	YFA,YDDD,YD,YDCOS,YJSURF,YJSERF,YDYDT,YTCOS
C--------|---------|---------|---------|---------|---------|---------|
        integer	ntet,n,jn,jt,jhb05,jn1,ji,in,in1,jh,jh1,jjh,j_jn
        integer	jjr,jr,jbb,jjn,jj,jt1,jt2,jii,jnn,jhp1,jeth,jtrap
	character*12 SHNAME
	DATA SHNAME/'dat/shth.dat'/
	DATA PLEJ2/
     .-5.000000E-01,-4.994000E-01,-4.976000E-01,-4.946000E-01,
     .-4.904000E-01,-4.850000E-01,-4.784000E-01,-4.706000E-01,
     .-4.616000E-01,-4.514000E-01,-4.400000E-01,-4.274000E-01,
     .-4.136000E-01,-3.986000E-01,-3.824000E-01,-3.650000E-01,
     .-3.464000E-01,-3.266000E-01,-3.056000E-01,-2.834000E-01,
     .-2.600000E-01,-2.354000E-01,-2.096000E-01,-1.826000E-01,
     .-1.544000E-01,-1.250000E-01,-9.440002E-02,-6.260002E-02,
     .-2.960002E-02, 4.599977E-03, 3.999998E-02, 7.659998E-02,
     . 1.144000E-01, 1.534000E-01, 1.936000E-01, 2.350000E-01,
     . 2.776000E-01, 3.214000E-01, 3.664000E-01, 4.126000E-01,
     . 4.599999E-01, 5.085999E-01, 5.584000E-01, 6.094000E-01,
     . 6.615999E-01, 7.150000E-01, 7.695999E-01, 8.253999E-01,
     . 8.823999E-01, 9.405999E-01, 9.999999E-01/

c=====Ripple

c*TRACE	write(*,*) 'BOWLA B nbtrag'
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
			YQSHTH=0.
cc			write(1,*) 'R[m]      Z[m]      Q/S[W/cmZ]'
    		NTET	=NTET1-1
		N	=N1-1
		Y0	=R/A
		DS	=ABS(RBMAX1-RBMIN1)*A/(N*IRB)
		DV	=DS*A
		YLOSS	=0.
Clean the sources
	do 1	JN	=1,N1
	do 10	JE	=JEB,IEB
		YANBA(JE,JN)=0.
		YACBA(JE,JN)=0.
		YATBA(JE,JN)=0.
		YAQBA(JE,JN)=0.
	do 10	JT	=1,NTET1
		YASBA(JE,JN,JT)=0.
 10	continue
		RC(JN)	=Y0+DX(JN)
 1	continue
	if(IHB.eq.0) return
c NBI out of plasma
		JHB05	=IHB/2
		YE3	=EB/PB
cc	include 'nbtrap.inc'
c for calculation of surface index II=JN(X(Rcrit)) in trapping analysis
c YR=(RJ-Rii/AB) normalised distance from ext. boundary RJ
c Trapping analysis with gyroradius and surface averaging
c Rlm = Gyrorad(in the midplain)*YDYH/sqrt(RCR)
 		RJ	=Y0+1.
		RJT	=Y0-1.
		RJ0	=Y0+DX(1)
		YDYH	=100.0
		YH	=0.01
		JN1	=N1
	do 80	JE	=JEB,IEB
		DRL(JE)=1./(0.0144*SQRT(EB*PB/(IEB-JE+1))/(BZ*A))
	do 80	JN	=1,N1
		ARD(JE,JN)=AR(JN)*DRL(JE)
 80	continue

	DO 8	JI	=1,200
		YR	=YH*(JI-1)
 81	IF(JN1.GT.1)					THEN
		JN	=JN1-1
		YRN1	=Y0+DX(JN1)+X(JN1)
		YRN	=Y0+DX(JN)+X(JN)
c*NEW*28MAY01vvvvvvvvvvvvvvvvvvvvvvvvvvvv
	 IF(YR.GE.(X(N1)+DX(N1)-DX(JN1)-X(JN1))
     .		.AND.YR.LT.(X(N1)+DX(N1)-DX(JN)-X(JN))) THEN

c*NEW*28MAY01^^^^^^^^^^^^^^^^^^^^^^^^^^
		II(JI)	=JN
						ELSE
	 	JN1	=JN

		go to	81
	 					ENDIF
C...For orbit averaging and gyrolosses \\\\\\\\\\\\\\\\\\\\\\\\\
		RJM(JI)	=RJ-YR
	do JE=JEB,IEB
		YPSI(JE,JI)=ARD(JE,JN)+(ARD(JE,JN1)-ARD(JE,JN))*
     *				(RJM(JI)-YRN)/(YRN1-YRN)
CCC		RLM(JE,JI)=0.
		RLM(JE,JI)=YDYH*sqrt(RJM(JI))/Y0/DRL(JE)
	if(CBMS1.ge.2.) RLM(JE,JI)=0.	! Finite Larmor radius off
	enddo
C...For orbit averaging and gyrolosses /////////////////////////
							ELSE
	 IF(JN1.EQ.1)	THEN
		JN1	=-1
		JCENTR=JI
				ENDIF
 82		IN1	=-JN1
		IN	=IN1+1

		YRN1	=Y0+DX(IN1)-X(IN1)
		YRN	=Y0+DX(IN)-X(IN)
c*NEW28MAY01	 IF(YR.GE.(RJ-YRN1).AND.YR.LT.(RJ-YRN)) THEN
	 IF(YR.GE.(X(N1)+DX(N1)+X(IN1)-DX(IN1))
     .		.AND.YR.LT.(X(N1)+DX(N1)+X(IN)-DX(IN))) THEN

c*NEW28MAY01^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		II(JI)	=IN1
	 					ELSE
	 	JN1	=JN1-1
		go to 82
	 					ENDIF
C...For orbit averaging and gyrolosses \\\\\\\\\\\\\\\\\\\\\\\\
		RJM(JI)	=RJ-YR
	do JE=JEB,IEB
		YPSI(JE,JI)=ARD(JE,IN)+(ARD(JE,IN1)-ARD(JE,IN))*
     *				(RJM(JI)-YRN)/(YRN1-YRN)
ccc		RLM(JE,JI)=0.
		RLM(JE,JI)=YDYH*sqrt(RJM(JI))/Y0/DRL(JE)
	if(CBMS1.ge.2.) RLM(JE,JI)=0.
	enddo
C...For orbit averaging and gyrolosses /////////////////////////
							ENDIF
 8	CONTINUE
		II(201) =N
C...For orbit averaging and gyrolosses \\\\\\\\\\\\\\\\\\\\\\\\
		RJM(201)=RJT
		RJM(1)=RJ
	do JE=JEB,IEB
		YPSI(JE,201)=ARD(JE,N1)
		YPSI(JE,1)=ARD(JE,N1)
CCC		RLM(JE,201)=0.
		RLM(JE,201)=YDYH*sqrt(RJM(201))/Y0/DRL(JE)
	if(CBMS1.ge.2.) RLM(JE,201)=0.
CCC		RLM(JE,1)=0.
		RLM(JE,1)=YDYH*sqrt(RJM(1))/Y0/DRL(JE)
	if(CBMS1.ge.2.) RLM(JE,1)=0.
	enddo
cc	open(39,file='dat\test.dat')
cc	write(39,*) ii
CC	write(39,*) (ypsi(3,jjjj),jjjj=1,201)
cc	close(39)

cc	include 'nbtrag.inc'

	do 9	JE	=JEB,IEB
		YE	=YE3/(IEB-JE+1)
  		YVEDE	=DS*SQRT(YE)/451.9
		YVE(JE)	=YVEDE*PB*YE*1.E-3	
 		YCU(JE)	=-CONTR*VNB(JE)*YVEDE*5.E-6
 9	CONTINUE

CHH 2 - loop beam`s height H = X(JN)+dX/2  HHHHHHHHHHHHHHHHH
	do 2	JJH	=1,JHB05
		JH	=JHB05-JJH+1
		JH1	=JH+1
		Y	=ELON1(JH1)*X(JH1)-ELON1(JH)*X(JH)
	if(JH.EQ.JHB05)	THEN
		YDH	=N*YDZ(JH)
		XJH	=X(JH)+0.5*YDH/(N*ELON1(JH))
		RE(JH)	=RC(JH)-(RC(JH1)-RC(JH))*(XJH-X(JH))/
     /			 (X(JH1)-X(JH))-TRIA1(JH)
			ELSE
		XJH	=X(JH)+0.5*Y/ELON1(JH)
		YDH	=N*YDZ(JH)
		RE(JH)	=0.5*(RC(JH1)+RC(JH)-TRIA1(JH)-TRIA1(JH1)) 
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
	if(jn.gt.1.and.RI(JN-1).lt.RI(JN))then
	 write(*,*) JN
	 write(*,*) JH
	 write(*,*) JHB05
	 write(*,*) N1
	endif
	enddo 
c 21	continue
	do 20	JN	=JH,N1
		DRE(JN)	=1./RE(JN)
		DRI(JN)	=1./RI(JN)
 20	continue
CZZ 22 - loop beam`s radius R = Z(JR)    ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
	do 22	JR	=1,IRB
	if(YAQB(3).le.0.)	goto 23
C---- no power in this pencil ------------------ go to 23 -------------
		JBB	=JEB
	do 221	JE	=JBB,IEB
		YS(JE)	=0.
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
		YR1	=0.
	ENDIF
	do 2200	JE	=JEB,IEB
		YCT2(JE)=0.
 2200	continue
	IF(CONTR.LT.0.)	THEN
c	include 'nbtcon.inc'
c*GG	include 'nbi/nbco_0.rip'
c
	include 'nbi/nbco_gg.inc'
ccc	include 'nbi/nbco_gg.inc'
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
cc			write(1,*) AZ(JR)*A*1.E-2,ZJH*A*1.E-2,YYYD
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
		YDEDJ	=1.E6/(1.6*EB)
	DO 4	JE=JEB,IEB
		Y	=(IEB-JE+1)*YDEDJ*YVE(JE)
	DO 4	JN=1,N1
		YSLEJ0(JE,JN)=0.0
		YSLEJ2(JE,JN)=0.0
	DO 4	JT=1,NTET1
	IF(YASBA(JE,JN,JT).GT.0.)	THEN
 		YASBA(JE,JN,JT)	=YASBA(JE,JN,JT)*Y
     	YSLEJ2(JE,JN)=YSLEJ2(JE,JN)+YASBA(JE,JN,JT)*2.5*PLEJ2(JT)
     	YSLEJ0(JE,JN)=YSLEJ0(JE,JN)+0.5*YASBA(JE,JN,JT)
					ENDIF
 4	CONTINUE
cc>>>>>Shinethrough
			
	write(38,*) YQSHTH,' Q shine through [MW] for source ',JSRNUM
			CLOSE(38)
cc>>>>>Shinethrough end
	return
	END
