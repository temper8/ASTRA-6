 		YQBP	=YAQBP(JE)*YDEX
		YF0	=ARD(JE,JN)-AZ(JR)
		YFA	=YF0-ARD(JE,N1)
	if(YFM0R.gt.YF0) 	then
			YFM0R=YF0
			jfig=0
				else
			jfig=1
				endif
		YDDD	=(YAQBP(JE)-YQBP)
	if(ARD(JE,JI).LT.YF0.AND.ARD(JE,N1).GT.YF0)	THEN
C...Trapping
		JTRAP	=0
ccc   	  if(YFA.LE.-YR2.AND.YFA.GE.YR2)	THEN
   	  if(YFA.LE.0..AND.YFA.GE.YR2)	THEN
c...Lost
   		YD	=0.
   		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
   		JNL	=JN1
   		JNR	=JN1
   		YJSURF	=1.
   	  ELSE
c...Kept
	call nbfoa(JTRAP,JLOSS)
	IF(JLOSS.EQ.0)	THEN
C...Kept
	if(jnr.lt.n1)jnr=jnr+1
	if(jnl.lt.n1)jnl=jnl+1
		YJSERF	=1./(JNR-JNL+1)
ctest
CCC	WRITE(*,*) 'TRAP',YJSERF,JNR,JNL
		YD	=YDDD*YJSERF
		YDCOS	=YD*YCT2(JE)
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
	    		ELSE
c...Lost
C	WRITE(1,*) 'trap_lost',JN,JR,JH,ixyu
		YD	=0.
		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
		JNL	=JN1
		JNR	=JN1
		YJSURF	=1.
			ENDIF
   	  ENDIF
	ELSE
c...Passing
		JTRAP	=1
cc   	  if(YFA.LE.YR1.AND.
		if(YFA.GE.YR2)	THEN
c...Lost
   		YD	=0.
   		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
   		JNL	=JN1
   		JNR	=JN1
   		YJSURF	=1.
             ELSE
c...Kept
	call nbfoa(JTRAP,JLOSS)
	IF(JLOSS.EQ.0)	THEN
C...Kept
	if(jnr.lt.n1)jnr=jnr+1
	if(jnl.lt.n1)jnl=jnl+1
		YJSERF	=1./(JNR-JNL+1)
ctest	WRITE(*,*) 'PASS',YJSERF,JNR,JNL
		YD	=YDDD*YJSERF
		YDCOS	=YD*(YC2+YCT2(JE))
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
	    		ELSE
c...Lost
C	WRITE(1,*) 'pass_lost',JN,JR,JH,ixyu
		YD	=0.
		YDCOS	=0.
		YTCOS	=YC2+YCT2(JE)
		YCT2(JE)=YC2
		JNL	=JN1
		JNR	=JN1
		YJSURF	=1.
			ENDIF
	  ENDIF
   	ENDIF
     		YANBA(JE,JN1)=YANBA(JE,JN1)+YDDD
     		YATBA(JE,JN1)=YATBA(JE,JN1)+YDDD*YTCOS
ctest		JNL=JN1		JNR=JN1
	do	J_JN=JNL,JNR
C    		YANBA(JE,JN1)=YANBA(JE,JN1)+YCDV*YAQBP(JE)
     		YAQBA(JE,J_JN)=YAQBA(JE,J_JN)+YD
     		YACBA(JE,J_JN)=YACBA(JE,J_JN)+YDCOS
     		YAQBP(JE)=YQBP
		YDYDT	=YD*YDT
	do  	JT	=JT1,JT2
		YASBA(JE,J_JN,JT)	=YASBA(JE,J_JN,JT)+YDYDT
	enddo
	enddo
	enddo
