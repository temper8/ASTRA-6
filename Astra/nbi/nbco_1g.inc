		YDT	=1./(JT2-JT1+1)
		YR2	=-SQRT(RJ*(RJ-RCR))
	if(RCR.LT.RJT)	THEN
		YR1	=-SQRT(RJT*(RJT-RCR))
		JI	=N1
	ELSE
		YR1	=0.
		JII	=(RJ-RCR)*YDYH
		JI	=II(JII)
	ENDIF
		YDY	=Y2-Y1
		YCDV	=DV*YDY
		RM0R	=Y0+DX(JN)+X(JN)
		JIM0R	=(RJ-RM0R)*YDYH+2
	IF(JIM0R.LE.0) JIM0R=1
		RM0L	=Y0+DX(JN)-X(JN)
		JIM0L	=(RJ-RM0L)*YDYH+2
	IF(JIM0L.GT.JII) THEN
		JIM0L	=JII
		RM0L	=RCR
				ENDIF
	IF(JIM0R.LE.0) JIM0R=1
CC		RM0R	=RJM(JIM0R)
CC		RM0L	=RJM(JIM0L)
		YSQR	=SQRT(RM0R*(RM0R-RCR))
		YSQL	=SQRT(RM0L*(RM0L-RCR))
		SQRCR	=SQRT(RCR)
        do 	JE	=JBB,IEB
		YFM0R	=ARD(JE,JN)-YSQR
		YFM0T	=ARD(JE,JN)+YSQR
CC		YFM0R	=YPSI(JE,JIM0R)-YSQR
		YFM0L	=ARD(JE,JN)-YSQL
CC		YFM0L	=YPSI(JE,JIM0L)-YSQL
		YDYS	=F(JE,JN)*YDY
     		YS(JE)	=YS(JE)+YDYS