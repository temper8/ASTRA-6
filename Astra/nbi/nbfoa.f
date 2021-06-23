	subroutine nbfoa(ITRAP,ILOSS)
c....	First orbit analysis taking account the finit 
c	gyroradius and the input source redistribution
c	due to the first orbit deviation
C input:
c	JIM0R,RM0R,YFM0R =	starting values of 
c	II index, Radius and F value in the midplane
c	RCR,RJ,RJT 	= critical, external,internal boundary radii 
C	YF0 		= F value at the birth position
c	ITRAP 	= 0/1 for banana/passing orbits
c output:
c	JNL,JNR	= minimum/maximum surface index of the orbit
C	ILOSS		= 0/1 if particle is kept/lost
	implicit none
	include	'for/parameter.inc'
	include 'nbi/nbicom.inc'
ccc	include 'nbi/nbicom.inc'
	include 'nbi/nbfoa.inc'
c	COMMON /CNBFOA/YPSI(3,201),RLM(3,201),RJM(201),II(201),
c     ,	RCR,RJ,RJT,JIM0R,JIM0L,RM0R,RM0L,YFM0R,YFM0L,YFM0T,YF0,
c     ,	YDYH,JNR,JNL,JE,jfig,ixyu,SQRCR,JCENTR
	double precision	Y,YY,RM1R,RM2R,YFM1R,YFM2R
	integer	ITRAP,ILOSS,JIM1R,JIM2R,JINR,JINL,JINL1,JNL1,N
	JIM2R	=JIM0R
	RM2R	=RM0R
	YFM2R	=YFM0R

 1	JIM1R	=JIM2R
	YFM1R	=YFM2R
	RM1R	=RM2R
C...motion to the external boundary RM2--->RJT
	JIM2R	=JIM2R-1
	IF(JIM2R.LT.1) 	THEN
	ILOSS=1
c			ixyu=11
		RETURN
				ENDIF
	RM2R	=RJM(JIM2R)
	IF(RM2R.LT.RCR) WRITE(*,*) 'PAS->TRAP_1',RM2R,RCR
     	YFM2R	=YPSI(JE,JIM2R) -SQRT(RM2R*(RM2R-RCR))
	Y	=(YF0-YFM2R)*(YF0-YFM1R)
	if(Y.GT.0.) GOTO 1
 	JINR	=JIM2R-SQRCR*RLM(JE,JIM2R)
	IF(JINR.LT.1) 	THEN
C...lost due to finit gyroradius
	ILOSS	=1 
c			ixyu=12

				RETURN		
				ELSE
	JNR	=II(JINR)
	ILOSS =0	
				ENDIF

	if(ITRAP.EQ.0) GOTO 33
	if(JFIG.eq.0) goto 44

 	JIM2R	=JIM0L
	YFM2R	=YFM0L
	RM2R	=RM0L

 2	JIM1R	=JIM2R
	YFM1R	=YFM2R
	RM1R	=RM2R
C...motion to the external boundary RM2--->RJ
	JIM2R	=JIM2R-1
	IF(JIM2R.LT.1) 	THEN
	ILOSS=1
c			ixyu=21
				RETURN
				ENDIF
	RM2R	=RJM(JIM2R)
	IF(RM2R.LT.RCR) WRITE(*,*) 'PAS->TRAP_1',RM2R,RCR
     	YFM2R	=YPSI(JE,JIM2R) -SQRT(RM2R*(RM2R-RCR))
	Y	=(YF0-YFM2R)*(YF0-YFM1R)
	if(Y.GT.0.) GOTO 2
	YY	=SQRCR*RLM(JE,JIM2R)
	JINL	=JIM2R-YY
	JINL1	=JIM2R+YY
	IF(JINL.LT.1.OR.JINL1.GE.201) 	THEN
C...lost due to finit gyroradius
	ILOSS	=1 		
c			ixyu=22
	RETURN
				ELSE
	JNL	=II(JINL)
	JNL1	=II(JINL1)
	IF(((JINL1-JCENTR)*(JINL-JCENTR)).LE.0)	JNL=1
	IF(JNL.GT.JNL1) JNL=JNL1
	ILOSS =0	
			ENDIF
				return
 44	JIM2R	=JIM0R
	RM2R	=RM0R
	YFM2R	=YFM0R

 4	JIM1R	=JIM2R
	YFM1R	=YFM2R
	RM1R	=RM2R
C...motion to the internal boundary RM2--->RJT
	JIM2R	=JIM2R+1
	IF(JIM2R.gT.201) 	THEN
	ILOSS=1
c			ixyu=41
	RETURN
				ENDIF
	RM2R	=RJM(JIM2R)
	IF(RM2R.LT.RCR) WRITE(*,*) 'PAS->TRAP_1',RM2R,RCR
     	YFM2R	=YPSI(JE,JIM2R) -SQRT(RM2R*(RM2R-RCR))
	Y	=(YF0-YFM2R)*(YF0-YFM1R)
	if(Y.GT.0.) GOTO 4
	YY	=SQRCR*RLM(JE,JIM2R)
	JINL	=JIM2R-YY
	JINL1	=JIM2R+YY
	IF(JINL.LT.1.OR.JINL1.GE.201)	THEN
C...lost due to finit gyroradius
	ILOSS	=1 
c			ixyu=42
				RETURN		
				ELSE
	JNL	=II(JINL)
	JNL1	=II(JINL1)
	IF(((JINL1-JCENTR)*(JINL-JCENTR)).LE.0)	JNL=1
	IF(JNL.GT.JNL1) JNL=JNL1
	ILOSS =0	
				ENDIF
				return
c...For trapped particles
 33	JIM2R	=JIM0R
ccc	YFM2R	=YFM0R
	RM2R	=RM0R
CC 33
	    	YFM2R	=YFM0T
 3	JIM1R	=JIM2R
	YFM1R	=YFM2R
	RM1R	=RM2R
C...motion to the Internal boundary RM2--->RJT
	JIM2R	=JIM2R+1
	IF(JIM2R.GT.201) 	THEN
	ILOSS=1
c			ixyu=1
				RETURN
				ENDIF
	RM2R	=RJM(JIM2R)
	IF(RM2R.LT.RCR) 	THEN
		ILOSS=1
c			ixyu=2
				RETURN
				ENDIF
     	YFM2R	=YPSI(JE,JIM2R) +SQRT(RM2R*(RM2R-RCR))
	Y	=(YF0-YFM2R)*(YF0-YFM1R)
	if(Y.GT.0.) GOTO 3
	YY	=SQRCR*RLM(JE,JIM2R)
	JINL	=JIM2R-YY
	JINL1	=JIM2R+YY
	IF(JINL.LT.1.OR.JINL1.GE.201)	THEN
C...lost due to finit gyroradius
	ILOSS	=1 		
c			ixyu=3
							return
							ELSE
	JNL	=II(JINL)
	JNL1	=II(JINL1)
	IF(((JINL1-JCENTR)*(JINL-JCENTR)).LE.0)	JNL=1
	IF(JNL.GT.JNL1) JNL=JNL1
	ILOSS	=0

							ENDIF
CTEST
 99	N=N1-1
	IF(JNL.LT.1) THEN
	WRITE(*,*) 'JNL<1' , JNL
	JNL=1
			ENDIF
	IF(JNR.LT.1) THEN
	WRITE(*,*) 'JNR<1' , JNR
	JNR=1
			ENDIF
	IF(JNR.GT.N1) THEN
	WRITE(*,*) 'JNR>N' , JNR,jiM2r,ITRAP
	JNR=N1
			ENDIF
	IF(JNL.GT.N1) THEN
	WRITE(*,*) 'JNL>N' , JNL,jiM2r,ITRAP
	JNL=N1
			ENDIF
	IF(JNL.GT.JNR) THEN
	WRITE(*,*) 'JNR<JNL' , JNR,JNL,ITRAP
	JNR=JNL
			ENDIF
	
	RETURN
	END
