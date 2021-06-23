C=======================================================================|
	SUBROUTINE	FOUREX(ARRIN,FREQ,NUMBER,ARROUT)
C-----------------------------------------------------------------------|
C Description:	    Fourier expansion in time for the function ARRIN(a,t)
C							of 2 arguments
C Input: ARRIN(*) - periodic radial dependent function of time to be expanded
C	 FREQ	  - fundamental harmonic frequency [Hz]
C	 NUMBER	  - number of harmonics in output is  0 <= 2*NUMBER+1 <= 21
C
C Output:ARROUT(ro,#) - a0(ro), a1(ro), fi1(ro), a2(ro), fi2(ro), ...
C	 		amplitudes and phases of Fourier harmonics
C			as functions of radius and #
C		if (NUMBER = 0) the output is time average only
C
C Example call from a model:	FOUREX(TE,CF2,5,CAR6):;
C	The Fourier harmonics ao(r), a1(r), fi1(r),..., fi5(r) of TE(r,t)
C	will be returned in CAR6, CAR7, CAR8,..., CAR16 respectively.
C
C Important: The subroutine can be called from a model only once, otherwise  
C				the second call spoils the previous. 
C	     To avoid the interference the same subroutine can be called
C				under another name (see sbr/four1.for)
C-----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	double precision ARRIN(*),ARROUT(NRD,*),YA(NRD,21)
	double precision FREQ,HARM,YSTART,YDT,YARG,YJARG,YYAA,YAMPL
	double precision YPHASE,YMAX1,YY,YGP,YTOLD
	integer	JNUM,j,NUMBER,JN,JCALL,LCALL,JTIME
	save	JCALL,YA,YTOLD
	data	JCALL/0/

	YARG	= GP2*FREQ*TIME
	if(JCALL.eq.0)	goto	10
	YDT	= TIME-YTOLD
	JTIME	= TIME*FREQ
C YSTART: beginning of a time period and lower limit of integral.
C NB!	- it is assumed that power oscillations have NO PHASE SHIFT
C	  			with respect to the moment TIME = 0
 	YSTART	= JTIME/FREQ
	if(TIME .ge. YSTART .and. YTOLD .lt. YSTART)	goto	10

	JNUM	= 0
	YTOLD	= TIME
	do	1	j = 1,NA1
		YA(j,1)	= YA(j,1)+ARRIN(j)*YDT
 1	continue
	if(NUMBER .eq. 0)	return

 2	JNUM	= JNUM+1
	JN	= 2*JNUM
	YJARG	= JNUM*YARG
	do	3	j = 1,NA1
	YYAA	= (ARRIN(j)-ARROUT(j,1))*YDT
	YA(j,JN)   = YA(j,JN)  +YYAA*cos(YJARG)
	YA(j,JN+1) = YA(j,JN+1)+YYAA*sin(YJARG)
 3	continue
	if(JNUM.lt.NUMBER)	goto	2
	return

 10	continue
	JCALL = 1
	YTOLD = TIME
	if(JCALL.eq.0)	YDT	= 0.
	JNUM	= 0
	do	4	j = 1,NA1
		ARROUT(j,1)	= FREQ*YA(j,1)
		YA(j,1)		= ARRIN(j)*YDT
 4	continue
	if(NUMBER .eq. 0)	return

	YMAX1 = 0.
 5	JNUM	= JNUM+1
	JN	= 2*JNUM
	YJARG	= JNUM*YARG
	do	6	j = 1,NA1
	YAMPL	= 2.*FREQ*sqrt(YA(j,JN)**2+YA(j,JN+1)**2)
	if (JNUM .eq. 1)	YMAX1 = max(YMAX1,abs(YAMPL))
	if (YA(j,JN) .ne. 0.)	then
		YPHASE	= atan(YA(j,JN+1)/YA(j,JN))
				else
		YPHASE	= 1.5707963
				endif
	if (YA(j,JN) .lt. 0.)	then
		if (YA(j,JN+1).gt.0.)	YPHASE = YPHASE+GP
		if (YA(j,JN+1).lt.0.)	YPHASE = YPHASE-GP
				endif
	ARROUT(j,JN)   = YAMPL
	ARROUT(j,JN+1) = YPHASE
	YYAA	  = (ARRIN(j)-ARROUT(j,1))*YDT
	YA(j,JN)  = YYAA*cos(YJARG)
	YA(j,JN+1)= YYAA*sin(YJARG)
 6	continue

C Remove jumps in phase:
	YGP=0.
	do	j=2,NA
		YY = ARROUT(j-1,JN+1)-ARROUT(j,JN+1)-YGP
		if( abs(YY) .gt. GP) YGP=YGP+sign(GP2,YY)
		ARROUT(j,JN+1) = ARROUT(j,JN+1)+YGP
	enddo
	ARROUT(NA1,JN+1) = ARROUT(NA,JN+1)

	if (JCALL .ne. 0)	then
C	if (JNUM  .eq. 1)	write(*,*)YMAX1
C Set phase to 0 when amplitude is too low
C	do	j=1,NA1
C		if( abs(ARROUT(j,JN)/YMAX1) .lt. .001) YGP=YGP+sign(GP2,YY)
C		ARROUT(j,JN+1) = ARROUT(
C	enddo
	endif
	if (JNUM .lt. NUMBER)	goto	5
	END
