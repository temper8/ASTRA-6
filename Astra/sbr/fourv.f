C=======================================================================|
	SUBROUTINE	FOURV(VARIN,FREQ,ID,NOHARM,VAROUT)
C---!---!---!---!---!---!---!---!---!---!---!---!---!---!---!---!---!---!
C Note!
C This is an old version of the subroutine FOURTR which does not use ID
C-----------------------------------------------------------------------|
C Description:	  Fourier expansion in time for the function VARIN(t)
C						of the single argument "t"
C Input: VARIN	- periodic function of time to be expanded
C	 FREQ	- fundamental harmonic frequency [Hz]
C	 ID	- subroutine call ID (1<=ID<=5)
C		  has to be given in order to avoid interference in stored data
C
C NB!		- it is implied that power oscillations have no phase
C		  shift with respect to the moment TIME = 0
C	 NOHARM	- 2*NOHARM+1 is a number of harmonics on output
C		  0 <= NOHARM <= 10
C
C Output:VAROUT(ro,#) - a0, a1, b1, a2, b2, ...
C	 - " -		Fourier harmonics
C
C Example call from a model:
C			CV1=betajb;	FOURV(CV1,CF2,1,1,CV6):;
C			CV2=WTOTB;	FOURV(CV2,CF2,2,5,CAR7):;
C			CV3=WETOTB;	FOURV(CV3,CF2,3,5,CAR8(12)):;
C	Here the first call will return ao, a1, b1 coef. of BETAJB(t)
C				in CV6, CV7, CV8, respectively.
C	The first 11 elements of the array CAR7 will be filled by
C			harmonics of the magnitude WTOTB, and so on.
C
C-----------------------------------------------------------------------|
	implicit none
	double precision VARIN,FREQ,VAROUT(*),YYA(5,21),YTOLD(5)
	double precision YSTART,YDT,YAMPL,YARG,YPHASE,YJARG,YYAA,YVAA
	integer	ID,NOHARM,JCALL(5),JNUM,JN,JTIME
	include	'for/parameter.inc'
	include  'for/const.inc'
	save	JCALL,YYA,YTOLD
	data	JCALL/5*0/
	YARG	= GP2*FREQ*TIME
	if(JCALL(ID) .eq. 0)	goto	10
	YDT	= TIME-YTOLD(ID)
	JTIME	= TIME*FREQ
	YSTART	= JTIME/FREQ
	if(TIME .gt. YSTART .and. YTOLD(ID) .le. YSTART)  goto	10

	JNUM	= 0
	YTOLD(ID)	= TIME
	YYA(ID,1)	= YYA(ID,1)+VARIN*YDT
	if(NOHARM .eq. 0)	return

 2	JNUM	= JNUM+1
	JN	= 2*JNUM
	YJARG	= JNUM*YARG
	YYAA	= (VARIN-VAROUT(1))*YDT
	YYA(ID,JN)  =YYA(ID,JN)  +YYAA*cos(YJARG)

	YYA(ID,JN+1)=YYA(ID,JN+1)+YYAA*sin(YJARG)
	if(JNUM.lt.NOHARM)	goto	2
	return

 10	if(JCALL(ID) .eq. 0)	YDT	= 0.
	JNUM	= 0
	JCALL(ID) = 1
	YTOLD(ID) = TIME
	VAROUT(1) = FREQ*YYA(ID,1)
	YYA(ID,1) = VARIN*YDT
	if(NOHARM .eq. 0)	return

 5	JNUM	= JNUM+1
	JN	= 2*JNUM
	YAMPL	= sqrt(YYA(ID,JN)**2+YYA(ID,JN+1)**2)

	if(YYA(ID,JN).ne.0.)	then
		YPHASE	= atan(YYA(ID,JN+1)/YYA(ID,JN))
				else
		YPHASE	= 1.5707963
				endif
	if( (YYA(ID,JN).lt.0.) )	then
		if (YYA(ID,JN+1).gt.0.)	YPHASE = YPHASE+GP
		if (YYA(ID,JN+1).lt.0.)	YPHASE = YPHASE-GP
				endif

	VAROUT(JN)	= 2.*FREQ*YAMPL
	VAROUT(JN+1)	= YPHASE
	YYAA	= (VARIN-VAROUT(1))*YDT
	YJARG	= JNUM*YARG
	YYA(ID,JN)	= YYAA*cos(YJARG)
	YYA(ID,JN+1)	= YYAA*sin(YJARG)
	if(JNUM.lt.NOHARM)	goto	5
	END
