C======================================================================|
	subroutine	FOURTR(VARIN,FREQ,HARM,VAROUT)
C----------------------------------------------------------------------|
C Description:	  Fourier expansion in time for the function VARIN(t)
C					of the single argument "t"
C Input: VARIN	- periodic function of time to be expanded
C	 FREQ	- fundamental harmonic frequency [Hz]
C
C	 HARM	- a number of harmonics on output is 2*HARM+1
C		  0 <= HARM <= 10
C
C Output: 
C	VAROUT(*) cointains amplitude and phase: A0, A1, ph1, A2, ph2,..
C 		  the phase is calculated as sin(w*t+phase)
C	   or Fourier harmonics a0, a1, b1, a2, b2, ...
C
C Example call from a model:
C		CV1=betajb;	FOURTR(CV1,CF2,1,CV6):;
C		    returns Fourier coef. of BETAJB(t) in CV6, CV7, CV8
C		CV2=WTOTB;	FOURTR(CV2,CF2,5,CAR7):;
C		    11 elements of the array CAR7 will be filled by
C		    harmonics of the magnitude WTOTB, and so on.
C----------------------------------------------------------------------|
	implicit none
	integer		NLOC
	parameter	(NLOC=220)
C  For DEC Alpha use in the next line:	integer*8
	integer LY(NLOC),L
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer	ID,JD,j,JN,JNUM,NOHARM
	double precision YYA(NLOC,21),YTOLD(NLOC),VARIN,VAROUT(*)
	double precision FREQ,HARM,YDT,YARG,YJARG,YVDT,YAMPL,YPHASE
	save	ID,LY,YYA,YTOLD
	data	ID/0/
	NOHARM = HARM
	L=loc(VARIN)
	if(ID.eq.0)	goto 2
	do 	1	J=1,ID
	   JD = J
	   if(L.eq.LY(J))	goto 3
 1	continue

 2	continue
C  1st call
	ID=ID+1
	if(ID.gt.NLOC)	then
		write(*,*)'Too many Fourier integrals: >',NLOC
		ID=NLOC
		return
	endif
C location of the 1st call
	LY(ID)=L
	YTOLD(ID) = TIME
	do	j = 1,2*NOHARM+1
	   YYA(ID,j) = 0.
	   VAROUT(j) = 0.
	enddo
	return

C  repeated call
 3	continue
	YDT	= TIME-YTOLD(JD)
	YARG	= GP2*FREQ*TIME
	JNUM	= 0
	YYA(JD,1)= YYA(JD,1)+VARIN*TAU
	if(NOHARM .eq. 0)	return
 4	JNUM	= JNUM+1
	JN	= 2*JNUM
	YJARG	= JNUM*YARG
	YVDT	= (VARIN-VAROUT(1))*TAU
	YYA(JD,JN)  =YYA(JD,JN)  +YVDT*cos(YJARG)
	YYA(JD,JN+1)=YYA(JD,JN+1)+YVDT*sin(YJARG)
	if(JNUM.lt.NOHARM)	goto	4

C  YDT*FREQ is a minimal number of periods for integration
	if (YDT*FREQ .lt. 1.)	return
	
	VAROUT(1)= YYA(JD,1)/YDT
	do	j=1,NOHARM
	   JN	= 2*j
C----------------------------------------------------------------------|
C Amplitude and phase: 
C	phase is determined as sin(omega*t+phase)
C
	   YAMPL = sqrt(YYA(JD,JN)**2+YYA(JD,JN+1)**2)
	   YPHASE = 0.5*GP-asin(YYA(JD,JN+1)/YAMPL)
	   if (YYA(JD,JN).lt.0.)   YPHASE = 2*GP-YPHASE
	   VAROUT(JN)	= 2.*YAMPL/YDT
	   VAROUT(JN+1)	= YPHASE
C----------------------------------------------------------------------|
C Re & Im
C
C	   VAROUT(JN)	= 2.*YYA(JD,JN)/YDT
C	   VAROUT(JN+1)	= 2.*YYA(JD,JN+1)/YDT
	enddo
C	write(*,*)YYA(JD,2),YYA(JD,3),VAROUT(2),VAROUT(3)/GP,CV14
C	write(*,*)JD,NOHARM
	YTOLD(JD) = TIME
	do	j=1,2*NOHARM+1
	   YYA(JD,j) = 0.
	enddo
	END
C======================================================================|
