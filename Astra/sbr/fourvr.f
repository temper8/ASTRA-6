C=======================================================================|
	SUBROUTINE	FOURVR(CYCLE,VARIN,FREQ,YNUM,VAROUT)
C-----------------------------------------------------------------------|
C Description:	    Fourier expansion in time for the function ARRIN(a,t)
C							of 2 arguments
C Input: CYCLE   - numbers of period for integration
C	 VARIN   - time dependent variable to be expanded
C	 FREQ	 - fundamental harmonic frequency [Hz]
C	 NUMBER	 - number of harmonics in output is 2*NUMBER+1 <= 21
C
C Output:VAROUT(#) - a0, a1, fi1, a2, fi2, ...
C	 		amplitudes and phases of Fourier harmonics
C			as functions of #
C
C Example call from a model:	FOURVR(3.,TE,CF2,5.,CAR31):;
C	The Fourier harmonics ao(r), a1(r), fi1(r),..., fi5(r) of TE(r,t)
C	will be returned in the first 5 places of CAR31.
C
C Important: The subroutine can be call only once, otherwise the 
C				second call will spoil the previous. 
C-----------------------------------------------------------------------|
        implicit none
	include	'for/parameter.inc'
	include  'for/const.inc'
	integer number,j,jnum,jtime,jtimeold,jn
	double precision VARIN,VAROUT(*),YAVAR(21),cycle
	double precision freq,ynum,yarg,ydt,ytold,yyaa,yphase,yjarg,yampl
	save	YAVAR
c
	number = ynum
	YARG	= GP2*FREQ*TIME !Phase of the main frequency
	YDT	= TIME-YTOLD
	jtime	= TIME*FREQ/cycle
c If the period has finished, define VAROUT and return it to ASTRA:
	if(jtime .gt. jtimeold) goto 10
c Integration for the 0th Fourier coefficient (time average):
	YAVAR(1) = YAVAR(1) + VARIN*YDT
c Integration for the Fourier coefficients:
 	do 4 JNUM = 1, number  !Harmonic number
           JN = 2*JNUM
	   YJARG = JNUM*YARG !Phase of the JNUMth harmonic
	   YYAA = (VARIN - VAROUT(1))*YDT
	   YAVAR(JN) = YAVAR(JN) + YYAA*cos(YJARG)
	   YAVAR(JN+1) = YAVAR(JN+1) + YYAA*sin(YJARG)
 4	continue
	jtimeold = jtime
	YTOLD	= TIME
	return
C End of the partial summation, go back to ASTRA
 10	continue
c Division by period, return Fourier coefficients to ASTRA:	
	VAROUT(1) = FREQ*YAVAR(1)/cycle !Division by period
	YAVAR(1) = VARIN*YDT !Begin for the next integration
c
 5	do 7 JNUM = 1, number
	  JN = 2*JNUM
	  YJARG = JNUM*YARG
	  YAMPL = 2.*FREQ*sqrt(YAVAR(JN)**2+YAVAR(JN+1)**2)/cycle
	  if(YAVAR(JN).ne.0.) then
	     YPHASE = atan(YAVAR(JN+1)/YAVAR(JN))
	  else
	     YPHASE = 1.5707963
	  endif
c
	  if (YAVAR(JN) .lt.0.) then
	     if (YAVAR(JN+1) .gt.0.) YPHASE = YPHASE + GP
	     if (YAVAR(JN+1) .lt.0.) YPHASE = YPHASE - GP
	  endif
c
	  VAROUT(JN) = YAMPL
	  VAROUT(JN+1) = YPHASE
c Begin for the next integration:
	  YYAA = (VARIN - VAROUT(1))*YDT
	  YAVAR(JN) = YYAA*cos(YJARG)
	  YAVAR(JN+1) = YYAA*sin(YJARG)
 7	continue
c
	YTOLD = TIME
	jtimeold = jtime
	return
	END
