C=======================================================================|
	SUBROUTINE	FOURARR(CYCLE,ARRIN,FREQ,YNUM,ARROUT)
C-----------------------------------------------------------------------|
C Description:	    Fourier expansion in time for the function ARRIN(a,t)
C							of 2 arguments
C Input: CYCLE    - numbers of period for integration
C	 ARRIN(*) - periodic radial dependent function of time to be expanded
C	 FREQ	  - fundamental harmonic frequency [Hz]
C	 NUMBER	  - number of harmonics in output is 2*NUMBER+1 <= 21
C
C Output:ARROUT(ro,#) - a0(ro), a1(ro), fi1(ro), a2(ro), fi2(ro), ...
C	 		amplitudes and phases of Fourier harmonics
C			as functions of radius and #
C
C Example call from a model:	FOUREX(TE,CF2,5,CAR6):;
C	The Fourier harmonics ao(r), a1(r), fi1(r),..., fi5(r) of TE(r,t)
C	will be returned in CAR6, CAR7, CAR8,..., CAR16 respectively.
C
C Important: The subroutine can be call only once, otherwise the 
C				second call will spoil the previous. 
C	     To avoid the interference the same subroutine can be called
C				under another name (see sbr/fourex.for)
C-----------------------------------------------------------------------|
        IMPLICIT NONE
	include	'for/parameter.inc'
	include  'for/const.inc'
	INTEGER number,j,jnum,jtime,jtimeold,jn
	double precision ARRIN(*),ARROUT(NRD,*),YA(NRD,21),cycle
	double precision freq,ynum,yarg,ydt,ytold,yyaa,yphase,yjarg,yampl
	double precision phasenear,diff
c
	save	YA
c
	number = ynum
	YARG	= GP2*FREQ*TIME !Phase of the main frequency
	YDT	= TIME-YTOLD
	jtime	= TIME*FREQ/cycle
c Integer number of completed periods 	
c If the period has finished, define ARROUT and return it to ASTRA:
	if(jtime .gt. jtimeold) goto 10
c Integration for the 0th Fourier coefficient (time average):
	do 1 j = 1,NA1
	   YA(j,1) = YA(j,1) + ARRIN(j)*YDT
 1	continue
c Integration for the Fourier coefficients:
 	do 4 JNUM = 1, number  !JNUM = armonic number
           JN = 2*JNUM
	   YJARG = JNUM*YARG !Phase of the JNUMth harmonic
	   do 3 j = 1,NA1
	      YYAA = (ARRIN(j)-ARROUT(j,1))*YDT
	      YA(j,JN) = YA(j,JN) + YYAA*cos(YJARG)
	      YA(j,JN+1) = YA(j,JN+1) + YYAA*sin(YJARG)
 3	   continue
 4	continue
	jtimeold = jtime
	YTOLD	= TIME
	return
C End of the partial summation, go back to ASTRA
 10	continue
c Division by period, return Fourier coefficients to ASTRA:	
	do 24 j = 1,NA1
	  ARROUT(j,1) = FREQ*YA(j,1)/cycle !Division by period
	  YA(j,1) = ARRIN(j)*YDT !Begin for the next integration
 24	continue
c
 5	do 7 JNUM = 1, number
	  JN = 2*JNUM
	  YJARG = JNUM*YARG
	  phasenear = 0.
	    do 6 j = 1,NA1
	       YAMPL = 2.*FREQ*sqrt(YA(j,JN)**2+YA(j,JN+1)**2)/cycle
	       if(YA(j,JN).ne.0.) then
	          YPHASE = atan(YA(j,JN+1)/YA(j,JN))
	       else
	          YPHASE = 0.5*GP
	       endif
c
	       if (YA(j,JN) .lt.0.) then
	          if (YA(j,JN+1) .gt.0.) YPHASE = YPHASE + GP
	          if (YA(j,JN+1) .lt.0.) YPHASE = YPHASE - GP
	       endif
c Avoid 2*pi jumps:
	       diff = YPHASE - phasenear
	       if (diff .gt. 0) then
		  if (diff .gt. GP) YPHASE = YPHASE - GP2
	       else
		  if (diff .lt. - GP) YPHASE = YPHASE + GP2
	       endif
c
	       ARROUT(j,JN) = YAMPL
	       ARROUT(j,JN+1) = YPHASE
c Begin for the next integration:
	       YYAA = (ARRIN(j)-ARROUT(j,1))*YDT
	       YA(j,JN) = YYAA*cos(YJARG)
	       YA(j,JN+1) = YYAA*sin(YJARG)
	       phasenear = YPHASE
 6	   continue
 7	continue
c
	YTOLD = TIME
	jtimeold = jtime
	return
	END
