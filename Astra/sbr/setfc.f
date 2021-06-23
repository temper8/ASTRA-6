C======================================================================|
C Adjust the control parameter YOUT in order to keep
C      the controlled quantity YIN at zero
C
C Algorithm:  The output signal is evaluated as
C                                         t
C       d(YOUT)               d(YIN)      /
C       ------- = YS*YIN + YD*------ + YI*|(YIN*dt)
C         dt                    dt        /
C                                         0
C   alternatively as
C	                                     t
C      d[ln(YOUT)]               d(YIN)      /
C      ----------- = YS*YIN + YD*------ + YI*|(YIN*dt)
C            dt                    dt        /
C                                            0
C Input:
C    YIN   Signal to be controlled
C    YS    Coefficient at the controlled signal
C    YD    Coefficient at its time derivative
C    YI    Coefficient at time integral
C Output:
C    YOUT  Control signal
C
C Example:
C    SETFC(CV5-CV1,CV4,CHE3,CHE4,0.)::.1<;
C       Here: CV5 is controlled quantity (to be equated to CV1)
C	      CV4 is actuator
C             (CHE3,CHE4,0.) are the control parameters (YS,YD,YI) 
C
C						(Pereverzev 04-APRIL-03)
C======================================================================|
	subroutine SETFC(YIN,YOUT,YS,YD,YI)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	double precision	YIN,YOUT,YS,YD,YI,YINT,YINO,YTO,YSTO,Y1
	save	YINO,YINT,YTO,YSTO
	data	YTO /-1.E10/ YINT/0./
C----------------------------------------------------------------------|
	if (YTO .lt. -1.E9)	then		! 1st entry only
	   YTO  = TIME
	   YINO = YIN
	   YSTO = YOUT
	   return
	endif
	if (YTO .ge. TIME)	return		! No time advance in between

	YINT = YINT+0.5*(YIN+YINO)*(TIME-YTO)
	Y1 = YD*(YIN-YINO)+(TIME-YTO)*(YS*YIN+YI*YINT)
C 1st option (linear law):	(YS,YD,YI) = (0.5,3.,0.)
	YOUT = YSTO+Y1
	YOUT = max(1.d-1,YOUT)
C 2nd option (exponential law):	(YS,YD,YI) = (0.05,.5,0.)
C	YOUT = YSTO+Y1*YOUT
C	write(*,'(1P,8E10.2)')TIME,YIN,YOUT
C     >		,(TIME-YTO)*YS*YIN,YD*(YIN-YINO)
	YTO  = TIME
	YINO = YIN
	YSTO = YOUT
	return
	end
