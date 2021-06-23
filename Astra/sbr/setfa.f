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
C----------------------------------------------------------------------|
C This subroutine differs from SETFC by restricting the minimum step of
C    change in the actuator and the minimum time delay between two
C    successive changes. Also the minimum YOUT limit is removed.
C
C Input:
C    YIN   Signal to be controlled
C    YS    Coefficient at the controlled signal
C    YD    Coefficient at its time derivative
C    YI    Coefficient at time integral
C
C Output:
C    YOUT  Control signal
C----------------------------------------------------------------------|
C Limitation:
C    YOUT+CIMP2*NE(NA1) > 0   or   YOUT+CIMP2 > 0
C Example:
C    SETFA(CV5-CV1,CV4,CHE3,CHE4,0.d0)::.1<;
C       Here: CV5 is controllable quantity (to be equated to CV1)
C	      CV4 is actuator
C             (CHE3,CHE4,0.) are the control parameters (YS,YD,YI) 
C
C						(Pereverzev 04-APRIL-03)
C======================================================================|
	subroutine SETFA(YIN,YOUT,YS,YD,YI)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	integer	J1
	double precision
     1		YIN,YOUT,YS,YD,YI,YINT,YINO,YTO,YSTO,Y1,YDT,YDOUT
	save	YINO,YINT,YTO,YSTO,YDT,YDOUT
	data	YTO /-1.d10/ YINT/0.d0/ YDT/0.d0/ YDOUT/1.d-4/
C----------------------------------------------------------------------|
	if (YTO .lt. -1.d9)	then		! 1st entry only
	   YTO  = TIME
	   YINO = YIN
	   YSTO = YOUT
	   return
	endif
	if (TIME .le. YTO+YDT)	return		! (t_1-t_0 < dt) condition

	YINT = YINT+0.5*(YIN+YINO)*(TIME-YTO)
	Y1 = YD*(YIN-YINO)+(TIME-YTO)*(YS*YIN+YI*YINT)
	J1 = Y1/YDOUT+0.5
	if (J1 .eq. 0)		return
C 1st option (linear law):
	YOUT = YSTO+YDOUT*J1
	YOUT = max(-CIMP2,YOUT)			! NIZ2=(F2+CIMP2)*NE
C	YOUT = max(-CIMP2*NE(NA1),YOUT)		! NIZ2= F2+CIMP2*NE
C 2nd option (exponential law):
C	YOUT = YSTO+Y1*YOUT
C	write(*,'(1P,8E10.2)')TIME,YIN,YOUT
C     >		,(TIME-YTO)*YS*YIN,YD*(YIN-YINO)
	YSTO = YOUT
	YTO  = TIME
	YINO = YIN
	return
	end
