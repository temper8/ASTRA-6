C The subroutine provides a total plasma current IPL ramp-up
C with a rate dI/DT limited as {.1MA/s < dI/dt < 1MA/s}
C			(Pereverzev July-07-89)
	subroutine	CONIPL(IST,IEN)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	J
	double precision IST,IEN,FRAMP,T,T1,T2,YRAMP,YQmin,YQedge
	FRAMP(T1,T2,T)	=MIN(1.d0,MAX(0.d0,(T-T1)/(T2-T1)))
C	FJUMP(T1,T)	=MIN(1.d0,MAX(0.d0,(T-T1)/1.D-6))
	if (IPL.gt.IEN)	goto	99
	YQmin	=0.
	do	1	J=1,NA1
	if (YQmin.le.1.d0/MU(J))	goto	1
	YQmin	=1.d0/MU(J)
1	continue
	YQedge	=1.d0/MU(NA)
	YRAMP	=1.d0-0.9*FRAMP(1.d0,1.05d0,YQedge-YQmin)
	IPL	=IPL+YRAMP
	if (IPL.gt.IEN)	IPL	=IEN
99	end
