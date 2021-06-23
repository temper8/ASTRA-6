C LINT []:	Internal inductance li(r) []
C			Pereverzev 27-DEC-90
C			/
C	      c*c	|
C     Li = ------------*|(Bpol**2)dV	  li = Li/(2*pi*R0)
C	   4*pi*Ipl*Ipl |
C			/
C			V
C
C   SI:	Li[H] = 2*Wi/Ipl**2 ;	li[dim.less] = Li[mkHn]/(0.2*pi*R0[m])
C   		     2*pi*R0[m]*li[dim.less] = 10*Li[mkHn]
C
	double precision function LINTR(YRO)
	implicit none
	double precision YRO,YK,YR
	integer  J,J1,JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	if (YRO .le. HRO)	then
	   YK = 1.
	elseif (YRO .gt. ROC)	then
	   YK = ROC/HRO
	else
	   YK = YRO/HRO
	endif
	JK = YK
	YR = JK+1.
	if (JK .ge. NA)	then
	   JK = NA
	   YR = JK+HROA/HRO
	endif
	LINTR = 0.
	do	J=1,JK
	   J1	 = J+J-1
	   LINTR = LINTR+J1*(2.*J*J-J1)*IPOL(J)*G22(J)/J*MU(J)**2
	enddo
	LINTR = LINTR+(YK*YK-JK*JK)*(YK*YK+JK*JK)*MU(JK+1)**2
     *		*IPOL(JK+1)*G22(JK+1)/YR
	LINTR = 0.5*LINTR*HRO*
     *		(YR/(YK*YK*MU(JK+1)*IPOL(JK+1)*G22(JK+1)))**2
	end
