C LI3 []:	Internal inductance li(r) [4]
C			Pereverzev 03-APR-08
C			/
C	      c*c	|
C     Li = ------------*|(Bpol**2)dV	  li = Li/(2*pi*R0)
C	   4*pi*Ipl*Ipl |
C			/
C			V
C
C   SI:	Li[H] = 2*Wi/Ipl**2 ;	li[dim.less] = Li[mkHn]/(0.2*pi*R0[m])
C   		     2*pi*R0[m]*li[dim.less] = 10*Li[mkHn]
C----------------------------------------------------------------------|
C The same as LINT except for calculation of the denominator 
C     as an integral of the current density
C----------------------------------------------------------------------|
	double precision function LI3R(YRO)
	implicit none
	double precision YRO,YK,YR,YIPL,YWBP
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
	YWBP = (YK*YK-JK*JK)*(YK*YK+JK*JK)*IPOL(JK+1)*MU(JK+1)**2
     *		*G22(JK+1)/(JK+1)
	YIPL = 0.
	do	J=1,JK
	   J1	 = J+J-1
	   YWBP = YWBP+J1*(2.*J*J-J1)*IPOL(J)*G22(J)/J*MU(J)**2
C	   YIPL = YIPL+0.25*(CU(j)+CU(j+1))*(VR(j)+VR(j+1))/IPOL(j)**2
	   YIPL = YIPL+0.5*(CU(j)*VR(j)+CU(j+1)*VR(j+1))/IPOL(j)**2
	enddo
C	YWBP = YWBP*HRO**3	! *1.25*GP*BTOR*BTOR/RTOR	! = WBPOL
C	YIPL = YIPL*HRO		! *IPOL(JK)/(GP2*RTOR)		! = IPL
	LI3R = 2.*HRO*YWBP*(5*GP*BTOR/IPOL(JK+1)/YIPL)**2
	end
