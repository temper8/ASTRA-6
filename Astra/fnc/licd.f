C LICD []: Internal inductance li(r) due to RF(NB) driven portion of a current.
C	   Assumption is made that the driven current radial distribution
C	   is the same as for parallel total current: <j||.B>/B0
C			Pereverzev 5-Nov-92
C			/
C	      c*c	|
C     Li = ------------*|(Bpol**2)dV	  li = Li/(2*pi*R0)
C	   4*pi*Ipl*Ipl |
C			/
C			V
C
C   SI:	Li[H] = 2*Wi/Ipl**2 ;	li[dim.less] = Li(mkHn)/(0.2*pi*R0(m))
C
	double precision FUNCTION LICDR(YR1)
	implicit none
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	integer	J,J1,JR
	double precision	YR1,YR,YCD1,YMUCD,YSINT
	IF(YR1.LE.HRO)	THEN
		YR = 1.
			ELSE
		YR = YR1/HRO
			ENDIF
	JR    = YR
	IF(JR.GE.NA)	JR = NA
	YCD1  = .1/(GP*BTOR)
	YMUCD = YCD1*CD(1)*VR(1)/IPOL(1)**2
	YSINT = 0.
	DO	1	J=1,JR
	J1    = J+J-1
	YSINT = YSINT+(J1*(YMUCD/J)**2*(2.*J*J-J1))*IPOL(J)/G22(J)/J
	YMUCD = YMUCD+YCD1*CD(J+1)*VR(J+1)/IPOL(J+1)**2
 1	continue
	YSINT = YSINT+(YR*YR-JR*JR)*(YR*YR+JR*JR)*(YMUCD/(JR+1))**2
     +		*IPOL(JR+1)/G22(JR+1)/(JR+1)
	if (YMUCD .lt. 0.001)	then
		LICDR = 0.
	else
	LICDR = .5*YSINT*HRO*((JR+1.)**2/(YR*YR*YMUCD*IPOL(JR+1)))**2
	endif
	END
