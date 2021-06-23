C=======================================================================
	double precision   function   WBPOLR(YR1)
C-----------------------------------------------------------------------
C WBPOL [MJ]:	Energy contents in the poloidal magnetic field
C	       /                               /
C         1    |                          1    | d\psi    IPOL
C Wi = -------*|(Bpol**2)dV (Gauss)  = -------*|(-----)^2 ---- G22 d\rho (SI)
C       8*pi   |                        2\mu_0 | d\rho    RTOR
C	       /                               /
C              V                               0
C
C  	Wi = 0.5*L_i*(I_pl)^2;		SI units: J, Hn, A
C
C			(G.Pereverzev 14-NOV-94)
C			   (corrected  6-JUL-99)
C-----------------------------------------------------------------------
	implicit none
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	double precision	YR,YR1
	integer	J,JK,J1
	YR	=YR1/HRO
	JK	=YR
	WBPOLR	=(YR*YR-JK*JK)*(YR*YR+JK*JK)*MU(JK+1)**2
     +		*IPOL(JK+1)*G22(JK+1)/(JK+1)
	DO	1	J=1,JK
	J1	=J+J-1
 1	WBPOLR	=WBPOLR+(J1*MU(J)**2*(2.*J*J-J1))*IPOL(J)*G22(J)/J
	WBPOLR	=1.25*GP*WBPOLR*HRO*HRO*HRO*BTOR*BTOR/RTOR
	end
C=======================================================================
