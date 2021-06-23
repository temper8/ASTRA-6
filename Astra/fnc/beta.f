C BETA %:	Beta toroidal  % 
C
C SI:	BETA=2\mu_0*p/B^2 = 12.8e-3*pi*(n_20)*(T_keV)/(B_T)^2 * 100[%]
C CGS:	BETA=8\pi*p/B^2   = 12.8e-3*pi*[0.1*NE*TE+...]/BTOR^2 * 100[%]
C
C    NOTE! NBI Pressure is included as in equilibrium
C					(Pereverzev 02-MAY-2006)
	double precision function BETAR(YR)
	implicit  none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision YV
	include 'for/yrjkdr.inc'
	YV = 0.
	BETAR = 0.
	do	J=1,JK
           YV = YV+VR(J)
           BETAR = BETAR+VR(j)*(TE(J)*NE(J)+TI(J)*NI(J)
     .		+.5*(PBLON(J)+PBPER(J)))		!*NB2EQL
        enddo
	YV = YV-YDR
	BETAR = BETAR-YDR*(TE(JK)*NE(JK)+TI(JK)*NI(JK)
     .		+.5*(PBLON(JK)+PBPER(JK)))
	BETAR = 0.402*BETAR/(YV*BTOR*BTOR)
	return
	end
