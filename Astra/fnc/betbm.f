C BETBMR [ ]:	Beta beam poloidal (r)
C
C   2*c*c		    2*c*c
C  ------*{Int(p*dS)-p} = - ------ * Int[S*dp]
C    J*J		     J*J
C
	double precision function BETBMR(YR)
	implicit none
	double precision YR,YB
	integer  J,JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	YB	=0.
	JK	=YR/HRO+0.5
	IF(JK.lt.1)	JK=2
	IF(JK.gt.NA)	JK=NA
	DO	1	J=2,JK
	YB	=YB+ELON(J)*AMETR(J)**2*(PBPER(J-1)-PBPER(J+1))
 1	CONTINUE
	BETBMR	=6.4E-4*GP2*YB*	0.5*
     .		(RTOR/(G22(JK)*IPOL(JK)*BTOR*JK*HRO*MU(JK)))**2
	END
