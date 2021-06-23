C QBTOT [MW]:	 Integral {0,R} (PBEAM) dV
C			(Polevoy 28.09.89)
	double precision FUNCTION QBTOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include  'for/status.inc'
	QBTOTR=VINT(PBEAM,YR)
	END
