C Q9TOT [[F9]/s]:	 Integral {0,R} (SF9TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q9TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q9TOTR=VINT(SF9TOT,YR)
	END
