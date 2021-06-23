C Q4TOT [[F4]/s]:	 Integral {0,R} (SF4TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q4TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q4TOTR=VINT(SF4TOT,YR)
	END
