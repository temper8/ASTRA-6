C Q5TOT [[F5]/s]:	 Integral {0,R} (SF5TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q5TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q5TOTR=VINT(SF5TOT,YR)
	END
