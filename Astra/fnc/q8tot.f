C Q8TOT [[F8]/s]:	 Integral {0,R} (SF8TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q8TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q8TOTR=VINT(SF8TOT,YR)
	END
