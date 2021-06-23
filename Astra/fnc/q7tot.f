C Q7TOT [[F7]/s]:	 Integral {0,R} (SF7TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q7TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q7TOTR=VINT(SF7TOT,YR)
	END
