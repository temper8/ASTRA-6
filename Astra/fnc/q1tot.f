C Q1TOT [[F1]/s]:	 Integral {0,R} (SF1TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q1TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q1TOTR=VINT(SF1TOT,YR)
	END
