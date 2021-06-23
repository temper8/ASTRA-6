C Q0TOT [[F0]/s]:	 Integral {0,R} (SF0TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q0TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q0TOTR=VINT(SF0TOT,YR)
	END
