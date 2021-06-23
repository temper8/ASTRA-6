C Q3TOT [[F3]/s]:	 Integral {0,R} (SF3TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q3TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q3TOTR=VINT(SF3TOT,YR)
	END
