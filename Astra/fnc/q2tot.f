C Q2TOT [[F2]/s]:	 Integral {0,R} (SF2TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q2TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q2TOTR=VINT(SF2TOT,YR)
	END
