C Q6TOT [[F6]/s]:	 Integral {0,R} (SF6TOT) dV
C			(Pereverzev 09-May-2008)
	double precision FUNCTION Q6TOTR(YR)
	implicit none
	double precision YR,VINT
	include 'for/parameter.inc'
	include 'for/status.inc'
	Q6TOTR=VINT(SF6TOT,YR)
	END
