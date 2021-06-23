C QOH [MW]:	 Integral {0,R} (POH) dV
C			(Pereverzev 12-FEB-90)
	double precision FUNCTION QOHR(YR)
	implicit none
	double precision POH
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QOHR=0.
	DO 1 J=1,JK
	include  'fml/poh'
 1	QOHR=QOHR+POH*VR(J)
	QOHR=HRO*(QOHR-POH*YDR)
	END
