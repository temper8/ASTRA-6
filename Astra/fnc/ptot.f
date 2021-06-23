C PTOT [MW/m#2]:	Total plasma heating power (r) [m]
C			(Yushmanov 12-MAY-87)
	double precision FUNCTION PTOTR(YR)
	implicit none
	double precision YR,YDH
	integer  JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	JK=YR/HRO+.5
	IF(JK.GE.NA1)		PTOTR=PETOT(NA1)+PITOT(NA1)
	IF(JK.LE.0)		PTOTR=PETOT(1)+PITOT(1)
	IF(JK.GE.NA1.OR.JK.LE.0)	RETURN
	YDH=YR/HRO+(.5-JK)
	PTOTR=
     >	YDH*(PETOT(JK+1)+PITOT(JK+1))+(1.-YDH)*(PETOT(JK)+PITOT(JK))
	END
