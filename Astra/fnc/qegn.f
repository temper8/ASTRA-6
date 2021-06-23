C QEGN [MW]:	 Integral {0,R} (PEGN) dV [=5/2*V'(R)*Gn(R)*Te(R)/625]
C			(Yushmanov 20-FEB-89)
	double precision FUNCTION QEGNR(YR)
	implicit none
	double precision YDR,YR
	integer  JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	JK=YR/HRO
		IF(JK.GE.NA)	THEN
	QEGNR=.001*(G11(NA)+G11(NA1))*GN(NA)*(TE(NA+1)+TE(NA))
	RETURN
		ELSEIF(JK.LT.0)	THEN
	QEGNR=0.
	RETURN
		ELSEIF(JK.EQ.0)	THEN
	QEGNR=YR*(G11(1)+G11(2))*GN(1)*(TE(2)+TE(1))/(HRO*1000.)
	RETURN
				ENDIF
	YDR=YR/HRO-JK
	QEGNR=((1.-YDR)*(G11(JK)+G11(JK+1))*GN(JK)*(TE(JK+1)+TE(JK))
     1	+YDR*(G11(JK+2)+G11(JK+1))*GN(JK+1)*(TE(JK+2)+TE(JK+1)))*.001
	END
