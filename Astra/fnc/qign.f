C QIGN [MW]: Integral {0,R} (PIGN) dV [=5/2*V'(R)*Gn(R)*Ti(R)*ni(R)/ne(R)/625]
C			(Yushmanov 20-FEB-89)
	double precision FUNCTION QIGNR(YR)
	implicit none
	double precision YDR,YR
	integer  JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	JK=YR/HRO
		IF(JK.GE.NA)	THEN
	QIGNR=.001*(G11(NA)+G11(NA1))*GN(NA)*(TI(NA+1)+TI(NA))
     1	*(NI(NA+1)+NI(NA))/(NE(NA+1)+NE(NA))
	RETURN
		ELSEIF(JK.LT.0)	THEN
	QIGNR=0.
	RETURN
		ELSEIF(JK.EQ.0)	THEN
	QIGNR=.001*YR*(G11(1)+G11(2))*GN(1)*(TI(2)+TI(1))/HRO
     1	*(NI(2)+NI(1))/(NE(2)+NE(1))
	RETURN
				ENDIF
	YDR=YR/HRO-JK
	QIGNR=((1.-YDR)*(G11(JK)+G11(JK+1))*GN(JK)*(TI(JK+1)+TI(JK))
     1	*(NI(JK+1)+NI(JK))/(NE(JK+1)+NE(JK))
     2	+YDR*(G11(JK+2)+G11(JK+1))*GN(JK+1)*(TI(JK+2)+TI(JK+1))
     3	*(NI(JK+1)+NI(JK+2))/(NE(JK+1)+NE(JK+2)))*.001
	END
