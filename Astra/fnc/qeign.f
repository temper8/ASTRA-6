C QEIGN [MW]:	 Integral {0,R} (PEIGN) dV
C			(Yushmanov 11-JAN-89)
C Geometry effect (VR -> G11) included by G.Pereverzev 21-JAN-2006
	double precision FUNCTION QEIGNR(YR)
	implicit none
	double precision YDR,YR,PEIGN
	integer  J,JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
C	INCLUDE 'for/yrjkdr.inc'
	QEIGNR=0.
	if(YR.le.0.)	return
	JK=YR/HRO+1
	if(JK.gt.NA)	JK=NA
	YDR=(JK-YR/HRO)*G11(JK)
	if(YR.ge.ROC)	YDR=(JK-ROC/HRO)*G11(JK)
	DO 1 J=1,JK
	include  'fml/peign'
 1	QEIGNR=QEIGNR+PEIGN*G11(J)
	QEIGNR=HRO*(QEIGNR-PEIGN*YDR)
	END
