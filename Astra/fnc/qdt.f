C QDT [MW]:	 Integral {0,R} ( PDT ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QDTR(YR)
	implicit none
	double precision PDT,SVDT
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDTR=0.
	DO 1 J=1,JK
	include  'fml/pdt'
 1	QDTR=QDTR+PDT*VR(J)
	QDTR=HRO*(QDTR-PDT*YDR)
	END
