C QDTF [MW]:	 Integral {0,R} ( PDT ) dV
C			(Yushmanov 11-JAN-89)
C       The function is obsolete (Pereverzev MAY-2006)
	double precision FUNCTION QDTFR(YR)
	implicit double precision (y)
C	implicit none
	double precision PDTF
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QDTFR=0.
	DO 1 J=1,JK
	include  'fml/pdtf'
 1	QDTFR=QDTFR+PDTF*VR(J)
	QDTFR=HRO*(QDTFR-PDTF*YDR)
	END
