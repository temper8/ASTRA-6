C QEDT [MW]:	 Integral {0,R} ( PEDT ) dV
C			(Yushmanov 11-JAN-89)
	double precision function QEDTR(YR)
	implicit double precision (y)
	double precision PEDT,SVDT,PAION,PDT
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QEDTR=0.
	DO 1 J=1,JK
	include  'fml/pedt'
 1	QEDTR=QEDTR+PEDT*VR(J)
	QEDTR=HRO*(QEDTR-PEDT*YDR)
	END
