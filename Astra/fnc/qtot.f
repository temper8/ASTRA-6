C QTOT [MW]:	 Integral {0,R} ( PETOT+PITOT ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION QTOTR(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	QTOTR=0.
	DO 1 J=1,JK
 1	QTOTR=QTOTR+(PETOT(J)+PITOT(J))*VR(J)
	QTOTR=HRO*(QTOTR-(PETOT(JK)+PITOT(JK))*YDR)
	END
