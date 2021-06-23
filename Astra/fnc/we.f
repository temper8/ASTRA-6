C WE [MJ]:	Integral {0:R} ( 3/2*NE*TE ) dV
C			(Yushmanov 11-JAN-89)
	double precision FUNCTION WER(YR)
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/yrjkdr.inc'
	WER=0.
	DO 1 J=1,JK
 1	WER=WER+NE(J)*TE(J)*VR(J)
	WER=HRO*(WER-NE(JK)*TE(JK)*YDR)*.0024
	END
