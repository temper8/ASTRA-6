C WALF [MJ]:	Integral {0,R} (Walf) dV
C			(Polevoy 06-JUN-95)
	double precision function WALFR(YR)
C	implicit none
	implicit double precision (y)
	double precision PDTF,YSALF,YDW,YX2,FNB2
	include	'for/parameter.inc'
	include	'for/const.inc'	
	include	'for/status.inc'
	include 'for/yrjkdr.inc'
	WALFR=0.
	DO 1 J=1,JK
	Y=(0.75*1.7725/1836.*24./17.*NI(J)/NE(J)/2.5)**0.666667
	YSALF	=2.*SQRT(TE(J))*TE(J)/(14.78+LOG(TE(J)))/NE(J)
	ECRIT	=4.*1836.*TE(J)*Y
	INCLUDE 'fml/pdtf'
	YX2	=3530./ECRIT
	YDW	=PDTF*YSALF*(1.-2.*FNB2(YX2)/YX2)*VR(J)
 1	WALFR=WALFR+YDW
	WALFR=HRO*(WALFR-YDW*YDR)*.5
	END
