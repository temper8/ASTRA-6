	SUBROUTINE	BSAL(BSTOK,BSTOT)
C	Jbs Alpha by C.S.CHANG /XIV Int.Conf. on Plasma Phys. & Contr.
C	Nucl. Fusion. Res., Wurzburg, Germany, 30 Sept.-7 Oct. 1992
C	IAEA-CN-56/E-3-8 
C			Polevoy 22.10.92,08-SEP-99
	include	'for/parameter.inc' 
	include 'for/const.inc'
	include 'for/status.inc'
	double precision BSTOK(*),YARR(NRD),IINT
	external IINT
C	Psi Hirsh = 2 PI Psi ASTRA THEN
	do	1	j=1,NA1
	include 'fml/pdtf'
C	YLIE	=15.
	YLIE	=14.7+log(TE(J)*SQRT(NE(J)*.1))
	YLII	=24.+0.5*log(TE(J)/NE(J))
	Y=(0.75*1.7725/1836.*YLII/YLIE*
     .	(ABS(NDEUT(J)/2.+NTRIT(J)/3.)/NE(J)))**0.666667
C	Y=(0.75*1.77245/1836.)**0.6666667
	TSALF	=2.*SQRT(TE(J))*TE(J)/NE(J)/YLIE
C	ECRIT	=4.*1836.*TE(J)*Y
	SB3	=0.1*ZEF(J)*(TE(J)*.05)*SQRT(TE(J)*.05)
	SC3	=SB3/0.75
	C1	=1.17*SQRT(SB3)/(1.+5.25*SC3)
	C2	=(1.+3.7*SC3-2.1*SC3*SC3)/(3.+30.*SC3)	
c	E0SATS	=3.52E3*NTRIT(J)*NDEUT(J)*SVDT*TSALF
 	E0SATS	=625.*PDTF*TSALF
     *		*(C1+(3.*C2-2.*C1+(C1-2.*C2)*SQEPS(J))*SQEPS(J))
 1	YARR(J)	=E0SATS
	DO	2	J=1,NA1-1
	YFP	=1.+(0.46*(1.+2.1/ZEF(J))*SQEPS(J)-
     -			1.46*(1.+.67/ZEF(J)))*SQEPS(J)
 2	BSTOK(J)=RTOR*GP2*1.6E-3*SQEPS(J)*(1.-2.*YFP/ZEF(J))
     *		*(YARR(J+1)-YARR(J))/(FP(J)-FP(J+1))
	BSTOK(NA1)=0.
	BSTOT	=IINT(BSTOK,ROC)
	END
