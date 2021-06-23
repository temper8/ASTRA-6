C FTLLMR [%]:	
C Effective trapped particle fraction Lin-Liu and Miller GA-A21820
C                                                        (Oct 94)
C Y.R.Lin-Liu and R.L.Miller, Phys.Plasmas 2(5), May 1995, pp.1666-1668
	double precision function FTLLMR(YR)
	implicit none
	double precision YR,YYR,YH,YFTUP,YFTLO
	integer JK
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
        FTLLMR=0.1
	IF(YR.LE.0.)	RETURN
	YYR=YR
	JK=YYR/HRO+1
	IF(JK.GT.NA) THEN
	   YYR=HRO*NA
	   JK=NA
	ENDIF
        IF (ABS(BMAXT(JK)).LT.0.0001 .OR. 
     +      ABS(BDB0(JK)).LT.0.0001) THEN
           FTLLMR=0.1
           GOTO 99
        ENDIF
        YH=min(.999999d0,BDB0(JK)/BMAXT(JK)*BTOR)
        YFTUP=1.-BDB02(JK)/BDB0(JK)**2*
     +        (1.-SQRT(1.-YH)*(1.+.5*YH))
        YFTLO=1.-BDB02(JK)*FOFB(JK)
	FTLLMR=.75*YFTUP+.25*YFTLO
99	RETURN
	END

