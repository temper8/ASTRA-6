C NECH [10#19/m#3]: Chord average density (r) [m]
C	Integral {0,r} ( NE ) dl / a
C			(Yushmanov 11-MAY-87)
	double precision function NECHR(YD)
	implicit none
	double precision YD,YR,YXO,YX,YDEL
	integer  J,JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	IF(YD.GE.HRO*NA)	THEN
		NECHR	=0.
		RETURN
			ENDIF
	JK=0
	YXO=0.
	DO 1 J=1,NA1
	YR=HRO*(J-.5)
	IF(J.EQ.NA1)	YR=HRO*NA
	IF(YR.LE.YD)	GO TO 1
	YX=SQRT(YR**2-YD**2)
	IF(JK.NE.0)	GO TO 2
	JK=1
	IF(J.NE.1)	GO TO 3
	NECHR=(YX-YXO)*2.*NE(1)
			GO TO 4
 3	YDEL=(YR-YD)/HRO
	NECHR=(YX-YXO)*(YDEL*NE(J-1)+(2.-YDEL)*NE(J))
			GO TO 4
 2	NECHR=NECHR+(YX-YXO)*(NE(J-1)+NE(J))
 4	YXO=YX
 1	CONTINUE
	NECHR=NECHR*.5/(HRO*NA1)
	end
