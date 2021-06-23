C NELA [10#19/m#3]: Density Line Average in the mid-plane (r) [m]
C	Integral {0,r} ( NE ) dl / a
C			(Pereverzev 22-12.97)
C 
	double precision FUNCTION NELAR(YD)
	implicit none
	double precision YD,YX,YXO,YDEL,YR
	integer  J,JK
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	IF(YD.GE.ABC)	THEN
		NELAR	=0.
		RETURN
			ENDIF
	JK=0
	YXO=0.
	DO 1 J=1,NA1
	YR=AMETR(j)
	IF(J.EQ.NA1)	YR=ABC
	IF(YR.LE.YD)	GO TO 1
	YX=SQRT(YR**2-YD**2)
	IF(JK.NE.0)	GO TO 2
	JK=1
	IF(J.NE.1)	GO TO 3
	NELAR=(YX-YXO)*2.*NE(1)
			GO TO 4
 3	YDEL=(YR-YD)/HRO
	NELAR=(YX-YXO)*(YDEL*NE(J-1)+(2.-YDEL)*NE(J))
			GO TO 4
 2	NELAR=NELAR+(YX-YXO)*(NE(J-1)+NE(J))
 4	YXO=YX
 1	CONTINUE
	NELAR=NELAR*.5/(ABC)
	END
