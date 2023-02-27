C Modules: >> ANLSTR,ARRNAM,ASKINT,ASXWIN,ASTWIN,WRFIGS <<
C          >> CMARK, CMARKT,CMARKP,CHECKU,CURV1  <<
C	   >> FMTF4, FMTF5, FMTF6, FMTXF4,FMTXF5,FMXF5, FMT456 <<
C	   >> KILLBL,length,NMARK, nextwd,noblan,PLOTGR,PLOTXY,PUTMAR <<
C	   >> RECTAN,ROUNDN,READF6,RZ2A,  RA2Z,  QUADIN,APPRID,nsymb  <<
C	   >> SCAL,  SCALA, SMOOTH,SORT,  SORTAB,TAMP,	UPCASE,VARNAM <<
C======================================================================|
	subroutine	SCAL(NOUT,SN,SO,OUT,NP,NDIM)
C----------------------------------------------------------------------|
C  Input:
C     NOUT	Total number of channels (NROUT or NTOUT)
C     SO(NOUT)	Current scales
C     NP	Upper index boundary for scale definition (= NDIM-3)
C     NDIM	Upper index boundary according to the array description
C     OUT(NDIM,NOUT) Array of curves for scale definition (ROUT or TOUT)
C  Output:
C     SN(NOUT)	Returned scales
C----------------------------------------------------------------------|
	implicit none
	integer NOUT,NP,J,JJ,NDIM
	double precision	SN(*),SO(*),OUT(NDIM,*),SC,SCALA
C----------------------------------------------------------------------|
C Exclude 3 edge/central points
	do	1	J=1,NOUT
	if(SO(J))	4,3,2
 2	SN(J)=SO(J)
				GO TO 1
 3	SN(J)=SCALA(OUT(1,J),NP)
				GO TO 1
 4	do	5	JJ=1,J-1
	if(SO(J).EQ.SO(JJ))	GO TO 7
 5	continue
	SC=0.
	do	6	JJ=J,NOUT
	if(SO(J).NE.SO(JJ))	GO TO 6
	SC=MAX(SC,SCALA(OUT(1,JJ),NP))
	if(ABS(SO(JJ)+JJ).LT..01) GO TO 8
 6	continue
	SN(J)=SC
				GO TO 1
 7	SN(J)=SN(JJ)
				GO TO 1
 8	SN(J)=SCALA(OUT(1,JJ),NP)
 1	continue
	end
C======================================================================|
	double precision	function SCALA(Y,NJ)
	implicit none
	integer NJ,J
	double precision	YS(6),Y(*),YMAX
	YMAX	=0.
	do	5	J=1,NJ
		if(YMAX.LT.ABS(Y(J)))	YMAX=ABS(Y(J))
 5	continue
	YS(1)	=1.0d-9
	YS(2)	=1.5d-9
	YS(3)	=2.0d-9
	YS(4)	=3.0d-9
	YS(5)	=5.0d-9
	YS(6)	=8.0d-9
 1	continue
	do 2 J	=1,6
		if (1.05*YMAX.LE.YS(J)) GO TO 4
 2	continue
	do 3 J	=1,6
		YS(J)	=10.d0*YS(J)
 3	continue
	go to 1
 4	SCALA	=YS(J)
	end
C======================================================================|
	subroutine	CMARK(NL,JPOS,SC,OS,NAME,STYL)
C Mark variable/scale in 1 & 2 modes
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	double precision	SC,OS,YS
	character*4	ST*10,NAME,F4,F5*5
	integer	NL,JY,JPOS,STYL,POINT(2)
	call	FMTF4(ST(1:4),SC)
	ST(5:5)=' '
	ST(6:10)=NAME
	if (OS)	1,4,2
 1	YS = -OS
	call	FMTF4(F4,YS)
	F5 = '-'//F4
	goto	3
 2	call	FMTF4(F4,OS)
	F5 = '+'//F4
 3	if (NL .lt. YSCMAX/2)	then
	    JY = DYLET+2
	else
	    JY = -DYLET-2
	endif
	call	textvm(JPOS+4*DXLET,NL+JY,F5,5)
 4	call	textvm(JPOS,NL,ST,10)
	if (STYL .le. 0)	return
	POINT(1) = JPOS+4*DXLET+4
	POINT(2) = NL-4
	call	NMARK(POINT,STYL)
	end
C======================================================================|
	subroutine	CMARKT(NL,JPOS,SC,OS,NAME,STYL)
C Mark variable/scale in 6th (time) mode
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	double precision  SC,OS,YS
	character*4 ST*10,NAME,F4
	integer NL,J,JPOS,JP,STYL,length,killbl,PT(2)
	ST(1:10)='          '
	ST(2:5)=NAME
	if (OS)	1,2,3
 1	YS = -OS
	call FMTF4(F4,YS)
	ST(6:10)='-'//F4
	JP=JPOS-DXLET
	if (ST(1:1).eq.' ')	JP=JP-DXLET
	j = killbl(ST(2:2),9)+1
	if (j .le. 6)		JP=JP+DXLET
	goto	4
 2	JP=JPOS
	J = 5
	goto	4
 3	call FMTF4(F4,OS)
	ST(6:10)='+'//F4
	JP=JPOS-DXLET
	if (ST(1:1).eq.' ')	JP=JP-DXLET
	j = killbl(ST(2:2),9)+1
	if (j .le. 6)		JP=JP+DXLET
 4	continue
	call textvm(JP,NL,ST,J)			! type name+yshift
	PT(1) = JP+3
	PT(2) = NL-5
	call FMTF4(F4,SC)
	J=NL+DYLET
	JP=JPOS+DXLET
	call textvm(JP,J,F4,4)			! type scale
	if (STYL .le. 0)	return
	call	NMARK(PT,STYL)
	end
C======================================================================|
	subroutine	CMARKP(NL,JPOS,NAME,STYL)
C Mark variable/scale in 4 & 5 (time-radial) modes
	implicit none
	character ST*5,NAME*4
	integer NL,JPOS,STYL
	integer	NBIT(6)
C diamond(1), o(111), +(43), *(42), x(120), #(35), $(36), 
	data NBIT /1,111,43,42,120,36/
	ST(2:5)=NAME
	ST(1:1)=' '
	if (STYL.ge.7 .and. STYL.le.12)  ST(1:1) = char(NBIT(STYL-6))
	call textvm(JPOS,NL,ST,5)
	end
C======================================================================|
	subroutine	PLOTXY(YARR,NP,JX,IX,IXO,IY,IYO,DMET,STYL,PT)
C-------------------------------------------------------------------
C   The subroutine displays NP points
C	of integer array IY vs IX to screen with the style=STYL
C	and puts marks with time interval equal to DMET(sec)
C   Entry:	YARR - time array (TTOUT)
C		NP,JX,IY,IX,DMET,STYL
C-------------------------------------------------------------------
	implicit none
	integer IX(*),IY(*),IXO(*),IYO(*),PT(*)
	integer JX,NP,EraseColor,STYL,J,J0,JMET,JJ,JPOINT,J1
	double precision    DMETO,DMET,YARR(NP)
	save EraseColor,DMETO
	data EraseColor/31/,DMETO/999999./
	call colovm(EraseColor)
	do	4	J0=0,1
	   do	1	J=1,NP
C	      write(*,*)j,j0,JX,IXO(J),IYO(J)
	      PT(2*J-1)	=JX+IXO(J)
	      PT(2*J)	=350-IYO(J)
 1	   continue
	   JMET	=0
	   do	2	JJ=1,NP-1+J0
	      if(YARR(JJ)-YARR(1).ge.JMET*DMETO.or.JJ.eq.1)	then
		 if(JJ.gt.1)	call	curvvm(0,JPOINT+1,PT(J1))
		 J1=2*JJ-1
		 call NMARK(PT(J1),STYL)
		 JMET	=JMET+1
		 JPOINT	=1
	      else
		 JPOINT	=JPOINT+1
	      endif
 2	   continue
	   if (JPOINT.ne.0)	call	curvvm(0,JPOINT,PT(J1))
	   if (J0.eq.1)	return			! J0=0 <- erasing
	   do	3	J=1,NP
	      IXO(J)=IX(J)
 3	      IYO(J)=IY(J)
	   DMETO	=DMET
	   call colovm(1)
 4	continue
	end
C======================================================================|
	subroutine	PLOTCR(NP,NPO,IX,IXOLD,IY,IYOLD,ICOLOR,STYL,PT)
C-----------------------------------------------------------------------
C   The subroutine displays NP points of the integer array IY
C	NP  is a number of points to plot
C	NPO is a number of points to erase, in addition,
C	NPO is a control parameter:
C	 NPO > 0  the old curve is erased, the drawn one is stored in IYOLD 
C	 NPO <= 0 a new curve IY(1:NP) is drawn, (IXOLD,IYOLD) are NOT used
C	 NPO = 0  no erasure, the drawn curve is stored in IYOLD
C	The points of the array IYOLD are used for erasing
C	curve of the previous call and are determined inside PLOTG1
C	STYL
C	PT is a working array 2*NB1
C   Input:	NP,IY,IYOLD,ICOLOR,STYL
C   Output:	IXOLD,IYOLD
C-----------------------------------------------------------------------
	implicit none
	integer IX(*),IXOLD(*),IY(*),IYOLD(*),PT(*)
	integer STYL,NP,EraseColor,NPO,J,ICOLOR
	save EraseColor
	data EraseColor /31/
	if (NPO .le. 0)	goto	2
C erase the old curve
	do	1	J=1,NPO
	PT(2*J-1)	=IXOLD(J)
	PT(2*J)		=IYOLD(J)
 1	continue
	call colovm(EraseColor)
	call CURV1(NPO,PT,STYL)

C draw a new curve
 2	do	3	J=1,NP
	PT(2*J-1)	=IX(J)
	PT(2*J)		=IY(J)
 3	continue
Colors: 1(Red) 2(Blue) 3(MeduimSeeGreen) 4(VioletRed) 5(Brown) 6(LightBlue)
C	7(Turquoise)
C STYL	1 2 3 4 5 6 7 	8 9 10 11 12 13 14	15 16 17 18
Color:	1 2 3 4 5 6 7   1 2 3  4  5  6  7 	 1  2  3  4  5  6 7
	call colovm(ICOLOR)
	call CURV1(NP,PT,STYL)
	if (NPO .lt. 0)	return
	do	5	J=1,NP
	IYOLD(J)	=IY(J)
	IXOLD(J)	=IX(J)
 5	continue
	end
C======================================================================|
	subroutine	CURV1(NP,PT,STYL)
C The subroutine has replaced older subroutine CURV
	implicit none
	integer PT(*),STYL,NP,LE,NF,J,j1,JJ,NM,LENG(6),jx,jy,PT1(2)
	equivalence	(PT1(1),jx),	(PT1(2),jy)
	save	LENG
	data	LENG /8,10,13,17,22,28/
	if(STYL) 1,4,2

 1	continue				! Draw dashed curves
	jj = -STYL
	if (jj .ge. 7)	jj = jj+1-jj/7*7
	if (jj .eq. 1)	goto 4
	LE = NP/LENG(jj-1)

	LE = 8
	do	j=1,jj
	   j1 = j+1
	   LE = LE+j1
	enddo
	j1 = jj
	if (jj .eq. 2)	LE = min(LE,16)
	if (jj .eq. 3)	LE = min(LE,8)
	if (jj .eq. 4)	LE = min(LE,4)
	NF = max(1,LE/4)
	do	j=1,NP,LE
	   j1 = min(NP-j+1,LE-NF)
	   call drcurv(0,j1,PT(2*j-1))
	enddo
	return

 2	continue				! Draw markers
	NM=NP/5
	NM=max(10,NP/5)
	LE=NM/5*STYL				! 1st marker position
	if (LE .ge. NM+2)	LE=LE-NM
	LE = max(1,LE)
C	call	getcolor(j1)
	do	jj=LE,NP,NM
	   J  = 2*jj
	   jx = PT(j-1)
	   jy = PT(j)/10
	   call NMARK(PT1,STYL)
C	   call	puto(jx,jy,j1,STYL)
	enddo
 4	continue				! Draw a solid curve
	call drcurv(0,NP,PT(1))
	end
C======================================================================|
C	subroutine	PUTMAR(TIME,YARR,NP,JX,IX,IY,STYL)
C Called from review only
C-------------------------------------------------------------------
C	implicit none
C	integer IX(*),IY(*),JX,PT(2),NP,STYL,JJ
C	double precision	TIME,YARR(*)
C	do	2	JJ=1,NP
C	   if (YARR(JJ).le.TIME .and. TIME.lt.YARR(JJ+1))	then
C		PT(1)	=JX+IX(JJ)
C		PT(2)	=350-IY(JJ)
C		call NMARK(PT(1),STYL)
C		return
C	   endif
C2	continue
C	end
C======================================================================|
	subroutine	NMARK_(POINT,STYL)
C used when NMARK_ is called from C
	implicit none
	integer POINT(2),STYL
	call	NMARK(POINT,STYL)
	end
C======================================================================|
	subroutine	NMARK(POINT,STYL)
	implicit none
	integer POINT(2),PT(32),DX(16,7),DY(16,7),STYL,J,N(7),JJ,IST
C IST definition shoud coincide with NBIT() in CMARK
C IST = 1-filled diamond, 2-o, 3-+, 4-*, 5-x(#), 6-<, 7-filled square
C STYL  <=7,              8    9    10   11      12   >=13
	save	N,DX,DY
	data	N /16,13,5,9,14,9,10/
	data	DX/
     1	 0, 3, 0,-3, 0, 0, 2, 0,-2, 0, 0, 1, 0,-1, 0, 0,	! 
     2	 3, 3, 2, 1,-1,-2,-3,-3,-2,-1, 1, 2, 3, 0, 0, 0,	! o
     3	 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	! +
     4	-3,-1, 2,-2, 1,-2, 2,-1, 3, 0, 0, 0, 0, 0, 0, 0,	! *
C     5	 1, 1,-1,-1, 3,-3, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0,	! #
     5	-2,-1, 0, 1, 2, 1, 0, 1, 2, 1, 0,-1,-2, 0, 0, 0,	! x
     6	 0, 0,-3, 0,-2, 0, 0,-2, 0, 0, 0, 0, 0, 0, 0, 0,	! <
     7	-2, 2, 2,-2,-2, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0/	! 
	data	DY/
     1	 3, 0,-3, 0, 3, 2, 0,-2, 0, 2, 1, 0,-1, 0, 1, 0,	! 
     2	 1,-1,-2,-3,-3,-2,-1, 1, 2, 3, 3, 2, 1, 0, 0, 0,	! o
     3	 0, 0, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,	! +
     4	 0, 0, 3,-3, 0, 3,-3, 0, 0, 0, 0, 0, 0, 0, 0, 0,	! *
C     5	 3,-3, 3,-3, 1, 1,-1,-1, 0, 0, 0, 0, 0, 0, 0, 0,	! #
     5	-2,-1, 0, 1, 2, 1, 0,-1,-2,-1, 0, 1, 2, 0, 0, 0,	! x
     6	 3,-3, 0, 3,-1,-1, 1, 1,-3, 0, 0, 0, 0, 0, 0, 0,	! <
     7	-2,-2, 2, 2,-1,-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0/	! 
	IST=max(1,min(STYL,7))
	do	1	JJ=1,N(IST)
	J=2*JJ
	PT(J-1)=POINT(1)+DX(JJ,IST)
1	PT(J)  =POINT(2)+DY(JJ,IST)
	call curvvm(0,N(IST),PT(1))
C return to the original point:
Call drawvm(0,....,POINT(1),POINT(2))
	end
C======================================================================|
Convert real positive R1>0 to a string F4*4
	subroutine	FMTF4(F4,R1)
	implicit none
	character F4*4,F6*6
	double precision	R1
	integer	J
	if(abs(R1).lt.1.e-9)	then
		F4=' 0. '
		return
			endif
	if(R1.lt.0)	then
		F4=' <0 '
		return
			endif
C	if(R1.ge.1.e11)	then
C	   call	FMT456(F6,R1,5)
C	   write(*,*)R1,'"',F6,'"'
C	else
	   call	FMT456(F6,R1,4)
C	endif
	do	5	J=3,6
	if(F6(J:J).eq.'.')	go to 4
 5	continue
	F4=F6(3:6)
	if(F4(2:3).eq.'0-')	then
	   read(F4(4:4),*)J
	   if(J.le.8)write(F4,'(1A3,1I1)')F4(1:1)//'e-',J-1
	endif
	if(F4(1:1).eq.' ' .and. F4(3:4).eq.'-9')
     >				write(F4,'(1A4)')F4(2:2)//'e-9'
	if(F4(2:4).eq.'0-4')	F4 ='.00'//F4(1:1)
	if(F4(2:3).eq.'0e')	then
	   read(F4(4:4),*)J
	   if(J.lt.9)write(F4,'(2A1,1I1,1A1)')F4(1:1),'e',J+1,' '
	   if(J.ge.9)write(F4,'(2A1,1I2)')F4(1:1),'e',J+1
	endif
C	if(F4(2:4).eq.'0e3')	F4 =F4(1:1)//'e4 '
C	if(F4(2:4).eq.'0e4')	F4 =F4(1:1)//'e5 '
C	if(F4(2:4).eq.'0e5')	F4 =F4(1:1)//'e6 '
C	if(F4(2:4).eq.'0e6')	F4 =F4(1:1)//'e7 '
C	if(F4(2:4).eq.'0e7')	F4 =F4(1:1)//'e8 '
C	if(F4(2:4).eq.'0e8')	F4 =F4(1:1)//'e9 '
C	if(F4(2:4).eq.'0e9')	F4 =F4(1:1)//'e10'
C	if(F4(2:4).eq.'0e4')	F4 =F4(1:1)//'e11'
			return
 4	do	3	J=6,1,-1
	if(F6(J:J).eq.'.')	go to 2
	if(F6(J:J).ne.'0')	go to 1
 3	continue
 2	J=J-1
 1	F4=F6(3:J)
	end
C======================================================================|
C The same as FMTF4 but does not suppress leading "0." if possible
	subroutine	FMTF40(F4,R1)
	implicit none
	character F4*4
	double precision	R1
	call	FMTF4(F4,R1)
	if (F4(1:1).eq.'.' .and. F4(4:4).eq.' ') F4='0.'//F4(2:3)
	end
C====================================================================
Convert real R1 to a string F5*5
	subroutine	FMTF5(F5,R1)
	implicit none
	character F4*4,F5*5,F6*6
	double precision	R1,R2
	if(abs(R1).lt.1.e-9)	then
		F5=' 0.  '
		return
			endif
	if(R1.lt.0)	then
		R2 = abs(R1)
		call FMTF4(F4,R2)
		F5='-'//F4
		return
			endif
	call	FMT456(F6,R1,5)
	F5=F6(2:6)
	if (F5(3:5).eq."000")	F5(3:5)='   '
	if (F5(4:5).eq."00")	F5(4:5)='  '
	if (F5(5:5).eq."0")	F5(5:5)=' '
	end
C====================================================================
Convert real R1 to a string F5*5
C The same as FMTF5 but without suppressing end zeros
	subroutine	FMTF50(F5,R1)
	implicit none
	character F4*4,F5*5,F6*6
	double precision	R1,R2
	if(abs(R1).lt.1.e-9)	then
		F5=' 0.  '
		return
			endif
	if(R1.lt.0)	then
		R2 = abs(R1)
		call FMTF4(F4,R2)
		F5='-'//F4
		return
			endif
	call	FMT456(F6,R1,5)
	F5=F6(2:6)
	end
C====================================================================
Convert real R1 to a string F6*6
	subroutine	FMTF6(F6,R1)
	implicit none
	character F4*4,F5*5,F6*6
	double precision	R1,R2
	if(abs(R1).lt.1.e-9)	then
		F6='    0.'
		return
			endif
	if(R1.lt.0)	then
		R2 = abs(R1)
		call FMTF5(F5,R2)
		F6='-'//F5
		goto	1
			endif
	call	FMT456(F6,R1,6)
 1	if (F6(3:6).eq."0000")	F6='    '//F6(1:2)
	if (F6(4:6).eq."000")	F6='   '//F6(1:3)
	if (F6(5:6).eq."00")	F6='  '//F6(1:4)
	if (F6(6:6).eq."0")	F6=' '//F6(1:5)
	if (F6(3:6).eq."    ")	F6='    '//F6(1:2)
	if (F6(4:6).eq."   ")	F6='   '//F6(1:3)
	if (F6(5:6).eq."  ")	F6='  '//F6(1:4)
	if (F6(6:6).eq." ")	F6=' '//F6(1:5)
C	if (F6(3:6).eq."0000")	F6(3:6)='    '
C	if (F6(4:6).eq."000")	F6(4:6)='   '
C	if (F6(5:6).eq."00")	F6(5:6)='  '
C	if (F6(6:6).eq."0")	F6(6:6)=' '
	end
C===================================================================
	subroutine	FMTXF4(F4,R1)
	implicit none
	character F4*5,F6*6
	double precision	R1,R2
	F4=' '
	if(R1.lt.0.)	F4='-'
	R2=abs(R1)
	call	FMT456(F6,R2,4)
	if (F6(4:6) .eq. '0-4')	then
		F6(3:6)='.00'//F6(3:3)
	endif
	F4(2:5)=F6(3:6)
	end
C====================================================================
	subroutine	FMTXF5(F4,R1)
	implicit none
	character F4*6,F6*6
	double precision	R1,R2
	F4=' '
	if(R1.lt.0.)	F4='-'
	R2=abs(R1)
	call	FMT456(F6,R2,5)
	F4(2:6)=F6(2:6)
	end
C====================================================================
	subroutine	FMTXF6(F4,R1)
	implicit none
	character F4*7,F6*6
	double precision	R1,R2
	F4=' '
	if(R1.lt.0.)	F4='-'
	R2=abs(R1)
	call	FMT456(F6,R2,6)
	F4(2:7)=F6(1:6)
	end
C====================================================================
	subroutine	FMXF5(F4,R1)
C the same as FMTXF5 but suppresses right "." or "0" when possible
	implicit none
	character F4*6,F6*6
	integer	J,JM
	double precision	R1,R2
	R2=abs(R1)
	if(R1.lt.0.)	then
		JM=1
		F4='-'
			else
		F4=' '
		JM=0
			endif
		call	FMT456(F6,R2,6-JM)
	do	5	J=1,6
	if(F6(J:J).eq.'.')	go to 4
 5	continue
	F4(1+JM:6)=F6(1+JM:6)
			return
 4	do	3	J=6,1,-1
	if(F6(J:J).eq.'.')	go to 2
	if(F6(J:J).ne.'0')	go to 1
 3	continue
 2	J=J-1
 1	F4(7+JM-J:6)=F6(1+JM:J)
	if(R1.ge.0.)	return
	F4(1:1)=' '
	F4(6+JM-J:6+JM-J)='-'
	end
C===================================================================
	subroutine	FMT456(F4,R2,N)
C R2>0 only
	implicit none
	character F4*6,T*25
	double precision	R,R2,ROUNDN
	integer	J,N,JM
	JM=6-N
	F4=' '
	if(R2.lt.1.e-9)	then
		F4=' 0.000'
		return
			endif
	R=ROUNDN(R2,5-JM)
	if(R.lt.1.e-4*10.d0**JM)	then
		write(T,11)R
 11		format(1F25.11)
		do	3	J=15,20+JM
		if(T(J:J).ne.'0')	go to 4
 3		continue
		J=20+JM
 4		F4(1+JM:5)=T(J:J+3-JM)//'-'
		write(F4(6:6),12)J-11-JM
		do	2	J=1,6
		if(F4(J:J).ne.'0'.and.F4(J:J).ne.' ')	return
		F4(J:J)=' '
 2		continue
 12		format(1I1)
		return
			endif
	R=ROUNDN(R2,6-JM)
	if(R.lt.1.e6/10.d0**JM)	then
	if(R.ge.1.e5/10.d0**JM)	R=ROUNDN(R2,7-JM)
		write(T,11)R
		do	5	J=8+JM,17-JM
		if(T(J:J).ne.' '.and.T(J:J).ne.'0')	go to 6
 5		continue
 6		F4(1+JM:6)=T(J:J+5-JM)
				return
			endif
	R=ROUNDN(R2,5-JM)
	if(R.lt.1.e13/10.d0**JM)	then
		write(T,11)R
		do	7	J=1,11
		if(T(J:J).ne.' '.and.T(J:J).ne.'0')	go to 8
 7		continue
 8		F4(1+JM:5)=T(J:J+3-JM)//'e'
		write(F4(6:6),12)10-J+JM
		return
			endif
	F4='  >e1 '
	write(F4(6:6),12)3-JM
	end
C=================================================================
	double precision function ROUNDN(R1,N)
	implicit none
	integer	N,J
	double precision    R,R1
	R=R1
	ROUNDN=R1
	if(R.ge.1.e13.or.R.lt.1.e-9)	return
	R=R*1.e9
	do	1	J=1,25
	if(R.lt.10.d0)	go to 2
	R=R/10.d0
1	continue
2	R=(R+50.d0/10.d0**N)*10.d0**(J-10)
	ROUNDN=R
	end
C=================================================================
	subroutine	READF6(LINE,F6,IERR)
C Read a number in LINE to 6-positinal field 
	implicit none
	integer	IERR,N,JE,JP,JM,J
	double precision	F6
	character*6	LINE,S
	IERR=0
Cpole	READ(LINE,10,ERR=1)N
10	FORMAT(1I6)
Cpole	F6=N
Cpole			return
1	JE=0
	JP=0
	JM=0
	do	5	J=1,6
	if(LINE(J:J).eq.'-')	JM=J
	if(LINE(J:J).eq.'e'.or.LINE(J:J).eq.'E') JE=J
	if(LINE(J:J).eq.'.')	JP=J
 5	continue
	if(JP.gt.0)	go to 4
	if(JE.gt.0)	go to 3
Cpole	if(JM.le.1)	go to 2
Cpole	S=LINE(1:JM-1)
	if(JM.gt.1)	go to 22
	READ(LINE,10,ERR=2)N
	F6=N
			return
 22	S=LINE(1:JM-1)

	READ(S,10,ERR=2)N
	F6=N
	S=LINE(JM:6)
	READ(S,10,ERR=2)N
	F6=F6*10.d0**N
			return
Cpole3	if(JE.le.1)	go to 3
 3	if(JE.le.1)	go to 2
	S=LINE(1:JE-1)
	READ(S,10,ERR=2)N
	F6=N
	S=LINE(JE+1:6)
	READ(S,10,ERR=2)N
	F6=F6*10.d0**N
			return
 4	READ(LINE,11,ERR=2)F6
 11	FORMAT(1F6.3)
			return
2	IERR=1
	write(*,*)'>>> ERROR: in "',LINE,'"'
	end
C=================================================================
	subroutine	UCASE(STRI)
	implicit none
	integer   KEY,J
	character STRI*(*)
	do	J=1,len(STRI)
	   KEY = ichar(STRI(J:J))
	   if(KEY.ge.97 .and. KEY.le.122)	KEY = KEY-32
	   STRI(J:J) = char(KEY)
	enddo
	end
C=================================================================
	subroutine	UPCASE(LENG,STRI)
	implicit    none
	integer     KEY,J,LENG
	character*1 STRI(LENG)
	do	1	J=1,LENG
	KEY = ichar(STRI(J))
	if(KEY.ge.97 .and. KEY.le.122)	KEY=KEY-32
	STRI(J) = char(KEY)
1	continue
	end
C======================================================================|
C	integer function IFNUM(STR,N)
C	implicit	none
C	character*1 STR(*),NU(16)
C	integer I,J,N,jj
C	data NU/'0','1','2','3','4','5','6','7','8','9',
C     +		'.','E','D','+','-',' '/
C        jj = 0
C	IFNUM=1					! True
C	do  1	J=1,N
C           do	I=1,16
C              if (STR(J).eq.NU(I))	goto 1
C           enddo
C           IFNUM = 0				! False
C           return
C 1	continue
C	end
C IFNUM ==============================================
	integer function IFNUM(STR,N)
	implicit	none
	character*1 STR(*),STRI*80
	integer	N,j
	double precision	R
	write(STRI,'(80A1)')(STR(j),j=1,N),char(0)
	IFNUM=1
	if (N .lt. 1)	return
	read(STRI(1:N),*,ERR=1)R
	return
 1	IFNUM=0
	end
C======================================================================|
	double precision	function	GETNUM(FIELD,ERCODE)
C----------------------------------------------------------------------|
        implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	character FIELD*(*),ZNUM*6
        integer	l,j,j1,ISHIFT,IFNUM,length,ERCODE,READR(48)
	save	ISHIFT
        data	ISHIFT/0/	READR/		! Total 48 elements
     >		    1,12,23,34,44,45,46,47,48,	! ZRD1X ->ZRD9X
     >		 2, 3, 4, 5, 6, 7, 8, 9,10,11,	! ZRD10X->ZRD19X
     >		13,14,15,16,17,18,19,20,21,22,	! ZRD20X->ZRD29X
     >		24,25,26,27,28,29,30,31,32,33,	! ZRD30X->ZRD39X
     >		35,36,37,38,39,40,41,42,43/	! ZRD40X->ZRD48X
C----------------------------------------------------------------------|
	if (ISHIFT .eq. 0)	then
	   do	1  j=1,NPRNAM
	      if(PRNAME(j).ne.'ZRD1  ')	goto	1
	      ISHIFT = j-1
	      goto	2
 1	   continue
	endif
 2	l  = len(FIELD)
        j  = index(FIELD(1:l),'ZRD')
        if (j .ne. 0)			goto	5
        j  = index(FIELD(1:l),'C')
        if (j .eq. 0)			goto	99
	j1 = min(length(FIELD(j:l)),l-j+1)
	ZNUM = FIELD(j:j+j1-1)
	do  3	j=1,NCFNAM
	   j1 = j
	   if (CFNAME(j) .eq. ZNUM)	goto	4
 3	continue
	goto	99
 4	ERCODE = 0
        GETNUM = CONSTF(j1)
        return

 5	j1 = index(FIELD(j+3:l),'X')
        if (j1 .ne. 0)	then
           if (j1 .gt. 3)		goto	99
           ZNUM = FIELD(j+3:j+j1+1)
        else
           ZNUM = FIELD(j+3:)
        endif
        if (IFNUM(ZNUM(1:),2) .eq. 0)	goto	99
        read(ZNUM(1:),*)j1
        if (j1.lt.1 .or. j1 .gt. 48)	goto	99
        GETNUM = DEVARX(ISHIFT+READR(j1))
        ERCODE = 0
        return
 99	ERCODE = 1
	end
C======================================================================|
	subroutine	STREAD(NCH,NFIELD,ARRAY,ERCODE)
C----------------------------------------------------------------------|
C	Reads one record group of the NBINP ("*.nbi") file 
C	and fills ARRAY(1:NFIELD) with data.
C	Numbers and references to ZRD??, ZRD??X and CONSTF_list
C	are allowed as records in the input file.
C ERCODE values:
C   0 -	Normal exit
C   1 - Unrecognized variable name
C   2 -	Read error
C   3 -	Wrong format
C   4 -	Array out of limits
C   5 -	Missing records
C----------------------------------------------------------------------|
        implicit none
        integer	NCH,NFIELD,j,IFNUM,ERCODE
	character    SFIELD(20)*12,STRING*132
        double precision	ARRAY(NFIELD),GETNUM
C----------------------------------------------------------------------|
        ERCODE = 4
	if (NFIELD .gt. 20)	return
        ERCODE = 5
 1      read(NCH,'(A)',err=99,end=3)STRING 		! Skip string
C        j = index(STRING,'//')				! Comment
        j = index(STRING,'!')				! Comment
        if (j .eq. 1)	goto	1
        if (j .eq. 0)	goto	2
        goto	97
 2      read(NCH,'(5A)',err=99)(SFIELD(j),j=1,NFIELD)
C	write(*, '(5A)')       (SFIELD(j),j=1,NFIELD)
        do	j=1,NFIELD
           call	UPCASE(12,SFIELD(j))
           if (IFNUM(SFIELD(j),12) .ne. 0)	then
              read(SFIELD(j),*,err=98)ARRAY(j)
           else
              ARRAY(j) = GETNUM(SFIELD(j),ERCODE)
              if (ERCODE .ne. 0)		goto	99
           endif
        enddo
        ERCODE = 0
 3      return
 97     continue
        ERCODE = 3
        return
 98     continue
        ERCODE = 2
        return
 99     continue
        ERCODE = 1
        end
C======================================================================|
	subroutine	ASKINT(NV,NVAR,NAME)
	implicit none
	include		'for/parameter.inc'
	integer     NVAR(*),NV,J
	character*4 NAME(*)
	double precision        VAR(NRW)
	do	1 J=1,NV
1	VAR(J)=NVAR(J)
	call ASKLIS(NV,VAR,NAME,4)
	do	2 J=1,NV
2	NVAR(J)=VAR(J)
	end
C======================================================================|
	subroutine	ASTWIN(NB,IBOX,NAME,SCALE,SHIFT,MOD10,YMODE)
	implicit	none
	double precision	SCALE(*),SHIFT(*)
	integer		NB,IBOX(*),YMODE,NRW16,NRW96,MOD10
	include		'for/parameter.inc'
	parameter	(NRW16=NRW+16,NRW96=NRW-128)
	integer		JMODE,JGR,j,j1,j2,jj,jn,jm,js,jb,jw,jsep,IB(NRW)
	integer		IP1(NRW),IP2(NRW),IP30(NRW),IP31(NRW),ASKTAB
	character*80	rows(NRW16),STR,TITLE,NAME(*)*4,KEY*1,TZER*1
	data	IP1/
     1		 1, 3, 5, 7,   9,11,13,15,   2, 4, 6, 8,  10,12,14,16,
     2		17,19,21,23,  25,27,29,31,  18,20,22,24,  26,28,30,32,
     3		33,35,37,39,  41,43,45,47,  34,36,38,40,  42,44,46,48,
     4		49,51,53,55,  57,59,61,63,  50,52,54,56,  58,60,62,64,
     5		65,67,69,71,  73,75,77,79,  66,68,70,72,  74,76,78,80,
     6		81,83,85,87,  89,91,93,95,  82,84,86,88,  90,92,94,96,
     7	 97,99,101,103, 105,107,109,111, 98,100,102,104,106,108,110,112,
     8	113,115,117,119,121,123,125,127,114,116,118,120,122,124,126,128/
C     9		NRW96*0/
	data	IP2/
     1		 1, 2, 9,10,   3, 4,11,12,   5, 6,13,14,   7, 8,15,16,
     2		17,18,25,26,  19,20,27,28,  21,22,29,30,  23,24,31,32,
     3		33,34,41,42,  35,36,43,44,  37,38,45,46,  39,40,47,48,
     4		49,50,57,58,  51,52,59,60,  53,54,61,62,  55,56,63,64,
     5		65,66,73,74,  67,68,75,76,  69,70,77,78,  71,72,79,80,
     6		81,82,89,90,  83,84,91,92,  85,86,93,94,  87,88,95,96,
     7	  97,98,105,106, 99,100,107,108,101,102,109,110,103,104,111,112,
     8	113,114,121,122,115,116,123,124,117,118,125,126,119,120,127,128/
C     9		NRW96*0/
	data	IP30/
     1		 1, 3, 5, 7,   2, 4, 6, 8,   9,11,13,15,  10,12,14,16,
     2		17,18,21,22,  25,26,29,30,  19,20,23,24,  27,28,31,32,
     3		33,34,37,38,  41,42,45,46,  35,36,39,40,  43,44,47,48,
     4		49,50,53,54,  57,58,61,62,  51,52,55,56,  59,60,63,64,
     5		65,66,69,70,  73,74,77,78,  67,68,71,72,  75,76,79,80,
     6		81,82,85,86,  89,90,93,94,  83,84,87,88,  91,92,95,96,
     7	  97,98,101,102,105,106,109,110, 99,100,103,104,107,108,111,112, 
     8	113,114,117,118,115,116,119,120,121,122,125,126,123,124,127,128/
C     9		NRW96*0/
	data	IP31/
     1		 1, 2, 5, 6,   3, 4, 7, 8,   9,10,13,14,  11,12,15,16,
     2		17,18,21,22,  19,20,23,24,  25,26,29,30,  27,28,31,32,
     3		33,34,37,38,  35,36,39,40,  41,42,45,46,  43,44,47,48,
     4		49,50,53,54,  51,52,55,56,  57,58,61,62,  59,60,63,64,
     5		65,66,69,70,  67,68,71,72,  73,74,77,78,  75,76,79,80,
     6		81,82,85,86,  83,84,87,88,  89,90,93,94,  91,92,95,96,
     7	  97,98,101,102,105,106,109,110, 99,100,103,104,107,108,111,112, 
     8	113,114,117,118,115,116,119,120,121,122,125,126,123,124,127,128/
C     9		NRW96*0/
	j = 0
	TZER= char(j)
	TITLE = "Presentation"//TZER
C              ----5----0----5----0----5----0----5----0----5----0----5
	STR = "Name|Box| Scale|Offset||Name|Box| Scale|Offset"//TZER
C	STR = "    |   |      |      ||    |   |      |      "//TZER
	jsep = index(STR,'||')+1
	JMODE = MOD10
	if (JMODE.eq.1)	then
	    JGR = 8
	elseif (JMODE.eq.2 .or. JMODE.eq.3)	then
	    JGR = 4
	elseif (JMODE.eq.6)	then
	    if (YMODE .eq. 1)	JGR = 4
	    if (YMODE .eq. 0)	JGR = 4
	    if (YMODE .eq.-1)	JGR = 2
	else
	    return
	endif
 10	continue
	do	j = 1,NRW
	   write(rows(j)(1:80),'(79X,1A1)')TZER
	enddo

C j  - ordinal box No.
C jb - box No. in the Astra nominations
C jw - position in the table
C jj - horizontal row
C js - position in the current row
C
	jn = 0
	jm = 0
	do  20	j = 1,NB
	    jb = IBOX(j)
	    if (jb .le. 0 )	then
		jm = jm+1
		goto	20
	    endif
	    if (JMODE .eq. 1)	jw = IP1(jb)
	    if (JMODE .eq. 2 .or. JMODE .eq. 3)	jw = IP2(jb)
	    if (JMODE .eq. 6)	then
		if (YMODE .eq. 1)	jw = jb
		if (YMODE .eq. 0)	jw = IP30(jb)
		if (YMODE .eq.-1)	jw = IP31(jb)
	    endif
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(rows(jj)(js+1:js+4),'(1A4)') NAME(j)
	    write(rows(jj)(js+6:js+8),'(1I3)') jb
	    call	num2str(SCALE(j),rows(jj)(js+10:js+15),6)
	    call	num2str(SHIFT(j),rows(jj)(js+17:js+22),6)
	    jn = max(jn,jj)
 20	continue
	if (jm .eq. 0 )	goto	22
	j2 = jn
	j1 = 2*jn+1
	do  21	j = 1,NB
	    if (IBOX(j) .gt. 0 )	goto	21
	    jw = j1
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(rows(jj)(js+1:js+4),'(1A4)') NAME(j)
	    write(rows(jj)(js+6:js+8),'(1I3)') IBOX(j)
	    call	num2str(SCALE(j),rows(jj)(js+10:js+15),6)
	    call	num2str(SHIFT(j),rows(jj)(js+17:js+22),6)
	    jn = max(jn,jj)
	    j1 = j1+1
 21	continue
 22	continue

	j1 =(jm+1)/2
	j = ASKTAB(TITLE,STR,rows,80,jn,JGR,j1)
	if (j .gt. 0)	goto	22
C	call	ASKCOL(TITLE,STR,rows,80,jn,JGR,j1)
C
C Warning: If two (or more) box # coincide, information will be 
C	   partly overrode and lost
C
C	write(*,'(1(A22,3X,A22))')
C     >		(rows(j)(1:22),rows(j)(25:46),j=1,jn)
C	write(*,'(10I5)')j1,NB,jn

	jn = 0
	do  30	j  = 1,NB
	    if (IBOX(j) .le. 0 )	goto	30
	    if (JMODE .eq. 1)	jw = IP1(IBOX(j))
	    if (JMODE .eq. 2 .or. JMODE .eq. 3)	jw = IP2(IBOX(j))
	    if (JMODE .eq. 6)	then
		if (YMODE .eq. 1)	jw = IBOX(j)
		if (YMODE .eq. 0)	jw = IP30(IBOX(j))
		if (YMODE .eq.-1)	jw = IP31(IBOX(j))
	    endif
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(NAME(j),'(1A4)',ERR=77) rows(jj)(js+1:js+4)
	    read(rows(jj)(js+6: js+8), *,ERR=77) IB(j)
	    read(rows(jj)(js+10:js+15),*,ERR=77) SCALE(j)
	    read(rows(jj)(js+17:js+22),*,ERR=77) SHIFT(j)
	    jn = max(jn,jj)
 30	continue
	if (jm .eq. 0 )	goto	32
	j1 = 2*jn+1
	do  31	j = 1,NB
	    if (IBOX(j) .gt. 0 )	goto	31
	    jw = j1
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(NAME(j),'(1A4)',ERR=77) rows(jj)(js+1:js+4)
	    read(rows(jj)(js+6: js+8), *,ERR=77) IB(j)
	    read(rows(jj)(js+10:js+15),*,ERR=77) SCALE(j)
	    read(rows(jj)(js+17:js+22),*,ERR=77) SHIFT(j)
	    jn = max(jn,jj)
	    j1 = j1+1
 31	continue
 32	continue
C	do	j = 1,NB
C	do	jj = 1,NB
C	    if (IB(j) .eq. IB(jj) .and. IB(j) .gt. 0 .and. j .ne. jj)
C     >			goto	78
C	enddo
C	enddo
 33	do	j = 1,NB
	    IBOX(j) = IB(j)
	enddo

	return
 77	write(*,*)
	write(*,'(2A)')'>>> INPUT ERROR encountered in ',
     &				'the dialog window "Presentation"'
	write(*,'(A23,A52,A1)')
     &	'                Line: "',rows(j)(1:52),'"'
	write(*,'(A52,$)')
     &	'     Enter "Y" to return, any other key to ignore > '
	KEY = 'X'
	read(*,'(:,A1)')KEY
	if (KEY .eq. 'Y' .or. KEY .eq. 'y')	goto	10
	return
C This warning being once ignored is repeated every next call
C 78	write(*,*)
C	write(*,'(A60,A7)')
C     &	">>> WARNING >>> You've put two (or more) curves in the same"
C     &	,' window'
C	write(*,'(A42)')'                 Information will be lost.'
C	write(*,'(A52,$)')
C     &	'     Enter "Y" to return, any other key to continue > '
C	KEY = 'X'
C	read(*,'(:,A1)')KEY
C	if (KEY .ne. 'Y' .or. KEY .ne. 'y')	goto	33
C	goto	10
	end
C======================================================================|
	subroutine   ASXWIN(NB,IBOX,NAME,SCALE,YSHIFT,XL,XR,MOD10,YMODE)
	implicit	none
	double precision	SCALE(*),YSHIFT(*),XL(*),XR(*),YY
	integer		NB,IBOX(*),YMODE,NRW16,NRW96,MOD10
	include		'for/parameter.inc'
	include		'for/const.inc'
	parameter	(NRW16=NRW+16, NRW96=NRW-128)
	integer		JMODE,JGR,j,j1,j2,jj,jn,jm,js,jb,jw,jsep,IB(NRW)
	integer		IP1(NRW),IP2(NRW),IP30(NRW),IP31(NRW),ASKTAB
	character*80	rows(NRW16),STR,TITLE
	character*1	NAME(*)*4,YKEY,TZER,ABNUM*6
	data	IP1/
     1		 1, 3, 5, 7,   9,11,13,15,   2, 4, 6, 8,  10,12,14,16,
     2		17,19,21,23,  25,27,29,31,  18,20,22,24,  26,28,30,32,
     3		33,35,37,39,  41,43,45,47,  34,36,38,40,  42,44,46,48,
     4		49,51,53,55,  57,59,61,63,  50,52,54,56,  58,60,62,64,
     5		65,67,69,71,  73,75,77,79,  66,68,70,72,  74,76,78,80,
     6		81,83,85,87,  89,91,93,95,  82,84,86,88,  90,92,94,96,
     7	 97,99,101,103, 105,107,109,111, 98,100,102,104,106,108,110,112,
     8	113,115,117,119,121,123,125,127,114,116,118,120,122,124,126,128/
C     9		NRW96*0/
	data	IP2/
     1		 1, 2, 9,10,   3, 4,11,12,   5, 6,13,14,   7, 8,15,16,
     2		17,18,25,26,  19,20,27,28,  21,22,29,30,  23,24,31,32,
     3		33,34,41,42,  35,36,43,44,  37,38,45,46,  39,40,47,48,
     4		49,50,57,58,  51,52,59,60,  53,54,61,62,  55,56,63,64,
     5		65,66,73,74,  67,68,75,76,  69,70,77,78,  71,72,79,80,
     6		81,82,89,90,  83,84,91,92,  85,86,93,94,  87,88,95,96,
     7	  97,98,105,106, 99,100,107,108,101,102,109,110,103,104,111,112,
     8	113,114,121,122,115,116,123,124,117,118,125,126,119,120,127,128/
C     9		NRW96*0/
	data	IP30/
     1		 1, 3, 5, 7,   2, 4, 6, 8,   9,11,13,15,  10,12,14,16,
     2		17,18,21,22,  25,26,29,30,  19,20,23,24,  27,28,31,32,
     3		33,34,37,38,  41,42,45,46,  35,36,39,40,  43,44,47,48,
     4		49,50,53,54,  57,58,61,62,  51,52,55,56,  59,60,63,64,
     5		65,66,69,70,  73,74,77,78,  67,68,71,72,  75,76,79,80,
     6		81,82,85,86,  89,90,93,94,  83,84,87,88,  91,92,95,96,
     7	  97,98,101,102,105,106,109,110, 99,100,103,104,107,108,111,112, 
     8	113,114,117,118,115,116,119,120,121,122,125,126,123,124,127,128/
C     9		NRW96*0/
	data	IP31/
     1		 1, 2, 5, 6,   3, 4, 7, 8,   9,10,13,14,  11,12,15,16,
     2		17,18,21,22,  19,20,23,24,  25,26,29,30,  27,28,31,32,
     3		33,34,37,38,  35,36,39,40,  41,42,45,46,  43,44,47,48,
     4		49,50,53,54,  51,52,55,56,  57,58,61,62,  59,60,63,64,
     5		65,66,69,70,  67,68,71,72,  73,74,77,78,  75,76,79,80,
     6		81,82,85,86,  83,84,87,88,  89,90,93,94,  91,92,95,96,
     7	  97,98,101,102,105,106,109,110, 99,100,103,104,107,108,111,112, 
     8	113,114,117,118,115,116,119,120,121,122,125,126,123,124,127,128/
C     9		NRW96*0/
C----------------------------------------------------------------------|
	j = 0
	TZER= char(j)
	TITLE = "Presentation"//TZER
C              ----5----0----5----0----5----0----5----0----5----0----5
	STR = "Name|Box| Scale|Offset| [a?,]| [,a?]||"//
     >	      "Name|Box| Scale|Offset| [a?,]| [,a?]"
     >		//TZER
	jsep = index(STR,'||')+1
C	STR = "    |   |      |      ||    |   |      |      "
C     >		//TZER
	JMODE = MOD10
	if (JMODE.eq.1)	then
	    JGR = 8
	elseif (JMODE.eq.2 .or. JMODE.eq.3)	then
	    JGR = 4
	elseif (JMODE.eq.6)	then
	    if (YMODE .eq. 1)	JGR = 4
	    if (YMODE .eq. 0)	JGR = 4
	    if (YMODE .eq.-1)	JGR = 2
	else
	    return
	endif
 10	continue
	do	j = 1,NRW
	   write(rows(j)(1:80),'(79X,1A1)')TZER
	enddo

C j  - ordinal box No.
C jb - box No. in the Astra nominations
C jw - position in the table
C jj - horizontal row
C js - position in the current row
C
	jn = 0
	jm = 0
	call	num2str(AB,ABNUM,6)
	do  20	j = 1,NB
	    jb = IBOX(j)
	    if (jb .le. 0 )	then
		jm = jm+1
		goto	20
	    endif
	    if (JMODE .eq. 1)	jw = IP1(jb)
	    if (JMODE .eq. 2 .or. JMODE .eq. 3)	jw = IP2(jb)
	    if (JMODE .eq. 6)	then
		if (YMODE .eq. 1)	jw = jb
		if (YMODE .eq. 0)	jw = IP30(jb)
		if (YMODE .eq.-1)	jw = IP31(jb)
	    endif
	    YY = XR(j)
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(rows(jj)(js+1:js+4),'(1A4)') NAME(j)
	    write(rows(jj)(js+6:js+8),'(1I3)') jb
	    call	num2str(SCALE(j),rows(jj)(js+10:js+15),6)
	    call	num2str(YSHIFT(j),rows(jj)(js+17:js+22),6)
	    call	num2str(XL(j),rows(jj)(js+24:js+29),6)
	    call	num2str(XR(j),rows(jj)(js+31:js+36),6)
	    if (YY .gt. AB) rows(jj)(js+31:js+36) = ABNUM
	    jn = max(jn,jj)
 20	continue
	if (jm .eq. 0 )	goto	22
	j2 = jn
	j1 = 2*jn+1
	do  21	j = 1,NB
	    if (IBOX(j) .gt. 0 )	goto	21
	    YY = XR(j)
	    jw = j1
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(rows(jj)(js+1:js+4),'(1A4)') NAME(j)
	    write(rows(jj)(js+6:js+8),'(1I3)') IBOX(j)
	    call	num2str(SCALE(j),rows(jj)(js+10:js+15),6)
	    call	num2str(YSHIFT(j),rows(jj)(js+17:js+22),6)
	    call	num2str(XL(j),rows(jj)(js+24:js+29),6)
	    call	num2str(XR(j),rows(jj)(js+31:js+36),6)
	    if (YY .gt. AB) rows(jj)(js+31:js+36) = ABNUM
	    jn = max(jn,jj)
	    j1 = j1+1
 21	continue
 22	continue
	j1 =(jm+1)/2
	j = ASKTAB(TITLE,STR,rows,80,jn,JGR,j1)
	if (j .gt. 0)	goto	22
C	call	ASKCOL(TITLE,STR,rows,80,jn,JGR,j1)
C
C Warning: If two (or more) box # coincide, information will be 
C	   partly overrode and lost
C
C	write(*,'(1(A22,3X,A22))')
C     >		(rows(j)(1:22),rows(j)(25:46),j=1,jn)
C	write(*,'(10I5)')j1,NB,jn

	jn = 0
	do  30	j  = 1,NB
	    if (IBOX(j) .le. 0 )	goto	30
	    if (JMODE .eq. 1)	jw = IP1(IBOX(j))
	    if (JMODE .eq. 2 .or. JMODE .eq. 3)	jw = IP2(IBOX(j))
	    if (JMODE .eq. 6)	then
		if (YMODE .eq. 1)	jw = IBOX(j)
		if (YMODE .eq. 0)	jw = IP30(IBOX(j))
		if (YMODE .eq.-1)	jw = IP31(IBOX(j))
	    endif
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(NAME(j),'(1A4)',ERR=77) rows(jj)(js+1:js+4)
	    read(rows(jj)(js+6: js+8), *,ERR=77) IB(j)
	    read(rows(jj)(js+10:js+15),*,ERR=77) SCALE(j)
	    read(rows(jj)(js+17:js+22),*,ERR=77) YSHIFT(j)
	    read(rows(jj)(js+24:js+29),*,ERR=77) XL(j)
	    read(rows(jj)(js+31:js+36),*,ERR=77) YY
	    if  (rows(jj)(js+31:js+36).ne.ABNUM)
     >	    read(rows(jj)(js+31:js+36),*,ERR=77) XR(j)
	    jn = max(jn,jj)
 30	continue
	if (jm .eq. 0 )	goto	32
	j1 = 2*jn+1
	do  31	j = 1,NB
	    if (IBOX(j) .gt. 0 )	goto	31
	    jw = j1
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(NAME(j),'(1A4)',ERR=77) rows(jj)(js+1:js+4)
	    read(rows(jj)(js+6: js+8), *,ERR=77) IB(j)
	    read(rows(jj)(js+10:js+15),*,ERR=77) SCALE(j)
	    read(rows(jj)(js+17:js+22),*,ERR=77) YSHIFT(j)
	    read(rows(jj)(js+24:js+29),*,ERR=77) XL(j)
	    if  (rows(jj)(js+31:js+36).ne.ABNUM)
     >	    read(rows(jj)(js+31:js+36),*,ERR=77) XR(j)
	    jn = max(jn,jj)
	    j1 = j1+1
 31	continue
 32	continue
C	do	j = 1,NB
C	do	jj = 1,NB
C	    if (IB(j) .eq. IB(jj) .and. IB(j) .gt. 0 .and. j .ne. jj)
C     >			goto	78
C	enddo
C	enddo
 33	do	j = 1,NB
	    IBOX(j) = IB(j)
	enddo

	return
 77	write(*,*)
	write(*,'(2A)')'>>> INPUT ERROR encountered in ',
     &				'the dialog window "Presentation"'
	write(*,'(A23,A52,A1)')
     &	'                Line: "',rows(j)(1:52),'"'
	write(*,'(A52,$)')
     &	'     Enter "Y" to return, any other key to ignore > '
	YKEY = 'X'
	read(*,'(:,A1)')YKEY
	if (YKEY .eq. 'Y' .or. YKEY .eq. 'y')	goto	10
	return
C This warning being once ignored is repeated every next call
C 78	write(*,*)
C	write(*,'(A60,A7)')
C     &	">>> WARNING >>> You've put two (or more) curves in the same"
C     &	,' window'
C	write(*,'(A42)')'                 Information will be lost.'
C	write(*,'(A52,$)')
C     &	'     Enter "Y" to return, any other key to continue > '
C	YKEY = 'X'
C	read(*,'(:,A1)')YKEY
C	if (YKEY .ne. 'Y' .or. YKEY .ne. 'y')	goto	33
C	goto	10
	end
C======================================================================|
	subroutine   ASKXGR(
     >		JCHAN,IBOX,NAME,YMODE,JNB,JXMODE,JGR,OUTFIG,OUTNAME
     >			   )
	implicit	none
	integer		JCHAN,JNB,IBOX(*),YMODE,ASKGRF,JMODE
	integer		NRW16,NRW96
	include		'for/parameter.inc'
	include		'for/const.inc'		! XOUT used
	include		'for/outcmn.inc'	! MOD10,TASK(1:3) used
	parameter	(NRW16=NRW+16, NRW96=NRW-128)
	integer		JXMODE,IFKEY,JDONE
	integer		JGR,j,j1,jj,jn,jm,js,jb,jc,jw,jsep,OUTFIG(*)
	integer		IP1(NRW),IP2(NRW),IP30(NRW),IP31(NRW),IB(NRW)
	character*80	rows(NRW16),STR,TITLE,OUTNAME(*)*8
	character*1	NAME(*)*4,YKEY,TZER
	data	IP1/
     1		 1, 3, 5, 7,   9,11,13,15,   2, 4, 6, 8,  10,12,14,16,
     2		17,19,21,23,  25,27,29,31,  18,20,22,24,  26,28,30,32,
     3		33,35,37,39,  41,43,45,47,  34,36,38,40,  42,44,46,48,
     4		49,51,53,55,  57,59,61,63,  50,52,54,56,  58,60,62,64,
     5		65,67,69,71,  73,75,77,79,  66,68,70,72,  74,76,78,80,
     6		81,83,85,87,  89,91,93,95,  82,84,86,88,  90,92,94,96,
     7	 97,99,101,103, 105,107,109,111, 98,100,102,104,106,108,110,112,
     8	113,115,117,119,121,123,125,127,114,116,118,120,122,124,126,128/
C     9		NRW96*0/
	data	IP2/
     1		 1, 2, 9,10,   3, 4,11,12,   5, 6,13,14,   7, 8,15,16,
     2		17,18,25,26,  19,20,27,28,  21,22,29,30,  23,24,31,32,
     3		33,34,41,42,  35,36,43,44,  37,38,45,46,  39,40,47,48,
     4		49,50,57,58,  51,52,59,60,  53,54,61,62,  55,56,63,64,
     5		65,66,73,74,  67,68,75,76,  69,70,77,78,  71,72,79,80,
     6		81,82,89,90,  83,84,91,92,  85,86,93,94,  87,88,95,96,
     7	  97,98,105,106, 99,100,107,108,101,102,109,110,103,104,111,112,
     8	113,114,121,122,115,116,123,124,117,118,125,126,119,120,127,128/
C     9		NRW96*0/
	data	IP30/
     1		 1, 3, 5, 7,   2, 4, 6, 8,   9,11,13,15,  10,12,14,16,
     2		17,18,21,22,  25,26,29,30,  19,20,23,24,  27,28,31,32,
     3		33,34,37,38,  41,42,45,46,  35,36,39,40,  43,44,47,48,
     4		49,50,53,54,  57,58,61,62,  51,52,55,56,  59,60,63,64,
     5		65,66,69,70,  73,74,77,78,  67,68,71,72,  75,76,79,80,
     6		81,82,85,86,  89,90,93,94,  83,84,87,88,  91,92,95,96,
     7	  97,98,101,102,105,106,109,110, 99,100,103,104,107,108,111,112, 
     8	113,114,117,118,115,116,119,120,121,122,125,126,123,124,127,128/
C     9		NRW96*0/
	data	IP31/
     1		 1, 2, 5, 6,   3, 4, 7, 8,   9,10,13,14,  11,12,15,16,
     2		17,18,21,22,  19,20,23,24,  25,26,29,30,  27,28,31,32,
     3		33,34,37,38,  35,36,39,40,  41,42,45,46,  43,44,47,48,
     4		49,50,53,54,  51,52,55,56,  57,58,61,62,  59,60,63,64,
     5		65,66,69,70,  67,68,71,72,  73,74,77,78,  75,76,79,80,
     6		81,82,85,86,  83,84,87,88,  89,90,93,94,  91,92,95,96,
     7	  97,98,101,102,105,106,109,110, 99,100,103,104,107,108,111,112, 
     8	113,114,117,118,115,116,119,120,121,122,125,126,123,124,127,128/
C     9		NRW96*0/
	j = 0
	TZER  = char(j)
	jdone = 1
	TITLE = "Select profiles for plotting"//TZER
C              ----5----0----5----0----5----0----5----0----5----0----5
	STR = "Box #|  Name  |Marker||"//
     >	      "Box #|  Name  |Marker"  //TZER
	jsep  = index(STR,'||')+1
	JMODE = MOD10
	JXMODE = XOUT+.49
	if (JMODE.eq.1)	then
	   JGR = 8
	elseif (JMODE.eq.2 .or. JMODE.eq.3)	then
	   JGR = 4
	elseif (JMODE.eq.6)	then
	   if (YMODE .eq. 1)	JGR = 4
	   if (YMODE .eq. 0)	JGR = 4
	   if (YMODE .eq.-1)	JGR = 2
	   JXMODE = -1
	else
	   return
	endif
 10	continue
	do	j = 1,NRW
	   OUTNAME(j) = '        '
	   write(rows(j)(1:80),'(79X,1A1)')TZER
	enddo
	JNB = JCHAN
	if (JMODE .ne. 1)	JNB = min(jchan,96)

C j  - ordinal box No.
C jb - box No. in the Astra nominations
C jm - number of empty boxes
C jw - position in the table
C jj - horizontal row
C js - position in the current row
	jn = 0
	jm = 0
	do  20	j = 1,JNB
	    jb = IBOX(j)
	    if (jb .le. 0 )	then
		jm = jm+1
		goto	20
	    endif
	    if (JMODE .eq. 1)	jw = IP1(jb)
	    if (JMODE .eq. 2 .or. JMODE .eq. 3)	jw = IP2(jb)
	    if (JMODE .eq. 6)	then
		if (YMODE .eq. 1)	jw = jb
		if (YMODE .eq. 0)	jw = IP30(jb)
		if (YMODE .eq.-1)	jw = IP31(jb)
	    endif
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(rows(jj)(js+1:js+5), '(1I4,1X)') jb
	    write(rows(jj)(js+7:js+14),'(2X,1A4,2X)') NAME(j)
	    write(rows(jj)(js+16:js+21),'(1I4,2X)') 0
C	    call	num2str(SCALE(j),rows(jj)(js+10:js+15),6)
	    jn = max(jn,jj)
 20	continue
	if (jm .eq. 0 )	goto	22
	jw = 2*jn+1
	do  21	j = 1,JNB
	    if (IBOX(j) .gt. 0 )	goto	21
	    js = jsep*(1-jw+jw/2*2)
	    jj = 1+(jw-1)/2
	    write(rows(jj)(js+1:js+5),'(1I4,1X)') IBOX(j)
	    write(rows(jj)(js+7:js+14),'(2X,1A4,2X)') NAME(j)
	    write(rows(jj)(js+16:js+21),'(1I4,2X)') 0
C	    call	num2str(SCALE(j),rows(jj)(js+10:js+15),6)
	    jn = max(jn,jj)
	    jw = jw+1
 21	continue
 22	continue
	j1 =(jm+1)/2		! jm number of switched off windows
	j = ASKGRF(TITLE,STR,rows,80,jn,JGR,j1,JXMODE)
	if (j .gt. 0)	goto	22
C	write(*,'(5(A))')('"',rows(j)(1:21),'","'
C     >			,rows(j)(24:44),'"',j=1,jn)
C	write(*,'(10I5)')j1,JNB,jn,jsep
C	write(*,*)(IBOX(j),j=1,JNB)
	if (j1 .eq. -1)	then
C	   write(*,*)"<Alt><Esc> pressed"
C	   jdone = 0
C	   if ( TASK(1:3) .eq. 'RUN' )	j1 = IFKEY(32)	! Set waiting mode
	endif

	jn = 0
	do  25	j  = 1,JNB
	   if (IBOX(j) .le. 0 )	goto	25
	   if (JMODE .eq. 1)	jw = IP1(IBOX(j))
	   if (JMODE .eq. 2 .or. JMODE .eq. 3)	jw = IP2(IBOX(j))
	   if (JMODE .eq. 6)	then
	      if (YMODE .eq. 1)	jw = IBOX(j)
	      if (YMODE .eq. 0)	jw = IP30(IBOX(j))
	      if (YMODE .eq.-1)	jw = IP31(IBOX(j))
	   endif
	   js = jsep*(1-jw+jw/2*2)
	   jj = 1+(jw-1)/2
	   write(OUTNAME(j),'(1A8)',ERR=77) rows(jj)(js+7:js+14)
	   if (OUTNAME(j) .ne. '        ')	then
 23	      if (OUTNAME(j)(1:1).ne.' ')	goto	24
	      OUTNAME(j)(1:) = OUTNAME(j)(2:)//'       '
	      goto	23
	   endif
 24	   j1 = index(rows(jj)(js+16:js+21),'.')
	   if (j1 .ge. 2)	then
	      read(rows(jj)(js+16:js+14+j1),*,ERR=77)OUTFIG(j)
	      read(rows(jj)(js+16+j1:js+21),*,ERR=77)jc
	   elseif (j1 .eq. 0)	then
	      read(rows(jj)(js+16:js+21),*,ERR=77)OUTFIG(j)
	      jc = 0
	   else
	      goto	77
	   endif
	   if (OUTFIG(j) .gt. 999)	goto	78
	   OUTFIG(j) = 1000*OUTFIG(j)+jc
	   jn = max(jn,jj)
 25	continue

	if (jm .eq. 0 )	goto	29
	jw = 2*jn+1
	do  28	j = 1,JNB
	   if (IBOX(j) .gt. 0 )	goto	28
	   js = jsep*(1-jw+jw/2*2)
	   jj = 1+(jw-1)/2
	   write(OUTNAME(j),'(1A8)',ERR=77) rows(jj)(js+7:js+14)
	   if (OUTNAME(j) .ne. '        ')	then
 26	      if (OUTNAME(j)(1:1).ne.' ')	goto	27
	      OUTNAME(j)(1:) = OUTNAME(j)(2:)//'       '
	      goto	26
	   endif
 27	   j1 = index(rows(jj)(js+16:js+21),'.')
	   if (j1 .ge. 2)	then
	      read(rows(jj)(js+16:js+14+j1),*,ERR=77)OUTFIG(j)
	      read(rows(jj)(js+16+j1:js+21),*,ERR=77)jc
	   elseif (j1 .eq. 0)	then
	      read(rows(jj)(js+16:js+21),*,ERR=77)OUTFIG(j)
	      jc = 0
	   else
	      goto	77
	   endif
	   if (OUTFIG(j) .gt. 999)	goto	78
	   OUTFIG(j) = 1000*OUTFIG(j)+jc
	   jn = max(jn,jj)
	   jw = jw+1
 28	continue
 29	continue
	jgr = 0
	do	j=1,JNB
	   jb = OUTFIG(j)/1000
	   jgr = max(jgr,jb)
	enddo
	if (jgr .eq. 0)	return

	if (jxmode.lt.-1 .or. jxmode.gt.5)	then
	   write(*,*)"Unknown data type"
	   return
	endif
	return

 77	write(*,*)
	write(*,'(2A)')'>>> INPUT ERROR encountered in ',
     &				'the dialog window "Output control"'
	write(*,'(A23,A52,A1)')
     &	'                Line: "',rows(j)(1:52),'"'
	write(*,'(A52,$)')
     &	'     Enter "Y" to return, any other key to ignore > '
	YKEY = 'X'
	read(*,'(:,A1)')YKEY
	if (YKEY .eq. 'Y' .or. YKEY .eq. 'y')	goto	10
	return
 78	write(*,*)
	write(*,'(2A)')'>>> INPUT ERROR encountered in ',
     &				'the dialog window "Output control"'
	write(*,'(A)')'    Fig.# cannot exceed 999'
	return
C This warning being once ignored is repeated every next call
C 79	write(*,*)
C	write(*,'(A60,A7)')
C     &	">>> WARNING >>> You've put two (or more) curves in the same"
C     &	,' window'
C	write(*,'(A42)')'                 Information will be lost.'
C	write(*,'(A52,$)')
C     &	'     Enter "Y" to return, any other key to continue > '
C	YKEY = 'X'
C	read(*,'(:,A1)')YKEY
C	if (YKEY .ne. 'Y' .or. YKEY .ne. 'y')	goto	33
C	goto	10
	end
C======================================================================|
	subroutine	WRFIGS(JNB,JXMODE,JGR,OUTFIG,OUTNAME,IBOX,
     >			ITIMES,TTOUT,TOUT)
C----------------------------------------------------------------------|
C The subroutine writes file in the directories AWD/out/ and AWD/xmg/ 
C     with curves selected in the dialog window ASKXGR (hot key "O")
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'		! NA1,NAB,ABC,ROC,VOLUME
	include 'for/status.inc'	! AMETR,RHO,FP,VOLUM
	include 'for/outcmn.inc'
	integer	JNB,JXMODE,JGR,OUTFIG(*),IBOX(*),IERR,JLR,JN0
	integer	ITIMES,j,jj,jr,jn,jm,js,jb,jc,jw,jl,NP1,MODEX,JNUM(NRW)
	integer		killbl,length,WarningColor,j0,j1,j2,j3,j4,j5
	character	OUTNAME(*)*8,FNAME*132,CVE*1,STRI*118
	logical		EXI
	double precision	YY,TTOUT(ITIMES),TOUT(ITIMES,NRW)
	save		JN0,JLR,WarningColor
	data		JN0/0/ JLR/426/ WarningColor/30/
C	write(*,*)jgr   
C	write(*,*)"Abscissa type",jxmode
C	write(*,*)(OUTFIG(j),j=1,JNB)
C Create output file name and open logical channel 7:
	j = length(RDNAME)
	j1 = length(EQNAME)
	j2 = length(AWD)
	if (j+j1+j2 .gt. 77)	then
	   write(*,*)">>> Warning the file name is too long. ",
     >		"The data may be lost"
	endif
	jr = index(RUNID,RDNAME(1:j))+j

C	call	CHECK_XMG_DIR
	write(*,*)
	FNAME=AWD(1:j2)//'xmg/'//char(0)
	inquire(file=FNAME,exist=EXI)
	if (.not.EXI)	then
	   write(STRI,'(2A)')'mkdir ',FNAME(1:j2+4)
	   call	system(STRI)
	endif
	FNAME=AWD(1:j2)//'xmg/'//EQNAME(1:j1)//'.'//RDNAME(1:j)//char(0)
	inquire(file=FNAME,exist=EXI)
	if (.not.EXI)	then
	   write(STRI,'(2A)')'mkdir ',FNAME(1:length(FNAME))
	   call	system(STRI)
	endif

	FNAME=AWD(1:j2)//'out/'//RDNAME(1:j)//'.'//EQNAME(1:j1)//char(0)
	jj = killbl(FNAME,132)
	call	SETFNA(FNAME,jj)
	call colovm(WarningColor)
	STRI='>>> Dataset is written in the file: '//FNAME(1:jj)//' '
	write(*,'(/,1X,A)')STRI(1:jj+40)
	JLR = XWH-124
	call textvm(JN0,JLR,STRI,jj+40)
	call	OPENWT(7,FNAME(1:jj),0,IERR)
	if(IERR.gt.0)	pause 'Data file error'
	FNAME=AWD(1:j2)//'xmg/'//EQNAME(1:j1)//'.'//RDNAME(1:j)
     +			//'/Fig'//char(0)
	j2 = killbl(FNAME,132)
C	write(*,*)FNAME(1:length(FNAME)),j2

	if     (jxmode .eq. 0)	then
	   write(STRI(1:10),'(A)')"  a, m    "
	   NP1 = NAB
	elseif (jxmode .eq. 1)	then
	   write(STRI(1:10),'(A)')"  a_N, d/l"
	   NP1 = NA1
	elseif (jxmode .eq. 2)	then
	   write(STRI(1:10),'(A)')"rho_N, d/l"
	   NP1 = NA1
	elseif (jxmode .eq. 3)	then
	   write(STRI(1:10),'(A)')" Psi, Vs  "
	   NP1 = NA1
	elseif (jxmode .eq. 4)	then
	   write(STRI(1:10),'(A)')" rho_V, m "
	   NP1 = NA1
	elseif (jxmode .eq. 5)	then
	   write(STRI(1:10),'(A)')"rho_pol, m"
	   NP1 = NA1
	elseif (jxmode .eq. -1)	then
	   write(STRI(1:10),'(A)')"  time, s "
	   NP1 = LTOUT-1
	endif

C	goto	11
C The next block selects all the curves with the same number
C and submits them to a Figure with this number.
C The figures are put in ascending order. The curve numbers are ignored.
	jl = 0					! jl - total amount of curves
	j3 = 1
	do	jj=1,jgr			! 999 - max Fig #
	   jc = 0				! Amount of curves in one Fig
	   do	10	j1=1,JNB			! NRW - max channel #
	      jb = OUTFIG(j1)/1000		! Fig #
	      if (jb .eq. 0)	goto	10
c	      write(*,*)j1,jj,jb
	      if (jb.ne.jj)	goto	10
	      if (jc .eq. 0)	then
		 FNAME = FNAME(1:j2)//char(0)	! FNAME = .../Fig
		 j5 = 0				! Check existance
 1		 j5 = j5+1
		 write(FNAME(j2+1:),'(1A1,1I3)')'.',j5
		 FNAME = FNAME(1:j2+4)//char(0)
		 jw = killbl(FNAME,j2+4)
		 inquire(FILE=FNAME(1:jw)//'.dat',EXIST=EXI)
		 if(EXI)	goto 1
		 write(*,'(/3A)')' >>> New Figure: "',FNAME(1:jw),'"'
		 call	OPENWT(8,FNAME(1:jw)//'.dat',0,IERR)
		 if (IERR .gt. 0)	goto	97
		 call	OPENWT(9,FNAME(1:jw)//'.par',0,IERR)
		 if (IERR .gt. 0)	goto	98
		 if (index(AWD(1:length(AWD)),'efda-itm') .eq. 0)
     >		 write(9,107)
		 if (TASK .eq. 'VIEW')	then
		    write(9,104)RSNAME(1:length(RSNAME))
		    write(9,105)
		 endif
		 write(9,106)FNAME(j2-2:jw),RUNID(1:jr),'"'	! Run & Fig ID
C		 write(*,'(1A80)')RUNID
C		 write(9,105)'Model file ',EQNAME(1:length(EQNAME))
C     >			,',    Data file ',RDNAME(1:length(RDNAME)),'"'
		 if (jxmode .eq. 0) write(9,101)"a [m]"
		 if (jxmode .eq. 1) write(9,101)"a_N"
		 if (jxmode .eq. 2) write(9,101)"rho_N"
		 if (jxmode .eq. 3) write(9,101)"Psi [Vs]"
		 if (jxmode .eq. 4) write(9,101)"rho_V [m]"
		 if (jxmode .eq. 5) write(9,101)"rho_pol [m]"
		 if (jxmode .eq. -1)write(9,101)"time [s]"
		 write(7,'(2A)')'Abscissa: ',STRI(1:10)
	      endif
	      STRI(11+10*jc:)=OUTNAME(j1)//'  '
	      if (jc .lt. 10) write(9,102)jc,OUTNAME(j1)
	      if (jc .ge. 10) write(9,103)jc,OUTNAME(j1)
	      jc = jc+1
	      CVE = char(96+jc)
	      write(*,'(3A,I3,A,I3,4A,I3))')
     >				 'Name  "',OUTNAME(j1)
     >				,'"     Box',IBOX(j1)
     >				,'"     Fig',jj,'(',CVE,')'
     >				,',     Curve',jc
	      write(7,'(2A,3(A,I3))')
     >				 'Name  ',OUTNAME(j1)
     >				,'	Box',IBOX(j1)
     >				,'	Fig',jj
     >				,'	Curve',jc
	      jl = jl+1
	      JNUM(jl) = IBOX(j1)
 10	   continue
	   if (jc .ne. 0)      then
	      j4 = jc+1
	      do	j0=1,NP1
		 if     (jxmode .eq. 0)	then
		    YY = AMETR(j0)
		 elseif (jxmode .eq. 1)	then
		    YY = AMETR(j0)/ABC
		 elseif (jxmode .eq. 2)	then
		    YY = RHO(j0)/ROC
		 elseif (jxmode .eq. 3)	then
		    YY = FP(j0)
		 elseif (jxmode .eq. 4)	then
		    YY = sqrt(VOLUM(j0)/VOLUME)
		 elseif (jxmode .eq. 5)	then
		    YY = sqrt((FP(j0)-FP(1))/(FP(NA1)-FP(1)))
		 elseif (jxmode .eq. -1)	then
		    YY = TTOUT(j0)
		 endif
		 j3 = jl-jc+1
		 if     (jxmode .ge. 0)	then
		    if (j0.eq.1)	then
C	write(*,'(3(I3))')(j,jnum(j),ibox(jnum(j)),j=jl-jc+1,jl)
		    endif
C		    write(8,'(1P,<j4>E12.4)')
		    write(8,'(1P,11E12.4)')
     >			YY,(ROUT(j0,ibox(jnum(j))),j=jl-jc+1,jl)
	 	 else
C		    write(8,'(1P,<j4>E12.4)')
		    write(8,'(1P,11E12.4)')
     >			YY,(TOUT(j0,ibox(jnum(j))),j=jl-jc+1,jl)
	         endif
	      enddo
	      close(8)
	      close(9)
	   endif
	enddo
 11	continue
C	return		! If enabled suppress writing file to out/ (unit 7)
	goto	13
C This block submits the curves to a Figure according to their numbering.
C Both the figures and curves within a figure are put in ascending order
C in exact correspondence with the pre-defined numbers.
C The curves with repeated numbers are ignored.
	j1 = 0
	jl = 0
	do	jj=1,jgr			! 999 - max Fig #
	   jw = -1
	   do	js=0,19				! 19  - max curve #
	      do	j=1,JNB	 		! NRW - max channel #
		 jb = OUTFIG(j)/1000 		! Fig #
		 if (jb.eq.jj .and. jj.ne.j1)	j1 = jj
		 CVE = '`'
		 do jm=1,JNB
		    jn = OUTFIG(jm)/1000		! Retrieve Fig #
		    if (jn .ne. jj)	goto	12
		    CVE = char(ichar(CVE)+1)
		    jc = OUTFIG(jm)-jn*1000 		! Curve #
		    if (js .ne. jc)	goto	12	! Ignore curve #
		    if (js .eq. jw)	goto	12	! Skip repeated curves
		    CVE = char(97+jc)
		    write(*,*)
     >				'Name  "',OUTNAME(jm)
     >				,'",    Box',IBOX(jm)
     >				,' ,    Curve',jc+1
     >				,'      Fig',jj,'(',CVE,')'
		    write(7,'(2A,3(A,I3))')
     >				 'Name  ',OUTNAME(jm)
     >				,'	Box',IBOX(jm)
     >				,'	Fig',jj
     >				,'	Curve',jc+1
		    jw = js
		    jl = jl+1
		    JNUM(jl) = IBOX(jm)
 12		    continue
		 enddo					! jm
	      enddo					! j
	   enddo					! js
	enddo						! jj
 13	continue

C	jw = min(5,jl)
	jj = 0
	js = 1
 14	continue
	do	j1=1,NP1
	   if     (jxmode .eq. 0)	then
	      YY = AMETR(j1)
	   elseif (jxmode .eq. 1)	then
	      YY = AMETR(j1)/ABC
	   elseif (jxmode .eq. 2)	then
	      YY = RHO(j1)/ROC
	   elseif (jxmode .eq. 3)	then
	      YY = FP(j1)
	   elseif (jxmode .eq. 4)	then
	      YY = sqrt(VOLUM(j1)/VOLUME)
	   elseif (jxmode .eq. 5)	then
	      YY = sqrt((FP(j1)-FP(1))/(FP(NA1)-FP(1)))
	   elseif (jxmode .eq. -1)	then
	      YY = TTOUT(j1)
	   endif
C	   if     (js .eq. 1)	then
	      if     (jxmode .ge. 0)	then
		 write(7,100)YY,(ROUT(j1,ibox(jnum(j))),j=js,jl)
	      else
		 write(7,100)YY,(TOUT(j1,ibox(jnum(j))),j=js,jl)
	      endif
C	   else
C	      if     (jxmode .ge. 0)	then
C		 write(7,100)(ROUT(j1,ibox(jnum(j))),j=js,jl)
C	      else
C		 write(7,100)(TOUT(j1,ibox(jnum(j))),j=js,jl)
C	      endif
C	   endif
	enddo
C	js = jw+1
C	write(7,*)"Endblock"
C	if (jw .eq. jl)	goto	15
C	jw = min(js+5,jl)
C	goto	14
 15	close(7)
	return
 100	format(1P,6(E11.3))
 101	format('    xaxis  label "',A,'"')
 102	format('    s',I1,' legend  "',A,'"')
 103	format('    s',I2,' legend  "',A,'"')
 104	format('    title "View file:   ',A,'"')
 105	format('    title size 1.000000')
 106	format('    subtitle "',3A)
 107	format('page size 421, 298'
C     .	/'with string'/'    string on'
C     .	/'string loctype view'/'string 0.35, 0.5'/'string color 14'
C     .	/'string rot 45'/'string font 0'/'string just 0'
C     .	/'string char size 1.0'/'string def "Any text"'
     .	)
 97	write(*,*)'Xmgrace file "',FNAME(1:j1),'.dat" creation error'
	return
 98	write(*,*)'Xmgrace file "',FNAME(1:j1),'.par" creation error'
	return
	end
C======================================================================|
C	The function returns the length of the string STRI
C	from the beginning till the first space, tab, '\0' or <CR>
C----------------------------------------------------------------------|
	integer	function length(STRI)
	implicit none
	character*1	STRI(*)
	character*1	TZER,TAB,TCR
	integer	j
	j = 0
	TZER=char(j)
	j = 9
	TAB=char(j)
	j = 13
	TCR=char(j)
	length	= 0
	do	j=1,132
C	write(*,*)'"',STRI(j),'"',j,ichar(STRI(j))
		if ( STRI(j).eq.' '
     &		.or. STRI(j).eq.TAB
     &		.or. STRI(j).eq.TZER
     &		.or. STRI(j).eq.TCR) return
	length	= j
	enddo
	end
C======================================================================|
C	The function returns the length of the string STRI
C	from the beginning till the first '\0'
C----------------------------------------------------------------------|
	integer	function lengt0(STRI)
	implicit none
	character*1	STRI(*),TZER
	integer	j
	j = 0
	TZER=char(j)
	lengt0	= 0
	do	j=1,132
C	write(*,*)'"',STRI(j),'"',j,ichar(STRI(j))
	   if ( STRI(j).eq.TZER )  return
           lengt0 = j
	enddo
	end
C======================================================================|
C	The function returns the length of the string STRI from 
C	the beginning till the last nonspace, followed by'\0' or by <CR>
C-----------------------------------------------------------------------
	integer	function lonlen(STRI)
	implicit none
	character*(*) STRI
	character*1   TZER,TAB,TCR
	integer	j
	j = 0
	TZER=char(j)
	j = 9
	TAB=char(j)
	j = 13
	TCR=char(j)
	do	j=1,len(STRI)
	   if (STRI(j:j).eq.TZER .or. STRI(j:j).eq.TCR)	return
	   if (STRI(j:j).ne.' ' .and. STRI(j:j).ne.TAB)	lonlen = j
	enddo
	end
C=======================================================================
C	The function returns the position of the 1st non-blanck (space 
C		or tabulation) in the string STRI
C		from the beginning till '\0' or <CR>
C	or  0	if no non-blancks found
C-----------------------------------------------------------------------
	integer	function noblan(STRI)
	implicit none
	character STRI*(*)
	character TZER,TCR,TAB
	integer	j
	j = 0
	TZER=char(j)
	j = 9
	TAB=char(j)
	j = 13
	TCR=char(j)
	noblan	= 0
	do	j=1,len(STRI)
	    if (STRI(j:j).eq.TZER .or. STRI(j:j).eq.TCR)	return
	    if (STRI(j:j).ne.' ' .and. STRI(j:j).ne.TAB)	then
		noblan = j
		return
	    endif
	enddo
	end
C=======================================================================
C   The function returns the position of the 1st character 
C		of 2nd word in the string STRI
C	or  0	if  '\0' or <CR> is encountered before
C-----------------------------------------------------------------------
	integer	function nextwd(STRI)
	implicit none
	character STRI*(*)
	integer	j1,j2,noblan,length
	nextwd = 0
C 1st nonblanck
	j1 = noblan(STRI)
	if (j1 .eq. 0)	return
C Position after the end of the 1st word
	j2 = length(STRI(j1:))+j1
C Position of the 2nd word
	j1 = noblan(STRI(j2:))
	if (j1 .ne. 0)	nextwd	= j2+j1-1
	end
C=======================================================================
C   The function returns the position of SYM in STRI(1:LENG) 
C	or  0  if not found
C	or '\0' or <CR> is encountered before
C				! similar to the intrinsic INDEX
C-----------------------------------------------------------------------
	integer function nsymb(STRI,SYM)
	implicit none
	character*1 STRI(*),SYM,TZER,TCR
	integer j
	j = 0
	TZER = char(j)
	j = 13
	TCR = char(j)
	nsymb = 0
	do	j=1,132
	    if (STRI(j).ne.SYM)	goto 1
	    nsymb = j
		return
 1	    if (STRI(j).eq.TZER .or. STRI(j).eq.TCR)	goto	2
	enddo
 2	continue
	nsymb = 0
	end
C=======================================================================
C   The subroutine re-arranges real array A(1:N)
C-----------------------------------------------------------------------
	subroutine	SORT(A,N)
	implicit none
	integer	j,jj,jl,jr,N
	double precision	A(N),Y
	jl = 1
 1	jr = jl
	Y = A(jl)
	do	jj=jl+1,N
	    if (A(jj) .lt. Y)	then
		Y = A(jj)
		jr = jj
	    endif
	enddo
	if (jr .ne. jl)	then
	    do	j=jr-1,jl,-1
		A(j+1) = A(j)
	    enddo
	    A(jl) = Y
	endif
	jl = jl+1
	if (jl .eq. N)	return
	goto	1
	end
C=======================================================================
C   The subroutine re-arranges real arrays A(1:N) and B(1:N)
C-----------------------------------------------------------------------
	subroutine	SORTAB(A,B,N)
	implicit none
	integer	j,jj,jl,jr,N
	double precision	A(N),B(N),YA,YB
	if (N .le. 1)	return
	jl = 1
 1	jr = jl
	YA = A(jl)
	YB = B(jl)
	do	jj=jl+1,N
	    if (A(jj) .lt. YA)	then
		YA = A(jj)
		YB = B(jj)
		jr = jj
	    endif
	enddo
	if (jr .ne. jl)	then
	    do	j=jr-1,jl,-1
		A(j+1) = A(j)
		B(j+1) = B(j)
	    enddo
	    A(jl) = YA
	    B(jl) = YB
	endif
	jl = jl+1
	if (jl .eq. N)	return
	goto	1
	end
C=======================================================================
C	The function removes blancks and tabs from the string STRI
C	and returns a length of the resultant string
C	The rest of the string is appended with char(0)
C-----------------------------------------------------------------------
	integer	function KILLBL(STRI,leng)
	implicit none
	integer     j,jj,leng,I4
	character*1 STRI(leng),TZER,TCR,TAB
	I4 = 0
	TZER=char(I4)
	I4 = 13
	TCR=char(I4)
	I4 = 9
	TAB=char(I4)
	jj = 1
	do	j=1,leng
		if (STRI(j).eq.TZER.or.STRI(j).eq.TCR)	goto  1
		if (STRI(j) .ne. ' ' .and. STRI(j) .ne. TAB)  then
			STRI(jj) = STRI(j)
			KILLBL = jj
			jj = jj+1
		endif
	enddo
 1	continue
	do	j=0,leng-jj
	   write(STRI(jj+j),'(1A1)')TZER
	enddo
C	write(*,*)'Length returned: ',KILLBL
	end
C=======================================================================
	subroutine	curvvm(id,npnts,array)
C The same as drcurv but without the factor 10 
C The chain: PLOTGR(obsolete) -> CURV(obsolete) -> CURVVM
C		 -> drawvm(PSADrawLine) is not used any more
C dimension array(2*npnts)
	implicit none
	integer	npnts,array(*),j,id
	if (npnts .eq. 1)
     >	   call	drawvm(id,array(1),array(2),array(1),array(2))
	if (npnts .le. 1) return

	call drawline(id,array,npnts)
C	do	j=1,2*npnts-3,2
C	   call drawvm(id,array(j),array(j+1),array(j+2),array(j+3))
C	enddo
	end
C=======================================================================
	subroutine	drcurv(id,npnts,array)
C The same as curvvm but the supplied integer array is multiplied 
C     by the factor 10 in order to enhance PS resolution
C     This factor is then removed in C function d1line
C       PLOTCR -> CURV1 -> DRCURV -> d1line
C dimension array(2*npnts)
	implicit none
	integer	npnts,array(*),j,id
	if (npnts .eq. 1)
     >	   call	d1line(id,array(1),array(2),array(1),array(2))
	if (npnts .le. 1) return

	call d1polyline(id,array,npnts)

C	do	j=1,2*npnts-3,2
C	   call d1line(id,array(j),array(j+1),array(j+2),array(j+3))
C	enddo
	end
C======================================================================
	subroutine	prntxt(string)
	implicit none
	character*1	TZER,TCR,ch,string(*)
	integer         I4,j,jj
	I4 = 0
	TZER=char(I4)
	I4 = 13
	TCR=char(I4)
	jj = 0
	do j=1,80
	   ch	=string(j)
	   if(ch.ne.TZER.and.ch.ne.' '.and.ch.ne.TCR)  jj = j
	enddo
	if (jj.ne.0) write(*,'(1X,80A1)') (string(j),j=1,jj)
	end
C=======================================================================
	character*6 function VARNAM(string,ierr)
C-----------------------------------------------------------------------
C The subroutine analizes a "string"
C	If the 1st position is tab or space the string 6*' ' is returned
C	If tabs are encountered on the end of the "string",
C		they are removed the "string" is appended with spaces
C		and VARNAM in the Astra standard is created,
C		in this case,	ierr=1 is returned
C	Otherwise, ierr=0 and VARNAM in the Astra standard are returned.
C	The last "X" if appears is replaced with " "
C-----------------------------------------------------------------------
	implicit none
	character*(*)	string
	character*1	TCR,TAB,symb,name*6
	integer         INT4,j,ierr,je
	INT4 = 13
	TCR=char(INT4)
	INT4 = 9
	TAB=char(INT4)
	INT4 = 0
	ierr = 0
	name = '      '
	je = 7
	do	j=1,6
	    symb = string(j:j)
	    if (INT4.eq.0 .and. (symb.eq.' ' .or. symb.eq.TAB))	then
C 		Ignore names starting with spaces and tabs
		VARNAM = '      '
		return
	    else
		INT4 = j
	    endif
	    if(symb.eq.TAB .or. symb.eq.' ' .or. symb.eq.TCR)	then
		je = j
		if(symb .ne. ' ')	ierr = 1
		goto	1
	    else
		name(j:j) = symb
	    endif
	enddo
 1	if (je.le.6)	then
	    do	j=je,6
		name(j:j) = ' '
	    enddo
	endif
	do	j=2,5
	    if (name(j:j+1) .eq. 'X ')	name(j:j+1) = '  '
	enddo
	if (name(6:6) .eq. 'X')	name(6:6) = ' '
	VARNAM = name
	end
C=======================================================================
	character*6 function ARRNAM(string)
C-----------------------------------------------------------------------
C The subroutine analizes a character*6 "string"
C	If the 1st position is tab or space the string 6*' ' is returned
C	If tabs are encountered on the end of the "string",
C		they are removed the "string" is appended with spaces
C	Trailing "X" is added when not present in string*6 
C	Finally ARRNAM in the Astra standard is created,
C-----------------------------------------------------------------------
	implicit none
	character*(*)	string
	character*1	TCR,TAB,symb,name*6
	integer         jj,j,je
	jj = 13
	TCR=char(jj)
	jj = 9
	TAB=char(jj)
	jj = 0
	name = '      '
	je = 7
	do	j=1,6
	    symb = string(j:j)
	    if (jj.eq.0 .and. (symb.eq.' ' .or. symb.eq.TAB))	then
C 		Ignore names starting with spaces and tabs
		ARRNAM = '      '
		return
	    else
		jj = j
	    endif
	    if(symb.eq.TAB .or. symb.eq.' ' .or. symb.eq.TCR)	then
		je = j
		goto	1
	    else
		name(j:j) = symb
	    endif
	enddo
 1	continue
C "je" - position of 1st blanck. Append with "X" (if absent) and spaces.
	if (je.le.6 .and. name(je-1:je-1) .ne. 'X')	then
	    name(je:je) = 'X'
	    je = je+1
	endif
	if (je.le.6)	then
	    do	j=je,6
		name(j:j) = ' '
	    enddo
	endif
	ARRNAM = name
	end
C=======================================================================
	subroutine	CHECKU
     >		(INTYPE,ABC,AB,XBDRY,YX,jrad,jbdry,STRING,FILENA)
C----------------------------------------------------------------------|
C Consistency check for grid array YX(1:jrad) and plasma boundary AB/ABC
C Input:
C	ABC	- 
C	AB	- 
C	YX(1:jrad) - array for a "radial" coordinate 
C	jrad	- YX array dimensionality
C	STRING	- U-file "Independent variable" description
C	FILENA	- U-file name
C Analyse array YX and returns proper values for XBDRY and jbdry
C Output:
C	INTYPE	- 
C	XBDRY	= ABC or AB depending on INTYPE selected
C	jbdry	- is determined from {YX(jbdry) <= ABC} or {YX(jbdry) <= AB}
C		  for	{INTYPE = 10} or {INTYPE = 11}, respectively
C		jbdry = jrad	if 	if YX(jrad) < ABC <= AB
C----------------------------------------------------------------------|
	implicit none
	integer		jrad,jbdry,j,jj,INTYPE,KILLBL,length
	double precision	ABC,AB,XBDRY,YX(jrad)
	character	FILENA*(*),STRING*132,STRAD*12

	jbdry = 0
	XBDRY = YX(jrad)
	STRAD = STRING(21:31)
	jj = KILLBL(STRING,80)
	call UPCASE(31,STRING)
	if	(STRING(1:8) .eq. 'MINORRAD')	then
	    INTYPE = 10
	elseif	(STRING(1:8) .eq. 'MAJORRAD')	then
	    INTYPE = 19
	elseif	(STRING(1:3) .eq. 'RHO')	then
C	elseif	(STRING(1:11) .eq. 'RHOTOROIDAL')	then
	    INTYPE = 12
	elseif	(STRING(1:12) .eq. 'POLOIDALFLUX')	then
	    INTYPE = 13
	else
	    write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'"',
     +		  '    Unrecognized "radial" variable. Input ignored.'
     +		//'    Allowed options are:'
     +		//' 	Minor Radius        m'
     +		//' 	Major Radius        m'
     +		//' 	Rho Toroidal, normalized'
     +		//' 	Poloidal Flux, normalized'
	    INTYPE = -1
	    return
	endif
	jj = KILLBL(STRAD,11)
	call UPCASE(1,STRAD)
	if ((INTYPE.eq.10 .or. INTYPE.eq.19) .and.
     +		STRAD(1:1).ne.'M')	then
	   write(*,*)">>> Warning: Inconsistency in the input data"
	   write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'": ',
     >	   '    radial grid is expected to be given in "m"'
	endif

C	if (XBDRY .lt. (ABC+AB)/2.d0)	INTYPE = 11

	if (INTYPE .eq. 10 .and. abs(XBDRY-AB) .gt. 0.3*AB/jrad)  then
C	   write(*,*)'>>> Inconsistency in the input data,  ',
C     >			' U-file "',FILENA(1:length(FILENA)),'"'
C	   write(*,'(1A26,1F5.3,1A22,1F5.3,1A1)')
C     >'              Edge radius ',XBDRY,'m does not match AB = ',AB,'m'
	   if(XBDRY.lt.AB)	then
C		write(*,'(1A51,1F5.3,1A1)')
C     >	  '     Warning: artificial data will be added at a = ',AB,'m'
		jbdry = jrad
	   else
		do	j=jrad,1,-1
		    if (YX(j) .gt. AB)	jbdry = j
		enddo
C		if (jbdry .lt. jrad)	write(*,'(1A31,1F5.3,1A13)')
C     >	'              Data beyond  a = ',AB,'m are ignored'
		XBDRY = AB
	   endif
	elseif (INTYPE.eq.11 .and. abs(XBDRY-ABC).gt.0.3*ABC/jrad) then
C	   write(*,*)'>>> Inconsistency in the input data,  ',
C     >			' U-file "',FILENA(1:length(FILENA)),'"'
C	   write(*,'(1X,1A25,1F5.3,1A22,1F6.3,1A1)')
C     >'             Core radius ',XBDRY,'m does not match ABC =',ABC,'m'
	   if(XBDRY.lt.ABC)	then
C		write(*,'(1A51,1F5.3,1A1)')
C     >	  '     Warning: Artificial data will be added at a = ',ABC,'m'
		jbdry = jrad
	   else
		do	j=jrad,1,-1
		    if (YX(j) .gt. ABC)	jbdry = j
		enddo
C		if (jbdry .lt. jrad)	write(*,'(1A31,1F5.3,1A13)')
C     >	'              Data beyond  a = ',ABC,'m are ignored'
		XBDRY = ABC
	   endif
	endif
	if (INTYPE .eq. 18)	then
	   write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'"',
     >			" Don't know a distance to the major axis.",
     >		     '           Set to RTOR'
	endif
	if (INTYPE .eq. 19)	then
C	   write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'"',
C     >			" Don't know distance to the major axis.",
C     >		     '           Set to RTOR'
	endif
	if (jbdry .eq. 0)	then
	    jbdry = jrad
	endif
	end
C=======================================================================
	subroutine	ANLSTR(string,fnam,lfnam,vnam,lvnam,factor,ierr)
C-----------------------------------------------------------------------
C The subroutine analizes a "string" pointing to an U-file in "exp" data
C Input:
C    string*132  Allowed format (case insensitive):
C
C    "A_name   U-file:U_file_name  <  Internal_name   [factor]:Number"
C
C Output:
C    fname  = U_file_name(1:lfnam)
C    vname  = Internal_name(1:lvnam)
C    factor = Number
C ierr = 1 is returned when a sign U-file (case unsensitive) is absent
C ierr = 2 is returned if a u-file name does not fit the required format
C ierr = 3 Internal_name includes more than one word
C-----------------------------------------------------------------------
	implicit none
	character*(*)	string,fnam,vnam
	character*1	TZER,TCR,TAB,symb,ufile*6,str*60
	integer	JU,JF,JN,j,js,jl,lfnam,lvnam,ierr
	integer	KILLBL,nextwd,noblan,length
	double precision	factor
	ierr = 0
	j = 0
	TZER=char(j)
	j = 13
	TCR=char(j)
	j = 9
	TAB=char(j)
C --------------------- Check for Ex/U-file command string
C Splitting "string" in two pieces
C Before:
C          "Name U-file:u_file_name < Internal_name factor:num_factor"
C After:
C  string: "Name U-file:u_file_name :num_factor"
C  str:    "Internal_name" -> STR40(1:jl)
	j  = len(string)
C	write(*,*)"Input string of length ",j
C	write(*,*)'"',string(1:j),'"'
	vnam = ' '
	jl = 0
	ju = index(string,'<')
	if (ju .eq. 0)	goto	3
	js = index(string(ju+1:),':')
	if (js .gt. 0)	then
	   str = string(ju+1:ju+js-1)
	   string = string(1:ju-1)//string(ju+js:)
	else
	   str = string(ju+1:)
	   js = len(str)
	   string = string(1:ju-1)
	endif
	ju = 1
	jl = 0
 1	jl = jl+1
	jf = nextwd(str(ju:))
	if (jf .eq. 0)	goto	2
	ju = ju+jf-1
	goto	1
 2	continue
C jl No. of words in str
C ju Position of the last word
C js Length of a string for analysis str
	if (jl.lt.1 .or. jl.gt.2)			  goto	99
	call   UCASE(str(1:js))
	if (jl.eq.2 .and. str(ju:ju+5).ne.'FACTOR')	  goto	99
	if (jl.eq.2)	     str(ju:ju+5)  = '      '
	jl = KILLBL(str,js)
	vnam = str(1:jl)
 3	continue
	lvnam = jl

C JU is a position of the 1st ":"
C JF is a position of the 2nd ":"
	JU = index(string,':')
	if (JU .eq. 0)	then
	   ierr = 1
	   return
	endif
	do   j=1,JU-1
	   fnam(j:j) = string(j:j)
	enddo

	fnam(JU:JU) = TZER
	j = KILLBL(fnam,len(fnam))
	ufile = fnam(j-5:j)
	call UPCASE(6,ufile)

C	write(*,*)'"',ufile(1:6),'"'
	if     (ufile .eq. 'U-FILE')	then
	   ierr = 0
	elseif (ufile .eq. 'J-FILE' .or.
     &		ufile .eq. 'X-FILE')	then
	   ierr = -1
C???	   fnam = str(1:jl)
C	   write(*,*)"Exit"
C	   return
	else
	   ierr = 1
	   return
	endif

	factor = 1.d0
	JN = len(string)
	JS = index(string(JU+1:),':')
	if (JS .ne. 0)	then
	   JS = JU+JS
	   JS = min(JS,JN)
	   read(string(JS+1:),*)factor
	else
	   JS = JN
	endif

	JL = index(string,'<')
	if (JL .ne. 0)	JS = min(JL,JS)
	JS = JS-1			! 2nd delimiter ':' or '<'
	j = noblan(string(JU+1:JS))
	str = string(JU+j:JS)
	j = min(JS-JU-j+1,len(str))
	jl = nextwd(str(1:j))
	if (jl .ne. 0)	then
	   j = jl-1
	   str(jl:) = ' '
	else
	   jl = j
	endif
	jl = length(str(1:jl))

C	write(*,*)'"',str(1:jl),'"',jl
	lfnam = jl
C Extension is not obligatory!
C	if (string(jl-3:jl-3) .ne. '.')	then
C		ierr = 2
C		return
C		endif
	if (  str(1:1+3) .ne.   'udb/'
     &	.and. str(1:1+5) .ne. './udb/' .and. ierr.ne.-1)	then
	   lfnam = lfnam+4
	   fnam='udb/'//str(1:jl)
	else
	   fnam = str(1:jl)
	endif
	fnam(lfnam+1:lfnam+1) = TZER
C	write(*,*)'"',fnam(1:lfnam+1),'"',lfnam,factor
	return
 99	ierr = 3
	end
C======================================================================|
C  Subroutine minimizes the value of functional
C  INTEGRAL(alfa*P(x)*(dU/dx)**2+(U-F)**2)*dx,
C  where FO(NO) is a function, given on the grid XOld(NOld)
C	P(x) is equal to unit now
C	ALFA=alfa<<0.01*XO(NO)**2 is regularizator
C	NO - number of old grid points
C	N=<NRD - number of new grid points
C	0<=XO(NO) - old grid |     both grids are arbitrary
C	0<=XN(N)  - new grid |     but XO(NO)=XN(N)
C	FO(NO) - origin function, given on the grid XO(NO)
C	FN(N) - smoothed function on grid XN(N)
C  The result is function FN(XN), given on the new grid
C	with additional conditions:
C	dFN/dx(x=0)=0 - cylindrical case and
C	FN(XN(N))=FO(XO(NO))
C----------------------------------------------------------------------|
	subroutine	SMOOTH(ALFA,NO,FO,XO,N,FN,XN)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	integer	NO,N,J,I
	double precision	ALFA,XO(*),FO(*),XN(*),FN(*),P(NRD)
	double precision	YF,YX,YP,YQ,YD,FJ
	if (N .gt. NRD .or. NO .le. 0)	then
		write(*,*)' >>> SMOOTH: array is out of limits'
		call	a_stop
	endif
	if (NO .eq. 1)	then
	   do	j=1,N
		FN(j) = FO(1)
	   enddo
	   return
	endif
	if (NO .eq. 2)	then
	  do	j=1,N
	   FN(j)=(FO(2)*(XN(j)-XO(1))-FO(1)*(XN(j)-XO(2)))/(XO(2)-XO(1))
	  enddo
	  return
	endif
	if (N .lt. 2)	then
		write(*,*)' >>> SMOOTH: no output grid is provided'
		call	a_stop
	endif
	if (abs(XO(NO)-XN(N)) .gt. XN(N)/N)	then
	    write(*,*)'>>> SMOOTH: grids are not aligned'
	    write(*,'(1A23,I4,F8.4)')'     Old grid size/edge',NO,XO(NO)
	    write(*,'(1A23,I4,F8.4)')'     New grid size/edge',N,XN(N)
	    call	a_stop
	endif
	do	1	j=2,N
	   YP = (XN(j)-XN(j-1))
	   if (YP .le. 0.d0)	then
	write(*,*)'>>> SMOOTH: new grid is not increasing monotonically'
	      write(*,'(A,I4,A,F8.4)')'Node ',j-1,'   Value',XN(j-1)
	      write(*,'(A,I4,A,F8.4)')'Node ',j,  '   Value',XN(j)
	      call	a_stop
	   endif
	   P(j)	=ALFA/YP/XO(NO)**2
 1	continue
	P(1)	=0.
	FN(1)	=0.
	I	=1
	YF	=(FO(2)-FO(1))/(XO(2)-XO(1))
	YX	=2.d0/(XN(2)+XN(1))
	YP	=0.
	YQ	=0.
	do	5	j=1,N-1
		if(XO(I) .gt. XN(j))	GO TO 4
 3		I	=I+1
		if(I .gt. NO)	I=NO
		if(I .ne. NO .and. XO(I) .lt. XN(j))	GOTO	3
		YF	=(FO(I)-FO(I-1))/(XO(I)-XO(I-1))
 4		FJ	=FO(I)+YF*(XN(j)-XO(I))
		YD=1.d0+YX*(YP+P(j+1))
		P(j)	=YX*P(j+1)/YD
		FN(j)	=(FJ+YX*YQ)/YD
		if (j .eq. N-1)	goto	5
		YX	=2.d0/(XN(j+2)-XN(j))
		YP	=(1.d0-P(j))*P(j+1)
		YQ	=FN(j)*P(j+1)
 5	continue
	FN(N)	=FO(NO)
	do	6	j=N-1,1,-1
		FN(j)	=P(j)*FN(j+1)+FN(j)
 6	continue
	end
C======================================================================|
	subroutine	SMOTH(ALFA,NO,FO,XO,N,FN,XN)
C----------------------------------------------------------------------|
C Same as SMOOTH but with a free boundary
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	integer	NO,N,J,I
	double precision	ALFA,XO(*),FO(*),XN(*),FN(*),P(NRD)
	double precision	YF,YX,YP,YQ,YD,FJ
	if (N .gt. NRD .or. NO .le. 0)	then
		write(*,*)' >>> SMOTH: array is out of limits'
		call	a_stop
	endif
	if (NO .eq. 1)	then
	   do	j=1,N
		FN(j) = FO(1)
	   enddo
	   return
	endif
	if (NO .eq. 2)	then
	  do	j=1,N
	   FN(j)=(FO(2)*(XN(j)-XO(1))-FO(1)*(XN(j)-XO(2)))/(XO(2)-XO(1))
	  enddo
	  return
	endif
	if (N .lt. 2)	then
		write(*,*)' >>> SMOTH: no output grid is provided'
		call	a_stop
	endif
	if (abs(XO(NO)-XN(N)) .gt. XN(N)/N)	then
	    write(*,*)'>>> SMOTH: grids are not aligned'
	    write(*,'(1A23,I4,F8.4)')'     Old grid size/edge',NO,XO(NO)
	    write(*,'(1A23,I4,F8.4)')'     New grid size/edge',N,XN(N)
	    call	a_stop
	endif
	do	1	j=2,N
	P(j)	=ALFA/(XN(j)-XN(j-1))/XO(NO)**2
 1	continue
	P(1)	=0.
	FN(1)	=0.
	I	=1
	YF	=(FO(2)-FO(1))/(XO(2)-XO(1))
	YX	=2.d0/(XN(2)+XN(1))
	YP	=0.
	YQ	=0.
	do	5	j=1,N-1
		if(XO(I) .gt. XN(j))	GO TO 4
 3		I	=I+1
		if(I .gt. NO)	I=NO
		if(I .ne. NO .and. XO(I) .lt. XN(j))	GOTO	3
		YF	=(FO(I)-FO(I-1))/(XO(I)-XO(I-1))
 4		FJ	=FO(I)+YF*(XN(j)-XO(I))
		YD=1.d0+YX*(YP+P(j+1))
		P(j)	=YX*P(j+1)/YD
		FN(j)	=(FJ+YX*YQ)/YD
		if (j .eq. N-1)	goto	5
		YX	=2.d0/(XN(j+2)-XN(j))
		YP	=(1.d0-P(j))*P(j+1)
		YQ	=FN(j)*P(j+1)
 5	continue
	FN(N) = ((2.d0-P(N-2))*FN(N-1)-FN(N-2))/(1.d0-(2.d0-P(N-2))*P(N-1))
	do	6	j=N-1,1,-1
		FN(j)	=P(j)*FN(j+1)+FN(j)
 6	continue
	end
C======================================================================|
	subroutine	TAMP(ALFA,X1,X2,NO,FO,XO,N,FN,XN)
C----------------------------------------------------------------------|
C Same as SMOOTH but P=P(ALFA)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	integer	NO,N,J,I
	double precision	ALFA,X1,X2,XO(*),FO(*),XN(*),FN(*)
	double precision	P(NRD),YF,YX,YP,YQ,YD,FJ
	if (N .gt. NRD .or. NO .le. 0)	then
		write(*,*)' >>> TAMP: array is out of limits'
		call	a_stop
	endif
	if (NO .eq. 1)	then
	   do	j=1,N
		FN(j) = FO(1)
	   enddo
	   return
	endif
	if (NO .eq. 2)	then
	  do	j=1,N
	   FN(j)=(FO(2)*(XN(j)-XO(1))-FO(1)*(XN(j)-XO(2)))/(XO(2)-XO(1))
	  enddo
	  return
	endif
	if (N .lt. 2)	then
		write(*,*)' >>> TAMP: no output grid is provided'
		call	a_stop
	endif
	if (abs(XO(NO)-XN(N)) .gt. XN(N)/N)	then
	    write(*,*)' >>> TAMP: input grids are not aligned'
	    write(*,'(1A23,I4,F8.4)')'     Old grid size/edge',NO,XO(NO)
	    write(*,'(1A23,I4,F8.4)')'     New grid size/edge',N,XN(N)
	    call	a_stop
	endif
	do	1	j=2,N
	P(j)	=ALFA/(XN(j)-XN(j-1))/XO(NO)**2
C	YX = max(0.,(4.*(XN(j)-X1)*(X2-XN(j))/(X2-X1)**2)**3)
C	if (XN(j).lt.X1 .or. XN(j).gt.X2)	P(j) = P(j)*YX
	P(j) = P(j)*exp(-((XN(j)-X1)/X2)**2)
 1	continue
	P(1)	=0.
	FN(1)	=0.
	I	=1
	YF	=(FO(2)-FO(1))/(XO(2)-XO(1))
	YX	=2.d0/(XN(2)+XN(1))
	YP	=0.
	YQ	=0.
	do	5	j=1,N-1
		if(XO(I) .gt. XN(j))	GO TO 4
 3		I	=I+1
		if(I .gt. NO)	I=NO
		if(I .ne. NO .and. XO(I) .lt. XN(j))	GOTO	3
		YF	=(FO(I)-FO(I-1))/(XO(I)-XO(I-1))
 4		FJ	=FO(I)+YF*(XN(j)-XO(I))
		YD=1.d0+YX*(YP+P(j+1))
		P(j)	=YX*P(j+1)/YD
		FN(j)	=(FJ+YX*YQ)/YD
		if (j .eq. N-1)	goto	5
		YX	=2.d0/(XN(j+2)-XN(j))
		YP	=(1.d0-P(j))*P(j+1)
		YQ	=FN(j)*P(j+1)
 5	continue
	FN(N)	=FO(NO)
	do	6	j=N-1,1,-1
		FN(j)	=P(j)*FN(j+1)+FN(j)
 6	continue
	end
C======================================================================|
	double precision  function	RA2Z(RTOR,YR,YA,YS,YE,YT)
C=========================================================09-FEB-97
C  Z [m] as a function of major and minor radii Z(R,a)
C  	above the plasma midplane 
C	YS = SHIF(A), YE=ELON(A), YT=TRIA(A), YV=SHIV(A)
C NB !  Up-down shift not allowed (G.P.)
C===========================================================Polevoi
	implicit none
	double precision    RTOR,YR,YA,YS,YE,YT,YF,YF1,YT2,YSIN2
		YF	=(YR-RTOR-YS)/YA
		YF1	=2.d0*YF*YT
		YT2	=(1.d0+YF1)
	if(YF.ge.1.d0.or.YF.lt.-1.d0.OR.YT2.LE.0.) 	goto 99
		YF1	=2.d0*YF*YT
		YSIN2	=(1.d0-YF)*(1.d0+YF)/YT2
	if(YSIN2.gt..001 .and. YT.gt..01)	then
		YT2	=2.d0*YT*YT
		YSIN2	=(SQRT(1.d0+2.d0*(YF1+YT2))-(1.d0+YF1))/YT2
			endif
		RA2Z	=YA*YE*SQRT(YSIN2)
	RETURN
 99		RA2Z	=0.
 	RETURN
	END
C======================================================================|
	double precision	function	RZ2A(YR,YZ,N)
C-----------------------------------------------------------------------
C	The function returns A(r,z) where
C
C	r = R0 + del(A) + A*[cos(theta)-tri(A)*sin^2(theta)]
C	z = UPD + A*elo(A)*sin(theta)
C
C	and del(j), elo(j), tri(j) are given as arrays[1:N] 
C					 on the grid A=AMETR(j)
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	integer	N,j,j1
	double precision YR,YZ,Y1,YAS,YAO,YHOR,YVER,YELO,YTRI,YA,QUADIN
	YAS = ((YZ-SHIV(N))/ELON(N))**2
	YA = sqrt(YAS+(YR-RTOR-SHIF(N))**2)
	do j=1,20
		YAO = YA
		YHOR = QUADIN(N,AMETR,SHIF,YA,YAS,j1)
		YELO = QUADIN(N,AMETR,ELON,YA,YAS,j1)
		YTRI = QUADIN(N,AMETR,TRIA,YA,YAS,j1)
		YVER = QUADIN(N,AMETR,SHIV,YA,YAS,j1)
		YAS = ((YZ-YVER)/YELO)**2
		Y1 = YAS
		if (YA .gt. 1.E-4)	Y1 = YAS/YA
		YA = sqrt(YAS+(YR-RTOR-YHOR+YTRI*Y1)**2)
C	write(*,'(1I5,2F11.7,10F7.3)')j,YA,YAO,YAS,Y1,YHOR,YELO,YTRI
		j1 = j
		if (abs(YA-YAO).lt.1.E-10)	goto	99
	enddo
 99	RZ2A = min(YA,AB)
C	write(*,'(1I5,5F15.10)')j,YR,YZ,YA
	end
C=======================================================================
	double precision function	QUADIN(NG,XGRID,FUNC,X,PQD,jj)
C-----------------------------------------------------------------------
C   Quadratic interpolation of a function FUNC(NG) 
C   		given on a grid XGRID(NG) to a position X
C Input		NG, XGRID(NG), FUNC(NG), X
C Output	QUADIN, PQD, jj
C-----------------------------------------------------------------------
	implicit none
	integer	NG,j,jj
	double precision
     1			XGRID(NG),FUNC(NG),X,YF1,YF2,YF3,
     2			PQD,Y,YY,YX,YD21,YD23,YD31,YDX1,YDX2,YDX3
C No extrapolation
	YX = max(X,0.d0)
	YX = min(YX,XGRID(NG))
	YY = 0.
	jj = 1
	do	4	j=1,NG
	    Y = XGRID(j)-YX
	    jj = j
	    if (Y) 1,5,3
 1	    YY = Y
	    goto	4
 3	    if (Y .gt. -YY)	jj = jj-1
	    goto	5	
 4	continue
 5	continue
	if (jj .le. 1)	jj = 2
	if (jj .ge. NG)	jj = NG-1

	YD21 = XGRID(jj)-XGRID(jj-1)
	YD23 = XGRID(jj)-XGRID(jj+1)
	YD31 = XGRID(jj+1)-XGRID(jj-1)
	YDX1 = YX-XGRID(jj-1)
	YDX2 = YX-XGRID(jj)
	YDX3 = YX-XGRID(jj+1)
	if (YD21 .le. 0. .or. YD23 .ge. 0.)	then
	    QUADIN = FUNC(jj)
	    PQD = 0.
	    return
	endif
	YF1 = FUNC(jj-1)/YD21/YD31
	YF2 = FUNC(jj)/YD21/YD23
	YF3 = FUNC(jj+1)/YD31/YD23
	QUADIN = YF1*YDX2*YDX3+YF2*YDX1*YDX3-YF3*YDX1*YDX2
	PQD = YF1*(YDX2+YDX3)+YF2*(YDX1+YDX3)-YF3*(YDX1+YDX2)
	end
C=======================================================================
	double precision function	RECTAN(NG,XGRID,FUNC,X)
C-----------------------------------------------------------------------
C	Step-function FUNC(X) is defined
C Input		NG 	    grid size
C		XGRID(1:NG) grid
C		FUNC(1:NG)  grid function
C		X	    function argument
C Output	RECTAN      function value
C		jj the grid point nearest to X (not used)
C-----------------------------------------------------------------------
	implicit none
	integer	NG,j,jj
	double precision	XGRID(NG),FUNC(NG),X,Y,YY,YX
C No extrapolation
	YX = max(X,0.d0)
	YX = min(YX,XGRID(NG))
	YY = 0.
	jj = 1
	do	4	j=1,NG
	    Y = XGRID(j)-YX
	    jj = j
	    if (Y) 1,5,3
 1	    YY = Y
	    goto	4
 3	    if (Y .gt. -YY)	jj = jj-1
	    goto	5	
 4	continue
 5	continue
	if (jj .le. 1)	jj = 1
	if (jj .ge. NG)	jj = NG

	RECTAN = FUNC(jj)
	end
C======================================================================|
	subroutine	APPRID(ID,IM,IY,IH,MI,IS)
C-----------------------------------------------------------------------
C The subroutine forms string RUNID and additionally returns 
C	date and time when those are not defined (calling from INIT)
C-----------------------------------------------------------------------
	implicit none
	character*60  STRI
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	integer	j,j1,jj,length,lonlen,nsymb,ID,IM,IY,IH,MI,IS,IERR
C	write(*,*)"ID =",ID
	if (ID .eq. 0)	call	TIMDAT(IY,IM,ID,IH,MI,IS)
	RUNID(1:17) = "=== ASTRA 5.3 ==="
	j = nsymb(VERSION,'Version')
	RUNID(11:13) = version(j+8:j+10)
	write(RUNID(18:),100)ID,IM,IY-100,IH,MI
	jj = length(EQNAME)
	j = min(14,jj)		! max mod file name length
	write(RUNID(33:),'(" === Model: ",A)')EQNAME(1:j)
	j = j+45
	jj = length(RDNAME)
	if (j+jj.le.61)	then
	   write(RUNID(j:),'(" === Data file: ",A," ===")')
     >							RDNAME(1:jj)
	elseif (j+jj.le.66)	then
	   write(RUNID(j:),'(" === Data: ",A," ===")')	RDNAME(1:jj)
	elseif (j+jj.le.70)	then
	   write(RUNID(j:),'(" < Data: ",A,10A1)')	RDNAME(1:jj)
     >		,' ',('=',j1=1,71-j-jj)
	elseif (j+jj.le.75)	then
	   write(RUNID(j:),'(" << ",A,10A1)')	 	RDNAME(1:jj)
     >		,' ',('=',j1=1,76-j-jj)
	else
	   write(RUNID(j:),'(" << ",A)')	RDNAME(1:min(77-j,jj))
	endif
	jj = lonlen(RUNID)
	if (jj.lt.79)	then
	   j1 = (80-jj)/2
	   do  j=80-j1,j1+1,-1
	      RUNID(j:j) = RUNID(j-j1:j-j1)
	   enddo
	   write(RUNID(1:j1),'(80A1)')(' ',j=1,j1)
	   write(RUNID(81-j1:80),'(80A1)')(' ',j=1,j1)
	endif
 100	format(1I3,2('-',1I2.2),1I3,':',1I2.2)
C	write(*,*)RUNID(1:jj)
	end
C======================================================================|
	subroutine makemovie(MOVIE,PNMNAME)
C----------------------------------------------------------------------|
C Input:
C        (=0,*) skip
C        (>0,x) continue
C        (=0,x) stop
C----------------------------------------------------------------------|
	implicit	none
	include	'for/parameter.inc'	! Remove after implementing
	include 'for/outcmn.inc'	! pop-up window
	character	STRI*32,PNMNAME*(*)
	integer		MOVIE,Magenta,jj,length
	external	length
	data		Magenta/14/
	save		Magenta
	if (PNMNAME(1:1) .eq. "*")	return		! Skip

	call	markloc(" calling px2pnm"//char(0))
	if (MOVIE .ne. 0)	then		! Save frame
	   jj = 1
	   if (MOVIE .gt. 9)	jj = 2
	   if (MOVIE .gt. 99)	jj = 3
	   if (MOVIE .gt. 999)	jj = 4
	   write(STRI,'(1I3.3)')MOVIE
C	   write(STRI,'(1I<jj>)')MOVIE
	   STRI = "Making Movie:  Frame # "//STRI(1:jj)//char(0)
	   call	colorb(Magenta)
	   call	textbf(20,XWH-124,STRI(1:jj),23+jj)
	   call redraw(0)
	   call	px2pnm(PNMNAME)
	   MOVIE = MOVIE+1
	else					! Save & Exit
	   call	colorb(Magenta)
	   call	textbf(20,XWH-124,"Saving movie in tmp/movie.mpeg  ",32)
	   call redraw(0)
	   jj = length(PNMNAME)
	   call system(".srv/esc2movie "//PNMNAME(1:jj)//char(0))
	   call system("rm tmp/*.ppm"//char(0))
	   call	textbf(20,XWH-124,"                                ",32)
	   PNMNAME='*'
	endif
	end
C======================================================================|
	subroutine getnames(YEQUNAME,YEXPNAME)
C----------------------------------------------------------------------|
C The subroutine can be called from C function, returns EQNAME, RDNAME
C	G.V.Pereverzev	16.02.2004
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/outcmn.inc'
	integer	j,j1,j2,length
	character	YEQUNAME*(*),YEXPNAME*(*)
	entry	getnames_(YEQUNAME,YEXPNAME)
C	write(*,*)len(EQNAME),len(RDNAME)
C	write(*,*)length(EQNAME),length(RDNAME)
	j1 = length(EQNAME)
	do   j	=1,j1
	   YEQUNAME(j:j) = EQNAME(j:j)
	enddo
	YEQUNAME(j1+1:j1+1) = char(0)
C	write(*,*)'"',(YEQUNAME(j:j),j=1,j1),'"',j1
	j2 = length(RDNAME)
	do   j	=1,j2
	   YEXPNAME(j:j) = RDNAME(j:j)
	enddo
	YEXPNAME(j2+1:j2+1) = char(0)
C	write(*,*)'"',(YEXPNAME(j:j),j=1,j2),'"',j2
	end
C======================================================================|
