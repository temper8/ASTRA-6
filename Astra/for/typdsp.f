C Modules >>> TYPDSP, PCOVA, UF1DWA, UF2DWA, OPENAP, OPENRD, OPENWT <<<
C >>> SETFRM, AFRAME, DNSTR, UPSTR,  TIMEDT, PUTXY,  NEGA,   SETFNA <<< 
C >>> SKIPM  <<< 
C======================================================================|
	subroutine	TYPDSP(NCH,YN,ITIMES,TTOUT,TOUT)
C NCH= 5 - terminal, 0 - file (old format), 1 - file (new format)
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	integer	ITIMES
	double precision	YN,YQ,TTOUT(ITIMES),TOUT(ITIMES,NRW)
	integer NCH,NP1,ITBE,ITEND,ITEN,IERR,MODEX,NLINSC,NVAR,NNJ,JN0
	integer JBE,JEND,J,JEN,JJ,J1,JLR,KILLBL,length,WarningColor
	character*6 CH6,STRI*118,STRMN*40,CONN(6),FNAME*132
	logical	EXI
	save	STRMN,CONN,JN0,NLINSC,JLR,WarningColor
	data	STRMN/' R=     a=     B=     I=     q=     <n>='/
     +	CONN/' CF   ',' CV   ',' CH   ',' CCD  ',' CBND ',' CRAD '/
     +	JN0/0/ NLINSC /50/ JLR/426/ WarningColor/30/
C NLINSC - maximum line number 
	MODEX = XOUT+.49
C MODEX = 0	[0,AB]  against "a"
C MODEX = 1	[0,ABC] against "a"
C MODEX = 2	[0,ROC] against "rho"
C MODEX = 3	[FP(1),FP(NA1)] against "psi"
C otherwise	Unknown option => MODEX=0
	NP1 = NAB
	if (MODEX .ge. 1 .and. MODEX .le. 3 .or. MOD10 .eq. 3)
     +	NP1 = NA1

	if(NCH.eq.5)	goto	300
	if(NCH.ne.0 .and. NCH.ne.1)	return
C Write output to a file:
	J = length(RDNAME)
	J1 = length(EQNAME)
	FNAME=AWD(1:length(AWD))//'dat/'//RDNAME(1:J)//'.'
     >			//EQNAME(1:J1)//char(0)
	JEND = KILLBL(FNAME,132)
	call	SETFNA(FNAME,JEND)
C	write(*,*)'Returned:  "',FNAME(1:JEND),'"',JEND
	call colovm(WarningColor)
	STRI='>>>  Data are written into the file: '//FNAME(1:JEND)//' '
	write(*,*)STRI(1:JEND+37)
C	call textvm(JN0,JLR,STRI,37+JEND)
	call	OPENWT(7,FNAME(1:JEND),0,IERR)
	if(IERR.gt.0)	pause 'Data file error'
C Creating UPSTRI
	STRI=XLINE1(2:17)
	STRI(17:)=STRMN
	call FMTF4(STRI(20:23),RTOR)
	call FMTF4(STRI(27:30),ABC)
	call FMTF4(STRI(34:37),BTOR)
	call FMTF4(STRI(41:44),IPL)
C Triangularity corrected MHD q (accoding to ITER guidelines)
C	YQ	=ELON(NA)**2
C	YD	=TRIA(NA)
C	YQ=(1.+YQ*(1.+YD**2*(2.-1.2*YD)))/(MU(NA)*(1.+YQ))
	YQ=1./MU(NA)
	call FMTF4(STRI(48:51),YQ)
	call FMTF4(STRI(57:60),YN)
	write(STRI(62:76),103)TIME
	call FMTF4(STRI(77:80),1000.*TAU)
 102	format('   Time',16(3X,1A4))
 103	format('Time=',1F6.3,' dt=')
 104	format(1X,1A120)
	write(7,104)STRI
	if(NCH.eq.1)	goto	400

C Writing radial data
	if(MOD10.le.5)			then
	    JBE=1
	    JEND=16
 3	    JEN=MIN0(NTOUT,JEND)
	    write(7,102)	(NAMET(J),J=JBE,JEN)
	    STRI=' '
	    call FMTXF4(STRI(1:5),TIME)
	    do	91	J=JBE,JEN
		JJ=7*(J-JBE)+8
 91		call FMTXF5(STRI(JJ:JJ+5),TOUT(LTOUT,J))
	    write(7,104)STRI
	    if(JEN.eq.NTOUT)	go to 4
	    JBE=JEN+1
	    JEND=JEN+16
				go to 3
 4	    JBE=1
	    JEND=16
 1	    JEN=MIN0(NROUT,JEND)
	    if     (MODEX .eq. 0 .or. MODEX .eq. 1)	then
		 write(7,'("     a  ",16(3X,1A4))')(NAMER(J),J=JBE,JEN)
	    elseif (MODEX .eq. 2)	then
		 write(7,'("     rho",16(3X,1A4))')(NAMER(J),J=JBE,JEN)
	    elseif (MODEX .eq. 3 .or. MOD10 .eq. 3)	then
		 write(7,'("     psi",16(3X,1A4))')(NAMER(J),J=JBE,JEN)
	    else
		 write(7,'("     ???",16(3X,1A4))')(NAMER(J),J=JBE,JEN)
	    endif
	    do 10 J=1,NP1
		STRI=' '
		do	11	JJ=JBE,JEN
		    J1=7*(JJ-JBE+1)+1
 11		    call FMTXF5(STRI(J1:J1+5),ROUT(J,JJ))
C Different options for a radial variable 
		if     (MODEX .eq. 0)	then
			call FMTXF5(STRI(1:5),AMETR(j))
		elseif (MODEX .eq. 1)	then
			call FMTXF5(STRI(1:5),AMETR(j))
		elseif (MODEX .eq. 2)	then
			call FMTXF5(STRI(1:5),RHO(j))
		elseif (MODEX .eq. 3 .or. MOD10 .eq. 3)	then
			call FMTXF5(STRI(1:5),FP(j))
				else
			call FMTXF5(STRI(1:5),AMETR(j))
				endif
		write(7,104)STRI
 10	    continue
	    if(JEN.eq.NROUT)	go to 2
	    JBE=JEN+1
	    JEND=JEN+16
				go to 1
	endif

C Writing time data
	if(MOD10.eq.6)	then
		JBE=1
		JEND=16
 6		JEN=MIN0(NTOUT,JEND)
		write(7,102)	(NAMET(J),J=JBE,JEN)
		do	15	J1=1,LTOUT-1
		STRI=' '
		call FMTXF5(STRI(1:5),TTOUT(J1))
		do	14	J=JBE,JEN
		JJ=7*(J-JBE)+8
 14		call FMTXF5(STRI(JJ:JJ+5),TOUT(J1,J))
 15		write(7,104)STRI
		if(JEN.eq.NTOUT)	go to 2
		JBE=JEN+1
		JEND=JEN+16
					go to 6
					endif
 400	continue
C Writing constants
 2	write(7,'(10X,1A80)')RUNID
	J1=0
	do	93 JEN=1,100
	   STRI=' '
	   do	J=1,16
	      J1=J1+1
	      if (J1.gt.NCFNAM)	go to 94
	      call FMTXF5(CH6,CONSTF(J1))
	      JJ=7*(J-1)+1
	      STRI(JJ:JJ+5)=CH6
	   enddo
 93	write(7,101) CONN(JEN),STRI
 94	write(7,101) CONN(JEN),STRI
 101	format(1X,1A6,1A111)
	if(NCH.eq.1)	goto	401
	close(7)
	return

 401	continue
C Writing radial data
	if(MOD10.le.5)			then
	JBE=1
	JEND=16
 402	JEN=MIN0(NTOUT,JEND)
	write(7,'(3X,"Time",16(3X,1A4))')(NAMET(J),J=JBE,JEN)
	STRI=' '
	call FMTXF4(STRI(1:5),TIME)
	do	J=JBE,JEN
	   JJ=7*(J-JBE)+8
	   call FMTXF5(STRI(JJ:JJ+5),TOUT(LTOUT,J))
	enddo
	write(7,104)STRI
	if(JEN.EQ.NTOUT)	go to 403
	JBE=JEN+1
	JEND=JEN+16
				go to 402
 403	JBE=1
	JEND=NRW
	JEN=MIN0(NROUT,JEND)
C Different options for radial variable 
	if (MODEX .eq. 2)	then
	   write(7,'(8X,"rho ",64(8X,1A4))')(NAMER(J),J=JBE,JEN)
	   do	j=1,NP1
		write(7,408)RHO(j),(ROUT(J,JJ),JJ=JBE,JEN)
	   enddo
	elseif (MODEX .eq. 3 .or. MOD10 .eq. 3)	then
	   write(7,'(8X,"psi ",64(8X,1A4))')(NAMER(J),J=JBE,JEN)
	   do	j=1,NP1
		write(7,408)FP(j),(ROUT(J,JJ),JJ=JBE,JEN)
	   enddo
	else
	   write(7,'(8X,"a   ",64(8X,1A4))')(NAMER(J),J=JBE,JEN)
	   do	j=1,NP1
		write(7,408)AMETR(j),(ROUT(J,JJ),JJ=JBE,JEN)
	   enddo
	endif
	endif

 408    format(1PE12.3,64(1PE12.3))

C Writing time data
	if(MOD10.ne.6)	goto	409
	STRI=' '
	write(7,104)STRI
	write(7,104)STRI
	write(7,104)STRI
	write(7,104)STRI
	JBE = 1
	JEND=MIN(NTOUT,NRW)
 404	continue
	JEN = JBE+7
	write(7,'(8X,"Time",64(8X,1A4))')(NAMET(J),J=JBE,JEN)
	do	J1=1,LTOUT-1
	   STRI=' '
	   call FMTXF5(STRI(1:5),TTOUT(J1))
	   write(7,408)TTOUT(J1),(TOUT(J1,J),J=JBE,min(JEN,JEND))
	enddo
	JBE = JBE+8
	if (JBE.lt.JEND)	goto	404
 409	close(7)
	return

 300	continue
C Output to the terminal
	if(MOD10.gt.5)	goto	305
	JBE=1
	JEND=16
 301	JEN=MIN0(NROUT,JEND)
	write(STRI,302)(NAMER(J),J=JBE,JEN)
 302	format(16(1X,1A4))
	call prntxt(STRI)
	do 304 J=1,NP1
	STRI=' '
	do	303	JJ=JBE,JEN
	J1=5*(JJ-JBE+1)-4
 303	call FMTXF4(STRI(J1:J1+4),ROUT(J,JJ))
	NNJ=9*J+3
	call prntxt(STRI)
 304	continue
	if(JEN.eq.NROUT)	return
	JBE=JEN+1
	JEND=JEN+16
	GO TO 301
 305	if(MOD10.gt.7)	return
	JBE=1
	JEND=15
 306	JEN=MIN(NTOUT,JEND)
	ITBE=1
	ITEND=NLINSC
 307	ITEN=MIN(LTOUT-1,ITEND)
	write(STRI,308)	(NAMET(J),J=JBE,JEN)
 308	format(1X,'Time',15(1X,1A4))
	call	prntxt(STRI)
	do	315	J1=ITBE,ITEN
	STRI=' '
	call FMTXF4(STRI(1:5),TTOUT(J1))
	do	314	J=JBE,JEN
	JJ=5*(J-JBE)+6
 314	call FMTXF4(STRI(JJ:JJ+4),TOUT(J1,J))
	NNJ=9*(J1-ITBE)+12
	call prntxt(STRI)
 315	continue
	if(ITEN.eq.LTOUT-1)	go to 316
	ITBE=ITEN
	ITEND=ITEN+NLINSC-1
				go to 307
 316	if(JEN.eq.NTOUT)	return
	JBE=JEN+1
	JEND=JEN+15
				GO TO 306
	END
C======================================================================|
	integer	function GETIME(TIME,TIMES,NNOUT)
C The function returns
C  if NNOUT=1  then	GETIME=1
C  otherwise
C     GETIME = an index of the array TIMES(1:NNOUT) element >= TIME
C
	implicit none
	double precision	TIME,TIMES(*),DT
	integer	NNOUT,j

	if (NNOUT .le. 1)	then
	   GETIME = 1
	   return
	endif

	if (TIME.le.TIMES(1))	then
	   GETIME = 1
	   TIME = TIMES(1)
	   return
	endif

	do	3	j=2,NNOUT
	   if (TIME.gt.TIMES(j))  goto	3
	   GETIME = j
	   return
 3	continue

	GETIME = NNOUT
	TIME = TIMES(NNOUT)
	end
C=======================================================================
	subroutine	PUTXY(IX,IY,ITIMES,TTOUT,TOUT)
C----------------------------------------------------------------------|
C Input:	IX
C		IY
C		MODEY
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer	ITIMES,IX,IY,IM,NX,NY,JX,JY,JLR,j,j1,jj,JC,JL
	integer IX0,IXM,MODEX,GETIME,JW,JW4,JN0,JN1,JN2,JNC
	double precision	TTOUT(ITIMES),TOUT(ITIMES,NRW)
	double precision
     1		DX,DY,YX,YX1,YY,YY1,YA,YA1,YD,YE,YT,YRHO,YFP,YFPC,
     2		RZ2A,QUADIN,YBR,YBZ
	character	STRI*80,STRIN*80,XF7*7
	JN0 = 0
	JLR = 426
	JLR = XWH-124
C	DLINER_ theGCA, 0, XWH-128, XWW-1, XWH-128);
C	DLINER_ theGCA, 0, XWH-110, XWW-1, XWH-110);
C	DLINER_ theGCA, 0, XWH-109, XWW-1, XWH-109);
	if (MOD10 .le. 0)	return
	if (MOD10 .eq. 7)	call NEGA(IM,ITIMES,TOUT)
	call	SETFRM(IM,IX0,IXM,NX,NY)
	write(STRI,'(20A4)')('    ',j=1,20)
	STRI(7:25) = '(x,y)=(     ,     )'
	JX = IX-10
	JY = IY-10
	if (IX0.gt.JX.or.JX.gt.IXM.or.IY0.gt.JY.or.JY.gt.IYM)  goto  6
	if (MOD10 .eq. 7)	NX = 2
	DX = 1./NX
	DY = 1./NY
	YX1 = (JX-IX0+0.)/(IXM-IX0)
	YY1 = 1.-(JY-IY0+0.)/(IYM-IY0)
	if (MOD10 .eq. 6)	then
C (window_width)/(step=IDX=23)/(n_labels)=592/23/25=1.0295652
	   YX = TINIT+1.029565*YX1*abs(TSCALE)
C	   STRI( 7:12) = "(t,y)="
C	   write(*,*)"From putxy",ix,iy,YX
C	   write(*,*)j,LTOUT,TTOUT(j),YX,TIME,TTOUT(LTOUT)
	   goto	2
	endif
	if (MOD10 .eq. 8)	then
	   YX = 5.*YX1*SCM
	   YY = (YY1-.5)*SCM*DYM/IDT/IDX
	   YA1= RZ2A(YX,YY,NAB)
C	   write(*,'(10F7.4)')YX1,YY1,YX,YY
	   YD = QUADIN(NAB,AMETR,SHIF,YA1,YY1,j1)
	   YE = QUADIN(NAB,AMETR,ELON,YA1,YY1,j1)
	   YT = QUADIN(NAB,AMETR,TRIA,YA1,YY1,j1)
C	   YV = QUADIN(NAB,AMETR,SHIV,YA1,YY1,j1)
C	   write(STRI1,'(i4,10F7.4)')j1,YX,YY,YZ1,YA1,
C     >		AMETR(j1-1),AMETR(j1),AMETR(j1+1)
C	   call	textvm(XWW-77*DXLET+2,JLR-50,STRI1(1:75),75)
	   STRI( 7:12) = "(r,z)="
	   STRI(32:37) = "(a,S)="
	   STRI(38:50) = '(     ,     )'
	   call  FMTF5(STRI(39:43),YA1)
	   call  FMTF5(STRI(45:49),YD)
	   STRI(57:62) = "(E,T)="
	   STRI(63:75) = '(     ,     )'
	   call  FMTF5(STRI(64:68), YE)
	   call  FMTF5(STRI(70:74),YT)
	   if (YA1 .lt. AB)	call	BRZ(yx,yy,ybr,ybz)
	   goto	4
	endif
	do	j=1,NX
	   YX1 = YX1-DX
	   if (YX1.lt.0)	goto	1
	enddo
 1	continue
	YX = (YX1+DX)/DX
	if (MOD10 .ge. 5)	goto	2
 	YA = 0.
	YFP = 0.
	YRHO = 0.
	YFPC = 1.125*FP(1)-0.125*FP(2)
	MODEX = XOUT+.49
	if (MOD10 .eq. 3)	MODEX = 3
	if (MOD10 .eq. 4)	MODEX = 0
	if	(MODEX .eq. 0)	then
		YA = YX*AB
		YRHO = QUADIN(NA1,AMETR,RHO,YA,YD,j1)
		YFP  = QUADIN(NA1,AMETR,FP, YA,YD,j1)
	elseif (MODEX .eq. 1)	then
		YA = YX*ABC
		YRHO = QUADIN(NA1,AMETR,RHO,YA,YD,j1)
		YFP  = QUADIN(NA1,AMETR,FP, YA,YD,j1)
	elseif (MODEX .eq. 2)	then
		YRHO = YX*ROC
		YA   = QUADIN(NA1,RHO,AMETR,YRHO,YD,j1)
		YFP  = QUADIN(NA1,RHO,FP,   YRHO,YD,j1)
	elseif (MODEX .eq. 3)	then
		YFP  = YFPC+(FP(NA1)-YFPC)*YX
		YA   = QUADIN(NA1,FP,AMETR,YFP,YD,j1)
		YRHO = QUADIN(NA1,FP,RHO  ,YFP,YD,j1)
	endif
	j = YRHO/HRO+1
	if (YRHO .gt. 0.5*(RHO(NA)+ROC))	j = NA1
C	write(*,'(4(2F7.4,2X))')YX,YA,YRHO,YRHO/ROC,YFP
	STRI(32:33) = "a="
	call  FMTF5(STRI(34:38),YA)
	YRHO = YRHO/ROC
	STRI(39:49) = 'm,   rho_t='
	call  FMTF5(STRI(50:54),YRHO)
	if (YFP .gt. YFPC)	then
	   YFP = sqrt((YFP-YFPC)/(FP(NA1)-YFPC))
	else
	   YFP = 0.
	endif
	STRI(55:64) = ",   rho_p="
	call  FMTF5(STRI(65:69),YFP)
	STRI(70:77) = ",  Node:"
	write(STRI(78:80),'(1I3)')j

 2	continue
	do	j=1,NY
		YY1 = YY1-DY
		if (YY1.lt.0)	goto	3
	enddo
 3	YY = (YY1+DY)/DY
	if (MOD10.le.1 .or. MOD10.ge.6 .or. MODEY .ne. -1)	goto  4
	if (JY-IY0 .gt. (IYM-IY0)/NY)	YY = YY-1.
 4	continue
	if (MOD10 .eq. 6)	then

	JN1 =-DXLET
	JN2 = YSCMAX+15+3*DYLET
	JC = 0
	JL = 0
C	write(STRIN,'(79X,1A1)')' '
C	write(STRI,'(79X,1A1)')' '

	   YY1 = max(TIME,TTOUT(LTOUT-1),TTOUT(LTOUT))
	   YY1 = min(YX,YY1)
	   YY1 = max(TTOUT(1),YY1)
	   j = GETIME(YX,TTOUT,LTOUT)
	   do 5	J1=1,NTOUT
	      JW=NWIND3(J1)-8*NSCR(MOD10)
	      if (NAMET(J1) .eq. '    ') JW = 0
	      if (JW.le.0 .or. JW.gt.8)	goto	5
C	      if (MODEY .ne. 0)	then		! Disable color numbers
		 call	DNSTR(j,ITIMES,TOUT)	! for the all modes
		 goto	5			! except MODEY=0
C	      endif
C Optionally numbers can be colored similar to curves (now suppressed):
	      if (MODEY .eq. 1)	JW4 = 2		! two-window mode
	      if (MODEY .eq.-1)	JW4 = 4		! four-window mode
	      if (MODEY .eq. 0)	JW4 = 1		! one-window mode
C	      JNC = (JW-1)/JW4
C	      call	colovm(JNC+2)
	      JC = JC+1		! Actual curve number in the mindow
	      jj = 8*JC-6
	      if (jj .gt. 74)	goto	5
	      if (jj .ge. 66)	JN1 = JN1-8*DXLET
	      if (jj .eq. 74)	JN1 = JN1-8*DXLET
C                   curve #, win #
C	      write(*,'(10I4)')jw,j1,jc,JN1,jj,JNC
	      call	FMTXF6(XF7,TOUT(j,j1))
	      STRI (jj:jj+6) = XF7
	      STRIN(jj:jj+6) = '  '//NAMET(J1)//' '
	      JL = max(JL,jj+6)
	      JN1 = JN1+8*DXLET
	      call	colovm(JC+2)
C	      write(*,*)'STRI:   "',STRI(jj:jj+6),'"'
C	      write(*,*)'STTRIN: "',STRIN(jj:jj+6),'"'
	      call	textvm(JN1,JN2,STRI(jj:),7)
	      call	textvm(JN1,JN2-DYLET+1,STRIN(jj:),7)
 5	   continue

	   call	colovm(2)
	   STRI(1:5)   = 'Time='
	   STRI(11:11) = 's'
	   call	FMTF50(STRI(6:10),YY1)
	   call	textvm(DXLET,JN2-3*DYLET+DYLET/2,STRI,11)
C	   write(*,*)STRI(1:11)
	   return
	endif
	call	FMTF5(STRI(14:18),YX)
	call	FMTF5(STRI(20:24),YY)
	call	colovm(3)
	goto	7
	entry	ERASXY
	JN0 = 0
	JN2 = YSCMAX+15+3*DYLET
 6	write(STRI,'(20A4)')('    ',j=1,20)
C	write(*,*)jn0,jlr,XWW-83*DXLET+2,DXLET,JLR-4*DYLET+DYLET/3
	call	colovm(0)
 7	call	textvm(JN0,JLR,"               ",15)
	call	textvm(XWW-83*DXLET+2,JLR,STRI(1:80),80)
C	if (MOD10 .eq. 6)
C     >	call	textvm(DXLET,JN2-3*DYLET+DYLET/2,"              ",14)
	end
C======================================================================|
	subroutine	BRZ(YR,YZ,YBR,YBZ)
C----------------------------------------------------------------------|
C Input:	{YR,YZ}={r,z} polar coordinates
C Output:	{YBR,YBZ}={B_r,B_z} magnetic field components
C Warning:	This subroutine cannot be used in the postview mode
C		because of too low accuracy of stored data
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	integer	jx,jy,j1,j3
	double precision
     1		YR,YZ,YA,YD,YD1,YE,YE1,YT,YT1,YBR,YBZ,YY,YS,YC,Y1,Y2,Y3,
     2		RZ2A,QUADIN
	YA = RZ2A(YR,YZ,NAB)
	if (YA .lt. 1.E-5 .or. YA.ge.ABC)	then
	   YBR = 0.
	   YBZ = 0.
	   return
	endif
	YD = QUADIN(NAB,AMETR,SHIF,YA,YD1,j1)
	YE = QUADIN(NAB,AMETR,ELON,YA,YE1,j1)
	YT = QUADIN(NAB,AMETR,TRIA,YA,YT1,j1)
	YS = (YZ-UPDWN)/(YA*YE)
	Y2 = YS*YS
	YC = sqrt(1.-Y2)			! YC & YY are two
	YY = (YR-RTOR-YD)/YA+YT*Y2		! equivalent definitions
	YC = sign(YC,YY)			! for cos(theta)
	YD =YE*(1.+YD1*YC)+YA*YE1*Y2+YC*Y2*(YE*YT+YA*(2.*YE1*YT-YE*YT1))
	Y1 = QUADIN(NAB,AMETR,MU,   YA,YY,j1)
	Y2 = QUADIN(NAB,AMETR,RHO,  YA,YY,j1)
	Y3 = QUADIN(NAB,AMETR,DRODA,YA,YY,j1)
	YBZ = BTOR*Y1*Y2*Y3/(YR*YD)
	YBR = YBZ*(1.+2.*YT*YC)*YS
	YBZ =-YBZ*YE*YC
C	Y2 = sqrt(YBR*YBR+YBZ*YBZ)
C	write(*,'(4(3F10.5,2X))')YA,YBR,YBZ,Y2,YS,YC
C	J1 = 75*YBR
C	J3 = 75*YBZ
C	call	drawvm(0,JX,JY,JX+J1,JY-J3)
	end
C======================================================================|
	subroutine	SETFRM(IM,XSC0,XSC,NX,NY)
C----------------------------------------------------------------------|
C Input:	MODEY
C Output:	XSC0	- x_left  of the graphic area
C		XSC	- x_right of the graphic area
C		NX	- for using in PUTXY
C		NY	- for using in PUTXY
C		IM	- for using in AFRAME
C		NST	- for using in AFRAME
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	integer	IM,XSC0,XSC, STRUP,STRDN,NX,NY
C IY0,IYM     - upper & lower plot boundaries
C DXM,DYM     - X & Y window dimension
C STRUP,STRDN - number of text strings up & down
C NX,NY       - number of windows in horisontal & vertical axises
C XSC0,XSC    - left & right plot boundaries
	if (NST .ge. 1)	goto	7
	if (IM.eq.11 .or. IM.eq.12 .or. IM.eq.21 .or. IM.eq.22)	goto 7
	if (MOD10 .le. 1 .or. MOD10 .ge. 7)	then
			      IM = MOD10
	elseif(MOD10.ge.2 .and. MOD10.le.5)	then
	   if (MODEY .eq. 1)  IM = 2
	   if (MODEY .eq. -1) IM = 3
	else
	   if (MODEY .eq. 1)  IM = 5
	   if (MODEY .eq. 0)  IM = 9
	   if (MODEY .eq. -1) IM = 6
	endif
	if (IM .le. 0 .or. IM .ge. 10)	then
		IM = 0
		return
	endif
			go to (1,2,3,2,5,6,2,8,9)IM
C mode # 1
 1	STRUP=2
	STRDN=5
	XSC0=0
	XSC=XSCMAX
	NX=4
	NY=2
	IYM=256
			goto 10
C modes 2,3,4,5 at y-mode = +1, dummy mode (IM=4 - not used)
 2	STRUP=2
	STRDN=4
	XSC0=0
	XSC=XSCMAX
	NX=2
	NY=1
	IYM=270
			goto 10
C modes 2,3,4,5 at y-mode = -1
 3	STRUP=2
	STRDN=4
	XSC0=0
	XSC=XSCMAX
	NX=2
	NY=2
	IYM=270
			goto 10
C mode # 6 (time) at y-mode=1 (2 windows)
 5	STRUP=1
	STRDN=1
	XSC0=6*DXLET
	XSC=XSCMAX
	NX=1
	NY=2
	IYM=320
			goto 10
C mode # 6 (time) at y-mode=-1 (4 windows)
 6	STRUP=1
	STRDN=1
	XSC0=6*DXLET
	XSC=XSCMAX
	NX=1
	NY=4
	IYM=320
			goto 10
 7	STRUP=2
	STRDN=4
	NST=NST+1
	XSC0=0
	XSC=XSCMAX
	NX=1
	NY=1
	if(NST.eq.1)THEN
	   XSC=XSCMAX/2.
	   if(IM.eq.21.or.IM.eq.22)THEN
	      NX=2
	      NY=2
	   endif
	endif
	if(NST.eq.2)THEN
	   XSC0=XSCMAX/2.
	   if(IM.eq.12.or.IM.eq.22)THEN
	      NX=2
	      NY=2
	   endif
	endif
	IYM=320
			goto 10
C mode 8 (equilibrium)
 8	STRUP=1
	STRDN=-2
	XSC=575
	XSC0=0
	NX=1
	NY=1
	IYM=360
			goto 10
C mode # 6 (time) at y-mode=0 (1 window), mode 9 (user's plot)
 9	STRUP=1
	STRDN=1
	XSC0=6*DXLET
	XSC=XSCMAX
	NX=1
	NY=1
	IYM=320

 10	IY0=STRUP*DYLET+1
	if (MOD10 .eq. 6)	IY0=IY0+1
	IYM=min(IY0+IYM,YSCMAX-STRDN*DYLET-1)
C	write(*,*)IY0,IYM,YSCMAX,YSCMAX-STRDN*DYLET-1
	DXM=(XSC-XSC0)/NX
	DYM=(IYM-IY0)/NY
	end
C======================================================================|
	subroutine	AFRAME(IM,XSC0,XSC,NX,NY)
C Subroutine draw frame for different modes
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	integer	IM, JJ, J, JN0, TIMWIN, lonlen, JX, JY,
     1		LENG, XP, XM, YP, YM, XSC, XSC0, NX, NY, JXSCM, LYM
	double precision	DY, YY, TIND
	character CH6*6,XF4*5,COMMENT*80
	data	LENG/3/JN0/0/
C IY0,IYM     - upper & lower grafic boundary
C DXM,DYM     - X & Y window dimension
C IDX,DY      - X & Y distance (in points) between X & Y axis labels
C IDT	      - distance (in labels) between longer labels in modes 6&8
C LENG        - label length
	NST = 0
	if (IM .le. 0)	return
	TIMWIN = 0
	if (IM.eq.5 .or. IM.eq.6 .or. IM.eq.8 .or. IM.eq.9)  TIMWIN = 1

 10	call 	colovm(1)
	call 	rectvm(0,JN0,JN0,XWW-1,XWH-1)
	if (NST.ne.2 .and. MOD10.eq.6)		then
	    call	colovm(1)
	    j = XWW-20*DXLET+1
	    jj = YSCMAX+DYLET
	    call	textvm(j,jj,'time, s',7)
	endif
C Vertical lines & Y-labels
	if (KPRI.ge.1 .and. KPRI.le.2)	then
	    write(COMMENT,'(A)')"Vertical lines"
	    j = lonlen(COMMENT)
	    call	pscom(COMMENT,j)
	endif
	call 	colovm(1)
C	write(*,*)
C	write(*,*)XSC0,XSC,DXM,IY0,IYM
	do	11	JJ=XSC0,XSC,DXM
	   JX=MIN0(XSC,JJ)
	   XP=MIN(XSC,JX+LENG)
	   XM=MAX(XSC0,JX-LENG)
	   call	drawvm(0,JX,IY0,JX,IYM)
C	   write(*,*)"X-coord",JX,"   V-lines",IY0," -> ",IYM
C Y-line labels
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(COMMENT,'(A)')"Y-line labels"
	      j = lonlen(COMMENT)
	      call	pscom(COMMENT,j)
	   endif
	   DY = (IYM-IY0)/20.
	   YY=dble(IYM)
 12	   continue
	      JY=YY
	      if (IM .ne. 8)	call 	drawvm(0,XM,JY,XP,JY)
C	      write(*,*)"x-label:",XM," ->",XP," @y=",JY
	      YY=YY-DY
	      if (sign(1.d0,DY)*(YY-dble(IY0)) .ge. 0.d0)	goto	12
C	   do	11  YY=dble(IYM),dble(IY0),-DY
 11	continue
C Horizontal lines
	if (KPRI.ge.1 .and. KPRI.le.2)	then
	    write(COMMENT,'(A)')"Horizontal lines"
	    j = lonlen(COMMENT)
	    call	pscom(COMMENT,j)
	endif
	do	14	JY=IYM,IY0,-DYM
	   YM=MAX(IY0,JY-LENG)
	   if (IM.eq.4 .or. IM.eq.5 .or. IM.eq.6)	then
		YP=JY
	   else
		YP=MIN0(IYM,JY+LENG)
	   endif
	   call	drawvm(0,XSC0,JY,XSC,JY)
C	write(*,*)"H-lines:",XSC0," ->",XSC,"   Y-coord",JY
C X-line labels
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(COMMENT,'(A)')"X-line labels"
	      j = lonlen(COMMENT)
	      call	pscom(COMMENT,j)
	   endif
	   JX = XSC0
	   IDX = 16
	   JXSCM = XSC-IDX
	   if (TIMWIN .eq. 1)	then
	      IDX = 23
	      JXSCM = XSC
	   endif
	   do	J=1,100
	      JX = JX+IDX
	      if (JX .gt. JXSCM)	goto	14
	      if (TIMWIN .eq. 1 .and. J/IDT*IDT .eq. J)	then
	         LYM = YM-2
	      else
	         LYM = YM
	      endif
	      if (LYM .gt. IY0+LENG)  call 	drawvm(0,JX,LYM,JX,YP)
C	      write(*,*)"y-label:",LYM," ->",YP," @x=",JX
	   enddo
 14	continue
	if (IM .eq. 8)	then
	   JY = (IY0+IYM)/2
	   do	j=0,10
	      jj = JY+IDX*j
	      if (jj .lt. IYM)	call 	drawvm(0,XM,JJ,XP,JJ)
	      jj = JY-IDX*j
	      if (jj .gt. IY0)	call 	drawvm(0,XM,JJ,XP,JJ)
	   enddo
	endif
	if (NST .eq. 1)	then
	   call	SETFRM(IM,XSC0,XSC,NX,NY)
	   goto	10
	endif
	if (MOD10 .ne. 6)	goto	43

C time-axis legend:
	call colovm(1)
C	XSC0 = 6*DXLET
	JJ = YSCMAX-DYLET+12
	do  42	J=0,XSCMAX,IDX
	    JX = (J-IDT)*IDT+XSC0
	    if (JX.gt.XSCMAX)	goto	42
C (right_label_pos)/(n_labels)=575/IDT=115
	    YY = abs(TSCALE)
	    TIND = TINIT+J*YY/115
	    if ( TINIT+YY.gt.10.0 .or.
     +	         TINIT+YY.gt.1.0 .and. YY.lt.0.1
     +	                          .or. YY.lt.0.01)	then
	       call FMTXF5(CH6,TIND)
	       call textvm(JX-2,JJ,CH6,6)
	    else
	       call FMTXF4(XF4,TIND)
	       call textvm(JX,JJ,XF4,5)
	    endif
 42	continue

C equilibrium plot
 43	if (MOD10 .ne. 8)	goto	45
C horizontal axis labels
	JJ = IYM+DYLET+2
	call colovm(1)
	do  44	J=1,5
	    JX = IDX*IDT*J-24
	    if (JX .gt. JXSCM)	goto	44
	    YY = J*SCM
	    call FMTXF4(XF4,YY)
	    call textvm(JX,JJ,XF4,5)
 44	continue
C vertical axis labels
	do	J=-1,1
	    JX=(IYM+IY0+DYLET)/2+IDT*IDX*J
	    YY=-J*SCM
	    call FMTXF4(XF4,YY)
	    call textvm(XSC+2,JX,XF4,5)
	enddo
C	call	colovm(1)			! Black
C	JJ = (RTOR+SHIF(1))*YSC8
C	JJ = (RTOR+SHIF(1))*YSC8
C	call	drawvm(0,JN0,YS0,JJ,YS0)	! Plot the mid-plane
 45	if (KPRI.ge.1 .and. KPRI.le.2)	then
	    write(COMMENT,'(A)')"Frame done"
	    j = lonlen(COMMENT)
	    call	pscom(COMMENT,j)
	endif
	end
C======================================================================|
	subroutine NEGA(J1,ITIMES,TOUT)
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	integer J1,ITIMES,NUM(4),JMIN,JMAX,NKL1,NKL2,MODK,IBEG,JJ,J
	double precision  TOUT(ITIMES,NRW)
	common	/AC_NEG1/ NUM,JMIN,JMAX,NKL1,NKL2,MODK(2)
	save IBEG
	data IBEG/0/
	IBEG=IBEG+1
	J1=11
	if (IBEG .lt. 3)	return
	do 	2	JJ=NKL1,NKL2
	   if(JJ.le.0)	goto 2
	   MODK(JJ)=0
	   do 	1	J=JMIN,JMAX
	      if (TOUT(J,NUM(2*JJ-1)) .ge. 0.0 .and.
     .	          TOUT(J,NUM(2*JJ))   .ge. 0.0) goto	1
	         MODK(JJ)=1
	         goto 2
 1	   continue
 2	continue
C	(0,0) - 11,	(-1,0) - 21,	(0,-1) - 12,	(-1,-1) - 22
	if(MODK(1).eq.0)	then
	   if(MODK(2).eq.0)	then
	      J1=11
	   else
	      J1=12
	   endif
	else
	   if(MODK(2).eq.0)	then
	      J1=21
	   else
	      J1=22
	   endif
	endif
	end
C======================================================================|
C Time dependences for radial output
 	subroutine	DNSTR(jt,ITIMES,TOUT)
C----------------------------------------------------------------------|
Curve to digit conversion
C  Input:  jt defines the current time
C
C    if mod10 != 6 or call from run then jt = LTOUT  
C
C  Data used from outcmn.inc
C	double precision	TOUT(1,1)
C	integer	LTOUT,YSCMAX,DXLET,DYLET,MOD10,NSCR(*),NTOUT
C	integer	NWIND3(*)
C	character*4 NAMET(*)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	integer jt,ITIMES,JN2,JN0,JEND,JB,FSHIFT,JL,JC,JW,JJ,J
	double precision	TOUT(ITIMES,*)
	character STRI*80,STRIN*80,XF4*5,XF7*7
	save	FSHIFT
	data	JN0/0/FSHIFT/10/
	if (jt .eq. 0)	then
	   jt = LTOUT
	   call	colovm(1)
	else
	   call	colovm(3)
	endif
	if (MOD10 .ne. 6)	goto	2
	JN0 = 6*DXLET
	JN2 = YSCMAX+FSHIFT+3*DYLET+5
	JC = 0
	JL = 0
	write(STRIN,'(79X,1A1)')' '
	write(STRI,'(79X,1A1)')' '
	do	1  j=1,NTOUT
	    JW=NWIND3(j)-8*NSCR(MOD10)
	    if (NAMET(j) .eq. '    ') JW = 0
	    if (JW.le.0 .or. JW.gt.8)	goto	1
	    JC = JC+1			! Actual curve number in the mindow
	    jj = 8*JC-6
	    if (jj .gt. 74)	goto	1
	    if (jj .ge. 66)	JN0 = 5*DXLET
	    if (jj .eq. 74)	JN0 = -DXLET
C                   curve #, win #, chan #, mode 6, screen #
C	write(*,'(10I4)')jw, j,  NWIND3(j), MOD10, NSCR(MOD10),jj
	    call FMTXF6(XF7,TOUT(jt,j))
	    STRI (jj:jj+6) = XF7
	    STRIN(jj:jj+6) = '  '//NAMET(J)//' '
	    JL = max(JL,jj+6)
 1	continue
	call textvm(JN0,JN2,STRI,JL)
	JN2=JN2-DYLET+1
	call textvm(JN0,JN2,STRIN,JL)
	return

 2	JB = 1
	JN0 = 0
	JN2 = YSCMAX-5*DYLET+FSHIFT+2
 3	JEND=MIN0(JB+15,NTOUT)
	do	4	J=JB,JEND
	call FMTXF4(XF4,TOUT(jt,J))
	if(NAMET(J).eq.' ')	XF4 = '    '
	JJ=5*(J-JB+1)-4
	STRI(JJ:JJ+4)=XF4
 4	continue
	JN2=JN2+2*DYLET+2
	JL = 5*(JEND-JB+1)
	call textvm(JN0,JN2,STRI,JL)
	write(STRI,'(16(1X,1A4))')	(NAMET(J),J=JB,JEND)
	JN2=JN2-DYLET+1
	call textvm(JN0,JN2,STRI,JL)
	JN2=JN2+DYLET-1
	if(JEND.eq.NTOUT.or.JEND.eq.NRW)	return
	JB=JB+16
! Temporary restriction to limit number of lines.
! Remove when the window size is increased.
C	if(JB.gt.64)	return
	goto	3
	end
C======================================================================|
C Upper string of a picture
 	subroutine	UPSTR(YN,YQ)
	implicit none
	double precision	YN,YQ
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer JN0,FSHIFT
	save	FSHIFT
	character STRMN*62,CHR*2
	data	JN0/0/FSHIFT/2/
	data 	STRMN/' R=     a=     B=     I=     q=     n=     '/
	call	FMTF40(STRMN(4:7),RTOR)
	call	FMTF40(STRMN(11:14),ABC)
	call	FMTF40(STRMN(18:21),BTOR)
	call	FMTF40(STRMN(25:28),IPL)
	call	FMTF40(STRMN(32:35),YQ)
	call	FMTF40(STRMN(39:42),YN)
	call	colovm(1)
	call	textvm(JN0,FSHIFT,XLINE1(2:20)//STRMN,61)
	call	textvm(JN0+55*DXLET,FSHIFT-DYLET+4,'_',1)
	call	colovm(3)
	write(CHR,'(1I2)')NSCR(MOD10)+1
	call	textvm(JN0+79*DXLET,FSHIFT+DYLET-2,CHR,2)	! Screen No.
C	call	textvm(JN0+79*DXLET,426-4*DYLET+DYLET/3,CHR,2)
	call 	colovm(1)
	call 	rectvm(0,JN0,JN0,XWW-1,XWH-1)
	end
C======================================================================|
 	subroutine	STROUT(STRI,LENG)
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	integer	 JLR,LENG,WarningColor
	character*1	STRI(LENG)
	save	WarningColor
	data	WarningColor/30/	JLR/426/
	if (TASK(1:3).eq.'BGD' .or. TASK(4:4).eq.'B')	return
	call	colovm(WarningColor)
	call	textvm(XWW-83*DXLET+2,JLR,STRI,LENG)
	end
C======================================================================|
 	subroutine	DSPMSG(STRI)
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	integer FSHIFT,J,lonlen,j1,WarningColor
	character	STRI*(*)
	save	FSHIFT,		WarningColor
	data	FSHIFT/2/	WarningColor/30/
	if (TASK(1:3).eq.'BGD' .or. TASK(4:4).eq.'B')	return
	J = min(19,lonlen(STRI))
	if (j.lt.19)write(STRI(J+1:),'(19A1)')(' ',j1=j+1,19)
	call	colovm(WarningColor)
	call	textvm(62*DXLET,FSHIFT,STRI,19)
	end
C======================================================================|
 	subroutine	TIMEDT(TIME,DT)
	implicit none
	double precision	TIME,DT
	integer FSHIFT,DXLET
	character*19	STRI
	DXLET = 8
	FSHIFT = 2
	STRI(1:5) = 'Time='
	STRI(11:16) = ' dt='
	call	FMTF50(STRI(6:10),TIME)
	call	FMTF50(STRI(15:19),DT)
	call	colovm(1)
	call	textvm(62*DXLET,FSHIFT,STRI,19)
	end
C=======================================================================
C Appending the list of constants to a PS file
	subroutine	PCOVA
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	character*6	CONN(22), STRI*80, CH6, TZER*1
	integer J,J1,J2,JJ,JDUM,JNY,JNX
	data	CONN/	 'CF1-> ','CF5-> ','CF9-> ','CF13->',
     d		'CV1-> ','CV5-> ','CV9-> ','CV13->',
     d		'CHE   ','CHI   ','CNB   ','CNBI  ',
     d		'CCD   ','CRF   ','CNEUT ','CPEL  ',
     d		'CBND  ','CFUS  ','CIMP  ','CMHD  ',
     d		'CRAD  ','CSOL  '/
C Writing constants
 	J1=0
	TZER = char(J1)
	JNY=540
	JNX=10
	do	93 J2=1,11
	JDUM=2*(J2-1)+1
	STRI=CONN(JDUM)
	do	91 J=1,4
	J1=J1+1
	if(J1.GT.NCFNAM)	go to 94
	call FMXF5(CH6,CONSTF(J1))
	JJ=7*(J-1)+8
 91	STRI(JJ:JJ+5)=CH6
	STRI(JJ+12:JJ+18)=CONN(JDUM+1)
	do	92 J=5,8
	J1=J1+1
	if(J1.GT.NCFNAM)	go to 94
	call FMXF5(CH6,CONSTF(J1))
	JJ=7*(J-1)+20
 92	STRI(JJ:JJ+5)=CH6
	JNY=JNY+17
 93	call textvm(JNX,JNY,STRI,75)
 94	continue
C Writing variables
	j1 = 1
	JNY=JNY+20
	JDUM = NPRNAM-48
	do	98	j=1,JDUM
	call FMXF5(CH6,DEVAR(j))
	STRI(j1:j1+19) = PRNAME(j)//'='//CH6//'     '
	j1 = j1+20
	if (j1.gt.70 .or. j.eq.JDUM)	then
		STRI(j1-5:) = TZER
C		write(*,*)STRI(1:j1-3)
		JNY=JNY+17
		call textvm(JNX,JNY,STRI,j1-5)
		j1 = 1
	endif
 98	continue
	end
C===================================================================
C	DEVID	- device ID e.g. "AUGD"
C	ARUNID	- Astra run ID (the lowest string of the graphic output)
C	UFNAME	- U-file name to write
C	TIME	- current time value
C	SIGNAM	- signal name
C	NPOINT	- number of points
C	NASC	- number of associated scalar quantities
C	PROCOD	- processing code
C	ARGUM	- argument values
C	SIGNAL	- function values
C	RTOR,AB,BTOR,IPL - standard Astra notations
C	YN	- average density
C	MODADD	- append the model file
C-----------------------------------------------------------------------
	subroutine	UF1DWA(DEVID,ARUNID,UFNAME,TIME,SIGNAM
     &		,NPOINT,NASC,PROCOD,ARGUM,SIGNAL
     &		,RTOR,AB,BTOR,IPL,YN,MODADD,MODEX)
	implicit none
	integer	NPOINT,NASC,PROCOD,MODEX,NSHOT,LUF,NCHU,length,j,LWRN
	double precision
     1		TIME,ARGUM(NPOINT),SIGNAL(NPOINT),YN,IPL,RTOR,AB,BTOR
	logical MODADD,THERE
	character*40  UFNAME
	character  DEVID*4,SIGNAM*4,ARUNID*80,STRING*132
	NSHOT = 0
	NCHU = 1
C	read(UFNAME(2:6),*)NSHOT

	LUF = length(UFNAME)
	open(NCHU,file='udb/'//UFNAME(1:LUF),status='UNKNOWN')
C	write(*,*)'"',UFNAME,'"',LUF
	write(NCHU,1001)NSHOT,DEVID,1,0
	write(NCHU,1002)
	write(NCHU,1003)NASC
	if (NASC .eq. 0)	write(NCHU,1005)
	if (NASC .eq. 1)	then
		write(NCHU,1004)TIME
		write(NCHU,1006)
		if (MODEX.eq.0 .or. MODEX.eq.1)	write(NCHU,1007)
		if (MODEX.eq.2)	write(NCHU,1010)
		if (MODEX.eq.3)	write(NCHU,1013)
	endif
	if (NASC .gt. 1 .or. MODEX .gt. 3)	then
	   write(NCHU,*)"Warning: don't know how to write 1D U-file"
	   write(*,*)"Warning: don't know how to write 1D U-file"
	endif
	write(NCHU,1008)SIGNAM
	write(NCHU,1009)PROCOD
	if (NASC .eq. 0)	write(NCHU,1011)NPOINT
	if (NASC .eq. 1)	write(NCHU,1012)NPOINT
	write(NCHU,'(1X,1P,6E13.5)')(ARGUM(J),J=1,NPOINT)
	write(NCHU,'(1X,1P,6E13.5)')(SIGNAL(J),J=1,NPOINT)
	write(NCHU,*)
     +	    ' ;----END-OF-DATA---------------COMMENTS:-----------'
	write(NCHU,'(1A80)')ARUNID
	write(NCHU,'(1A4,1F6.2,3(1A8,1F6.2),1A13,1F6.2,1A7)')
     +	    ' R =',RTOR,'m,   a =',AB,'m,   B =',BTOR,
     +	    'T,   I =',IPL,'MA,   <n_e> =',.1*YN,'E20m^-3'
	if(.TRUE.)	goto 219
C	if(.not.MODADD)	goto 219
C The model listing can be optionally appended to the U-file:
	LWRN	= length(ARUNID(69:77))
	write(NCHU,*)'>>>  Model  "',ARUNID(69:68+LWRN),'"  listing:'
	inquire(FILE='equ/'//ARUNID(69:68+LWRN),EXIST=THERE)
	if(.not.THERE)   then
		write(*,*) '>>> U-file write error:  Model "equ/',
     +		ARUNID(69:68+LWRN),'" not found'
		goto   219
	endif
	STRING="cat "//"equ/"//ARUNID(69:68+LWRN)//
     >					" >> udb/"//UFNAME(1:LUF)
	call	system(STRING)
 219	write(*,*)'>>>  Dataset "',SIGNAM,'" is written in the 1D '
     &	    //'U-file "udb/'//UFNAME(1:LUF),'"'
	close(NCHU)
	return
 1001	format(1I7.5,1A4,2I2,16X,'; Shot #  ')
 1002	format(1X,9HDD-MMM-YY,21X,
     +		'; Shot date -  Ufiles ASCII file system')
 1003	format(1I4,27X,'; # of associated scalars')
 1004	format(1P,1E13.6,18X,'; Scalar value, LABEL FOLLOWS:')
 1005	format(' Time:',15X,'sec',7X,'; Independent variable: time')
 1006	format(' Time:',15X,'sec',7X,';')
 1007	format(' MINOR RAD',11X,'m',9X,'; Independent variable: a')
 1010	format
     +	(' Rho toroidal, normalized',6X,'; Independent variable: Rho')
 1013	format
     +	(' Poloidal flux, normalized',5X,'; Independent variable: Psi')
 1008	format( ' Function            ??        ;',
     +		' Function name:   "',1A4,'"')
 1009	format(1I2,29X,'; Processing code: ..., 4 - Astra output')
 1011	format(1I11,20X,';-# of  time  pts  t, F(t) (data follow)')
 1012	format(1I11,20X,';-# of radial pts  X, F(X) (data follow)')
	end
C=======================================================================
	subroutine	UF2DWA(DEVID,UFNAME,SIGNAM,JN,NRP,NASC,PROCOD,
     &		PRMARK,TIMOD4,RTOR,AB,BTOR,IPL,YN,MODADD,YWA,YWB,YWC)
C-----------------------------------------------------------------------
C Note:	2D U-file is written with the radial number of points NRP
C	as mapped from the full grid size NB1 because NA1 can vary in
C	time and is not suitable for the U-file fixed grid
C  YWA(*), YWB(*), YWC(*) working arrays 
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	integer	j,j1,j2,JN,NTP,NRP,NASC,PROCOD,IERR,JNT,NSHOT,NCHU,NCHR
	integer	SKIPM,JAB,I,jj,LWRN,length
	double precision
     1		PRMARK(*),TIMOD4(*),YN,IPL,YWA(*),YWB(*),YWC(*),
     2		RTOR,AB,BTOR,TEMPR,YTIME,SCL,DOWN,ALFA,SIGNAL(NRD)
C	byte	JNT2(NRD)
	integer*2	JNT2(NRD)
	logical MODADD,THERE
	character*4  CHAR4,DEVID,SIGNAM,UFNAME*10,STRING*132
C NTP - total number of time slices: 
	NTP = 0
	do	j = 1,IPOUT-1
	   if (PRMARK(j) .ge. 0)	NTP = NTP+1
	enddo
	NSHOT = 0
	NCHU = 1
	NCHR = 2
C	read(UFNAME(2:6),*)NSHOT
	open(NCHU,file='udb/'//UFNAME,status='UNKNOWN')
	write(NCHU,1001)NSHOT,DEVID,2,0
	write(NCHU,1002)
	write(NCHU,1003)NASC
	if (NASC .ne. 0)	then
	   write(NCHU,*)"Warning: don't know how to write 2D U-file"
	   write(*,*)"Warning: don't know how to write 2D U-file"
	   write(*,*)"         associated scalar is not defined"
	endif
C NOTE! a radially contiguous u-file is written
	write(NCHU,1007)
	write(NCHU,1005)
	write(NCHU,1008)SIGNAM
	write(NCHU,1009)PROCOD
	write(NCHU,1012)NRP
	write(NCHU,1011)NTP
	do	j=1,NRP
	   YWC(j)	=(j-1.)/(NRP-1.)
	enddo
	write(NCHU,'(1X,1P,6E13.5)')(YWC(J)*AB,J=1,NRP)
	j1 = 0
	j2 = 0
	do	30	j = 1,NTP
	   if (PRMARK(j) .lt. 0)	goto	30
	   j1 = j1+1
	   j2 = j2+1
	   if (j2 .eq. 1)	then
	      write(NCHU,'(1X,1P,E13.5,$)')TIMOD4(j1)
	   elseif (j2 .eq. 6)	then
	      write(NCHU,'(1P,E13.5)')TIMOD4(j1)
	      j2 = 0
	   else
	      write(NCHU,'(1P,E13.5,$)')TIMOD4(j1)
	   endif
 30	continue
	if (j2 .ne. 0)write(NCHU,*)

C	write(NCHU,'(1X,1P,6E13.5)')(ARGT(J),J=1,NTP)

	call 	OPENRD(NCHR,RSNAME,1,IERR)
	if(IERR.gt.0)	pause ' >>> Profiles file error'
	j = 0
	j = SKIPM(NCHR,j)		! Returned value is not used
	read(NCHR,ERR=38,END=39)CHAR4
	if (NXOUT.gt.0 .and. NGR.gt.0) read(NCHR,ERR=38,END=39)JNT
C	read(NCHR,ERR=38,END=39)(CHAR4,j=1,61),(JNT,j=1,8)
C     .		,(CHAR4,J=1,NROUT),(TEMPR,I=1,NROUT),JNT
C     .		,(CHAR4,I=1,NTOUT),(TEMPR,I=1,NTOUT+1)
C     .		,(JNT,I=1,11)
C		if (NXOUT .gt. 0 .and. NGR .gt. 0)	then
C		read(NCHR,ERR=38,END=39)(JNT,j=1,3*NGR)
C     .	,	(TEMPR,j=1,3*NGR+GDEY(NGR)+NGRIDX(NGR)-1)
C     .	,	(JNT,j=1,3*NARRX)
C		endif
	do  215  J2 = 1,IPOUT-1
	    read(NCHR,ERR=38,END=39)JNT
	    if(JNT.ne.0) read(NCHR,ERR=38) TEMPR
	    read(NCHR,END=39)YTIME
C  Skip  CONSTF, DEVAR, LINEAV, ROC
	    read(NCHR,END=39)TEMPR
C	    if(JNT.eq.0)	goto 32
C	    read(NCHR,ERR=38)((TEMPR,I=1,NTOUT+1),J=1,JNT)
C 32	    read(NCHR,END=39)YTIME
C  Skip  CONSTF, DEVAR, LINEAV, ROC
C	    read(NCHR,END=39)(TEMPR,I=1,NCFNAM+NPRNAM+3)
	    read(NCHR,END=39)JAB,JAB,(JNT,I=1,10),(TEMPR,I=1,10)
C  Retrieve  AMETR
	    read(NCHR)SCL,DOWN,(JNT2(jj),jj=1,JAB)
	    do	j=1,JAB
C		YWA(J) = (DOWN+SCL*(JNT2(J)+128)/255.)/AB
		YWA(J) = (DOWN+SCL*(JNT2(J)+32768)/65535.)/AB
	     enddo
	    YWA(JAB) = 1.
C  Skip  SHIF, ELON, TRIA, RHS1, RHS2, FP
	    do j=1,6
		read(NCHR)TEMPR,TEMPR,(JNT2(jj),jj=1,JAB)
	    enddo
	    do  37	J1=1,NROUT
		read(NCHR,ERR=38)SCL,DOWN,(JNT2(jj),jj=1,JAB)
		if(j1.ne.JN .or. PRMARK(j2).lt.0)  goto 37
		do	j=1,JAB
C		    SIGNAL(J)=DOWN+SCL*(JNT2(J)+128)/255.
		    SIGNAL(J)=DOWN+SCL*(JNT2(J)+32768)/65535.
		enddo
		ALFA = .0001
		CALL	SMOOTH(ALFA,JAB,SIGNAL,YWA,NRP,YWB,YWC)
		write(NCHU,'(1X,1P,6E13.5)')(YWB(J),J=1,NRP)
 37	    continue
 215    continue
 39	close(NCHR)

	write(NCHU,*)
     +		' ;----END-OF-DATA---------------COMMENTS:-----------'
	write(NCHU,'(1A80)')RUNID
	write(NCHU,'(1A4,1F6.2,3(1A8,1F6.2),1A13,1F6.2,1A7)')
     +	    ' R =',RTOR,'m,   a =',AB,'m,   B =',BTOR,
     +	    'T,   I =',IPL,'MA,   <n_e> =',.1*YN,'E20m^-3'
	if(.TRUE.)	goto 219
C	if(.not.MODADD)	goto 219
C The model listing can be optionally appended to the U-file:
	LWRN	= length(RUNID(69:77))
	write(NCHU,*)'>>>  Model  "',RUNID(69:68+LWRN),'"  listing:'
	inquire(FILE='equ/'//RUNID(69:68+LWRN),EXIST=THERE)
	if(.not.THERE)   then
		write(*,*) '>>> U-file write error:  Model "equ/',
     +		RUNID(69:68+LWRN),'" not found'
		goto   219
	endif
	STRING="cat "//"equ/"//RUNID(69:68+LWRN)//" >> udb/"//UFNAME
	call	system(STRING)
 219	write(*,*)'>>>  Dataset "',SIGNAM,
     &		'" is written in the 2D U-file "udb/',UFNAME,'"'
	close(NCHU)
	return
 1001	format(1I7.5,1A4,2I2,16X,'; Shot #  ')
 1002	format(1X,9HDD-MMM-YY,21X,
     +		'; Shot date -  Ufiles ASCII file system')
 1003	format(1I4,27X,'; # of associated scalars')
 1004	format(1P,1E13.6,18X,'; Scalar value, LABEL FOLLOWS:')
 1005	format(' Time:',15X,'sec',7X,'; Independent variable: time')
 1006	format(' Time:',15X,'sec',7X,';')
 1007	format(' MINOR RAD',11X,'m',9X,'; Independent variable: Rho')
 1008	format( ' Function            ??        ;',
     +		' Function name:   "',1A4,'"')
 1009	format(1I2,29X,'; Processing code: 4 - Astra output')
 1011	format(1I11,20X,';-# of  time  pts  t, F(t) (data follow)')
 1012	format(1I11,20X,';-# of radial pts  X, F(X) (data follow)')
 38	write(*,*)'Read file PROFIL.DAT error'
	pause
	end
C=======================================================================
C	This is Sun version only
C-----------------------------------------------------------------------
C	For IBM workstation changes are required to the subroutines
C	OPENAP (search for "IBM" below)
C=======================================================================
C	Subroutines	OPENRD,OPENWT 
C   Input:	NCHAN - channel number
C		FILNAM - file name (character variable or string CALL 
C		operator). Extention (name.ext) is necessary for string.
C		ITYPE - format status (0 - default, > 0 - 'UNFORMATTED')
C   Output:	IERR - output status (0 - okay, 1 - warning, 2 - error)
C----------------------------------------------------------------------|
	subroutine	OPENAP(NCHAN,FILNAM,ITYPE,IERR)
	implicit 	none
	integer		NCHAN,ITYPE,IERR
	character*(*)	FILNAM
	if(ITYPE.eq.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1,
C For Sun, Alpha:  access=		For IBM: position=
     .		ACCESS='APPEND')
C     .		POSITION='APPEND')

	if(ITYPE.gt.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1,
C For Sun, Alpha:  access=		For IBM: position=
     .		ACCESS='APPEND',FORM='UNFORMATTED')
C     .		POSITION='APPEND',FORM='UNFORMATTED')
	IERR=0
	return
1	IERR=2
	end
C=======================================================================
	subroutine	OPENRD(NCHAN,FILNAM,ITYPE,IERR)
	implicit 	none
	integer		NCHAN,ITYPE,IERR,length
	character*(*)	FILNAM
	if(ITYPE.eq.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1)
	if(ITYPE.gt.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1,
     *			FORM='UNFORMATTED')
	rewind(NCHAN)
	IERR=0
	return
1	IERR=2
        write(*,*)'>>> ERROR:  Cannot open file "',
     *		FILNAM(1:length(FILNAM)),'"'
	end
C=======================================================================
	subroutine	OPENWT(NCHAN,FILNAM,ITYPE,IERR)
	implicit 	none
	integer		NCHAN,ITYPE,IERR
	character*(*)	FILNAM
	if(ITYPE.eq.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1)
	if(ITYPE.gt.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1,
     *			FORM='UNFORMATTED')
	rewind(NCHAN)
	IERR=1
	return
1	continue
	if(ITYPE.eq.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='NEW',ERR=2)
	if(ITYPE.gt.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='NEW',ERR=2,
     *			FORM='UNFORMATTED')
	rewind(NCHAN)
	IERR=0
	return
2	IERR=2
	end
C=======================================================================
	integer function SKIPM(NCHR,NCHL)
C Skip model & model.log records
! A retrieving routine is in the file .srv/listres.f
C If NCHL =/= 0 then SCRATCH file containing model.log is writen
C Returns	(a returned value is used in the postviewer only)
C   0 for the old format (version before 5.1)
C   1 for the new format (version 5.2 and later)
	implicit none
	integer NCHR,NCHL,j,n
	character STRI*132,STR*32,CH1*1
	read(NCHR) STR
	if (STR.ne."^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")	then
	   rewind(NCHR)
	   SKIPM = 0
	   return	! -> Old format file
	endif
	SKIPM = 1	! -> New format file
	j = 0
	n = 1
	STR = "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
 10	read(NCHR) CH1,STRI(1:ichar(CH1))
	if (ichar(CH1).eq.32 .and. STR.eq.STRI(1:32)) then
	   j = j+1	! No. of adjacent strings 32*"^"
	   n = n+1	! Total No. of strings 32*"^" already encountered
	   if (j .eq. 2)	then
	      goto (11,12,13,14),n
 11	      continue	! n=1 and n=2 never happens
 12	      write(*,*)"Error in Post-view file: wrong format"
 13	      if (NCHL .eq. 0)	return
	      close(NCHL)
	      NCHL = 0
 14	      return
	   endif
	   if (n.eq.2 .and. NCHL.ne.0) open (NCHL,STATUS='SCRATCH')
	else
	   j = 0
	endif
	if (j.ne.0 .or. n.ne.2 .or. NCHL.eq.0)	 goto	10
	write(NCHL,'(A)')STRI(1:ichar(CH1))
C	write( *  ,'(A)')STRI(1:ichar(CH1))
	goto	10
	end
C======================================================================|
	subroutine	SETFNA(FNAME,NLEN)
C FNAME - input name (without blanks) is appended with an extension.
C	  The extension is the ordinal number of the file
C NLEN  - returns length of the new name
C
	implicit none
	integer	killbl,nvar,nlen,length,j
	character FNAME*(*)
	logical	EXI
	NLEN = length(FNAME)
	NVAR = 0
 1	NVAR = NVAR+1
	write(FNAME(NLEN+1:),'(1A1,1I3)')'.',NVAR
	j = killbl(FNAME,NLEN+4)
C	write(*,*)'Inquired:  "',FNAME(1:j),'"',j
	inquire(FILE=FNAME,EXIST=EXI)
	if(EXI)		goto 1
	NLEN = j
	return
	end
C======================================================================|
