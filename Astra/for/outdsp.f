C======================================================================|
C Modules   >>  OUTDSP, DRFOOT, SETARX <<
C----------------------------------------------------------------------|
 	subroutine	OUTDSP(MARK,JIFNEW,PT,IYO,ITIMES,TTOUT,TOUT)
C----------------------------------------------------------------------|
C Drawing options:
C    X - axis
C	0 <= rho <= ROC=RHO(NA1)
C	0 <= a <= ABC
C	0 <= a <= AB
C
C MARK = 1	Put marks
C MARK = 0	Solid lines
C MARK =-1	Dashed lines
C JIFNEW  =  0	Re-draw (erase) the previous curves
C JIFNEW =/= 0	New curves only
C JIFNEW  <  0	Don't mark resonances q=m/n
C JIFNEW > 10	Call from Review. (JIFNEW-10) is used to control erasing
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'

	integer ICVMX
	parameter (ICVMX=32)		! <= 32 curves in a screen
	integer MARK,JIFNEW,PT(*),ITIMES,IYO(ITIMES,*),
     1	        IXO(NRD,ICVMX),IX(2*NRD),PTM(2),PTMO(2,NRDX,NRW),
     2		IY(2*NRD),IWN(16),JWC(16),EQOLD(260,11),NUM(4),jk(5),
     3		IST,IQ1,IQ2,YS0,jnl,JS,NKL1,NKL2,MODEX,IMODE,SKIPM,
     4		IYM0,LTOUT1,LTOUT2,JFNEW,FSHIFT,STYL,JX,JY,jxout,
     5		JW,JNC,JNW,JW4,JWX4,JW8,JUP,MODK,INT2,jzero,jcv,JBOX,
     6		IY0GR,IYMN,IYMX,JDSP,JPOS,NPTMO(NRW),jlx(8),IWX(8),
     7		JMIN,JMAX,NP1O,NP1,NP,lonlen,length,jwx,ierr,i,jab,is,j,
     8		j1,jj,jsco,jn,jc,jpnt,jsc,jpnto,jposy,jarr,jtyp,
     9		WarningColor,EraseColor	
	double precision	TTOUT(ITIMES),TOUT(ITIMES,NRW),SC(NRW),
     1	 	YX,YROUT,YA,YL,YR,YQ1,YQ2,
     1          TEMPR,SCL,DOWN,YS,YXR,YXL,XROUT,YSC8,YFI,YZ,QUADIN,ABSC
	character	STRI*80,XF4*5,ST*9,CHAR4*4,CHAR6*6
	save	EQOLD,IXO,NP1O,PTMO,NPTMO,YQ1,YQ2,IQ1,IQ2,IWN,
     1		WarningColor,EraseColor
	common	/AC_NEG1/ NUM,JMIN,JMAX,NKL1,NKL2,MODK(2)
	data 	jzero/0/ FSHIFT/10/ NP1O/0/
!		 WarningColor/30/   EraseColor/31/
C----------------------------------------------------------------------|
	JFNEW = JIFNEW
	if (JFNEW .ge. 10)	JFNEW = JFNEW-10
	if (JFNEW .ne. 0)	then
	   do	J=1,16
	      IWN(J) = 0
	   enddo
	endif

C Suppress previous markloc called from RADOUT
	call	markloc("OUTDSP called from IFKEY"//char(0))

C Different options for abscissa. 
C Implemented for modes 0,1,2,3
	MODEX = XOUT+.49
	if     (MODEX .eq. 0)	then
C Draw up to AB against "a"
	   NP1 = NAB
	elseif (MODEX .eq. 1)	then
C Draw up to ABC against "a"
	   NP1 = NA1
	elseif (MODEX .eq. 2)	then
C Draw up to ROC against "rho"
	   NP1 = NA1
	elseif (MODEX .eq. 3)	then
C Draw up to ROC against "Psi"
	   NP1 = NA1
	else
C Unknown option
	   return
	endif
	if (MOD10 .eq. 3)	then
	   NP1 = NA1
	endif

	if(MARK .eq. 0)	IST=1
	if(MARK .eq. 1)	IST=8
	if(MARK .eq.-1)	IST=15
	IYMN=YSCMAX-IYM
	IYM0=IYM-DYM
	IYMX=YSCMAX-IY0
	JY = 3500
C----------------------------------------------------------------------|
C Mode 1:
	if (MOD10 .ne. 1)	goto	2
	call	add2loc("Drawing mode 1"//char(0))
	do	JJ=1,NROUT
	   do	J=1,NP1
		ROUT(j,jj) = ROUT(j,jj)+OSHIFR(jj)
	   enddo
	enddo
	call	add2loc("Calling SCAL from OUTDSP"//char(0))
	call	SCAL(NROUT,SC,SCALER,ROUT(4,1),NP1-3,NRD)
C	write(*,*)NROUT,ROUT(4,1),NP1-3,NRD
C	write(*,*)(SC(j),SCALER(j),j=1,NROUT)
	do	J=1,8
	   IWX(J) = 0
	enddo

C Plot curves
	call	add2loc("Loop 10 in OUTDSP"//char(0))
	jcv = 0
	do	10	jj=1,NROUT
	   JW = NWIND1(jj)-16*NSCR(MOD10)
	   if (NAMER(jj) .eq. '    ') JW = 0
	   if (JW.le.0 .or. JW.gt.16)	goto 10
	   JW4 = (JW-1)/4
	   JWX4 = JW-4*JW4-1
	   JW8 = (JW-1)/8
	   JUP = JW4-2*JW8
	   JX = DXM*JWX4
C  JW4 =0,0,0,0,1,1,1,1,2...  JW8=0,0,0,0,0,0,0,0,1,1...
C  JWX4=0,1,2,3,0,1,2,3,0...  JUP=0,0,0,0,1,1,1,1,0,0...
	   JBOX = JW-(JW-1)/8*8
	   YL = max(0.d0,ABSC(GRAL(jj)))
	   YR = ABSC(GRAP(jj))
	   if (YL.ge.YR)	YL = 0.d0
	   if (YL.lt.YR .and. 
     +		(YL.gt.AMETR(1)/ABC .or. YR.lt.0.999) )	then
C Draw a bar for an extended axis:
	      PT(1) = JX+160*YL
	      PT(2) = 350-(IYMN+128*(1-JUP))-2*IWX(JBOX)
	      PT(3) = JX+160*YR
	      call	colovm(JW8+2)
	      call	drawvm(0,PT(1),PT(2)-1,PT(3),PT(2)-1)
	      call	drawvm(0,PT(1),PT(2)-2,PT(3),PT(2)-2)
	      IWX(JBOX) = IWX(JBOX)+1
	   endif
	   jxout = 0
	   do	j=1,NP1
	      YX = ABSC(AMETR(j)) 
	      if (YX.lt.YL .or. YX.gt.YR)	goto	11
	      jxout = jxout+1
	      if (jxout.eq.1 .and. j.gt.1)	then ! left edge interpolation
		 IX(jxout) = JX
		 YA=ROUT(J,jj)+(ROUT(J-1,jj)-ROUT(J,jj))*(YL-YX)/(YA-YX)
		 YROUT = min(max(YA/SC(jj),-7.d0),7.d0)
		 JDSP  = 10*(DYM*YROUT+IYMN+128*(1-JUP))
		 IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
		 jxout = jxout+1
	      endif
	      IX(jxout) = JX+160*(YX-YL)/(YR-YL)
	      YROUT = min(max(ROUT(J,jj)/SC(jj),-7.d0),7.d0)
	      JDSP  = 10*(DYM*YROUT+IYMN+128*(1-JUP))
	      IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
 11	      continue
	      if (YX.le.YL .and. YX.gt.YR)	then ! right edge interpolation
		 jxout = jxout+1
		 IX(jxout) = JX+160
		 YA=ROUT(J,jj)+(ROUT(J-1,jj)-ROUT(J,jj))*(YR-YX)/(YA-YX)
		 YROUT = min(max(YA/SC(jj),-7.d0),7.d0)
		 JDSP  = 10*(DYM*YROUT+IYMN+128*(1-JUP))
		 IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
	      endif
	      YA = YX
	   enddo
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A6,1A4,1A1)')'Plot "',NAMER(jj),'"'
	      j = lonlen(STRI)
	      call	pscom(STRI,j)
	   endif
	   jc = JW8+2
	   STYL = (jc-1)*MARK
	   jcv = jcv+1
	   if (jcv .le. ICVMX)	then
	      call	PLOTCR
     >		(jxout,IWN(JW),IX,IXO(1,jcv),IY,IYO(1,jcv),jc,STYL,PT)
	   endif
	   IWN(JW) = jxout
	   do	J=1,NP1
	      ROUT(j,jj) = ROUT(j,jj)-OSHIFR(jj)
	   enddo
 10	continue

C Plot dots
	do	j=1,8
	    jlx(j) = 0
	enddo
	jsco = 0
	jn = 0
	do	14	jxout=1,NXOUT
	   if (NWINDX(jxout) .eq. 0)	goto	14
	   CHAR6 = NAMEX(jxout)
	   if (CHAR6(1:1) .eq. ' ')	goto	14
	   do	jj=1,NARRX
	      if (EXARNM(jj) .eq. CHAR6)	jn = jj	
	   enddo
	   if (jn .eq. 0)	then
	      jj = length(EQNAME)
	      write(*,*)
	      write(*,*)'>>> Error in the model "',EQNAME(1:jj),'"'
	      write(*,*)'>>> Unknown data set "',CHAR6,'" is requested'
	      call	a_stop
	   endif
	   jpnt = NPTM(jn)
	   if (jpnt .le. 0)	goto	14
	   jsc = NWINDX(jxout)
	   if (abs(SC(jsc)) .lt. 1.1E-7)	call
     >	        SCAL(1,SC(jsc),SCALER(jsc),DATAX(1,jn),jpnt,NRDX)
	   JW = NWIND1(jsc)-16*NSCR(MOD10)
	   if (JW.LE.0 .or. JW.GT.16)	goto	14
C	   write(*,*)EXARNM(jn),' is drawn in box',NWIND1(jsc),SC(jsc)
	   JW4 = (JW-1)/4
	   JWX4 = JW-4*JW4-1
	   JW8 = (JW-1)/8
	   JUP = JW4-2*JW8
	   JX = DXM*JWX4
C  JW4 =0,0,0,0,1,1,1,1,2...  JW8=0,0,0,0,0,0,0,0,1,1...
C  JWX4=0,1,2,3,0,1,2,3,0...  JUP=0,0,0,0,1,1,1,1,0,0...
	   call colovm(JW8+2)
	   if (jsc .ne. jsco)	then
	      jc = 5		! 13
	   else
	      jc = jc+1
	   endif

	   if (JFNEW .eq. 0)	then
	      call	colovm(EraseColor)
	      do	j=1,NPTMO(jxout)
		 call	NMARK(PTMO(1,j,jxout),jc)
	      enddo
	   endif

	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A6,1A6,1A1)')'Dots "',EXARNM(jn),'"'
	      j = lonlen(STRI)
	      call	pscom(STRI,j)
	   endif
	   YL = max(0.d0,ABSC(GRAL(jsc)))
	   YR = ABSC(GRAP(jsc)) 
	   if (YL.ge.YR)	YL = 0.d0
	   j1 = 0
	   do	  13	j=1,jpnt
	      YX = XAXES(j,jn)
	      if (MODEX .ge. 1 .and. YX .gt. ABC)	goto	13
	      YA = ABSC(YX)
	      if (YA .gt. 1.)	  goto	13
	      if (YA.lt.YL .or. YA.gt.YR)	goto	13
	      j1 = j1+1
	      PTM(1) = JX+160*(YA-YL)/(YR-YL)
	      YROUT=max((DATAX(j,jn)+OSHIFR(jsc))/SC(jsc),-7.d0)
	      YROUT=min(YROUT,7.d0)
	      JDSP	=DYM*YROUT+IYMN+128*(1-JUP)
	      PTM(2)	=350-min(max(JDSP,IYMN),IYMX)
	      call	colovm(JW8+2)
	      call	NMARK(PTM,jc)
	      PTMO(1,j1,jxout) = PTM(1)
	      PTMO(2,j1,jxout) = PTM(2)
 13	   continue
	   NPTMO(jxout) = j1

C (JPOS,jnl)=(ix,iy)
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A)')'Mark & time for dots'
	      j = lonlen(STRI)
	      call	pscom(STRI,j)
	   endif
	   JWX = JW-8*JW8
	   jlx(JWX) = jlx(JWX)+1
	   JPOS=JX+DXM-6*DXLET
	   jnl=(1+jlx(JWX))*DYLET+FSHIFT+(IYM-IY0-DYM)*JUP+3
	   call FMTXF4(XF4,TOUTX(jn))
	   call textvm(JPOS,jnl,XF4,5)
C	if (EXARNM(jn).eq.'CAR4X ')write(*,*)EXARNM(jn),TIME,
C     *		TOUTX(jn),jn,KOGDA(jn),NPTM(jn),DATAX(1,jn)
	   PTM(1) = JPOS+5.5*DXLET
	   PTM(2) = jnl-DYLET/2+1
	   call	NMARK(PTM,jc)

	   jsco = jsc
	   if (jc .eq. 7)	jc = 1
 14	continue

	do	15	jj=1,NROUT
	   JW = NWIND1(jj)-16*NSCR(MOD10)
	   if (NAMER(jj) .eq. '    ') JW = 0
	   if (JW.LE.0 .or. JW.GT.16)	goto 15
	   JW4 = (JW-1)/4
	   JWX4 = JW-4*JW4-1
	   JW8 = (JW-1)/8
	   JUP = JW4-2*JW8
	   JX = DXM*JWX4
	   STYL = (JW8+1)*MARK
	   call colovm(JW8+2)
	   JPOS=JX+JW8*DXM/2
	   jnl=DYLET+FSHIFT+(IYM-IY0+2+DYLET)*JUP+2
	   call CMARK(jnl,JPOS,SC(jj),OSHIFR(jj),NAMER(jj),STYL)
 15	continue
	if (JIFNEW .ge. 10)	return

C Erase/put resonance radii
	call colovm(EraseColor)
	do	17	J=1,640,160
		if(IQ1 .eq. 0)		GO TO 16
		JPOS = J+160*YQ1-1
		call drawvm(0,JPOS,IYM0,JPOS,IYM0-4)
 16		if(IQ2 .eq. 0)		GO TO 17
		JPOS = J+160*YQ2-1
		call drawvm(0,JPOS,IYM0,JPOS,IYM0-4)
		JPOS = JPOS+2
		call drawvm(0,JPOS,IYM0,JPOS,IYM0-4)
 17	continue
	call colovm(WarningColor)
	IQ1 = 0
	IQ2 = 0
	do	J=1,NAB-1
	    if (MU(J).gt.1. .and. MU(J+1).le.1.)	IQ1 = J
	    if (MU(J).gt..5 .and. MU(J+1).le..5)	IQ2 = J
	enddo
	if (KPRI.ge.1 .and. KPRI.le.2)	then
	    write(STRI,'(1A)')'Mark resonance radii'
	    j = lonlen(STRI)
	    call	pscom(STRI,j)
	endif
	if (IQ1 .ne. 0)	YQ1 = ABSC(AMETR(IQ1)) 
	if (IQ2 .ne. 0)	YQ2 = ABSC(AMETR(IQ2)) 
	do	19	J=1,640,160
		if(IQ1.eq.0)		GO TO 18
		JPOS=J+160*YQ1-1
		call drawvm(0,JPOS,IYM0,JPOS,IYM0-4)
 18		if(IQ2.eq.0)		GO TO 19
		JPOS=J+160*YQ2-1
		call drawvm(0,JPOS,IYM0,JPOS,IYM0-4)
		JPOS=JPOS+2
		call drawvm(0,JPOS,IYM0,JPOS,IYM0-4)
 19	continue
	return
C----------------------------------------------------------------------|
C Mode 2 & 3:
 2	if (MOD10.ne.2 .and. MOD10.ne.3)	goto	30
	call	add2loc("Drawing mode 2"//char(0))
	IY0GR=IYMN
	if (MODEY .eq. -1) IY0GR=IYMN+DYM
	do	jj=1,NROUT
	   do	J=1,NP1
		ROUT(j,jj) = ROUT(j,jj)+OSHIFR(jj)
	   enddo
	enddo
	call	SCAL(NROUT,SC,SCALER,ROUT(4,1),NP1-3,NRD)
	do	J=1,8
	   IWX(J) = 0
	enddo

	jcv = 0
	do	20	jj=1,NROUT
	   JW = NWIND2(jj)-16*(NSCR(MOD10)/2)
	   if (NAMER(jj) .eq. '    ') JW = 0
	   if (JW.LE.0 .or. JW.GT.16)	go to 20
	   JW4	= (JW-1)/4
	   JWX4 = JW-4*JW4-1
	   IMODE = NSCR(MOD10)-(NSCR(MOD10)/2)*2
	   if(JWX4/2.ne.IMODE)	go to 20
C  JW4 =0,0,0,0,1,1,1,1,2...  JWX4=0,1,2,3,0,1,2,3,0...
	   JNC = (JW-1)/2
	   JNW = JW-2*JNC-1
C  JNC =0,0,1,1,2,...  JNW=0,1,0,1,0,...
	   JX = DXM*JNW
	   JBOX = JW-(JW-1)/2*2		! 1,2,1,2,1,2,1,2,...
	   YL = ABSC(GRAL(jj)) 
	   YR = ABSC(GRAP(jj)) 
	   if (YL.ge.YR)	YL = 0.d0
	   if (YL.lt.YR .and. 
     +		(YL.gt.AMETR(1)/ABC .or. YR.lt.0.999) )	then
	      PT(1) = JX+320*YL
	      PT(2) = 350-IYMN-2*IWX(JBOX)-2
	      PT(3) = JX+320*YR
	      call	colovm(JW4+2)
	      call	drawvm(0,PT(1),PT(2),PT(3),PT(2))
	      call	drawvm(0,PT(1),PT(2)+1,PT(3),PT(2)+1)
	      IWX(JBOX) = IWX(JBOX)+1
	   endif

	   jxout = 0
	   do	j=1,NP1
	      YX = ABSC(AMETR(j)) 
	      if (YX.lt.YL .or. YX.gt.YR)	goto	21
	      jxout = jxout+1
	      if (jxout.eq.1 .and. j.gt.1)	then ! left edge interpolation
		 IX(jxout) = JX
		 YA=ROUT(J,jj)+(ROUT(J-1,jj)-ROUT(J,jj))*(YL-YX)/(YA-YX)
		 YROUT = min(max(YA/SC(jj),-7.d0),7.d0)
		 JDSP  = 10*(DYM*YROUT+IY0GR)
		 IY(jxout) = min(max(JDSP,10*IYMN),10*IYMX)
		 jxout = jxout+1
	      endif
	      IX(jxout) = JX+320*(YX-YL)/(YR-YL)
	      YROUT = min(max(ROUT(J,jj)/SC(jj),-7.d0),7.d0)
	      JDSP  = 10*(DYM*YROUT+IY0GR)
	      IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
 21	      continue
	      if (YA.le.YR .and. YX.gt.YR)	then ! right edge interpolation
		 jxout = jxout+1
		 IX(jxout) = JX+320
		 YA=ROUT(J,jj)+(ROUT(J-1,jj)-ROUT(J,jj))*(YR-YX)/(YA-YX)
		 YROUT = min(max(YA/SC(jj),-7.d0),7.d0)
		 JDSP  = 10*(DYM*YROUT+IY0GR)
		 IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
	      endif
	      YA = YX
	   enddo
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A6,1A4,1A1)')'Plot "',NAMER(jj),'"'
	      j = lonlen(STRI)
	      call	pscom(STRI,j)
	   endif
	   jc = JW4+2
	   STYL = (jc-1)*MARK
	   jcv = jcv+1
	   if (jcv .le. ICVMX)	then
	      call	PLOTCR
     >		(jxout,IWN(JW),IX,IXO(1,jcv),IY,IYO(1,jcv),jc,STYL,PT)
	   endif
	   IWN(JW) = jxout
	   do	J=1,NP1
	      ROUT(j,jj) = ROUT(j,jj)-OSHIFR(jj)
	   enddo
 20	continue

	do	j=1,8
	    jlx(j) = 0
	enddo
	jsco = 0
	do	23	jxout=1,NXOUT
	   if (NWINDX(jxout) .eq. 0)	goto	23
	   CHAR6 = NAMEX(jxout)
	   if (CHAR6(1:1) .eq. ' ')	goto	23
	   do	jj=1,NARRX
	      if (EXARNM(jj) .eq. CHAR6)	jn = jj	
	   enddo
	   jpnt = NPTM(jn)
	   if (jpnt .le. 0)	goto	23
	   jsc = NWINDX(jxout)
	   if (abs(SC(jsc)) .lt. 1.1E-7)	call
     >	        SCAL(1,SC(jsc),SCALER(jsc),DATAX(1,jn),jpnt,NRDX)
	   JW = NWIND2(jsc)-16*(NSCR(MOD10)/2)
	   if (JW.LE.0 .or. JW.GT.16)	goto	23
	   JW4 = (JW-1)/4
	   JWX4 = JW-4*JW4-1
	   IMODE = NSCR(MOD10)-(NSCR(MOD10)/2)*2
	   if (JWX4/2 .ne. IMODE)	goto	23
	   JNC	= (JW-1)/2
	   JNW = JW-2*JNC-1
	   JX = DXM*JNW
	   STYL = (JW4+1)*MARK
	   call colovm(JW4+2)
	   if (jsc .ne. jsco)	then
	      jc = 5
	   else
	      jc = jc+1
	   endif
	   if (JFNEW .eq. 0)	then
	      call	colovm(EraseColor)
	      do	j=1,NPTMO(jxout)
		 call	NMARK(PTMO(1,j,jxout),jc)
	      enddo
	   endif

	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A6,1A6,1A1)')'Dots "',EXARNM(jn),'"'
	      j = lonlen(STRI)
	      call	pscom(STRI,j)
	   endif
	   YL = ABSC(GRAL(jsc)) 
	   YR = ABSC(GRAP(jsc)) 
	   if (YL.ge.YR)	YL = 0.d0
	   j1 = 0
	   do	  22	j=1,jpnt
	      YX = XAXES(j,jn)
	      if (MODEX.ne.0 .and. YX.gt.ABC)	goto	22
	      YA = ABSC(YX)
	      if (YA .gt. 1.)			goto	22
	      if (YA.lt.YL .or. YA.gt.YR)	goto	22
	      j1 = j1+1
	      PTM(1) = JX+320*(YA-YL)/(YR-YL)
	      YROUT=max((DATAX(j,jn)+OSHIFR(jsc))/SC(jsc),-7.d0)
	      YROUT=min(YROUT,7.d0)
	      JDSP	=DYM*YROUT+IY0GR
	      PTM(2)	=350-min(max(JDSP,IYMN),IYMX)
	      call	colovm(JW4+2)
	      call	NMARK(PTM,jc)
	      PTMO(1,j1,jxout) = PTM(1)
	      PTMO(2,j1,jxout) = PTM(2)
 22	   continue
	   NPTMO(jxout) = j1

C (JPOS,jnl)=(ix,iy)
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A)')'Mark & time for dots'
	      j = lonlen(STRI)
	      call	pscom(STRI,j)
	   endif
	   JWX = JW-2*JNC
	   jlx(JWX) = jlx(JWX)+1
	   JPOS=JX+DXM-6*DXLET
	   jnl=(1+jlx(JWX))*DYLET+FSHIFT+3
	   call FMTXF4(XF4,TOUTX(jn))
	   call textvm(JPOS,jnl,XF4,5)
C	if (EXARNM(jn).eq.'TEX   ')write(*,*)EXARNM(jn),TIME,
C     *		TOUTX(jn),jn,KOGDA(jn),NPTM(jn),DATAX(1,jn)
	   PTM(1) = JPOS+5.5*DXLET
	   PTM(2) = jnl-DYLET/2+1
	   call	NMARK(PTM,jc)

	   jsco = jsc
	   if (jc .eq. 7)	jc = 1
 23	continue

	do	201	jj=1,NROUT
	    JW = NWIND2(jj)-16*(NSCR(MOD10)/2)
	    if (NAMER(jj) .eq. '    ') JW = 0
	    if (JW.LE.0 .or. JW.GT.16)	goto	201
	    JW4	= (JW-1)/4
	    JWX4 = JW-4*JW4-1
	    IMODE = NSCR(MOD10)-(NSCR(MOD10)/2)*2
	    if(JWX4/2.ne.IMODE)		goto	201
	    JNC	= (JW-1)/2
	    JNW = JW-2*JNC-1
	    JX = DXM*JNW
	    STYL = (JW4+1)*MARK
	    call colovm(JW4+2)
	    JPOS=JX+JW4*DXM/4
	    jnl=DYLET+FSHIFT
C	write(*,*)jj,jnl,JPOS,SC(jj),OSHIFR(jj),NAMER(jj),STYL
	    call CMARK(jnl,JPOS,SC(jj),OSHIFR(jj),NAMER(jj),STYL)
C	write(*,*)jj,NROUT
 201	continue
	NP1O = NP1

C	write(*,*)MOD10,MODEY
	return
C----------------------------------------------------------------------|
C Mode 4 & 5:
 30	if (MOD10.ne.4 .and. MOD10.ne.5)	goto	4
C----------------------------------------------------------------------|
C Mode 6:
 4	if(MOD10.ne.6)	goto	5
	call	add2loc("Drawing mode 6"//char(0))
	LTOUT1 = 1
	LTOUT2 = LTOUT
	do	J=1,16
	   JWC(J) = 0		! Ordinal curve number in the current mindow
	enddo
C right_label_position=JDX*JDMX=23*5*5=575 (see typdsp.f)
	do	J=1,LTOUT-1
	   YROUT    = (TTOUT(J)-TINIT)*575/abs(TSCALE)
	   IYO(J,ICVMX+1) = 6*DXLET+YROUT
	   if(YROUT.lt.0)		LTOUT1 = J+1
	   if(IYO(J,ICVMX+1) .le. XSCMAX-1)	LTOUT2 = J-1
	enddo
	LTOUT2 = LTOUT2-LTOUT1+1
	if (LTOUT2.lt.3)	then ! No plot for small time
	   return
	endif
	do	jj=1,NTOUT
	   do	j=LTOUT1,LTOUT1+LTOUT2
	      TOUT(j,jj) = TOUT(j,jj)+OSHIFT(jj)
	   enddo
	enddo
	call	SCAL(NTOUT,SC,SCALET,TOUT(LTOUT1,1),LTOUT2,ITIMES)

	jcv = 0
	do	43	JJ=1,NTOUT
	   JW=NWIND3(JJ)-8*NSCR(MOD10)
	   if (NAMET(JJ) .eq. '    ') JW = 0
	   if (JW.le.0 .or. JW.gt.8)	go to 43
	   if (MODEY .eq. 1)	JW4 = 2		! two-window mode
	   if (MODEY .eq.-1)	JW4 = 4		! four-window mode
	   if (MODEY .eq. 0)	JW4 = 1		! one-window mode
	   JNC = (JW-1)/JW4			! current_window_No.-1
	   JNW = JW-JW4*JNC-1			! current_curve_No.-1
	   JWC(JNW+1) = JWC(JNW+1)+1
	   do	J=1,LTOUT
	      YROUT = max(TOUT(J,JJ)/SC(JJ),-7.d0)
	      YROUT = min(YROUT,7.d0)
	      JDSP  = 10*(DYM*YROUT+IYMN+(JW4-1-JNW)*DYM)
	      JDSP  = max(JDSP,10*IYMN)
	      IYO(J,ICVMX+2) = JY-min(JDSP,10*IYMX)
	   enddo
	   jnl = 0
	   if (JFNEW .eq. 0) jnl = LTOUT2
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A6,1A4,1A1)')'Plot "',NAMET(jj),'"'
	      j = lonlen(STRI)
	      call	pscom(STRI,j)
	   endif
	   jc = JWC(JNW+1)		! Ordinal curve number in its window
	   STYL = jc*MARK
	   jc = jc+1
	   jcv = jcv+1
	   if (jcv .le. ICVMX)	then
	      call PLOTCR(LTOUT2+1,jnl,IYO(LTOUT1,ICVMX+1),IXO(1,1),
     >	               IYO(LTOUT1,ICVMX+2),IYO(LTOUT1,jcv),jc,STYL,PT)
	   endif
	   JPOS=0
	   if (MODEY .eq. 0)	then
	      if (JWC(JNW+1) .gt. 13)	goto	42		! max 13 curves
	      jnl = FSHIFT+DYLET*(2*JWC(JNW+1)-1)
	   endif
	   if (MODEY .eq. 1)	then
	      if (JWC(JNW+1) .gt. 5)	goto	42		! max 5 curves
	      jnl = FSHIFT+DYLET*(2*JWC(JNW+1)-1)+JNW*DYM
	   endif
	   if (MODEY .eq.-1)	then
	      if (JWC(JNW+1) .gt. 2)	goto	42		! max 2 curves
	      jnl = FSHIFT+DYLET*(2*JWC(JNW+1)-1)+JNW*DYM
	   endif

C	   write(*,*)JW,JNW+1,JNC+1,JWC(JNW+1),STYL,"  ",NAMET(JJ),jnl
	   call CMARKT(jnl,JPOS,SC(JJ),OSHIFT(jj),NAMET(JJ),STYL)
 42	   continue
	   do	j=LTOUT1,LTOUT1+LTOUT2
		TOUT(j,jj) = TOUT(j,jj)-OSHIFT(jj)
	   enddo
 43	continue
C----------------------------------------------------------------------|
C Mode 7:
 5	if(MOD10.ne.7)	go to 6
	call	add2loc("Drawing mode 7"//char(0))
	if(LTOUT.lt.2)		return
	I=0
	do 52 JJ=1,4
	   NUM(JJ)=0
	   do 51 J=1,NTOUT
	      if (NWIND7(J)-4*NSCR(MOD10) .eq. JJ)	then
		 NUM(JJ) = J
		 goto 52
	      endif
 51	   continue
 52	continue
	NKL1=2
	NKL2=1
	if (NUM(1).ne.0.and.NUM(2).ne.0)NKL1=1
	if (NUM(3).ne.0.and.NUM(4).ne.0)NKL2=2
	if (NKL1.eq.2.and.NKL2.eq.1) goto 6
C  Time interval : Tmin7 < t < Tmax7, Mark interval - Tmet
	JMIN	=1
	JMAX	=0
	do 53 J=1,LTOUT-1
	if (TTOUT(J) .le. TIM7(1))	JMIN=J
	if (TTOUT(J) .le. TIM7(2))	JMAX=J
 53	continue
	if (JMIN .gt. JMAX)	goto 6
	call SCAL(NTOUT,SC,SCALET,TOUT,ITIMES,ITIMES)
	do	57	JJ=nkl1*2,nkl2*2,2
	   JPOSY=MODK(JJ/2)*(IYM-IY0)/2	
	   JX=JJ/4*XSCMAX/2+MODK(JJ/2)*XSCMAX/4
	   do	56	J=JMIN,JMAX
	      YROUT=max(TOUT(J,NUM(JJ-1))/SC(NUM(JJ-1)),-1.d0)
	      YROUT=min(YROUT,1.d0)
	      JDSP	=(2-MODK(JJ/2))*(IYM-IY0)/2*YROUT+IYMN+JPOSY
	      JDSP	=max(JDSP,-MODK(JJ/2)*(IYM-IY0)/2+IYMN)
	      j1	=J-JMIN+1
	      IYO(j1,ICVMX+2)=min(JDSP,IYMX)
	      XROUT=max(TOUT(J,NUM(JJ))/SC(NUM(JJ)),-1.d0)
	      XROUT=min(XROUT,1.d0)
	      JDSP	=(2-MODK(JJ/2))*XSCMAX/4*XROUT
	      JDSP	=max(JDSP,XSCMAX/4*(-MODK(JJ/2)))
	      IYO(j1,ICVMX+1)=min(JDSP,XSCMAX/2)
 56	   continue
	   STYL	=TIM7(4)+6.5
	   call PLOTXY(TTOUT(JMIN),JMAX-JMIN+1,JX,IYO(1,ICVMX+1),
     .		IYO(1,JJ),IYO(1,ICVMX+2),IYO(1,JJ-1),TIM7(3),STYL,PT)
C	   write(*,*)"After PLOTXY"
	   call FMTF4(ST(1:4),SC(NUM(JJ-1)))
	   ST(5:9)=NAMET(NUM(JJ-1))
	   call colovm(1)
	   call textvm(JJ/4*XSCMAX/2,DYLET+FSHIFT,ST,9)
	   call FMTF4(ST(1:4),SC(NUM(JJ)))
	   ST(5:9)=NAMET(NUM(JJ))
	   call textvm(JJ/2*XSCMAX/2-DXLET*10,IYM-2*DYLET+FSHIFT,ST,9)
 57	continue
C----------------------------------------------------------------------|
C Mode 8:
 6	if(MOD10.ne.8)	go to 7
C	write(*,*)"Outdsp 8",jifnew,KPRI,LEQ(5)
	call	add2loc("Drawing mode 8"//char(0))
	YS0 = (IY0+IYM)/2			! Mid-plane position
	YSC8 = IDX*IDT/SCM
	j = (RTOR+SHIF(1))*YSC8
	call	colovm(1)
	call	drawvm(0,jzero,YS0,j,YS0)	! Plot the mid-plane

C Plot the complete wall structure (Pixmap # 1)
	if (jifnew .ne. 0)	then
C	   write(*,*)"Outdsp 8",jifnew,KPRI,CNFILE
	   call	drconf(YS0,YSC8)
	   call	redraw(1)
	endif
	if (NBFLAG .eq. 0 .or. NBFILE(1:1) .eq. '*')	goto	67
	if (jifnew .ne. 0)	then	! Create/update NBI Pixmap # 2
	   call	drawfoot(jifnew)
	   call	redraw(2)
	   NBFLAG = 0
	endif
 67	continue
C	call	makebnd(YS0,YSC8)

C or plot the AWALL boundary in blue instead
	if (CNFILE(1:1) .eq. '*')	then
	   call	colovm(3)
	   call	drawvm(0,350,35,400,35)
	   call	textvm(410,35+DYLET/2,'Wall',4)
	   do	j=1,130
	      YFI=GP*(j-1)/64.
	      YZ=AB*ELONM*SIN(YFI)
	      YR=RTOR+AB*(COS(YFI)+0.5*TRICH*(COS(2.*YFI)-1.))
	      PT(2*j-1)=10.*YR*YSC8
	      PT(2*j)=10.*(YS0-YZ*YSC8)
	   enddo
	   j=0
	   jn=130
	   call d2polyline(j,PT,jn)
	endif

C if (data file includes NAMEXP BND) then (NBND > 0);
C or (NBND == 8) after calling ESC/SPIDER with no boundary points provided;
C     NBND == 0 otherwise

	if (LEQ(5).le.1 .or. TASK(1:4).eq.'VIEW')	then
	   call	DRAW3M (jifnew,DYLET,YS0,YSC8,PT,EQOLD,
     >			NA1,NAB,RTOR,AMETR,SHIF,UPDWN,ELON,TRIA)
	elseif (LEQ(5).eq.2)	then
	   call	BNDRAW(JIFNEW,YS0,YSC8,IYO,NBND,NBNT,TIME,BNDARR)
	   call colovm(EraseColor)
	   do  JJ=1,10
	      do	J=1,130
		 PT(2*J-1) = EQOLD(2*J-1,JJ)
		 PT(2*J)   = EQOLD(2*J,JJ)
	      enddo
	      j = 0
	      jn=130
	      call d2polyline(j,PT,jn)
	   enddo
	   call	drawconfig(0,YS0,YSC8,EQOLD)
	elseif (LEQ(5).eq.3)	then
	   call	BNDRAW(JIFNEW,YS0,YSC8,IYO,NBND,NBNT,TIME,BNDARR)
	   call	DRAWSPFLUX (jifnew,DYLET,YS0,YSC8)
C	   call	DRAW3M (jifnew,DYLET,YS0,YSC8,PT,EQOLD,
C     >			NA1,NAB,RTOR,AMETR,SHIF,UPDWN,ELON,TRIA)
	endif

 60	continue
	jc = -1
	do	69	jxout=1,NXOUT
	    if (NWINDX(jxout) .eq. 0)		goto	69
	    CHAR6 = NAMEX(jxout)
	    if (CHAR6(1:1) .eq. ' ')		goto	69
	    if (jxout .gt. 1)	then
		do	jj=jxout-1,1,-1
		    if (NAMEX(jj) .eq. CHAR6)	goto	69
		enddo
	    endif
	    do	jj=1,NARRX
		if (EXARNM(jj) .eq. CHAR6)	jn = jj
	    enddo
C	    jpnt = NPTM(jn)
	    jarr = IFDFAX(jn)
	    if (jarr .le. 0)			goto	69
	    jtyp = NTYPEX(jarr)
	    if (jtyp .lt. 18)			goto	69
	    jpnt = NGRIDX(jarr)
	    if (jpnt .le. 0)			goto	69
	    jc = jc+1
	    js = GDEX(jarr)

C		if (jtyp .eq. 18 .or. jtyp .eq. 19)	then
C		    write(*,*)EXARNM(jn),jn,jpnt,js,DATARR(js)
C		    write(*,*)(DATARR(jj),jj=js+1,js+jpnt)
C		elseif (jtyp .eq. 20)	then
C		    write(*,*)'Points:',jpnt,'  Address:',js
C		    write(*,*)(DATARR(jj),jj=js,js+jpnt-1)
C		    write(*,*)(DATARR(jj),jj=js+jpnt,js+2*jpnt-1)
C		    write(*,*)(DATARR(jj),jj=js+2*jpnt,js+3*jpnt-1)
C		    write(*,'(A,10F6.2)')"CAR6X  ",(CAR6X(j),j=10,1,-1)
C		endif
	    if (JFNEW .eq. 0)	then
		call	colovm(EraseColor)
		do	j=1,NPTMO(jxout)
		   call	NMARK(PTMO(1,j,jxout),7)
		enddo
	    endif
	    call colovm(jc+2)
	    do	j=1,jpnt
		if (jtyp .eq. 18)		then
		    YR = DATARR(js)
		    YZ = DATARR(js+j)
		elseif (jtyp .eq. 19)		then
		    YR = DATARR(js+j)
		    YZ = DATARR(js)
		elseif (jtyp .eq. 20)	then
		    YR = DATARR(js-1+j)
		    YZ = DATARR(jpnt+js-1+j)
		else
		    write(*,*)'Unknown input-grid type'
		endif
		PTM(1)=YR*YSC8
		PTM(2)=YS0-YZ*YSC8
		call	NMARK(PTM,7)
		PTMO(1,j,jxout) = PTM(1)
		PTMO(2,j,jxout) = PTM(2)
	    enddo
	    call textvm(510,100+2*DYLET*jc,NAMEX(jxout),6)
	    NPTMO(jxout) = jpnt
 69	continue
C----------------------------------------------------------------------|
C Mode 9:				User's drawing
 7	if (MOD10 .ne. 9)	go to 8
	call	add2loc("User drawing mode"//char(0))
C	call	MYDRAW			! Example of user's drawing
C----------------------------------------------------------------------|
C Mode 0:				No graphic output
 8	if (MOD10 .ne. 0)		return
	end
C======================================================================|
C	call		BNDRAW(YS0,YSC8,IYO,NBND,NBNT,TIME,BNDARR)
	subroutine	BNDRAW(ifnew,YS0,SC8,IYO,NBND,NBNT,TIME,BNDARR)
C----------------------------------------------------------------------|
C IFNEW  =  0	Re-draw (erase) the previous curves
C IFNEW =/= 0	New curves only
C IFNEW  <  0	Don't mark resonances q=m/n
C IFNEW  > 10	Call from Review. (JIFNEW-10) is used to control erasing
	implicit	none
	integer	 	ifnew,j,j1,j2,jj,NBND,NBNT,YS0,PTM(2),IYO(2,*)
	double precision SC8,YS,YX,YXL,YXR,YZ,TIME,BNDARR(*)
C----------------------------------------------------------------------|
	j2 = 1+NBND/32
	if (IFNEW .ne. 0)	goto	60	! 
	j = 0			! Erase the old points
C	j = 31			! Draw the old points with shadow color
	call colovm(j)
	do	j=1,NBND,j2
	   PTM(1) = IYO(1,j)
	   PTM(2) = IYO(2,j)
	   call	NMARK(PTM,4)
C	   call	CNMARK(PTM,3,9)
	enddo

 60	continue
C Boundary points: (red)
	call	colovm(2)
	if (NBNT .gt. 1)	goto	62
 61	continue
	jj = max(1,NBNT)
C	write(*,*)"Drawing boundary points"
	do	j=1,NBND,j2
	   j1 = NBNT+(2*j-1)*jj
C	   write(*,'(1I5,2F8.3)')j,BNDARR(j1),BNDARR(j1+jj)
	   PTM(1) = BNDARR(j1)*SC8
	   PTM(2) = YS0-BNDARR(j1+jj)*SC8
C	   call	  NMARK(PTM,5)			!Use (PTM,5) for crosses
	   call	  NMARK(PTM,4)			!Use (PTM,4) for *
C	   call	  CNMARK(PTM,3,9)		! Put bitmap rectangle 4x4
	   IYO(1,j) = PTM(1)
	   IYO(2,j) = PTM(2)
	enddo
	return

 62	continue
	if (TIME .le. BNDARR(1))	goto	61
	if (TIME .ge. BNDARR(NBNT))	goto	61

	do	j=1,NBNT			! Find current time
	   if (TIME .gt. BNDARR(j)) jj = j
	enddo
	YS  = BNDARR(jj+1)-BNDARR(jj)
	YXL = (TIME-BNDARR(jj))/YS
	YXR = (TIME-BNDARR(jj+1))/YS
	do	j=1,NBND,j2
	   j1 = jj+(2*j-1)*NBNT
	   YX = YXL*BNDARR(j1+1)-YXR*BNDARR(j1)
	   j1 = jj+2*j*NBNT
	   YZ = YXL*BNDARR(j1+1)-YXR*BNDARR(j1)
	   PTM(1) = YX*SC8
	   PTM(2) = YS0-YZ*SC8
	   call	  NMARK(PTM,4)
C	   call	CNMARK(PTM,3,9)
	   IYO(1,j) = PTM(1)
	   IYO(2,j) = PTM(2)
	enddo
	end
C======================================================================|
	subroutine	DRAW3M(jifnew,DYLET,YS0,YSC8,PT,EQOLD,
     >				NA1,NAB,RTOR,AMETR,SHIF,UPDWN,ELON,TRIA)
C----------------------------------------------------------------------|
	implicit	none
	integer		j,jj,jn,JX,jifnew,BGColor,NA1,NAB,YS0,DYLET
	integer		PT(*),EQOLD(260,11)
	double precision YSC8,GP,RTOR,AMETR(*),ELON(*),TRIA(*)
	double precision YR,YZ,YFI,UPDWN,SHIF(*)
	save	BGColor,GP
	data	BGColor/31/ GP/3.1415926536/
C----------------------------------------------------------------------|
C Redraw magnetic surfaces:
	do  65	J=1,10
C Erase:
	   if (jifnew .eq. 0)	then
	      call colovm(BGColor)
	      do	JJ=1,130
		 PT(2*JJ-1)=EQOLD(2*JJ-1,J)
		 PT(2*JJ)=EQOLD(2*JJ,J)
	      enddo
	      jj=0
	      jn=130
	      call d2polyline(jj,PT,jn)
	   endif

C New configuration:
	   JX = max(1.,0.1*NA1*J-1)
	   JX = min(NA1,JX)
	   if (J.le.10)	call colovm(14)		! Outermost flux surface
	   if (J.eq.10)	then
	      call colovm(2)
	      JX = NA1
	      if (NA1.ne.NAB)	then
		 call drawvm(0,350,55,400,55)
		 call textvm(410,55+DYLET/2,'Transport boundary',18)
	      endif
	   endif
	   if (J .eq. 11)	then ! Never happens, reserved for boundary
	      call colovm(13)			! Pink
	      JX = NAB
	      call drawvm(0,350,75,400,75)
	      call textvm(410,75+DYLET/2,'Plasma boundary',15)
	   endif
	   do	JJ=1,130
	      YFI=GP*(JJ-1)/64.
	      YZ=UPDWN+AMETR(JX)*ELON(JX)*sin(YFI)
	      YR=RTOR+SHIF(JX)+AMETR(JX)*(cos(YFI)
     +		+0.5*TRIA(JX)*(cos(2.*YFI)-1.))
	      PT(2*JJ-1)=10.*YR*YSC8
	      PT(2*JJ)=10.*(YS0-YZ*YSC8)
C Save picture:
	      EQOLD(2*JJ-1,J)=PT(2*JJ-1)
	      EQOLD(2*JJ,J)=PT(2*JJ)
	   enddo
	   jj=0
	   jn=130
	   call d2polyline(jj,PT,jn)
 65	continue
	PT(1)=(RTOR+SHIF(1))*YSC8
	PT(2)=(YS0-UPDWN*YSC8)
C	call	cnmark(PT,1,2)
C	write(*,*)"CNMARK called"
	end
C======================================================================|
	subroutine drconf(YS0,YSC8)
C----------------------------------------------------------------------|
	implicit	none
	include 'for/parameter.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	character 	STRI*64
	integer		j,j1,ix1,ix2,iy1,iy2,YS0,length
	integer		ierr, nSHOT, Ndim, NdGC, NGC
	integer		ixbeg(40), lenix(40), valix(40)
        double precision	YSC8,xyGC(750,2),YR,YZ
C	equivalence	(WORK1(1,1),xyGC(1,1))
	save		ierr
	data		ierr/1/
C----------------------------------------------------------------------|
        if (CNFILE(1:1) .eq. '*')	return
        if (ierr .eq. 1)	goto	5
        if (KPRI.ge.1 .and. KPRI.le.2)	goto	5
	return
 5	continue
	j1 = length(CNFILE)
	call	OPENRD(7,CNFILE(1:j1),0,ierr)
	read(7,*)STRI
	read(7,*)nSHOT
	read(7,*)Ndim
	if (Ndim .gt. 750)	goto	11
	read(7,*)NGC
	if (NGC .gt. 40)	goto	11
	read(7,*)((xyGC(j,j1),j1=1,2),j=1,Ndim)
	read(7,*)STRI
	read(7,*)(ixbeg(j),j=1,NGC)
	read(7,*)STRI
	read(7,*)(lenix(j),j=1,NGC)
	read(7,*)STRI
	read(7,*)(valix(j),j=1,NGC)
	close(7)
	do  10	j1=1,NGC
	   if (valix(j1) .eq. 0)	goto	10
	   call	colovm(valix(j1))
	   do	j=ixbeg(j1),ixbeg(j1)+lenix(j1)-1
	      ix1 = ix2
	      iy1 = iy2
	      ix2 = 10.*YSC8*xyGC(j,1)
	      iy2 = 10.*(YS0-YSC8*xyGC(j,2))
	      if (j .gt. ixbeg(j1))	call  d2line(1,ix1,iy1,ix2,iy2)
C	      if (j .eq. ixbeg(j1))	then
C		 write(*,*)j1,ixbeg(j1),ixbeg(j1)+lenix(j1)-1
C		 write(*,'(4F9.4)')xyGC(j,1),xyGC(j,2),
C     >	    xyGC(ixbeg(j1)+lenix(j1)-1,1),xyGC(ixbeg(j1)+lenix(j1)-1,2)
C	      endif
	   enddo
 10	continue
	return
 11	continue
	close(7)
	write(*,*)"Configuration file is too long"
	end
C======================================================================|
	subroutine makebnd(YS0,YSC8)
C----------------------------------------------------------------------|
C This subroutine reads AUG file YGC and translates it to the Astra
C input format, i.e. file AWD/exp/cnf/*
C----------------------------------------------------------------------|
	implicit	none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/outcmn.inc'
	integer		j,j1,nch,PTM(2),YS0
	integer	ierr, nSHOT, nEDIT, Ndxy, Nxy, itopbt, indxpt
	real*4	xysp(750,2),tSHOT,psip,tSHf,Rgeo,zgeo,Rmag,zmag
	double precision	YSC8
	data	Ndxy/750/ nEDIT/0/ itopbt/1/
C----------------------------------------------------------------------|
	tSHOT = TIME
	nSHOT = 18091
C	write(*,*)iERR,"AUGD","FPG",nSHOT,nEDIT,tSHOT,itopbt,Ndxy
C        call	kkfpsp   (iERR  ,"AUGD","FPG",nSHOT,nEDIT  ,tSHOT,
C     <                 itopbt,Ndxy,xysp,Nxy,indxpt,tSHf ,
C     <                 Rgeo,zgeo,Rmag,zmag)

C Boundary points: (red)
	call	colovm(2)
	do	j=1,Nxy
	   PTM(1) = xysp(j,1)*YSC8
	   PTM(2) = YS0-xysp(j,2)*YSC8
	   call	  NMARK(PTM,5)			! Crosses
	enddo
	if (indxpt .gt. 0)	then
	   PTM(1) = xysp(indxpt,1)*YSC8
	   PTM(2) = YS0-xysp(indxpt,2)*YSC8
	   call	  NMARK(PTM,2)
	   call	  puto(PTM(1),PTM(2),3,5)		! cross
	endif
	PTM(1) = Rgeo*YSC8
	PTM(2) = YS0-zgeo*YSC8
	call	  puto(PTM(1),PTM(2),2,2)		! circle
	PTM(1) = Rmag*YSC8
	PTM(2) = YS0-zmag*YSC8
	call	  puto(PTM(1),PTM(2),3,7)		! fcircle
C	call	  NMARK(PTM,2)
	return

C        call	kkEQpsp  (iERR  ,expnam,dianam,nSHOT,nEDIT ,tSHOT,
C    >                 psip,
C    <                 Ndxy,xysp,Nxy,work,        tSHf )

C Output:
C   xysp(1:Ndxy,2) - array of (r,z) coordinates
C   Nxy   - number of points stored in xysp
C   tSHf  - time point [s] found in shot
C     &			xysp, Nxy, work, tSHf)
 100	format(1P,2(2E12.4,3X))
 101	format(1P,2(E12.4,3X),A)
 102	format(I7,A)
 103	format(1P,E12.4,A)
	write(*,*)Nxy
	NCH = 6
C	call	OPENWT(NCH,"exp/aug.bnd",0,j1)
	write(NCH,*)"AUG configuration data file"
	write(NCH,102)nSHOT,"              Shot number"
	write(NCH,103)tSHOT,     "         Time requested [s]"
	write(NCH,103)tSHf,      "         Time found [s]"
	write(NCH,102)Nxy,"              number of points stored in xysp" 
	write(NCH,101)Rgeo,zgeo,"              (r,z) geometric axis [m]"
	write(NCH,101)Rmag,zmag,"              (r,z) magnetic axis [m]"
	write(NCH,102)indxpt,"              index of x_point within xysp"
	write(NCH,100)xysp(indxpt,1),xysp(indxpt,2)
	write(NCH,100)(xysp(j,1),j=1,Nxy)
	write(NCH,100)(xysp(j,2),j=1,Nxy)
	close(NCH)
	call	a_stop
	end
C======================================================================|
	subroutine makeAUGcnf
C----------------------------------------------------------------------|
C This subroutine reads AUG file YGC and translates it to the Astra
C input format, i.e. file AWD/exp/cnf/*
C----------------------------------------------------------------------|
	implicit	none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/outcmn.inc'
	integer		j,j1
	integer*4	ierr, nSHOT, nEDIT, Ndim, NdGC, NGC
	integer*4	ixbeg(40), lenix(40), valix(40)
	character*8	GCnam(40)
	real*4	xyGC(750,2)
	data	Ndim/750/ NdGC/40/ nEDIT/0/ nSHOT/16324/
C----------------------------------------------------------------------|
	nSHOT = 14513
C	call	kkGCd0 ( ierr, "AUGD", "YGC", nSHOT, nEDIT, Ndim, xyGC,
C     &			 NdGC, NGC, ixbeg, lenix, valix, GCnam )
 100	format(1P,2(2E12.4,3X))
 101	format(10(I7))
 102	format(I7,A)
	call	OPENWT(7,"exp/cnf/aug",0,j1)		! IERR=j1
	write(7,*)"AUG configuration data file"
	write(7,102)nSHOT,"              Shot number"
	write(7,102)Ndim,"              1st dimension of xyGC file" 
	write(7,102)NGC,"              Dimension of GC-element arrays" 
	write(7,100)((xyGC(j,j1),j1=1,2),j=1,Ndim)
	write(7,*)"ixbeg"
	write(7,101)(ixbeg(j),j=1,NGC)
	write(7,*)"lenix"
	write(7,101)(lenix(j),j=1,NGC)
	write(7,*)"valix"
	write(7,101)(valix(j),j=1,NGC)
	close(7)
	return
	end
C======================================================================|
       	double precision function	ABSC(YIN)
C----------------------------------------------------------------------|
C Input: MODEX,YIN,FP
C Output: Value a=YIN mapped to the current abscissa
C----------------------------------------------------------------------|
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'			! XOUT,
	include 'for/status.inc'
	include 'for/outcmn.inc'			! MOD10 is used
	double precision	YAB,YIN,QUADIN,RFA,Y
	integer	MODEX,j
C----------------------------------------------------------------------|
	MODEX = XOUT+.49
	if     (MODEX .eq. 0)	then
C Draw up to AB against "a"
	   YAB = AB
	elseif (MODEX .eq. 1)	then
C Draw up to ABC against "a"
	   YAB = ABC
	elseif (MODEX .eq. 2)	then
C Draw up to ROC against "rho"
	   YAB = ROC
	elseif (MODEX .eq. 3)	then
C Draw up to ROC against "Psi"
	   YAB = 1.
	else
C Unknown option
	   return
	endif

	if (MOD10.eq.3 .or. MODEX.eq.3)	then
	   ABSC = QUADIN(NA1,AMETR,FP,YIN,Y,j)
	   ABSC = (ABSC-FP(1))/(FP(NA1)-FP(1))
	   return
	endif

	if     (MODEX .eq. 0 .or. MODEX .eq. 1)	then
	   ABSC = min(1.d0,max(0.d0,YIN/YAB))
	elseif (MODEX .eq. 2)	then
	   ABSC = RFA(YIN)/YAB
	else
	   write(*,*)"It cannot be"
	endif
	end
C======================================================================|
       subroutine	DRAWFOOT(jifnew)
C----------------------------------------------------------------------|
C YRBMN maximum radius of the footprint 
C YRBMX minimum radius of the footprint 
C YHBM  the upshift of the beam footprint
C YASP  the aspect ratio of the beam footprint
C YQ    the beam power
C NBFILE,XWW,XWH,DYLET,IY0,IYM,IDX,IDT,SCM are taken from "outcmn.inc"
C----------------------------------------------------------------------|
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'			! Parameter CNB1 is used
	include 'for/status.inc'		! Array WORK1 is used
	include 'for/outcmn.inc'
        double precision  YS0,YSC8,YRBMN,YRBMX,YHBM,YASP,YH,YQ
	character STRI*16
	logical	EXI
	integer	j,jj,jifnew,JL,JN,PT(10),length,ERCODE
C----------------------------------------------------------------------|
C	write(*,*)"Drawing footprint",jifnew
	call	createpixmap(2)	! All calls except the first are ignored
	jn = length(NBFILE)
	open(2,file=NBFILE(1:jn),status='OLD')
	JN = 0
 1	continue
	if(JN .eq. anint(CNB1))	goto	2
	JN = JN+1
C Read a record for one beam source
	call	STREAD(2,20,WORK1(1,JN),ERCODE)
C	write(*,*)">>> NB drawing: After STREAD: ercode =",ERCODE
C	write(*,*)JN
C	write(*,'(1P,5(E12.4))')(WORK1(j),j=1,20)
	goto	(1,94,95,96,97,98),ERCODE+1
 2	close(2)
	JN = 0
	JL = 0
	if (jifnew .eq. 0)	call	redraw(2)	! Erase previous
	call	cleare(2,0,0,XWW-1,XWH-1)
	call	colovm(14)			! Footprint (magenta)
C	write(*,*)jifnew
 3	continue
	JN = JN+1				! The ordinal beam number
	YQ = WORK1(1,JN)
	if (YQ .lt. 1.d-2)	goto	4
	YHBM = WORK1(11,JN)
	YASP = WORK1(15,JN)
	YRBMX = WORK1(12,JN)
	YRBMN = WORK1(13,JN)
C Beam footprint drawing:
	JL = JL+1
	if (JN .le. 9)	then
	   write(STRI(1:2),'(1X,I1)')JN
	else
	   write(STRI(1:2),'(I2)')JN
	endif
	STRI(3:4) = '  '
	write(STRI(5:10),'(1F6.3)')YQ
C	write(*,*)STRI(1:10)
C	write(*,*)"Beam",JN,"   Line",JL,"   Power",YQ,"   Time",TIME
	call	textnb(16,35-3+JL*DYLET,STRI(1:10),10)
	call	textnb(10,35-3,'Beam  Power',11)

	YS0 = (IY0+IYM)/2
	YSC8 = IDX*IDT/SCM
	YH = 0.5*YASP*(YRBMX-YRBMN)
	PT(1) = max(YRBMN,0.d0)*YSC8
	PT(2) = YS0+(-YHBM-YH)*YSC8
	PT(3) = YRBMX*YSC8
	PT(4) = PT(2)
	PT(5) = PT(3)
	PT(6) = YS0+(-YHBM+YH)*YSC8
	PT(7) = PT(1)
	PT(8) = PT(6)
	PT(9) = PT(1)
	PT(10) = PT(2)

	j=2
	jj=5
	call drawline(j,PT,jj)
C Use if footprints are not rectangles: (then replace PT(j) -> 10.*PT(j))
C	do	j=1,2*5-3,2
C	   call d2line(2,PT(j),PT(j+1),PT(j+2),PT(j+3))
C	enddo
 4	continue
	if(JN .lt. anint(CNB1))	goto	3
	return
 94	j = length(NBFILE)
	write(*,*)'>>> NB drawing >>> Error in file "',NBFILE(1:j),'"'
     >		 , ': unrecognized variable name'
        call	a_stop
 95	j = length(NBFILE)
	write(*,*)'>>> NB file "',NBFILE(1:j),'" read error'
        call	a_stop
 96	j = length(NBFILE)
	write(*,*)'>>> NB file "',NBFILE(1:j),'" wrong format:'
        write(*,*)'            More records expected than available.'
        call	a_stop
 97	write(*,*)'>>> NB drawing STREAD: array out of limits'
        call	a_stop
 98	j = length(NBFILE)
	write(*,*)'>>> NB file "',NBFILE(1:j),'" wrong format:'
        write(*,'(A,I2,2A,I2,A)')'     ',JN-1,' beam records ',
     >		'available, ',int(CNB1),' records required.'
        call	a_stop
	end
C======================================================================|
       subroutine	DRFOOT1(JN,YRBMN,YRBMX,YHBM,YASP,YQ)
C----------------------------------------------------------------------|
C JN is the ordinal beam number,
C YRBMN maximum radius of the footprint 
C YRBMX minimum radius of the footprint 
C YHBM  the upshift of the beam footprint
C YASP  the aspect ratio of the beam footprint
C YQ    the beam power
C----------------------------------------------------------------------|
	implicit none
	include 'for/parameter.inc'
	include 'for/const.inc'
	include 'for/outcmn.inc'
	character STRI*16
        double precision	YS0,YSC8,YRBMN,YRBMX,YHBM,YASP,YH,YQ
	integer	JL,JNB,JN,PT(10),j,jj
	save	JL,JNB
	data	JNB/0/
C----------------------------------------------------------------------|
        if (MOD10 .ne. 8)	return
	if (JNB .eq. 0)		then
C	   write(*,*)"1st entry"
	   call	createpixmap(2)
	   call	createpixmap(5)
	   JNB = anint(CNB1)
	   goto	1
	endif
	if (JN .eq. JNB .and. JL .eq. 0 .and. YQ .lt. 0.01)	then
C	   write(*,*)"All erased"
	   call	redraw(2)			! Erase everything when
	   call	cleare(2,0,0,XWW-1,XWH-1)	! the total power is 0
	   return				! 
	endif
 1	continue
	if (JN .eq. 1)	JL = 0
	if (YQ .lt. 0.01)	return

C	if (JL .eq. 0)	write(*,*)
C	write(*,*)JN,JL,YQ,TIME
	call	savepm(2,5)
	call	redraw(5)			! Erase previous drawing
	call	colovm(14)			! Footprint (magenta)
	if (JL .eq. 0)	then
	   call	cleare(2,0,0,XWW-1,XWH-1)	! Create pixmap from scratch
	   call textnb(10,35-3,'Beam  Power',11)
	endif
	if (JN .le. 9)	then
	   write(STRI(1:2),'(1X,I1)')JN
	else
	   write(STRI(1:2),'(I2)')JN
	endif
	STRI(3:4) = '  '
	write(STRI(5:10),'(1F6.3)')YQ
	JL = JL+1
	call textnb(16,35-3+JL*DYLET,STRI(1:10),10)
	YS0 = (IY0+IYM)/2
	YSC8 = IDX*IDT/SCM
	YH = 0.5*YASP*(YRBMX-YRBMN)
	PT(1) = max(YRBMN,0.d0)*YSC8
	PT(2) = YS0+(-YHBM-YH)*YSC8
	PT(3) = YRBMX*YSC8
	PT(4) = PT(2)
	PT(5) = PT(3)
	PT(6) = YS0+(-YHBM+YH)*YSC8
	PT(7) = PT(1)
	PT(8) = PT(6)
	PT(9) = PT(1)
	PT(10) = PT(2)

	j=2
	jj=5
	call drawline(j,PT,jj)
	do	j=1,10-3,2
	   write(*,*)j
	enddo
C Use if footprints are not rectangles: (then replace PT(j) -> 10.*PT(j))
C	do	j=1,2*5-3,2
C	   call d2line(2,PT(j),PT(j+1),PT(j+2),PT(j+3))
C	enddo
	call	redraw(2)			! Draw new
	call	cleare(5,0,0,XWW-1,XWH-1)	! Erase auxiliary pixmap[5]
	call	redraw(0)
	end
C======================================================================|
	subroutine SETARX(ICALL)
C----------------------------------------------------------------------|
!
!	Add treatment for NGRIDX()=1
!
C All arrays are mapped to the WHOLE radial grid [1,NB1]
C This can cause an inconsistency when NA1 varies in time.
C The time evolution of the input data is taken from
C	 - data array		|  for arrays
C
C Then it is stored for the current time in the arrays
C 	EXT(NRD,NARRX)	- (description in the file for/status.inc)
C----------------------------------------------------------------------|
	implicit  none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer jj,N1,N2,ICALL,jar,jstim,jto,jentim,kn,j3,jtn,jt,jt0
	integer jx,jy,NP1,j1,N11
	double precision	REX(NRDX),QEX(NRDX),XA(NRD),DA(NRD)
	double precision	QUADIN,RORZ,RZ2A,YDT,YDTA,YDTB,Y,Y1
C----------------------------------------------------------------------|
C  NARRX    maximal number of arrays readable from a data file 
C  NTARR    maximal number of time slices for all arrays (total)
C  NGR	    number of actually defined groups (grid + data)
C----------------------------------------------------------------------|
C Input
C	ICALL	= 0 - call from REVIEW (no transfer to EXT(,) is needed)
C		> 0 - call from STEPON
C		= 1 - time interpolation off
C		= 2 - time interpolation on
C	DATARR(NRDX*NTARR) - data array
C   Let  1 <= j <= NGR is an ordinal number of a group in DATARR
C	TIMEX(j)  - time for this group
C	NTYPEX(j) - type of grid for the group j
C	NGRIDX(j) - number of grid points for the group j
C	GDEX(j)   - pointer (in DATARR) to the grid for the group j
C	GDEY(j)   - pointer (in DATARR) to the data for the group j
C	KTO(j)	  - pointer (ordinal number) in the array EXARNM
C		so that EXARNM(KTO(j)) gives the name of the quantity j
C   Let   1 <= kn <= NARRX is an ordinal number of q-ty EXARNM(kn)
C	KOGDA(kn)  - pointer to a position in the array TIMEX
C	EXARNM(kn) - name*6 of the quantity kn
C Output
C	IFDFAX(kn)   - current pointer to data set in DATARR
C	NPTM(kn)     - number of data points within a<=AB
C	XAXES(jj,kn) - "radial" grid for displayed data
C	DATAX(jj,kn) - array for displayed data
C	EXT(jj,kn)   - smoothed curve
C----------------------------------------------------------------------|
C Pointer is returned to the root window after calling ESC
	call	markloc("SETARX"//char(0))
C	write(*,*)"SETARX Entry",TE(1)
	if (NGR .eq. 0)	return
	do  50	jar=1,NGR
	   jj = 0
	   if (jar .lt. NGR)	then
	      if (KTO(jar+1) .ne. KTO(jar))	then
		 jj = jar	! jj -> group end
		 jentim = KOGDA(KTO(jar+1))-1 ! jentim -> last time
	      endif
	   else
	      jj = jar
	      jentim = NGR
	   endif
	   if (jj .eq. 0)	goto	50
	   KN = KTO(jj)
	   jstim = KOGDA(KN)
	   jto = jstim
	   do	j3=jstim+1,jentim
	      if (TIMEX(j3) .lt. TIMEX(j3-1))	  goto	99
	      if (TIMEX(j3) .eq. TIMEX(j3-1))	  goto	98
	      if (TIMEX(j3) .le. time)	  jto = j3
	   enddo
	   IFDFAX(KN) = jto
	   jtn = min(jto+1,jentim)
	   if (time .le. TIMEX(jstim))	jtn = jstim
	   jt  = jtn
	   if (2.*TIME .gt. TIMEX(jtn)+TIMEX(jto))	jt  = jto
	   jt0 = 0
	   if (ICALL .gt. 1 .or. jto .eq. jtn)	goto	1
C	only one run needed
	   jt = jto
	   if (2.*TIME .gt. TIMEX(jtn)+TIMEX(jto))	jt  = jtn

 1	   continue
	   N1 = NGRIDX(jt)
	   N2 = NTYPEX(jt)
	   jx = GDEX(jt)
	   jy = GDEY(jt)

C	write(*,*)N2,KN,EXARNM(KN)
C	if (EXARNM(KN).eq.'CAR1X ') write(*,*)N2
	if (N2.gt.31 .and. EXARNM(KN).eq.'CAR1X   ') then
C	if (EXARNM(KN).eq.'TIX ' .or. EXARNM(KN).eq.'NEX ') then
C	if (EXARNM(KN).eq.'CAR1X') then
C	if (EXARNM(KN).eq.'TIX ') then
	   write(*,'(1A12,1A6,1A1,1A21,1I4,1A3,1I4,1A20,1I3,1A3,1I3)')
     >		'Quantity: "',EXARNM(KN),'"'
     >		,',  Records_in_TIMEX:',jstim,' ->',jentim
C          write(*,*)'All time slices:'
C	   write(*,'(6F9.3)')(TIMEX(j3),j3=jstim,jentim)
	   write(*,'("Times: ",3F9.5,",  Selected",F9.5,3X,3I5)')
     >		TIMEX(jto),time,TIMEX(jtn),TIMEX(jt),jt,jto,jtn
	   write(*,*)'Selected group & address',jt,jx,
     >		',   Grid_type:',N2,',   Grid_size:',N1
	   if (N2 .ge. 10)	then
	      write(*,*)'Grid:'
	      write(*,'(1P,6E13.5)')(DATARR(j3),j3=jx,jx+N1-1+N2/18)
	   endif
	   write(*,*)'Data:'
	   write(*,'(1P,6E13.5)')(DATARR(j3),j3=jy,jy+N1-1)
	endif
C----------------------------------------------------------------------|
C The following is done below:
C	(1) The grid in "a", XA(NP1), and the data DA(NP1) on this grid
C	    are defined by 
C		(i)  mapping the original grid to the "a" grid REX(N11)
C		(ii) transfer (SMOOTH) from {REX(N11),QEX(N11)} to {XA,DA} 
C	(2) XAXES(N1,KN) is defined which is as "a" grid for exp-dot plots
C	    DATAX(N1,KN) data on this grid
C	(3) EXT(NRD,KN) smoothed input arrays interpolated in time
C----------------------------------------------------------------------|
	if (N2 .ge. 21)	then
	    write(*,*)'>>> Warning:    Quantity  ',EXARNM(KN),
     >		' Unknown input type =',N2,',  data ignored'
	    goto	49
	endif

	if (N2 .ge. 10)	goto	9
C Set equidistant X-grid
	do	j3=1,N1
	    REX(j3) = (j3-1.)/(N1-1.)
	    QEX(j3) = DATARR(jy+j3-1)
	    DATAX(j3,KN) = QEX(j3)
	enddo
	NP1 = NA1
C----------------------------------------------------------------------|
	if (N2 .eq. 0)	then
C	    write(*,*)
C     >		'Input data for "',EXARNM(KN),'" ',
C     >		' use equidistant grid in "a": 0 <= a <= AB',
C     >	'                                        1 <= j <= NAB'
	    NP1 = NAB
	    do	j3 = 1,NP1
		XA(j3) = AMETR(j3)/AB
	    enddo
	    do	j3=1,N1
		XAXES(j3,KN) = AB*REX(j3)
	    enddo
	    goto	8
	endif
C------------------------------------------------------ INTYPE = 0,1,2,3
	if (N2 .eq. 1)	then
C	    write(*,*)
C     >		'Input data for "',EXARNM(KN),'" ',
C     >		' use equidistant grid in "a": 0 <= a <= ABC',
C     >	'                                        1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = AMETR(j3)/ABC
	    enddo
	    do	j3=1,N1
		XAXES(j3,KN) = ABC*REX(j3)
	    enddo
	    goto	8
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 2)	then
C	    write(*,*)
C     >		'Input data for "',EXARNM(KN),'" ',
C     >		' use equidistant grid in "rho_tor": 0 <= RHO <= ROC',
C     >	'                                           1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = RHO(j3)/ROC
	    enddo
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 3)	then
	    write(*,*)
C     >		'Input data for "',EXARNM(KN),'" ',
C     >		' use equidistant grid in "rho_pol":',
C     >	'      sqrt{[Psi(j)-Psi(1)]/[Psi(NA1)-Psi(1)]}    1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = sqrt((FP(j3)-FP(1))/(FP(NP1)-FP(1)))
	    enddo
C AMETR(XA(1:NP1)) is given;	QUADIN=AMETR(REX(j)) is returned;
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 4)	then
	    write(*,*)
     >		'Input data for "',EXARNM(KN),'" ',
     >		'use equidistant grid in "rho_vol":',
     >	'          sqrt{[V(j)-V(1)]/[V(NA1)-V(1)]}    1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = sqrt((FP(j3)-FP(1))/(FP(NP1)-FP(1)))
	    enddo
	    write(*,*)	'Option is not implemented'
	    call	a_stop
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 5)	then
	    write(*,*)
     >		'Input data for "',EXARNM(KN),'" ',
     >		' use equidistant grid in "Psi":',
     >	'       [Psi(j)-Psi(1)]/[Psi(NA1)-Psi(1)]    1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = (FP(j3)-FP(1))/(FP(NP1)-FP(1))
	    enddo
	    write(*,*)	'Option is not implemented'
	    call	a_stop
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 6)	then
	    write(*,*)
     >		'Input data for "',EXARNM(KN),'" ',
     >		' use equidistant grid in "Phi":',
     >	'       [Phi(j)-Phi(1)]/[Phi(NA1)-Phi(1)]    1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = (FP(j3)-FP(1))/(FP(NP1)-FP(1))
	    enddo
	    write(*,*)	'Option is not implemented'
	    call	a_stop
	endif
C----------------------------------------------------------------------|
	if (N2 .gt. 6)	then
	    write(*,*)'GRIDTYPE =',N2
	    write(*,*)	'Option is not implemented'
	    call	a_stop
	endif
C----------------------------------------------------------------------|
	do	j3 = 1,N1
	   XAXES(j3,KN) = QUADIN(NP1,XA,AMETR,REX(j3),Y,j1)
	enddo
 8	continue
	N11 = N1
	goto	30

 9	continue
	if (N2 .ge. 20)	goto	20
C--------------------- INTYPE = 10,11,12,13,19 (u-file input) ---------|
C Input grid
	RORZ = DATARR(jx)
	if (N2 .eq. 18 .or. N2 .eq. 19)	jx = jx+1
	do	j3=1,N1
	    XAXES(j3,KN) = DATARR(jx+j3-1)
	    DATAX(j3,KN) = DATARR(jy+j3-1)
	    REX(j3) = XAXES(j3,KN)
	    QEX(j3) = DATARR(jy+j3-1)
	enddo
	N11 = N1
	NP1 = NA1
C----------------------------------------------------------------------|
	if (N2 .eq. 10)	then
C	    write(*,*)
C     >	    'Input data for "',EXARNM(KN),'" ',
C     >	    ' use arbitrary grid in "a": 0 <= a <= AB',
C     >	    '                                      1 <= j <= NAB'
	    NP1 = NAB
	    do	j3 = 1,NP1
		XA(j3) = AMETR(j3)/AB
	    enddo
 10	    if (REX(N11) .gt. (1.+.5/N1)*AB)	then
		N11 = N11-1
		if (N11 .gt. 1)    goto	10
	    endif
	    if (REX(N11)  .lt. (1.-.5/N1)*AB)	N11 = N1+1
	    REX(N11) = AB
	    do	j3=1,N11-1
		REX(j3) = REX(j3)/REX(N11)
	    enddo
	    goto	18
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 11)	then
C	    write(*,*)
C     >	    'Input data for "',EXARNM(KN),'" ',
C     >	    ' use arbitrary grid in "a": 0 <= a <= ABC',
C     >	    '                                      1 <= j <= NP1 = NA1'
	    do	j3 = 1,NP1
		XA(j3) = AMETR(j3)/ABC
	    enddo
 11	    if (REX(N11) .gt. (1.+.5/N1)*ABC)	then
		N11 = N11-1
		if (N11 .gt. 1)    goto	11
	    endif
	    if (REX(N11)  .lt. (1.-.5/N1)*ABC)	N11 = N1+1
	    REX(N11) = ABC
	    do	j3=1,N11-1
		REX(j3) = REX(j3)/REX(N11)
	    enddo
	    goto	18
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 12)	then
C	    write(*,*)
C     >	    'Input data for "',EXARNM(KN),'" ',
C     >	    ' use arbitrary grid in "rho": 0 <= RHO <= ROC',
C     >	    '                                      1 <= j <= NP1 = NA1'
	    do	j3 = 1,NP1
		XA(j3) = RHO(j3)/ROC
	    enddo
 12	    if (REX(N11) .gt. (1.+.5/N1))	then
		N11 = N11-1
		if (N11 .gt. 1)    goto	12
	    endif
	    if (REX(N11) .lt. (1.-.5/N1))	N11 = N1+1
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 13)	then
C	    write(*,*)
C     >	    'Input data use for "',EXARNM(KN),'" ',
C     >	    ' arbitrary normalized grid in "rho_pol":',
C     >	'      sqrt{[Psi(j)-Psi(1)]/[Psi(NA1)-Psi(1)]},   1 <= j <= NA1'
C Note! No check if the boundary grid value is = 1 !
	    do	j3 = 1,NP1
		XA(j3) = sqrt((FP(j3)-FP(1))/(FP(NP1)-FP(1)))
	    enddo
 13	    if (REX(N11) .gt. (1.+.5/N1))	then
		N11 = N11-1
		if (N11 .gt. 1)    goto	13
	    endif
	    if (REX(N11) .lt. (1.-.5/N1))	N11 = N1+1
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 14)	then
C	    write(*,*)
C     >	    'Input data use for "',EXARNM(KN),'" ',
C     >	    ' arbitrary normalized grid in "rho_Vol":',
C     >	'           sqrt{[V(j)-V(1)]/[V(NA1)-V(1)]},   1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = sqrt((VOLUM(j3)-VOLUM(1))/VOLUME)
	    enddo
 14	    if (REX(N11) .gt. (1.+.5/N1))	then
		N11 = N11-1
		if (N11 .gt. 1)    goto	14
	    endif
	    if (REX(N11) .lt. (1.-.5/N1))	N11 = N1+1
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 15)	then
C	    write(*,*)
C     >	    'Input data use for "',EXARNM(KN),'" ',
C     >	    ' arbitrary normalized grid in "Psi":',
C     >	'          [Psi(j)-Psi(1)]/[Psi(NA1)-Psi(1)],   1 <= j <= NA1'
	    do	j3 = 1,NP1
		XA(j3) = (FP(j3)-FP(1))/(FP(j3)-FP(NP1))
	    enddo
 15	    if (REX(N11) .gt. (1.+.5/N1))	then
		N11 = N11-1
		if (N11 .gt. 1)    goto	15
	    endif
	    if (REX(N11) .lt. (1.-.5/N1))	N11 = N1+1
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 16)	then
C	    write(*,*)
C     >	    'Input data for "',EXARNM(KN),'"  use arbitrary ',
C     >	    'dimensional grid in "Phi": 0 <= Phi <= \pi*BTOR*ROC^2',
C     >	    '                                      1 <= j <= NP1 = NA1'
C Note!  Usage of the normalized grid in Phi is undesirable for a SOL
C	 configuration because an uncertainty in Phi(edge) may result
C	 in a strong deformation of the normalized grid.
C        Therefore, below a mapping is done assuming that 
C	 the dimensional Phi coordinate is used!
	    do	j3 = 1,NP1
		XA(j3) = RHO(j3)/ROC
	    enddo
	    Y = 1./(GP*BTOR)
	    do	j3=1,N1
		REX(j3) = sqrt(REX(j3)*Y)
	    enddo
 16	    if (REX(N11) .gt. ROC*(1.+.5/N1))	then
		N11 = N11-1
		if (N11 .gt. 1)    goto	16
	    endif
	    if (REX(N11) .lt. ROC*(1.-.5/N1))	N11 = N1+1
	    do	j3=1,N11-1
		REX(j3) = REX(j3)/ROC
	    enddo
	endif
C----------------------------------------------------------------------|
	if (N2 .eq. 17)	then
	    write(*,*)'GRIDTYPE =',N2
	    write(*,*)	'Option is not implemented'
	    call	a_stop
	endif
	if (N2 .ge. 18)	goto	19

	do	j3 = 1,N1
	   XAXES(j3,KN) = QUADIN(NP1,XA,AMETR,REX(j3),Y,j1)
	enddo
 18	continue
	REX(N11) = 1.
	QEX(N11) = DATAX(min(N1,N11),KN)

	goto	30
 19	continue
C------------------------------------------------------- INTYPE = 18 --|
	if (N2 .eq. 18)	then
C	    write(*,*)
C     >'Data set is given at a fixed distance from the major torus axis'
C     >		,'     {r=',RTOR,', z=z_i}, {z_i}-array follows'
C	    write(*,'(1P,6E13.5)')(REX(j3),j3=1,N1)
	    do	j3=1,N1
C		XAXES(j3,KN) = RZ2A1(RORZ,REX(j3))
		XAXES(j3,KN) = RZ2A(RORZ,REX(j3),NP1)
		REX(j3) = XAXES(j3,KN)
	    enddo
	endif
C------------------------------------------------------- INTYPE = 19 --|
	if (N2 .eq. 19)	then
C	    write(*,*)
C     >'Input data are mapped to the mid-plane: {r, z=0}, r-data follow'
C	    write(*,'(1P,6E13.5)')RORZ,ABC,AB,RTOR
C	    write(*,'(1P,6E13.5)')(REX(j3),j3=1,N1)
	    do	j3=1,N1
CC		XAXES(j3,KN) = RZ2A1(REX(j3),RORZ)
		XAXES(j3,KN) = RZ2A(REX(j3),RORZ,NP1)
		REX(j3) = XAXES(j3,KN)
	    enddo
C	    write(*,'(1P,6E13.5)')(REX(j3),j3=1,N1)
C	    read(*,*)
	endif
C	if (N2 .ne. 18 .and. N2 .ne. 19)	goto	30	!<= 17
	goto	21

 20	continue
C------------------------------------------------------- INTYPE = 20 --|
	if (N2 .eq. 20)	then
C Input data are given on {r,z} plane
	   N11 = N1
	   NP1 = NAB
C	   write(*,*)'Input data (dots to plot):'
C	   write(*,'(1P,6E13.5)')(DATARR(jx-1+j3),j3=1,N1)
C	   write(*,'(1P,6E13.5)')(DATARR(jx+N1-1+j3),j3=1,N1)
C	   write(*,'(1P,6E13.5)')(DATARR(jy-1+j3),j3=1,N1)
	   do	j3=1,N1
	      Y = DATARR(jx+j3-1)
	      Y1 = DATARR(jx+N1+j3-1)
	      XAXES(j3,KN) = RZ2A(Y,Y1,NP1)
	      DATAX(j3,KN) = DATARR(jy+j3-1)
	      REX(j3) = XAXES(j3,KN)
	      QEX(j3) = DATARR(jy+j3-1)
	   enddo
C	   write(*,*)'Input data are mapped to the mid-plane:',
C     >			' {r, z=0}, r-data follow'
C	   write(*,'(1P,6E13.5)')(REX(j3),j3=1,N1)
C	   write(*,'(1P,6E13.5)')(XAXES(j3,KN),j3=1,N1)
	endif
 21	call	SORTAB(REX,QEX,N1)
C	call	SORTAB(XAXES(1,KN),DATAX(1,KN),N1)
 22	if (REX(N11) .gt. (1.+.5/N1)*AB)	then
	    N11 = N11-1
	    if (N11 .gt. 1)    goto	22	
	endif
	if (REX(N11)  .lt. (1.-.5/N1)*AB)	N11 = N1+1
	QEX(N11) = QEX(min(N1,N11))
	do	j3=1,N11-1
	    REX(j3) = REX(j3)/AB
	enddo
	REX(N11) = 1.
	NP1 = NAB
	do	j3 = 1,NP1
	    XA(j3) = AMETR(j3)/AB
	enddo
C----------------------------------------------------------------------|
 30	NPTM(KN) = min(N1,N11)
	TOUTX(KN) = TIMEX(jt)
	if (ICALL .eq. 0)	goto	49

C RECTAN takes QEX values from the nearest grid point
C Effectively the grid function is replaced with a step-function
C The procedure is not suitable when a derivative of QEX is needed
C	do   j3=1,NP1
C	   YD(j3) = RECTAN(N11,REX,QEX,XA(j3))
C	enddo
C       ALFA = 0.001
C	call	SMOOTH(ALFA,NP1,YD,XA,NP1,DA,XA)

C This is added to avoid too long extrapolation to the magnetic axis
	if ( N11 .gt. 1 )	then
	   if ( REX(2)-REX(1) .lt. REX(1) )	REX(1)=0.
	endif

C	if (EXARNM(KN).eq.'NEX   ') then
	if (N2.eq.31 .and. EXARNM(KN).eq.'NEX   ')  then
	   write(*,*)
	   write(*,*)TIMEX(jt),TIMEX(jtn),TIMEX(jto),jt,jtn,jto,KN
	   write(*,*)'Selected original data for ',EXARNM(KN),
     >		',  Type =',N2,',  Grid =',N1,'  ->',N11,jx
	   write(*,*)'Grid:'
	   write(*,'(1P,6E13.5)')(DATARR(jx+j3),j3=0,N1-1)
	   write(*,*)'Data:',NEX(1)
	   write(*,'(1P,6E13.5)')(DATARR(j3),j3=jy,jy+N1-1)
	   write(*,*)'Input data (dots to plot):'
	   write(*,'(1P,6E13.5)')(XAXES(j3,KN),j3=1,N1)
	   write(*,'(1P,6E13.5)')(DATAX(j3,KN),j3=1,N1)
	   write(*,*)'Re-ordered data:'
	   write(*,'(1P,6E13.5)')(AB*REX(j3),j3=1,N1)
	   write(*,'(1P,6E13.5)')(QEX(j3),j3=1,N1)
	   write(*,*)'Selected profile: ',EXARNM(KN),jy,jstim,jt,N1,N11
	   write(*,'(1P,6E13.5)')(QEX(j3),j3=1,N11)
	   write(*,*)'From grid of',N11,' nodes'
	   write(*,'(1P,6E13.5)')(REX(j3),j3=1,N11)
	   write(*,*)'To grid of',NP1,' nodes',ABC,AB,AMETR(NP1),NAB
	   write(*,'(1P,6E13.5)')(XA(j3),j3=1,NP1)
C	   write(*,'(1P,6E13.5)')(EXT(j3,KN),j3=1,NP1)
	endif

C All input data are mapped to the grid XA(1:NP1) in the variable "a"
	call	SMOOTH(FILTER(jt),N11,QEX,REX,NP1,DA,XA)

C Print the result of mapping:
	if (N2.gt.30 .and. EXARNM(KN).eq.'TIX   ')  then
C	if (EXARNM(KN).eq.'MUX   ')  then
	    write(*,*)EXARNM(KN)
	    write(*,'(1P,6E13.5)')(XA(j3),j3=1,NP1)
	    write(*,'(1P,6E13.5)')(DA(j3),j3=1,NP1)
C	    write(*,'(1P,6E13.5)')(EXT(j3,KN),j3=1,NP1)
	endif
C----------------------------------------------------------------------|
C			 data interpolation
C jto - pointer to the previous time
C jtn - pointer to the subsequent time
C jt  - pointer to the current time
C	if jto=/=jtn two runs are accomplished: with jt=jtn and jt=jto
C	the order of runs depends on the time selected for exp output
C	namely, the last run determines current XAXES and DATAX

	if (jto .eq. jtn)	then	! no time dependence
	    do	j3=1,NRD
		EXT(j3,KN) = DA(j3)
	    enddo
	    goto	49
	endif

	if (jt0 .eq. 0)	then
	    jt0 = 1
	    if (jt .ne. jto)	then
		jt = jto
	    else
		jt = jtn
	    endif
	    do	j3=1,NRD
		EXT(j3,KN) = DA(j3)
	    enddo
	    goto 1
	endif
	ydt  = (TIMEX(jtn)-TIMEX(jto))
	ydta = (TIMEX(jtn)-time)/ydt
	ydtb = (time-TIMEX(jto))/ydt
	if (jto .eq. jt)	then
	    do	j3=1,NRD
		EXT(j3,KN) = EXT(j3,KN)*ydtb+DA(j3)*ydta
	    enddo
	else
	    do	j3=1,NRD
		EXT(j3,KN) = EXT(j3,KN)*ydta+DA(j3)*ydtb
	    enddo
	endif
 49	continue
 50	continue
C	write(*,*)"NEX:"
C	write(*,'(1P,6E13.5)')(NEX(j3),j3=1,6)
	return

 98	write(*,*)'Quantity  ',
     >		EXARNM(KN),' Input times repeated'
	call	a_stop 
 99	write(*,*)'Quantity  ',
     >		EXARNM(KN),' Input times out of order'
	call	a_stop 
	end
C======================================================================|
