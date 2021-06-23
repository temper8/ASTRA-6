C=======================================================================
C Modules >> IFKEY, STUFF, ADDMOD, LINEAV << 
C===================================================10.03.95 G.P.=======
	integer function	IFKEY(IFKL)
C----------------------------------------------------------------------|
C IFKL = 257	floating point exception event. IFKEY(257) 
C			is called from signal_handler for/serv.c
C IFKL = 256	call from initial iteration loop,
C       		move to KEY analysis skipping D[PRT]OUT checks,
C			no drawing unless in "DSP" mode
C IFKL = 0	standard call
C
C IFKL = -1	set PAUSE mode
C
C IFKL = ichar('G')	equivalent to pressing "G" in interactive mode
C
C The following logical units are used
C 1 - for files equ/log/model (once on entry and on request "I")
C		equ/txt/model (on request "L")
C 12,13 - for equ/model.log file (once on entry)
C 3 - for post-viewer file (first on entry, then periodically)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	include	'for/timeoutput.inc'
C  ICVMX    max No. of curves in one screen (no more than 32)
	integer ICVMX
	parameter(ICVMX=32)
	double precision PRMARK(NTIMES)
C PT dimension: 4*NRD(Mode 5,8) 320(7) 2*NTIMES(Mode 6) 2*NRD(Modes 1-4)
	integer PT(4*NRD+2*NTIMES),ITO(NTIMES,ICVMX+2)
	character	NAMEP(NTIMES)*6,OUTNAME(NRW)*8
	character*10	UNAMES(NRW),DEFUNA,CNSFIL*40,STRI*132,PNMNAME*80
	character*80	HELP(28),PSNAME,STR,STRB,ANAMES(NRW),BNAMES(NRW)
	integer	IFKL,IFKEY_,POLLEVENT,OUTFIG(NRW)
	integer LENGTH,KILLBL,INDEX,WAITEVENT,KIBM,KASCII,getatr,ASKTAB
	integer MARK,J,JJ,NNN,LTOUTO,JTOUT,LRJJ,IDSP,IERR,MOVIE,Magenta,
     1		LWRN,LRDN,IFLAG,INT4,IRET,NTRUN,IM,XSC0,XSC,NX,NY,
     2		NN80,NN0,NN1,FLEN,MODEX,IX,IY,NU1,j2,J1,JIT,WarningColor
	logical MODADD,IFPPF,ISTORY
	integer*2 YWD(NRD)
	double precision
     1		DEVARO(NCONST),LINEAV,TIMOD4(NTIMES),CHORDN,Y,ABD,ALFA,
     2		YWA(NRD),YWB(NRD),YWC(NRD),TIMEB,TROUT,TPOUT,SWATCH
	save	ITO,IFLAG,TROUT,TPOUT,MARK,JIT,LTOUTO,IDSP,WarningColor
	save	NAMEP,PNMNAME,MOVIE,Magenta
!	save	OUTFIG,OUTNAME
	data	PRMARK/NTIMES*0./  TROUT/-99999./     TPOUT/-99999./ 
     1		IFPPF/.FALSE./	   IFLAG/0/  DEFUNA/'      .tmp'/
     2		NN0/0/   NN1/1/    NN80/80/       LRJJ/426/
     3		JTOUT/0/ LTOUTO/0/ MARK /0/       IDSP/0/	JIT/0/
     4		WarningColor/30/   OUTFIG /NRW*0/ PNMNAME/'***'/
     4		Magenta/14/
C ASCII codes: ^C 3  <Esc>27 <Space>32  % 37  * 42  . 46  / 47  ? 63
C	  0 48  1 49  2 50  3 51  4 52  5 53  6 54  7 55  8 56  9 57
C	  A 65  B 66  C 67  D 68  E 69  F 70  G 71  H 72  I 73  J 74
C	  K 75	L 76  M 77  N 78  O 79  P 80  Q 81  R 82  S 83  T 84
C	  U 85  V 86  W 87  X 88  Y 89  Z 90
	data	(HELP(j),j=1,10)/
     1	'Esc or / - STOP',
     2	'<space>  - Set Pause mode, one-step advance',
     3	' <CR>    - Return to Run mode',
     4	'  0      - No graphic output',
     5	'  1,2,3  - Radial profile graphs',
     6	'  4,5    - Time evolution of radial profiles',
     7	'  6      - Time evolution of local/global quantities',
     8	'  7      - Trace of the discharge in a phase space',
     9	'  8      - Magnetic flux surfaces',
     &	'  .      - Curve style'/
	data	(HELP(j),j=11,20)/
     1	'  A      - Adjust colors (for a color monitor only)',
     2	'  B & N  - Backward & forward screen scan',
     3	'  C & V  - Constants & main Variables control',
     4	'  D      - Time step, output times & subroutine call control',
     5	'  F & P  - Store displayed curves in a file',
     6	'  G & Q  - Save the screen as a PostScript file',
     7	'  H & ?  - Show operative keys',
     8	'  I      - Save the model constants for the next run',
     9	'  L      - Model listing',
     &	'  M      - Mark time slice[s] in the modes 4, 5, 7'/
	data	(HELP(j),j=21,28)/
     1	'  R      - Refresh screen',
     2	'  S      - New scales',
     3	'  T      - Type numerical values of the current curves',
     4	'  U      - Write U-file; (1D in modes 1-3,6; 2D in modes 4,5)',
     5	'  W      - Re-arrange windows',
     6	"  X      - What's X-axis/grid?",
     7	'  Y      - Shift curve up/down',
     8	' '/
C----------------------------------------------------------------------|
	entry	IFKEY_(IFKL)
	NTRUN = NTIMES
	if (IFKL .eq. 257)	goto	97
	call	markloc("IFKEY"//char(0))
	CHORDN = LINEAV(1)
	IFKEY = 0
	IFKEY_ = 0
	if (IFKL .eq. -1)	then
C This sets "pause" mode each time when IFKEY(j=-1) is called
	   TASK(1:3)='DSP'
	   IFKL = 0
	endif
	if (IFKL .eq. 258)	then		! Call from INIT
	   TTOUT(1) = -1.d10
	   call	SETTSC(NTRUN)
	   return
	endif
	if (IFKL.gt.257 .or. IFKL.lt.0)	then
	   write(*,*)'>>> IFKEY: wrong input parameter. Call ignored.'
	   return
	endif
C	write(*,*)LTOUTO,TIME,DTOUT,int(TIME/DTOUT)
	if	(IFKL.gt.0 .and. IFKL.lt.256)	then
	   if (TASK(1:3) .eq. 'BGD')	then
	      write(*,*)"Illegal IFKEY parameter"
	      return
	   endif
	   JIT = 0
	   KEY = IFKL
	   goto	10
	elseif	(IFKL.eq.256)			then
	   if (TASK(1:3) .eq. 'BGD')	return
	   JIT = JIT+1
C	   write(*,*)'"',TASK(1:3),'"',JIT
	   call	PUTJIT(JIT)
C	   goto	1
	   goto	431	! Enable drawing during initial iiterations
	endif
	if (IFKL.ne.256 .and. TASK(1:3).ne.'BGD' .and. TASK(4:4).ne.'B')
     &		call TIMEDT(TIME,1000.*TAU)
	if (MOD10.gt.3.and.MOD10.lt.8)			goto 432
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
	if (TIME-TROUT+.5*TAU.lt.DROUT)			goto 432
 431	call	markloc("IFKEY (profile drawing)"//char(0))
	TROUT=TIME
	if (TASK(1:3).eq.'BGD' .or. TASK(4:4).eq.'B')	goto 432
C	call	ADDTIME(CPT)
	call	markloc("TIMOUT 1st call from IFKEY"//char(0))
	call	TIMOUT
	call	markloc("RADOUT 1st call from IFKEY"//char(0))
	call	RADOUT
	if (MOD10.eq.4 .or. MOD10.eq.5)	then
	   call	markloc("SMODE5 1st call from IFKEY"//char(0))
	   call SMODE5(MARK,PT,PRMARK,NAMEP,NTIMES)
	else
	   call	markloc("OUTDSP 1st call from IFKEY"//char(0))
	   call OUTDSP(MARK,0,PT,ITO,NTRUN,TTOUT,TOUT)
	endif
	call	markloc(" UPSTR 1st call from IFKEY"//char(0))
	call UPSTR(CHORDN,1./MU(NA))
	j = 0
	call	markloc(" DNSTR 1st call from IFKEY"//char(0))
	if (MOD10.le.7)	call DNSTR(j,NTRUN,TOUT)
	call	redraw(0)
	call	makemovie(MOVIE,PNMNAME)
	call	markloc(" "//char(0))
C	call	ADDTIME(CPTGRA)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
 432	continue
	if (LTOUT.eq.1)					goto 433
	if (TIME-TTOUT(LTOUT-1)+.5*TAU .lt. DTOUT)	goto 437
	call	markloc("IFKEY (saving time traces)"//char(0))
	if (LTOUT.lt.NTRUN)				goto 433
	do 	J	=1,NTRUN-1
	   do	JJ	=1,NTOUT
	      TOUT(J,JJ)=TOUT(J+1,JJ)
	   enddo
	   TTOUT(J)= TTOUT(J+1)
	enddo
	LTOUT	= NTRUN-1
 433	JJ = 0
 434	continue
C Moves the time window when the current time goes beyond the right border
	if (TIME.gt.TINIT+1.025*abs(TSCALE) .and. TSCALE.lt.0.)	then
	   TINIT = TINIT+.2*abs(TSCALE)
	   JJ = 1
	   goto	434
	endif
C Moves the time window when the current time stays beyond the left border
 435	continue
	if (TIME.lt.TINIT .and. TSCALE.lt.0.)	then
	   TINIT = TINIT-abs(TSCALE)
	   JJ = 1
	   goto	435
	endif
	call	markloc("TIMOUT 2nd call from IFKEY"//char(0))
	call	TIMOUT
	if (IFPPF) IFPPF = ISTORY(1,'TIMOUT',NTRUN,TTOUT,TOUT)
	TTOUT(LTOUT)=TIME
	LTOUT	= LTOUT+1
	JTOUT	= JTOUT+1
	if (TASK(1:3).eq.'BGD')	goto	437
C	if (TASK(1:3).eq.'BGD' .or. TASK(4:4).eq.'B')	goto	437
	if (MOD10.lt.6.or.MOD10.gt.7)	goto	437
C	call	ADDTIME(CPT)
	if (JJ .ne. 0)	goto	201 ! Renew plot if scale has been changed
	call	markloc("RADOUT 2nd call from IFKEY"//char(0))
	call	RADOUT
	call	markloc("TIMOUT 2nd call from IFKEY"//char(0))
	call	TIMOUT
	if (MOD10.eq.4 .or. MOD10.eq.5)	then
	call	markloc("TIMOUT 3rd call from IFKEY"//char(0))
	   call	markloc("SMODE5 2nd call from IFKEY"//char(0))
	   call SMODE5(MARK,PT,PRMARK,NAMEP,NTIMES)
	else
	   call	markloc("OUTDSP 2nd call from IFKEY"//char(0))
	   call OUTDSP(MARK,0,PT,ITO,NTRUN,TTOUT,TOUT)
	endif
	call	markloc(" UPSTR 2nd call from IFKEY"//char(0))
	call	UPSTR(CHORDN,1./MU(NA))
	j = 0
	call	markloc(" DNSTR 2nd call from IFKEY"//char(0))
	call	DNSTR(j,NTRUN,TOUT)
	call	redraw(0)
C 436	continue
C	call	ADDTIME(CPTGRA)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
 437	continue
C The next line suppresses writing a view file during the iteration loop
	if (IFKL .eq. 256)	goto	1
	if (TIME-TPOUT+.5*TAU.lt.DPOUT)	goto	1
	call	markloc("RADOUT 3rd call from IFKEY"//char(0))
	call	RADOUT
	if (LTOUTO.ne.0)		then
				!!! Logical units 1 & 2 are used
	   if (IFPPF) IFPPF = ISTORY(1,'RADOUT',NTRUN,TTOUT,TOUT)
	   j = length(FILEX)
	   if (j .ne. 0)	then
C Write ASCII JSP file
	      call OPENAP(3,FILEX(1:length(FILEX))//'a4jse',1,IERR)
	      if(IERR.gt.0)write(*,*)'>>> IFKEY: jse file append error'
	      call	AJSE(3,1)		! Append pre-jst file
	      close(3)
	      call OPENAP(3,FILEX(1:length(FILEX))//'a4jst',1,IERR)
	      if(IERR.gt.0)write(*,*)'>>> IFKEY: jst file append error'
	      call	AJST(3,1)		! Append pre-jst file
	      close(3)
	      call OPENAP(3,FILEX(1:length(FILEX))//'a4jsp',1,IERR)
	      if(IERR.gt.0)write(*,*)'>>> IFKEY: jsp file append error'
	      call	AJSP(3,1)		! Append pre-jsp file
	      close(3)
	   endif
	   call 	OPENAP(3,RSNAME,1,IERR)
	   if(IERR.gt.0)write(*,*)'>>> IFKEY: Review file append error'
	else
! Check of user's request for PPF !!! Logical units 1 & 2 are used
	   IFPPF = ISTORY(1,'INITAL',NTRUN,TTOUT,TOUT)
	   j = length(FILEX)
	   if (j .ne. 0)	then
C Write ASCII JSP file
	      call OPENWT(3,FILEX(1:length(FILEX))//'a4jse',1,IERR)
	      if(IERR.gt.1)write(*,*)'>>> IFKEY: jse file open error'
	      call	AJSE(3,0)		! open pre-jst file
	      close(3)
	      call OPENWT(3,FILEX(1:length(FILEX))//'a4jst',1,IERR)
	      if(IERR.gt.1)write(*,*)'>>> IFKEY: jst file open error'
	      call	AJST(3,0)		! open pre-jst file
	      close(3)
	      call OPENWT(3,FILEX(1:length(FILEX))//'a4jsp',1,IERR)
	      if(IERR.gt.1)write(*,*)'>>> IFKEY: jsp file open error'
	      call	AJSP(3,0)		! open pre-jsp file
	      close(3)
	   endif
C The file RSNAME='.res/profil.dat' (default name) is used in 7 places: 
C  here, outdsp (modes 4,5), typdsp (writing 2D U-file ),
C  in review*3, in listres
	   LWRN = length(EQNAME)
	   LRDN = length(RDNAME)
	   do	j=1,40
	      if (j.gt.LWRN) EQNAME(j:j)=' '
	      if (j.gt.LRDN) RDNAME(j:j)=' '
	   enddo
	   call		OPENRD(12,'equ/'//EQNAME(1:LWRN),0,IERR)
	   if(IERR.gt.1) write(*,*)'>>> IFKEY: Model file "equ/',
     >		EQNAME(1:LWRN),'" open error'
	   call 	OPENWT(3,RSNAME(1:length(RSNAME)),1,IERR)
	   if(IERR.gt.1) write(*,*)'>>> IFKEY: Review file open error'
	   CNSFIL='equ/log/'//EQNAME(1:LWRN)//char(0)
	   inquire(FILE=CNSFIL,EXIST=MODADD)
	   if (MODADD)	then
	      call	OPENRD(1,CNSFIL,0,IERR)
	      if (IERR.gt.1)  write(*,*)'>>> IFKEY: File "',
     >				CNSFIL(1:LWRN+8),'" open error'
	      call	ADDMOD(3,12,1)
	      close(1)
	   else
	      call	ADDMOD(3,12,0)
	   endif
	   close(12)
	   write(3)RDNAME,EQNAME,VERSION,XLINE1
     .	,	IYEAR,IMONTH,IDAY,IHOUR,IMINUT,NCFNAM,NPRNAM
     .	,	NROUT,(NAMER(J),J=1,NROUT),(SCALER(J),J=1,NROUT)
     .	,	NTOUT,(NAMET(J),J=1,NTOUT),(SCALET(J),J=1,NTOUT)
     .	,	HRO,NB1,NSBR,NGR,NXOUT,(LEQ(j),j=1,7)
C Note Change the cycle in NEQNS,(LEQ(j),j=1,NEQNS)
C Presently LEQ is not used by review.f and need not be stored
C     .	,	HRO,NB1,NSBR,NGR,NXOUT
	   if (NXOUT .gt. 0 .and. NGR .gt. 0)	then
C Total length: 3*NGR*int+(3*NGR+GDEY(NGR)+NGRIDX(NGR)-1)*real+3*NARRX*int
	      write(3)
     .		(KTO(j),j=1,NGR),(NGRIDX(j),j=1,NGR),(NTYPEX(j),j=1,NGR)
     .	,	(TIMEX(j),j=1,NGR),(GDEX(j),j=1,NGR),(GDEY(j),j=1,NGR)
     .	,	(DATARR(j),j=1,GDEY(NGR)+NGRIDX(NGR)-1)
     .	,	(NAMEX(j),j=1,NARRX),(NWINDX(j),j=1,NARRX)
     .	,	(KOGDA(j),j=1,NARRX)
	   endif
	endif
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
	LTOUTO = LTOUT-JTOUT
	write(3)JTOUT
	if (JTOUT.eq.0)		goto 438
	write(3)(TTOUT(J),(TOUT(J,JJ),JJ=1,NTOUT),J=LTOUTO,LTOUT-1)
	JTOUT	=0
 438	write(3)TIME
	write(3)(CONSTF(J),J=1,NCFNAM),(DEVAR(J),J=1,NPRNAM)
     .	,		ABC,ROC,CHORDN,1./MU(NA)
	write(3)NA1,NAB,(0,j=1,10),(0.d0,j=1,10)

	call	RHSEQ
	call	STUFF(3,AMETR,NAB,YWD)
	call	STUFF(3,SHIF,NAB,YWD)
	call	STUFF(3,ELON,NAB,YWD)
	call	STUFF(3,TRIA,NAB,YWD)
	call	STUFF(3,EQFF,NAB,YWD)
	call	STUFF(3,EQPF,NAB,YWD)
	call	STUFF(3,FP,NAB,YWD)
	do	J=1,NROUT
	   call	STUFF(3,ROUT(1,J),NAB,YWD)
	enddo
	close(3)
	TPOUT=TIME
	TIMOD4(IPOUT) = TPOUT
	call FMTF6(NAMEP(IPOUT),TPOUT)
	if (IPOUT.lt.NTRUN)  IPOUT=IPOUT+1
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++|
 1	continue
	call	markloc("IFKEY (loop)"//char(0))
	KEY = 0
	if (TASK(1:3).ne.'BGD' .and. TASK(4:4).ne.'B')	call  redraw(0)
C	write(*,*)TASK(1:3),JIT
C__________ Check Pause time condition
	if (TIME.ge.TPAUSE .and. IDSP.eq.0)	then
	   IDSP = 1
	   KEY = 32
	   goto	10
	endif
C	if (TIME.lt.TPAUSE)	IDSP = 0
C	call	ADDTIME(CPTGRA)
C__________ Check EXIT condition
	if (TIME-TEND+.1E-7 .ge. DPOUT+TAU)	goto	97
	if (TASK(1:3) .eq. 'BGD')   return
C	call	ADDTIME(CPT)
C__________ Polling events
	KEY = 0
	if (TASK(1:3) .eq. 'RUN')	then
	   KIBM = pollevent(KEY)
C KIBM = 1 - <Ctrl> was pressed
C KIBM = 2 - <Alt>  was pressed
C KEY=318 - the root window was closed
C KEY=322 and then KEY=319 - the root window is opened
C	   if (KEY.gt.317)	write(6,*)KEY
C	   if  (KEY  .eq. 0)	write(*,*)"Exit 0"
	   if  (KEY  .eq. 0)	return
	   if  (KIBM.eq.1.and.(KEY.eq.99.or.KEY.eq.67))	goto	97 ! <Ctrl>+C
	   call	markloc("IFKEY (key analysis)"//char(0))
C	   write(*,*)TASK,KIBM,KEY,TIME
C	   j=getatr(jj)
	   if  (KIBM.eq.65006)	then
	      TASK = 'RUN '
	      KIBM  = 0
	      return
	   endif
	   if  (KIBM.eq.65005)	then
	      TASK = 'RUNB'
C	      write(*,*)"UnmapNotify ",TASK
	      return
	   else
	      TASK = 'RUN '
C             write(*,*)"MapNotify   ",TASK,KIBM
	   endif
	   if (KIBM.eq.1 .or. KIBM.eq.2)	goto	49
	   if (KEY.ge.97 .and.KEY.le.122)	KEY=KEY-32
C	   write(*,*)KIBM,KEY,TIME
	endif
C__________ Waiting events
	if (TASK(1:3) .eq. 'DSP')	then
	   if (IFLAG.eq.1)	then
C IFLAG is equal 1 if
C	(1) <Space> is pressed 
C	(2) when a subroutine-call KEY is pressed.
C	One time step and re-drawing is done
	      KEY = 0
	      IFLAG = 0
	      call	markloc("RADOUT 4th call from IFKEY"//char(0))
	      call	TIMOUT
	      call	RADOUT
	      if (MOD10.eq.4 .or. MOD10.eq.5)	then
		 call	markloc("TIMOUT 4th call from IFKEY"//char(0))
		 call	markloc("SMODE5 3rd call from IFKEY"//char(0))
		 call	SMODE5(MARK,PT,PRMARK,NAMEP,NTIMES)
	      else
		 call	markloc("OUTDSP 3rd call from IFKEY"//char(0))
		 call	OUTDSP(MARK,0,PT,ITO,NTRUN,TTOUT,TOUT)
	      endif
	      call	markloc(" UPSTR 3rd call from IFKEY"//char(0))
	      call	UPSTR(CHORDN,1./MU(NA))
	      call	markloc("TIMEDT called from IFKEY"//char(0))
	      if (IFKL .ne. 256)	call	TIMEDT(TIME,1000.*TAU)
	      j = 0
	      call	markloc(" DNSTR 3rd call from IFKEY"//char(0))
	      if (MOD10 .le. 7)		call	DNSTR(j,NTRUN,TOUT)
	   endif
	   call	PUTXY(IX,IY,NTRUN,TTOUT,TOUT)
	   KASCII = 0
	   KIBM = waitevent(KASCII,ix,iy)
	   KEY  = KASCII
	   if (KEY .eq. 0)	goto	1
C	   if (KEY .eq. 27)	KEY = 32
	   if (KEY .lt.127 .and. (KIBM.eq.1 .or. KIBM.eq.2)) goto 49
	   if (KIBM .lt. 65000) goto	1
	   KIBM   = KIBM-65000
C	   if (KIBM .eq. 3)	write(*,*)TASK,KIBM	! Unmap
C	   if (KIBM .eq. 4)	write(*,*)TASK,KIBM	! Map
C	   if(KIBM.gt.1) write(*,*)'Integer ',KIBM+65000,' was returned'
C	   if(KIBM.eq.1) write(*,*)'A mouse button was pressed'
C	   if(KIBM.eq.0) write(*,*)'The key [',char(KASCII),'] was pressed'
C	   if(KIBM.eq.4) write(*,*)'The rootwindow was exposed',KASCII
C	   if(KIBM.eq.5) write(*,*)'The configuration was changed',KASCII
	   if(KIBM.ge.361.and.KIBM.le.364) then
	      call	mvcursor(KIBM,ix,iy)
	      goto	1
	   endif
C	   call	PUTXY(IX,IY,NTRUN,TTOUT,TOUT)
	   if (  KEY .ge .97) KEY = KEY-32
C	   write(*,*)'DSP   ',KASCII,KIBM,IFLAG,KEY
	endif
C 2	continue
C----------------------------------------------------------------------|
C__________'CR'
 10	if (KEY .eq. 13)	then
	   TASK='RUN '
	   call rcurso
	   call ERASXY(IX,IY)
			goto 1
	endif
C__________'%'
	if (KEY .eq. 37)	then
	   write(6,*)
	   call	wrtime(6,'  >>> Astra run time  ',22,swatch(Y),-1.d0)
	   call	CPUSE(6)
	   goto	1
	endif
C__________'P'
 	if (KEY.ne.80)	goto 12
C  Optional output format (G.W.Pacher)
	   call	markloc("TIMOUT 5th call from IFKEY"//char(0))
	   call	TIMOUT
	   call	markloc("TYPDSP called from IFKEY"//char(0))
	   call TYPDSP(1,CHORDN,NTRUN,TTOUT,TOUT)
			goto 1
C__________'A'
 12	if (KEY.ne.65)	goto 121
C		call SETCOL(COLTAB)
C		call colovm(11)
C			goto 201
C__________'N'
 121	if (KEY.ne.78)	goto 13
	   if (MOD10.ge.8 .or. MOD10.lt.0)	goto	201
		NSCR(MOD10) = NSCR(MOD10)+1
		J  = NROUT
		JJ = 8				! MOD10 = 2,3,6
		if (MOD10.eq.6 .or. MOD10.eq.7)	J = NTOUT
		if (MOD10 .eq. 1)	JJ=16
		if (MOD10 .eq. 4)	JJ=2
		if (MOD10 .eq. 5)	JJ=2
		if (MOD10 .eq. 7)	JJ=4
		if (J .gt. JJ*NSCR(MOD10))	goto	201
			NSCR(MOD10) = 0
			goto 201
C__________'B'
 13	if (KEY.ne.66)	goto 14
	   if (MOD10.ge.8 .or. MOD10.lt.0)	goto	201
		NSCR(MOD10) = NSCR(MOD10)-1
		if (NSCR(MOD10) .ge. 0)		goto	201
		J  = NROUT
		JJ = 8				! MOD10 = 2,3,6
		if (MOD10.eq.6 .or. MOD10.eq.7)	J = NTOUT
		if (MOD10 .eq. 1)	JJ=16
		if (MOD10 .eq. 4)	JJ=2
		if (MOD10 .eq. 5)	JJ=2
		if (MOD10 .eq. 7)	JJ=4
		NSCR(MOD10) = (J-1)/JJ
			goto 201
C__________'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'
 14	if (KEY.lt.48 .or. KEY.gt.57)	goto 16
		if (MOD10 .ge. 2 .and. MOD10 .le. 5)	then
		   if (MOD10 .eq. KEY-48)	then
		      MODEY = -MODEY
		   else
		      MODEY = 1
		   endif
		endif
		if (MOD10 .eq. 1 .and. KEY .eq. 49)	NSCR(MOD10) = 0
		if (MOD10 .eq. 6)	then
		   if (KEY .eq. 54)	then
		      MODEY = MODEY+1
		      if (MODEY.eq.2)	MODEY = -1
		   else
		      MODEY = 1
		   endif
		endif
		if (MOD10 .ne. KEY-48)	then
		   MOD10=KEY-48
C		   NSCR=0
		   call	erasrw
		   IM = 1
		   NST = 0
		   if (MOD10 .eq. 7)	call NEGA(IM,NTRUN,TOUT)
		   call	SETFRM(IM,XSC0,XSC,NX,NY)
		   call	AFRAME(IM,XSC0,XSC,NX,NY)
C		   call	AFRAME
		   j = XOUT+.49
		   call	ASRUMN(j)			! Task menu
		   call	textbf(0,XWH-104,RUNID,80)	! Task ID
		endif
		goto 201
C__________'F'
 16	if (KEY.ne.70)	goto 17
	   call	markloc("TIMOUT 6th call from IFKEY"//char(0))
	   call	TIMOUT
	   call TYPDSP(0,CHORDN,NTRUN,TTOUT,TOUT)
			goto 1
C__________'I'
 17	if (KEY.ne.73)	goto 18
	   LWRN	= length(EQNAME)
	   CNSFIL='equ/log/'//EQNAME(1:LWRN)//char(0)
	   call 	OPENWT(1,CNSFIL,0,IERR)
	   if(IERR.gt.1)	then
	   write(*,*)'>>> IFKEY: file "',CNSFIL(1:LWRN+8),'" open error'
	      goto	18
	   endif
C	goto 171
C MODEL.log file format for versions 5.2 (not used any more)
C	   write(1,*)'The description of the data see in the file ',
C     .				'"for/const.inc"'
C	   write(1,*)'Variables:',NPRNAM
C	   write(1,'(1P,8E11.3)')(DEVAR(J),J=1,NPRNAM)
C	   write(1,*)'Constants:',NCFNAM
C	   write(1,'(1P,8E11.3)')(CONSTF(J),J=1,NCFNAM)
C	   write(1,*)'Time control (common block A_OUTPUT):',24
C	   write(1,'(1P,4E11.3)')(DELOUT(j),j=1,24)
C	   goto	173
C 171	continue
C       New format of the equ/MODEL.log file for versions => 5.3
	   write(1,'(3(1A,1I1))')
     &		' Start file for version ',AVERS,'.',ARLEAS,'.',AEDIT
	   write(1,*)
     &		'Variables as described in the file "for/const.inc"'
	   do	J=1,NPRNAM
	      if (PRNAME(J) .eq. 'ZRD1  ')	goto	172
	      write(1,'(1A6,1A2,1P,8E11.3)')PRNAME(J),' =',DEVAR(J)
	   enddo
 172	   continue
	   write(1,'(A)')' Constants:'
	   do	J=1,NCFNAM
	      write(1,'(1A6,1A2,1P,8E11.3)')CFNAME(J),' =',CONSTF(J)
	   enddo
	   write(1,'(A,I2)')
     &		' Control parameters (common block A_OUTPUT):',22
	   do	J=1,22			! Don't save TPAUSE and TEND
	      write(1,'(1A6,1A2,1P,8E11.3)')SRNAME(J),' =',DELOUT(J)
	   enddo
	   write(1,*)'Color table (description in for/outcmn.inc):',32
	   write(1,'(4(2I4,3X))')(COLTAB(j),j=1,64)
	   close (1)
	   write(*,*)"Default start file is modified"
			goto 1
C__________'E'free
 18	continue
	if (KEY.ne.69)	goto 19
	call caution
			goto 1
C__________'.'
 19	if (KEY.ne.46)	goto 20
		MARK = MARK+1
		if (MARK .eq. 2) MARK = -1
			goto 201
C__________'R'
 20	if (KEY.ne.82)	goto 21
 201		continue
		if (IFKL.eq.256 .and. TASK(1:3).ne.'DSP')	goto 1
		call	erasrw
		IM = 1
		NST = 0
		if (MOD10 .eq. 7)	call NEGA(IM,NTRUN,TOUT)
		call	SETFRM(IM,XSC0,XSC,NX,NY)
		call	AFRAME(IM,XSC0,XSC,NX,NY)
C		call	AFRAME
		j = XOUT+.49
		call	ASRUMN(j)			! Task menu
		call	textbf(0,XWH-104,RUNID,80)	! Task ID
		call	markloc("RADOUT 5th call from IFKEY"//char(0))
		call	RADOUT
		call	markloc("TIMOUT 7th call from IFKEY"//char(0))
		call	TIMOUT
		if (MOD10.eq.4 .or. MOD10.eq.5)	then
		   call	markloc("SMODE5 4th call from IFKEY"//char(0))
		   call SMODE5(MARK,PT,PRMARK,NAMEP,NTIMES)
		else
		   call	markloc("OUTDSP 4th call from IFKEY"//char(0))
		   call OUTDSP(MARK,1,PT,ITO,NTRUN,TTOUT,TOUT)
		endif
		call	markloc(" UPSTR 4th call from IFKEY"//char(0))
		call	UPSTR(CHORDN,1./MU(NA))
		if (IFKL.ne.256)call	TIMEDT(TIME,1000.*TAU)
		call	markloc(" DNSTR 4th call from IFKEY"//char(0))
		j = 0
		if (MOD10.le.5 .or. MOD10.eq.7)	call DNSTR(j,NTRUN,TOUT)
		if (MOD10.eq.6 .and. KPRI.eq.0)	call DNSTR(j,NTRUN,TOUT)
		if (KPRI.ge.1 .and. KPRI.le.2)	then
		   STRI='>>>  The figure is stored in the file: '
     /					//PSNAME(1:FLEN)
		   write(*,*)STRI(1:39+FLEN)
		   call	colovm(WarningColor)
		   if (KPRI .eq. 1)	call	PCOVA
		   STRI='The figure is stored in the file: '
     /					//PSNAME(1:FLEN)
		   call	textvm(40,980,STRI,34+FLEN)
		   call	PSCLOS
		   KPRI = 0
		endif
		if (IFKL .eq. KEY)	return
		goto 1
C__________'U'
 21	if (KEY.ne.85)	goto 219
	if (MOD10.ge.7)	goto 1
	if (NUF.gt.NRD)  goto	91
C----------------------------------------------------------------------|
C Different options for radial array 
	if (MOD10.lt.1 .or. MOD10.gt.4)   goto   215
	call ASKUNA(NROUT,UNAMES,NAMER,DEFUNA)
	jj = 0
	do	j=1,NROUT
	   if (UNAMES(j).ne.DEFUNA)	jj=1
	enddo
	if (jj .eq. 0)	goto   1

	call	markloc("RADOUT 6th call from IFKEY"//char(0))
	call	RADOUT
	MODEX = XOUT+.49
	if (MOD10 .eq. 4)	goto   211
	if	(MODEX .eq. 1)	then
C			 Write up to ABC against "a"
	   NU1 = NA1
	   ABD = ABC
	   do	j = 1,NA1
		YWC(j) = AMETR(j)/ABC
	   enddo
	elseif	(MODEX .eq. 2)	then
C			 Write up to ROC against "rho"
	   NU1 = NA1
	   ABD = 1.
	   do	j = 1,NA1
		YWC(j) = RHO(j)/ROC
	   enddo
	else
C			 Write up to AB against "a"	(default)
	   if (MODEX .ne. 0)	write(*,*)
     ,	   ">>>  Warning: Unknown radial mode.  Writing anyway [0,AB]"
	   NU1 = NAB
	   ABD = AB
	   do	j = 1,NAB
		YWC(j) = AMETR(j)/AB
	   enddo
	endif
C Radial coordinate in a U-file in [m] (presently)
	if (MOD10 .eq. 1)	MODADD = .FALSE.
	if (MOD10 .ne. 1)	MODADD = .TRUE.
	do	210	jj=1,NROUT
	   if (UNAMES(jj).eq.DEFUNA)	goto 210
	   ALFA	=1.d-4
C Transfer to an equidistant radial grid
	   do	j=1,NUF
		YWA(j)	=(j-1.)/(NUF-1.)
	   enddo
	   call	SMOOTH(ALFA,NU1,ROUT(1,jj),YWC,NUF,YWB,YWA)
	   do	j=1,NUF
		YWA(j)	=YWA(j)*ABD
	   enddo
	   call	UF1DWA('AUGD',RUNID,UNAMES(jj),TIME,NAMER(jj),NUF,1,4,
     &		YWA,YWB,RTOR,AB,BTOR,IPL,CHORDN,MODADD,MODEX)
 210	continue
	goto   1
 211	continue
C-----------------------------------------------------------------------|
      if (MOD10 .ne. 4)   goto   215
C  Unknown option
	if (MODEX .ne. 0)	write(*,*)">>>  Warning: ",
     ,	   "2D U-file incompatible radial mode.  Writing anyway [0,AB]"
C  Write up to AB against "a"	(default)
	do	212	jj=1,NROUT
	if (UNAMES(jj).eq.DEFUNA)	goto 212
	MODADD = .TRUE.
	call  UF2DWA('AUGD',UNAMES(jj),NAMER(jj),jj,NUF,0,4,PRMARK,
     &	      TIMOD4,RTOR,AB,BTOR,IPL,CHORDN,MODADD,YWA,YWB,YWC)
 212	continue
	goto   1
C-----------------------------------------------------------------------|
 215  if (MOD10 .ne. 5)   goto   216
	write(*,*)
     >	'>>>  WARNING: U-file writing is not implemented in this mode'
	goto   1
C-----------------------------------------------------------------------|
 216	continue
      if (MOD10 .ne. 6)   goto   1
	call	markloc("TIMOUT 8th call from IFKEY"//char(0))
	call	TIMOUT
	call ASKUNA(NTOUT,UNAMES,NAMET,DEFUNA)
	do	218	jj=1,NTOUT
	if (UNAMES(jj).eq.DEFUNA)	goto 218
	if (MODEY.eq.1)	then
C Do not append the model
		MODADD = .FALSE.
	else
C The model is written to the end of U-file.
		MODADD = .TRUE.
	endif
	call	UF1DWA('AUGD',RUNID,UNAMES(jj),TIME,NAMET(jj),LTOUT-1,0,
     &	     4,TTOUT(1),TOUT(1,jj),RTOR,AB,BTOR,IPL,CHORDN,MODADD,MODEX)
 218	continue
			goto 1
C__________'W'
 219	if (KEY.ne.87)	goto 22
	if (MOD10.eq.1)			call ASKINT(NROUT,NWIND1,NAMER)
	if (MOD10.eq.2.or.MOD10.eq.3)	call ASKINT(NROUT,NWIND2,NAMER)
	if (MOD10.eq.4.or.MOD10.eq.5)	call ASKINT(NROUT,NWIND4,NAMER)
	if (MOD10.eq.6)			call ASKINT(NTOUT,NWIND3,NAMET)
	if (MOD10.eq.7)			call ASKINT(NTOUT,NWIND7,NAMET)
			goto 201
C__________'L'
 22	if (KEY.ne.76)	goto 23
	   LWRN	= length(EQNAME)
	   j1 = 0
	   do	j2=1,LWRN
	      if (EQNAME(j2:j2) .eq. '/')	j1 = j2
	   enddo
	   if (j1 .eq. 0)	then
	      CNSFIL = 'equ/txt/'//EQNAME(1:LWRN)//char(0)
	   else
	      CNSFIL = 'equ/'//EQNAME(1:j1)//'txt/'
     &			//EQNAME(j1+1:LWRN)//char(0)
	   endif
	   call 	OPENRD(1,CNSFIL,0,IERR)
	   if (IERR .gt. 0) write(*,*)
     >		'>>> IFKEY: "equ/txt/',EQNAME(1:LWRN),'" file error'
	   if (IERR .gt. 0)	goto 70
	   j2=0
 224	   if (j2 .lt. 0)	goto 226
	   do 221 J1=1,int(1.333*YSCMAX/DYLET)-1
	      if (j2 .lt. 0)	goto 225
	      read(1,'(1A80)',END=222)STR
	      goto 225
 222	      j2=-1
 225	      NNN=(J1-1)*DYLET+1
	      if (j2.lt.0)		then
		 STRB='                                  '
		 call prntxt(STRB)
	      else
		 call prntxt(STR)
	      endif
 221	   continue
	   goto	224
 226	   close (1)
	   goto 1
C__________'S'
 23	if (KEY.ne.83)	goto 231
		if (MOD10.le.5)	then
			call ASKLIS(NROUT,SCALER,NAMER,5)
		else
			call ASKLIS(NTOUT,SCALET,NAMET,5)
		endif
			goto 201
C__________'Y'
 231	if (KEY.ne.89)	goto 232
		if (MOD10.le.5)	then
			call ASKLIS(NROUT,OSHIFR,NAMER,6)
		else
			call ASKLIS(NTOUT,OSHIFT,NAMET,6)
		endif
			goto 201
C__________'X'
 232	if (KEY.ne.88)	goto 233
C	STR = "Test window"//char(0)
C	call	wintitle(STR)
	MODEX = XOUT+.49
	if	(MOD10 .eq. 0)	then
	   write(*,*)'X-axis:   none'
	elseif  (MOD10 .eq. 3)	then
	   write(*,112)'poloidal flux,  1 < j < NA1 =',NA1
	elseif  (MOD10 .eq. 4)	then
	   write(*,114)'0 < a < AB =',AB,'m     1 < j < NAB =',NAB
	elseif  (MOD10 .eq. 5)	then
	   write(*,112)'major radius in the mid-plane'
	elseif  (MOD10 .eq. 6)	then
	   write(*,*)'X-axis:   time [s]'
	elseif  (MOD10 .eq. 7)	then
	   write(*,*)'X-axis:   phase space'
	elseif  (MOD10 .eq. 8)	then
	   write(*,*)'X-axis:   major radius [m]'
	elseif  (MOD10 .eq. 9)	then
	   write(*,*)"User's plot"
	elseif  (MODEX .eq. 0)	then
	   write(*,114)'0 < a < AB =',AB,'m,     1 < j < NAB =',NAB
	elseif  (MODEX .eq. 1)	then
	   write(*,114)'0 < a < ABC =',ABC,'m,    1 < j < NA1 =',NA1
	elseif  (MODEX .eq. 2)	then
	   write(*,114)'0 < rho < ROC =',ROC,'m,    1 < j < NA1 =',NA1
	elseif  (MODEX .eq. 3)	then
	   write(*,112)'0 < Psi < FP(NA1),  1 < j < NA1 =',NA1
	else
	   write(*,*)'X-axis:   Unknown option'
	endif
 112	format('X-axis:   ',1A29,1I3)
 114	format('X-axis:   ',1A,F5.2,1A,1I3)
			goto 201
C__________'O'_________ Write file for figure production
 233	if (KEY.ne.79)	goto 234
	if (MOD10 .eq. 6)  then
	   call	ASKXGR(NTOUT,NWIND3,NAMET,MODEY,jj,j1,j2,OUTFIG,OUTNAME)
	   call	markloc("TIMOUT 9th call from IFKEY"//char(0))
	   call	TIMOUT
	   call	WRFIGS(jj,j1,j2,OUTFIG,OUTNAME,NWIND3,NTRUN,TTOUT,TOUT)
	elseif (MOD10 .le. 3)	then
	   call	ASKXGR(NROUT,NWIND1,NAMER,MODEY,jj,j1,j2,OUTFIG,OUTNAME)
C	   write(*,'(1A,1I6)')(outname(j),outfig(j),j=1,NROUT)
	   call	WRFIGS(jj,j1,j2,OUTFIG,OUTNAME,NWIND1,NTRUN,TTOUT,TOUT)
	endif
			goto 201
C__________'J'_________ Test field
 234	if (KEY.ne.74)	goto 235
	call	system("ipcs -s")	! Report active semaphore sets
	call	system("ipcs -m")	! Report active shared memory segments
			goto 201
C__________'K'_________ Test field
 235	if (KEY.ne.75)	goto 236
			goto 201
C__________'Z'_________ Test field
 236	if (KEY.ne.90)	goto 24
C	write(*,*)"Z"
C	write(*,*)"Runtime ="!,etime(tarray)
C     &		,"   user time =",tarray(1)
C     &		,"   system time =",tarray(2)
C	write(*,*)"Runtime =",dtime(tarray)
C     &		,"   user time =",tarray(1)
C     &		,"   system time =",tarray(2)
			goto 201
C Positions 10,20,30,40,50 are marked:
C               ---------1---------2---------3---------4---------5---
C	STRI = "ASTRA_name|     U-file_name  | shot_No | Signal |"
C	do	j = 1,NROUT
C	   write(ANAMES(j)(1:11),'(1A4,7X)')NAMER(j)
C	   write(ANAMES(j)(12:30),'(1A9,10X)')"text/text"
C	   write(ANAMES(j)(31:40),'(1A10)')'xxxxx     '
C	   write(ANAMES(j)(41:),'(1A4)')NAMER(j)
C	enddo
C               ---------1---------2---------3---------4---------5---
	STR = "Test window"//char(0)
	STRI = "    |    |    |    ||    "//
     ."|    |    |    "//char(0)
	jj = 1
	do	j = 1,NROUT/8
	   write(ANAMES(j)(1:4),'(1A4)')NAMER(jj)
	   write(ANAMES(j)(6:9),'(1A4)')NAMER(jj+1)
	   write(ANAMES(j)(11:14),'(1A4)')NAMER(jj+2)
	   write(ANAMES(j)(16:19),'(1A4)')NAMER(jj+3)
	   write(ANAMES(j)(22:25),'(1A4)')NAMER(jj+4)
	   write(ANAMES(j)(27:30),'(1A4)')NAMER(jj+5)
	   write(ANAMES(j)(32:35),'(1A4)')NAMER(jj+6)
	   write(ANAMES(j)(37:40),'(1A4)')NAMER(jj+7)
	   jj = jj+8
	enddo
	jj = ASKTAB(STR,STRI,ANAMES,NN80,NROUT/8,0,0)
	jj = 1
	do	j = 1,NROUT/8
	   NAMER(jj)   = ANAMES(j)(1:4)
	   NAMER(jj+1) = ANAMES(j)(6:9)
	   NAMER(jj+2) = ANAMES(j)(11:14)
	   NAMER(jj+3) = ANAMES(j)(16:19)
	   NAMER(jj+4) = ANAMES(j)(22:25)
	   NAMER(jj+5) = ANAMES(j)(27:30)
	   NAMER(jj+6) = ANAMES(j)(32:35)
	   NAMER(jj+7) = ANAMES(j)(37:40)
	   jj = jj+8
	enddo
	do	j = 1,NRW
	   BNAMES(j) = NAMER(j)
	enddo
	do	j = 1,NRW
	   NAMER(j) = BNAMES(j)
	enddo
			goto 201
C__________'D'
 24	if (KEY.ne.68)	goto 25
		NDTNAM=24+4*NSBR
		TIMEB=TIME
		MODEX = XOUT+.49
		call ASKLIS(NDTNAM,DELOUT,DTNAME,3)
		if (int(DELOUT(13)) .ne. NB1)	then
		   write(*,*)">>> NB1 re-definition ignored"
		endif
		DELOUT(13) = NB1
		NUF = DELOUT(14)
		NBND = DELOUT(19)
		XFLAG = DELOUT(20)
		j = XOUT+.49
		if (j.lt.0 .or. j.gt.3)	then
		   write(*,*)">>> Unknown X-axis. Redefinition ignored"
		   j = MODEX
		   XOUT = MODEX
		endif
		if (j .ne. MODEX)	call	xaxis(j)
		IF(TIME.GE.TIMEB)	goto	201
		TROUT=TIME
		TTOUT(LTOUT-1)=TIME
		TPOUT=TIME
		do 241 J=1,NSBR
 241		TEQ(J)=TIME
			goto 1
C__________'C'
 25	if (KEY.ne.67)	goto 26
	call ASKLIS(NCFNAM,CONSTF,CFNAME,2)
			goto 1
C__________'V'
 26	if (KEY.ne.86)	goto 27
	do	J=1,NPRNAM
	   DEVARO(J)=DEVAR(J)
	enddo
	INT4 = NPRNAM-48		! INT4 = NPRNAM - No. of ZRDs
	call ASKLIS(INT4,DEVAR,PRNAME,1)
	do	262	J=1,NPRNAM
	if (IFDFVX(J) .gt. 3)	DEVAR(J)=DEVARO(J)
	if (ABS(DEVAR(J)-DEVARO(J)).gt.1.d-6*ABS(DEVAR(J))) IFDFVX(J)=3
 262	continue
			goto 1
C__________'M'
 27	if (KEY.ne.77)	goto 30
	if (MOD10.eq.1)			call
     &	  ASXWIN(NROUT,NWIND1,NAMER,SCALER,OSHIFR,GRAL,GRAP,MOD10,MODEY)
	if (MOD10.eq.2.or.MOD10.eq.3)	call
     &	  ASXWIN(NROUT,NWIND2,NAMER,SCALER,OSHIFR,GRAL,GRAP,MOD10,MODEY)
	if (MOD10.eq.6)			call
     &	  ASTWIN(NTOUT,NWIND3,NAMET,SCALET,OSHIFT,MOD10,MODEY)
	if (MOD10.le.3 .or. MOD10.eq.6)		goto 201
	if (MOD10.ne.7)	goto 271
	call ASKLIS(4,TIM7,NAM7,7)
			goto 201
 271	if (MOD10.ne.4 .and. MOD10.ne.5)	goto 30
	INT4 = -MAX(4,IPOUT-1)
	call ASKLIS(INT4,PRMARK,NAMEP,8)
C	if (MOD10.eq.4.or.MOD10.eq.5)	call
C     &	  ASXWIN(NROUT,NWIND4,NAMER,SCALER,OSHIFR,GRAL,GRAP,MOD10,MODEY)
C----------------------------------------------------------------------|
			goto 201
C__________'H, ?'
 30	if (KEY.ne.72.and.KEY.ne.63)	goto	31
		write(*,*)
		write(*,*)"The following keys are operable in ",
     >				"this mode:"
		do 301 J=1,28
 301		call prntxt(HELP(J))
		goto 1
C__________'/', <ESC>
C 31	if (KEY.ne.47.and.KEY.ne.3.and.KEY.ne.27)	goto	32
C 31	if (KEY.ne.47.and.KEY.ne.27)			goto	32
 31	if (KEY.ne.47)		goto	32		! "/"
	goto	99
C__________'T'
 32	if (KEY.ne.84)	goto	42
	   call	markloc("TIMOUT 10th call from IFKEY"//char(0))
	   call	TIMOUT
	   call	TYPDSP(5,CHORDN,NTRUN,TTOUT,TOUT)
	   goto 1
C__________'space'
 42	if (KEY.ne.32)	goto	43
		KEY = 0
		if (TASK(1:3) .eq. 'DSP')	goto	76
		if (TASK(1:3) .eq. 'RUN')	then
			TASK='DSP '
			ix = 0
			iy = 0
			call pcurso
			goto	1
		endif
C__________'G', 'Q'
 43	if (KEY.ne.71 .and. KEY.ne.81)	goto	70
	   LRDN = length(RDNAME)
	   LWRN = length(EQNAME)
	   PSNAME='dat/'//RDNAME(1:LRDN)//'-'//EQNAME(1:LWRN)//char(0)
	   FLEN = KILLBL(PSNAME,NN80)
 	   if (KEY.eq.71)	INT4 = NN0
 	   if (KEY.eq.81)	INT4 = NN1
	   call PSOPEN (PSNAME, INT4, IRET)
	   FLEN = KILLBL(PSNAME,NN80)
	   if (IRET.eq.0)	then
	      if (KEY .eq. 71)	KPRI = 1
	      if (KEY .eq. 81)	KPRI = 2
	      goto	201
	   endif
	   if (IRET.eq.1)	then
	      STRI='>>>  Can not open file: '//PSNAME(1:FLEN)
	      call colovm(WarningColor)
	      call textvm(NN0,LRJJ,STRI,24+FLEN)
	   endif
	   KEY = 0
	   goto	70
 49	continue
C----------------------------------------------------------------------|
	if (KIBM .ne. 2)	goto	60
C-------- <Alt>'M' or  <Alt>'m'
	if (KEY.ne.77 .and. KEY.ne.109)	goto	50
	if (PNMNAME(1:1) .eq. "*")	then	! odd pressing <AltM>
	   LRDN = length(RDNAME)
	   LWRN = length(EQNAME)
	   PNMNAME='tmp/'//RDNAME(1:LRDN)//'-'//EQNAME(1:LWRN)//char(0)
	   MOVIE = 1
	else
	   MOVIE = 0
	endif
	goto	70
 50	continue
	if  (KEY.eq.47)	goto	97 ! <Alt>+/
	if (KIBM.eq.2 .and. (KEY.ge.32 .and. KEY.le.126) )	then
	   write(*,*)'  "<Alt>+<',char(KEY),'>"  pressed'
	endif
 60	continue
	if (KIBM.eq.1 .and. (KEY.ge.32 .and. KEY.le.126) )	then
	   write(*,*)'  "<Ctrl>+<',char(KEY),'>"  pressed'
	endif
	jj = 0
	if (KEY .gt. 90)	KEY = KEY-32
	KEY = KEY-64
C	write(*,*)KIBM,KEY,char(KEY)
	do	j=1,NSBR
C          write(*,*)j,DTEQ(4,J)
	   if (ABS(KEY-DTEQ(4,j)) .lt. 0.1)	jj = 1
	enddo
	if (TASK(1:3).eq.'DSP' .and. jj.eq.1)	goto	76
	if (TASK(1:3).eq.'DSP')			goto	1
	if ( jj .eq. 1)				goto	77
 70	continue
	if (KEY.eq.27)		then
	   write(*,'(/2A)')'Use key "/" for exit',char(7)	! Beep
	elseif (KEY.ne.0 .and. KIBM.eq.0)	then
	   write(*,*)'Unrecognized key: "',char(KEY),'"',KEY,char(7)
	endif
	KEY = 0
	goto	1
 76	IFLAG = 1				! for DSP mode only
 77	IFKEY = 0
	IFKEY_= 0
C	call	ADDTIME(CPTGRA)
C	write(*,*)"Exit"
	return
 91	write(*,*)'>>> No. of radial points in a U-file is too large:'
	write(*,*)'           Decrease NUF or increase NB1'
	goto	1
 97	continue
	if (PNMNAME(1:1) .ne. "*")	then	! even pressing <AltM>
	   call	makemovie(0,PNMNAME)
	endif
	write(STR,'(1I2.2,1H:,1I2.2,1H:,1I2.2)')IHOUR,IMINUT,ISEC
	call	TIMDAT(IYEAR,IMONTH,IDAY,IHOUR,IMINUT,ISEC)
	if (TASK(1:3).ne.'BGD' .and. TASK(4:4).ne.'B')	call	endvm
	if (IFKL .eq. 257)	then
	   write(6,'(A)')' >>> ASTRA error >>>'
	   write(6,'(A,F11.6,A)')"    Floating point exception at  t =",
     >			TIME,' sec'
	   call	wrtime(6,'    Run time',12,swatch(Y),-1.d0)
	else
	   TIMEB = swatch(Y)
	   j1 = TIMEB
	   j2 = j1/3600
	   jj = (j1-3600*j2)/60
	   j1 = timeb-60*jj-3600*j2
	   write(6,'(A,I4.2,2(A1,I2.2))')
     >		'>>> ASTRA normal exit >>>  Run time',j2,':',jj,':',j1
	endif
	call	CPUSE(6)
 99	continue
	call	a_stop
C	call	freeshm
C	call	system("rm tmp/*.ppm >& /dev/null"//char(0))
C	call	system ("if ( -e ${HOME}/bin/cln ) ${HOME}/bin/cln")
C	write(*,*)"Exit"
C	stop
	end
C======================================================================|
	subroutine	SETTSC(JTIMS)
C----------------------------------------------------------------------|
C Define time scales (former piece of READAT
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'		! TAU,TAUMIN,TAUMAX,D?OUT,TSCALE
	include	'for/outcmn.inc'	! EQNAME is used
	integer	JTIMS,j,length
	double	precision YS(10),YSMAX
	logical EXILOG
C Read file equ/MODEL.log
	j = length(EQNAME)
C	write(*,*)'"',EQNAME(1:j),'"',j
	inquire(file='equ/log/'//EQNAME(1:j),exist=EXILOG)
C	if (EXILOG)write(*,*)'"equ/log/'//EQNAME(1:j),'"'
	if ( EXILOG )	then
	   TAU = TAUMIN
	else
	   TAUMAX = .01*VOLUME
	   YS(1) = 0.00000010
	   YS(2) = 0.00000015
	   YS(3) = 0.00000020
	   YS(4) = 0.00000025
	   YS(5) = 0.00000030
	   YS(6) = 0.00000040
	   YS(7) = 0.00000050
	   YS(8) = 0.00000075
 41	   do	J=1,8
	      YSMAX = YS(J)
	      if (1.1*TAUMAX .le. YS(J))	goto 42
	   enddo
	   do	J=1,8
	      YS(J) = 10.*YS(J)
	   enddo
	   goto 41
 42	   TSCALE = JTIMS*YSMAX*.1
	   TAUMAX = YSMAX
 43	   continue
	   do	J=1,8
C 115=(right_label_position)/IDT=575/5
	      if (TSCALE*115/JTIMS .le. YS(J))	goto 44
	   enddo
	   do	J = 1,8
	      YS(J) = 10.*YS(J)
	   enddo
	   goto 43
 44	   TSCALE = YS(J)

	   TAUMAX = .01*VOLUME
	   YS(1) = 0.00000010
	   YS(2) = 0.00000015
	   YS(3) = 0.00000020
	   YS(4) = 0.00000025
	   YS(5) = 0.00000050
	   YS(6) = 0.00000075
 45	   do	J=1,6
	      YSMAX = YS(J)
	      if (1.1*TAUMAX .le. YS(J))	goto 46
	   enddo
	   do	J=1,6
	      YS(J) = 10.*YS(J)
	   enddo
	   goto 45
 46	   TAUMAX = YSMAX
	   TAUMIN = .001*TAUMAX
	   TAU    = TAUMIN
	   DTOUT  = .01*TAUMAX
	   DROUT  = .02*TAUMAX
	   DPOUT  = .2*TAUMAX
	endif
C	write(*,*)TAU,TAUMIN,TAUMAX,TSCALE,DROUT,DTOUT,DPOUT,YSMAX
	end
C======================================================================|
	subroutine	SMODE5(MARK,PT,PRMARK,NAMEP,ITIMES)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'		! NAB,AB are used
	include	'for/outcmn.inc'
	integer ICVMX
	parameter (ICVMX=32)		! <= 32 curves in a screen
	integer	MARK,ITIMES,PT(*),IX(2*NRD),IY(2*NRD),JTIM,
     1		j,jj,ierr,int2,jab,i,is,NP,NP1,jxout,jx,jx1,jy,jy1,jnl,
     2		jc,JDSP,IYMN,IYMX,STYL,jpos,FSHIFT,jk(5),SKIPM,lonlen
	double precision	PRMARK(*),
     1	 	SC(NRW),TEMPR,SCL,DOWN,YWA(NRD),YWB(NRD),YS,YL,YR,YX
     2		,YROUT,YA,YQ1,YQ2,YXR,YXL
	character*6	NAMEP(*),STRI*80,CHAR4*4,ST*9
	integer*2	INTY(NRD)
	save	YQ1,YQ2
	data 	FSHIFT/10/
C----------------------------------------------------------------------|
	IYMN=YSCMAX-IYM
	IYMX=YSCMAX-IY0
	JY = 3500
C----------------------------------------------------------------------|
C Mode 4 & 5:
	call	add2loc("Drawing mode 4"//char(0))
C Plots against (a/AB) (scale does not change with time)
C Scale is determined by the current radial distribution
	call	SCAL(NROUT,SC,SCALER,ROUT(4,1),NAB-3,NRD)
	JTIM=0
	call 	OPENRD(3,RSNAME,1,IERR)
	if(IERR.gt.0)	pause ' >>> Profiles file error'
	j = 0
	j = SKIPM(3,j)			! Returned value is not used
	read(3,ERR=38,END=39)CHAR4
	if (NXOUT.gt.0 .and. NGR.gt.0) read(3,ERR=38,END=39)INT2
 31	read(3,ERR=38,END=39)INT2
	if (INT2.ne.0)	read(3,ERR=38)TEMPR
	read(3,END=39)TEMPR
	read(3,END=39)TEMPR
	if (JTIM .gt. ITIMES)	goto	381
	JTIM=JTIM+1
	read(3,END=39)JAB,JAB,(INT2,I=1,10),(TEMPR,I=1,10)
C	write(*,*)"AMETR:",JAB
	read(3)SCL,DOWN,(INTY(j),j=1,JAB)
	do	j=1,JAB
	   YWA(j) = DOWN+(32768+INTY(j))*SCL/65535.
	enddo
C	write(*,*)(A(j),j=1,jab)
C	write(*,*)"SHIF:"
	read(3)SCL,DOWN,(INTY(j),j=1,JAB)
	do	j=1,JAB
	   YWB(j) = DOWN+(32768+INTY(j))*SCL/65535.
	enddo
C	write(*,*)"ELON, TRIA, A, B, FP:"
	do jj=1,5
	   read(3)SCL,DOWN,(INTY(j),j=1,JAB)
	enddo
	
	do	37	jj=1,NROUT
	read(3,ERR=38)SCL,DOWN,(INTY(j),j=1,JAB)
	if(NAMER(jj).eq.'    ' .or. PRMARK(JTIM).lt.0) goto 37
	NP=NWIND4(jj)-2*NSCR(MOD10)
	JX=(NP-1)*DXM
	if (NP.ne.1 .and. NP.ne.2)	goto	37
	DOWN = DOWN+OSHIFR(jj)
	jnl=DYLET+FSHIFT
	do IS=1,5
	   jk(IS)=0
	enddo

	call FMTF4(ST(1:4),SC(jj))
	ST(5:9)=NAMER(jj)
	call colovm(1)
	call textvm((NP-1)*DXM,jnl,ST,9)
	YS = OSHIFR(jj)
	if (YS)	341,344,342
 341	YS = -YS
	call FMTF4(CHAR4,YS)
	ST = '-'//CHAR4
	goto	343
 342	call FMTF4(CHAR4,YS)
	ST = '+'//CHAR4
 343	call textvm((NP-1)*DXM+4*DXLET,jnl+DYLET+2,ST,5)
 344	continue

	jxout = NP1
	if (MOD10 .eq. 4)	then
	   NP1 = JAB
	   YL = GRAL(jj)/YWA(NP1)
	   YR = GRAP(jj)/YWA(NP1)
	   if (YL.ge.YR)	YL = 0.d0
	   if (YL.lt.YR .and. (YL.gt.0.001 .or. YR.lt.0.999))	then
	      PT(1) = NP+320*YL
	      PT(2) = 350-IYMN
	      PT(3) = NP+320*YR
	      call	colovm(2)
	      call	drawvm(0,PT(1),PT(2),PT(3),PT(2))
	      call	drawvm(0,PT(1),PT(2)+1,PT(3),PT(2)+1)
	      call	drawvm(0,PT(1),PT(2)+2,PT(3),PT(2)+2)
	   endif
	   jxout = 0
	   do	j=1,NP1
	      YX = YWA(j)/YWA(NP1)
	      if (YX.lt.YL .or. YX.gt.YR)	goto	33
	      jxout = jxout+1
	      if (jxout.eq.1 .and. j.gt.1)	then ! left edge interpolation
		 IX(jxout) = JX
		 YQ1 = DOWN+(32768+INTY(j-1))*SCL/65535.
		 YQ2 = DOWN+(32768+INTY(j))*SCL/65535.
		 YA = YQ2+(YQ1-YQ2)*(YL-YX)/(YA-YX)
		 YROUT = min(max(YA/SC(jj),-7.d0),7.d0)
		 JDSP  = 10*(DYM*YROUT+IYMN)
		 if (MODEY .eq. -1)	JDSP = JDSP+10*DYM
		 IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
		 jxout = jxout+1
	      endif
	      IX(jxout) = JX+320*(YX-YL)/(YR-YL)
	      YROUT=(DOWN+(32768+INTY(j))*SCL/65535.)/SC(jj)
	      YROUT=min(max(YROUT,-7.d0),7.d0)
	      JDSP	=10*(DYM*YROUT+IYMN)
	      if (MODEY .eq. -1)	JDSP = JDSP+10*DYM
	      IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
 33	      continue
	      if (YA.le.YR .and. YX.gt.YR)	then ! right edge interpolation
		 jxout = jxout+1
		 IX(jxout) = JX+320
		 YQ1 = DOWN+(32768+INTY(j-1))*SCL/65535.
		 YQ2 = DOWN+(32768+INTY(j))*SCL/65535.
		 YA = YQ2+(YQ1-YQ2)*(YR-YX)/(YA-YX)
		 YROUT = min(max(YA/SC(jj),-7.d0),7.d0)
		 JDSP  = 10*(DYM*YROUT+IYMN)
		 if (MODEY .eq. -1)	JDSP = JDSP+10*DYM
		 IY(jxout) = JY-min(max(JDSP,10*IYMN),10*IYMX)
	      endif
	      YA = YX
	   enddo
	   NP1 = jxout
	else
C Mode 5
C Take AMETR(x,t) and SHIF(x,t) from "profile.dat"
	   NP1 = 2*JAB
	   do	j=1,JAB
C	      YROUT=(DOWN+(128+INTY(j))*SCL/255.)/SC(jj)
	      YROUT=(DOWN+(32768+INTY(j))*SCL/65535.)/SC(jj)
	      YROUT=min(max(YROUT,-7.d0),7.d0)
	      JDSP	=10*(DYM*YROUT+IYMN)
	      if (MODEY .eq. -1)	JDSP = JDSP+10*DYM
	      IY(j) = JY-min(max(JDSP,10*IYMN),10*IYMX)
	      IY(JAB+j) = IY(j)
	   enddo
	   do j = 1,JAB
		YXR = (YWB(j)+YWA(j))/AB
		YXL = (YWB(j)-YWA(j))/AB
		IX(JAB+j) = JX+160*(1.+min(1.d0,YXR))
		IX(JAB+1-j) = JX+160*(1.+max(-1.d0,YXL))
		IY(JAB+1-j) = IY(JAB+j)
	   enddo
	   IX(1) = min(IX(1),JX)
	   IX(NP1) = max(IX(NP1),JX+320)
	endif
	jc = 31		! non-marked profiles (shadow color)
	do	IS=1,5
	   if (PRMARK(JTIM) .eq. IS)	then
	      if (jk(IS) .ne. 0)	goto 36
	      jk(IS) = 1
	      jc = IS+1
	   endif
	enddo
 36	continue
	STYL = (jc-1)*MARK
	if (PRMARK(JTIM).eq.0)	STYL=0
	if (KPRI.ge.1 .and. KPRI.le.2)	then
	    write(STRI,'(1A6,1A6,1A1)')'Plot "',NAMEP(JTIM),'"'
	    j = lonlen(STRI)
C	    write(*,'(1A6,1A6,1A1,I3)')'Plot "',NAMEP(JTIM),'"',j
	    call	pscom(STRI,j)
	endif
	call	PLOTCR(NP1,-1,IX,jx1,IY,jy1,jc,STYL,PT)
	if (PRMARK(JTIM) .eq. 0)	goto	37
	JPOS=3*DXLET+PRMARK(JTIM)*DXM/6.5+DXM*(NP-1)
	call CMARKP(jnl,JPOS,NAMEP(JTIM),STYL)
 37	continue
	go to 31
 38	pause	'Read file PROFIL.DAT error'
		return
 381	write(*,*)">>> Too many time slices. Data skipped"
	write(*,*)"    Use Astra post-viewer or reduce DPOUT"
 39	close(3)
	return
	end
C======================================================================|
	subroutine	wrtime(nch,string,len,tim,time)
	implicit none
	integer	 nch,jh,jm,js,len
	character*(*)	string
	double precision	 tim,time,t1
	js = tim
	jh = js/3600
	jm = (js-3600*jh)/60
	t1 = tim-60*jm-3600*jh
	js = t1
	if (time .lt. 0.)	goto	1
	write(nch,'(2A,I4.2,2(A1,I2.2),F8.1,A1)')string(1:len),
     >		char(9),jh,':',jm,':',js,100.*tim/time,'%'
	return
 1	continue
	write(nch,'(2A,I4.2,2(A1,I2.2))')string(1:len),
     >		char(9),jh,':',jm,':',js
	end
C===================================================10.03.95 G.P.=======
	subroutine	CPUSE(nch)
	implicit none
	character*32	STRI
	double precision	Y
	integer	nch,j,j1,j2,length
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	CPTOT = CPT+CPTEQL
	do j=1,NSBR
	   CPTOT = CPTOT+CPTSBR(j)
	enddo
	write(nch,'(A,I8)')    "    Total time steps  ",NSTEPS
	if (NSTEPS .eq. 0)	return
c	write(nch,*)TIME,TSTART,CPTOT,CPT,CPTEQL
c	write(nch,*)CPTSBR
	Y = (TIME-TSTART)/NSTEPS
	if (Y .lt. 1.d-1)	then
	   Y = 1.d3*Y
	   write(nch,'(A,F6.3,A)')"    Average time step   ",Y," msec"
	else
	   write(nch,'(A,F6.3,A)')"    Average time step   ",Y," sec"
	endif
	write(nch,'(A,F6.3,A)')"    CPU per time step   "
     >			,CPTOT/NSTEPS," sec"
	Y = CPTOT/(TIME-TSTART)
	if (Y .lt. 60.)	then
	   write(nch,'(A,F6.3,A)')"    CPU per 1 sec       ",Y," sec"
	else
	   call	wrtime(nch,'    CPU per 1 sec ',18,Y,-1.d0)
	endif
	call	wrtime(nch,'    Total CPU time',18,CPTOT,CPTOT)
	call	wrtime(nch,'    Transport core',18,CPT,CPTOT)
	call	wrtime(nch,'    Equilibrium   ',18,CPTEQL,CPTOT)
	j2 = 1
	do j1=1,NSBR
	   j = min(6,length(DTNAME(24+4*j1)))
	   if (j1 .eq. IFSBX(j2))	then
	      call	wrtime(nch,'    Xroutine   "'//
     >			DTNAME(24+4*j1)(1:j)//'"',17+j,CPTSBR(j1),CPTOT)
	      j2 = j2+1
	   else
	      call	wrtime(nch,'    Subroutine "'//
     >			DTNAME(24+4*j1)(1:j)//'"',17+j,CPTSBR(j1),CPTOT)
	   endif
	enddo
	write(nch,*)
	end
C======================================================================|
	subroutine ADDMOD(NCHW,NCHM,NCHL)
C Write model & model.log records in the header of a post-view file
C 	 Both are preceded by one line 32*"^" 
C             and two such lines are written after model.log
C Unit NCHW (profile data file) must be open
C Unit NCHM (model) must be open
C Unit NCHL (model.log) must be open if nonzero
C Note:	  1) NCHL =/= NCHM;	2) Empty lines are skipped.
	implicit none
	integer NCHW,NCHM,NCHL,NCH,NCHI,LSTRI,j,LSTR,leng0,killbl,ierr
	character STRI*133,STR*133,INCNAM*80
	write(NCHW)"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
	NCH = NCHM
	NCHI = 0
 1	continue
	write(STR,'(133A1)')(' ',j=1,132),char(0)
	read(NCH,FMT='(1A132)',ERR=13,END=11)STR
	LSTR = leng0(STR)
C '#' - pre-processor command
	if (STR(1:8) .ne. "#include")	goto	10
C	write(*,'(3A,I5)')'"',STR(1:LSTR),'"',LSTR
	LSTRI = index(STR(1:),';')-1
	if (LSTRI .eq. -1)	LSTRI=LSTR
	INCNAM = STR(9:LSTRI)//char(0)
	j = killbl(INCNAM,LSTRI-8)
C	write(*,'(3A,3I5)')'"',INCNAM(1:j),'"',j,LSTRI,LSTR
	NCHI = NCHM+1
	call	OPENRD(NCHI,"equ/"//INCNAM(1:j)//char(0),0,IERR)
	if (IERR.gt.1)	goto	97
 2	read(NCHI,FMT='(1A132)',ERR=98,END=3)STRI
	if (STRI(1:8) .eq. "#include")	goto	96
	ierr = leng0(STRI)
	write(NCHW)char(ierr),STRI(1:ierr)
	goto	2
 3	continue
	close(NCHI)
 4	continue
	NCHI = 0
	if (LSTRI .lt. LSTR)
     >		write(NCHW)char(LSTR-LSTRI),STR(LSTRI+1:LSTR)
	goto	1
 10	write(NCHW)char(LSTR),STR(1:LSTR)
	goto	1
 11	continue
	write(NCHW)char(32),"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
	if (NCHL.eq.0 .or. NCHL.eq.NCH) goto	12
	NCH = NCHL
	goto	1
 12	write(NCHW)char(32),"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
	return
 13	write(*,*)'>>> Error reading model file'
	return
 96	write(*,*)'>>> Error in include file: "',INCNAM(1:j),'"'
	write(*,*)'    Recurrent "#include" is not allowed. '
     >			,'Saving model failed'
	write(*,*)
	goto	4
 97	write(*,*)'>>> Cannot open include file ',
     >			'"equ/'//INCNAM(1:j),'". Saving model failed'
	write(*,*)
	goto	4
 98	write(*,*)'>>> Cannot read include file ',
     >			'"equ/'//INCNAM(1:j),'". Saving model failed'
	write(*,*)
	goto	3
	end
C=======================================================================
	integer	function leng0(STRI)
C Function returns a position of the last non-blanck symbol in STRI
	implicit none
	character*1 STRI(132),HT
	integer  j
	j = 9
	HT = char(j)
	do	j=132,1,-1
	   leng0 = j
	   if (STRI(j).ne.' ' .and. STRI(j).ne.HT)	return
	enddo
	end
C======================================================================|
	subroutine	STUFF(NUNIT,ARRAY,NA,BITARR)
	implicit none
	integer	NUNIT,NA,JJ
	double precision	ARRAY(NA),YUP,YDN,YSC
	integer*2	BITARR(NA)
	YUP =-1.E37
	YDN = 1.E37
	do	JJ=1,NA
	   YUP = max(ARRAY(JJ),YUP)
	   YDN = min(ARRAY(JJ),YDN)
	enddo
	YSC = YUP-YDN
	if (YSC .eq. 0.)	goto	10
	do	JJ=1,NA
	   YUP = (ARRAY(JJ)-YDN)/YSC
	   BITARR(JJ)=65535*YUP-32768
c	   BITARR(JJ)=255*YUP-128
	enddo
 10	write(NUNIT)YSC,YDN,BITARR
	end
C======================================================================|
C LINEAV [10#19/m#3]: Horizontal chord average density (r) [m]
C	Integral {0,r} ( NE ) dl / a
	double precision function	LINEAV(jj)
	implicit none
	include  'for/parameter.inc'
	include  'for/const.inc'
	include  'for/status.inc'
	integer	j,jj
	LINEAV=2.*AMETR(1)*NE(1)
	do	j=2,NA
	   LINEAV=LINEAV+(AMETR(j)-AMETR(j-1))*(NE(j)+NE(j-1))
	enddo
	LINEAV=0.5*(LINEAV+(ABC-AMETR(NA))*(NE(NA1)+NE(NA)))/ABC
	return
	END
C======================================================================|
	subroutine PUTJIT(jit)
	implicit none
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	character	STRI*80
	integer	jit,lonlen,FSHIFT
	data	FSHIFT/2/
	call	colovm(2)			! Iteration # in blue
	write(STRI,'(a,i3,2x)')"  Iteration  #",jit
	call	textvm(64*DXLET,FSHIFT,STRI,lonlen(STRI))
	end
C======================================================================|
	subroutine ADDTIME(anyTime)
	implicit none
	double precision	anyTime, cpuTime, swatch
	cpuTime = swatch(anyTime)
	end
C======================================================================|
	subroutine	REPORT
C----------------------------------------------------------------------|
C When in background mode this subroutine allows 
C----------------------------------------------------------------------|
	implicit 	none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	character*80	string
	double	precision	swatch,Y
	integer	j,j1,jb,jl,length
C	if (TASK(1:3) .ne. 'BGD')	return
C	write(*,*)'AWD = "',AWD(1:length(AWD)),'"'
	jl = length(RSNAME)
C	write(*,*)'RSNAME = "',RSNAME(1:jl),'"'
	j1 = index(RSNAME(1:jl),'profil.dat')
	if (j1 .eq. 0)	then
	   jb = 1
 1	   j = index(RSNAME(jb:jl),'/')
	   if (j .ne. 0)	then
	      jb = jb+j
	      goto	1
	   endif
	   string = './tmp/status.'//RSNAME(jb:jl)//char(0)
	else
	   string = './tmp/status'//char(0)
	endif
C	write(*,*)'Status file name = "',string(1:length(string)),'"'
	open(7,file=string(1:length(string)),status="UNKNOWN",err=2)
	if (TIME .lt. 1.d1)	then
	   write(7,'(A/A,F9.6)')RUNID,'     Time =',TIME
	else
	   write(7,'(A/A,F8.3)')RUNID,'     Time =',TIME
	endif
	call	wrtime(7,'  >>> Astra run time  ',22,swatch(Y),-1.d0)
	call	CPUSE(7)
	close(7)
	return
 2	write(*,*)'>>> ERROR opening status file ',
     &		'"',string(1:length(string)),'". Call ignored.'
	end
C======================================================================|
	subroutine	A_POSTMORTEM
C----------------------------------------------------------------------|
C Print time table before exit 
C----------------------------------------------------------------------|
	implicit 	none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	character*80	string
	double	precision	swatch,Y
	write(6,'(A/A,F8.3)')RUNID,'     Time =',TIME
	call	wrtime(6,'  >>> Astra run time  ',22,swatch(Y),-1.d0)
	call	CPUSE(6)
	return
	end
C======================================================================|
