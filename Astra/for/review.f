C======================================================================|
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
C		NCURVS	maximal # of radial profiles stored
C		NRMAX	maximal # of time slices to be read from a file
C		IROUT	current # of the radial picture
C		ITOUT	current # of the time record
C		NROUT	# of radial channels in the main run (outcmn.inc)
C		NROUT+7	total # of radial profiles stored for each time slice
C		NNROUT	total # of radial records
C		NRTIMS	=min(NRMAX,NCURVS/(NROUT+7)) max # of radial records
C		NTMAX	# of times stored for time curves (see "parameter.inc")
C		NNTOUT	total # of time records
	integer j,IROUT,ITOUT,NNROUT,NNTOUT,INA1,INAB,IRTYPE,IFKL,IFKEY
	double precision AVNE, TROUT, YX
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     >			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
	double precision TTOUT(NTMAX),TOUT(NTMAX,NRW)
	equivalence (WORK1D(1),TTOUT(1)),(WORK1D(NTMAX+1),TOUT(1,1))
	integer	iargc,KASCII,NEVENT,KIBM,waitevent,length,	! PTN(2),
     >		ix,iy,jl,jr,JM,JX,JX0,JXM,JNX,JNY,LOCTIM
	character STRING*132
C----------------------------------------------------------------------|
	do j=0,iargc()
	   call	getarg(j,STRING)
	   jl = length(STRING)
C	   if (j .eq. 0) write(*,'(3a,1i6)')
C     &		'Task: "',STRING(1:jl)	!,'",   PID =',getpid()
C	   if (j .gt. 0) write(*,'(a,i2,3a)')
C     &		'      Argument #',j,'  "', STRING(1:jl),'"'
	   if (j .eq. 1) write(RSNAME,'(a)')STRING(1:jl)
C	   write(*,*)'RSNAME = "',RSNAME(1:jl),'"'
	enddo
C----------------------------------------------------------------------|
	do	j=1,NRW
	   NWIND1(j) = j
	   NWIND2(j) = j
	   NWIND3(j) = j
	   NWIND7(j) = j
	enddo
C	data	NWIND4/
C     &		1,5,9,13,3,7,11,15,2,6,10,14,4,8,12,16,
C     &		17,0,0,0,19,0,0,0,18,0,0,0,20,67*0./
C     7		1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,44*0./
	do	j=1,4
	   NWIND4(j) = 4*j-3
	   NWIND4(j+4) = NWIND4(j)+2
	   NWIND4(j+8) = NWIND4(j)+1
	   NWIND4(j+12) = NWIND4(j)+3
	enddo
	do	j=17,NRW
	   NWIND4(j) = NWIND4(j-16)+16
	enddo
	call	RDRES(0)
	call	RADOUT
	IFKL	=82
1	continue
	if (IFKL.eq.32 .or. IFKL.eq.82)  call	ERASXY(IX,IY)
C    IFKL = 32 (change current time), 82 (load & redraw)
C or IFKL = upper_case(ASCII_code)
	j = IFKEY(IFKL)
	NEVENT	=0
	KASCII	=0
2	NEVENT	=waitevent(KASCII,ix,iy)

C	write(*,*)char(KASCII),KASCII,NEVENT
C	write(*,*)'Review: NEVENT =',NEVENT,'   KASCII =', KASCII
	call PUTXY(ix,iy,NTMAX,TTOUT,TOUT)
C	if (NEVENT.eq.65001 .and. KASCII.eq.0)	call testmarkers(ix,iy)
C	PTN(1) = ix
C	PTN(2) = iy
C	if (NEVENT.eq.65001 .and. KASCII.eq.0)	then
C	   call	NMARK(PTN,1)
C	   PTN(1) = PTN(1)+10
C	   call	NMARK(PTN,2)
C	   PTN(1) = PTN(1)+10
C	   call	NMARK(PTN,3)
C	   PTN(1) = PTN(1)+10
C	   call	NMARK(PTN,4)
C	   PTN(1) = PTN(1)+10
C	   call	NMARK(PTN,5)
C	   PTN(1) = PTN(1)+10
C	   call	NMARK(PTN,6)
C	   PTN(1) = PTN(1)+10
C	   call	NMARK(PTN,7)
C	endif

C Button press-right_movement-release event
	if (NEVENT.gt.65000 .and. KASCII.eq.0 .and. MOD10.eq.6)	then
	   call	SETFRM(JM,JX0,JXM,JNX,JNY)
	   JX = IX-10
C	   if (JX0.gt.JX)  write(*,*)"Beyond left"
C	   if (JX.gt.JXM)  write(*,*)"Beyond right"
	   YX  = TINIT+1.029565*(JX-JX0)*abs(TSCALE)/(JXM-JX0)
	   if (NEVENT.eq.65001)	then
	      jl = LOCTIM(YX,TTOUT,LTOUT)
C	      write(*,'(/1A,2I7,1A,1F10.3)')"Button  pressed",
C     >		ix,jl,"   t =",TTOUT(jl)
	   endif
	   if (NEVENT.eq.65002)	then
	      jr = LOCTIM(YX,TTOUT,LTOUT)
C	      write(*,'( 1A,2I7,1A,1F10.3)')"Button released",
C     >		ix,jr,"   t =",TTOUT(jr)
C	      write(*,*)'Press "Z" to restore original scale'
	      if (jr-jl .lt. 10)	then
		 IFKL = 32
		 goto	3
	      endif
	      call	SETTSC(TTOUT(jl),TTOUT(jr))
	      IFKL = 82
	      goto	1
	   endif
	endif
 3	continue
	if (NEVENT.eq.0     .and. KASCII.eq.0)	goto	2
	if (KASCII.lt.127 .and. NEVENT.eq.1)	then
C	   write(*,*)'  "<Ctrl>+<',char(KASCII),'>"  pressed'
	   goto 2
	endif
	if (KASCII.lt.127 .and. NEVENT.eq.2)	then
C	   write(*,*)'  "<Alt>+<',char(KASCII),'>"  pressed'
	   goto 2
	endif
	if (NEVENT.eq.0   .and. KASCII.ge.65000)	then
	   KIBM = KASCII-65000
	else
	   KIBM	=NEVENT-65000
	endif
	if (KIBM.eq.0 .and. KASCII.eq.32)	KIBM = 363
	IFKL	=KASCII
	if (KASCII .eq. 32)	IFKL=KASCII-32
	if (KASCII .GE. 97)	IFKL=KASCII-32
C <PgUp>
	if(KIBM.eq.498 .or. KIBM.eq.365)	then
	   IFKL	=82
	   call	RDRES(-1)		
C	   TINIT = TSTART
	   goto	1
	endif
C <PgDn>
	if(KIBM.eq.504 .or. KIBM.eq.366)	then
	   call	RDRES(1)		
	   IFKL	=82
	   goto	1
	endif
C <-
	if(KIBM.eq.361)	then
	   if(MOD10.ge.6.and.MOD10.le.7)	then
	      ITOUT	=max(1,ITOUT-1)
	      TIME	=TTOUT(ITOUT)
	   else
	      IROUT	=max(1,IROUT-1)
	      TIME	=TROUT(IROUT)
	   endif
	   IFKL	=32
	endif
C ->
	if(KIBM.eq.363)	then
	   if(MOD10.ge.6.and.MOD10.le.7)	then
	      ITOUT	=min(NNTOUT,ITOUT+1)
	      TIME	=TTOUT(ITOUT)
	   else
	      IROUT	=min(NNROUT,IROUT+1)
	      TIME	=TROUT(IROUT)
	   endif
	   IFKL	=32
	endif
C <arrow down>
	if(KIBM.eq.364)	then
	   if(MOD10.ge.6.and.MOD10.le.7)	then
	      ITOUT	=min(NNTOUT,ITOUT+10)
	      TIME	=TTOUT(ITOUT)
	   else
	      IROUT	=min(NNROUT,IROUT+10)
	      TIME	=TROUT(IROUT)
	   endif
	   IFKL	=32
	endif
C <arrow up>
	if(KIBM.eq.362)	then
	   if(MOD10.ge.6.and.MOD10.le.7)	then
	      ITOUT	=max(1,ITOUT-10)
	      TIME	=TTOUT(ITOUT)
	   else
	      IROUT	=max(1,IROUT-10)
	      TIME	=TROUT(IROUT)
	   endif
	   IFKL	=32
	endif
C <Home>
	if(KIBM.eq.496 .or. KIBM.eq.360)	then
	   if(MOD10.ge.6 .and. MOD10.le.7)	then
	      ITOUT	=1
	      TIME	=TTOUT(ITOUT)
	   else
	      IROUT	=1
	      TIME	=TROUT(IROUT)
	   endif
	   IFKL	=32
	endif
C <End>
	if(KIBM.eq.502 .or. KIBM.eq.367)	then
	   if(MOD10.ge.6 .and. MOD10.le.7)	then
	      ITOUT	=NNTOUT
	      TIME	=TTOUT(ITOUT)
	   else
	      IROUT	=NNROUT
	      TIME	=TROUT(IROUT)
	   endif
	   IFKL	=32
	endif
	goto	1
	end
C======================================================================|
	subroutine	RDRES(NLOAD)
C----------------------------------------------------------------------|
C NLOAD - 0 - 1st call; >0 - load up; <0 - load down;
C NROUT - # of radial channels in use
C NTOUT - # of time channels in use
C KROUT - total # of radial records (profile-time slices) in the review file
C KTOUT - total # of time records in the review file
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer	KROUT,KTOUT,MTOUT,NRTIMS,NLOAD,LOCTIM
	integer	LTOUT1,LTOUTO,IRTYPE,SKIPM,index,nsymb,length
	integer	J,J0,J1,J4,jt,jr,jst,jsr,JJ,JJT,JJJ,IDENT,JTOUT
	integer	IERR,JINC,INAB,INA1,NNTOUT,NNROUT,ITOUT,IROUT,NCHL
	double precision TTOUT(NTMAX),TOUT(NTMAX,NRW)
	equivalence (WORK1D(1),TTOUT(1)),(WORK1D(NTMAX+1),TOUT(1,1))
	double precision YY,YT,YTIME,SCM8(20),RSCAL,RDOWN,QBND,	! tst,ten
     1		AVNE,ROCT,DDVAR,ABCT,CONS,TROUT
	character	RUNAST*79,VNAM*6,STR4*4,STRI*132,VVERS*5
	character*40	FILENA,FILNAD
	logical		EXIRES
	integer*2	RROUT(NRD*NCURVS)
C--------------------------------------	 5*NCURVS ~(NRD/NB1)*NCURVS
	common	/AVIEW_RDAT/ RROUT,RSCAL(5*NCURVS),RDOWN(5*NCURVS)
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     >			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
	common	/AVIEW_CV/ CONS(NCONST,NRMAX),DDVAR(NCONST,NRMAX)
     .			,ABCT(NRMAX),ROCT(NRMAX),QBND(NRMAX)
	save	LTOUTO,NRTIMS,KROUT,KTOUT,jst,jsr
	data	SCM8/.1,.15,.2,.25,.3,.4,.5,.6,.8,1.,1.25,1.5,
     &			2.,2.5,3.,4.,5.,6.,8.,10./	LTOUTO/1/
C----------------------------------------------------------------------|
	if(NLOAD.ne.0)	goto	30
C First entry:---------------------------------------------------------|
C Read model features: model, exp data; trunk values, names, scales, ...
C-------
C Read file "tmp/astra.log".  Define RSNAME, FILNAD
C    RSNAME	View-file name for the last run
C    FILNAD	View-file name if explicitly specified
	jj = length(RSNAME)
	if (jj .gt. 0) inquire(file=RSNAME(1:jj),exist=EXIRES)
	if (EXIRES) then
C	   write(*,*)jj,'RSNAME = "',RSNAME(1:jj),'"'
	   goto	24
	endif
C	write(*,*)jj,' ? RSNAME = "',RSNAME,'"'
	FILENA='tmp/astra.log'
	call	OPENRD(1,FILENA,0,IERR)
	if(IERR.eq.0)	goto	22
 21	pause	'>>>> Start file reading error <<<<'
 22	read(1,'(1A79)',ERR=21,END=23)RUNAST
	VNAM=RUNAST(1:6)
	if(VNAM.eq.'PROFT:')	RSNAME	=RUNAST(7:)
	if(VNAM.eq.'RSNAM:')	FILNAD	=RUNAST(7:)
	goto 22
 23	close (1)
	if(length(RSNAME) .eq. 0)	then
	   RSNAME = FILNAD(1:length(FILNAD))
	else
	   RSNAME = '.res/'//RSNAME(1:length(RSNAME))
	endif
C	write(*,*)'FILNAD = "',FILNAD,'"'
C	write(*,*)'RSNAME = "',RSNAME,'"'
 24	continue
C-------
C Initial information, profile number, names and scales
C    RDNAME(1:8) - data  file name
C    EQNAME(1:8) - model file name
	XLINE1(1:1) = ' '
	call 	OPENRD(3,RSNAME,1,IERR)
	if(IERR.gt.1)	then
	    write(*,*)'>>> Post-view file does not exist'
	    call	a_stop
	endif
	NCHL = 2 ! Scratch file with model.log to be connected to unit 2
	IRTYPE = SKIPM(3,NCHL)
C IRTYPE == 0 - old format (version < 5.2)
C	write(*,*)"Done",IRTYPE
	if (IRTYPE.ne.0)	then
	   read(3) RDNAME,EQNAME,VERSION,XLINE1
     ,	,	IYEAR,IMONTH,IDAY,IHOUR,IMINUT,NCFNAM,NPRNAM
     ,	,	NROUT,(NAMER(J),J=1,NROUT),(SCALER(J),J=1,NROUT)
     ,	,	NTOUT,(NAMET(J),J=1,NTOUT),(SCALET(J),J=1,NTOUT)
     ,	,	HRO,NB1,NSBR,NGR,NXOUT		!,(J4,j=1,7)
C	   write(*,*)RDNAME,EQNAME,VERSION,XLINE1
C     ,	,	IYEAR,IMONTH,IDAY,IHOUR,IMINUT,NCFNAM,NPRNAM
C     ,	,	NROUT,(NAMER(J),J=1,NROUT),(SCALER(J),J=1,NROUT)
C     ,	,	NROUT,(NAMER(J),J=1,NROUT),(SCALER(J),J=1,NROUT)
	else
	   read(3) RDNAME(1:8),EQNAME(1:8),XLINE1(2:17)
     ,	,	IYEAR,IMONTH,IDAY,IHOUR,IMINUT,NCFNAM,NPRNAM
     ,	,	NROUT,(NAMER(J),J=1,NROUT),(SCALER(J),J=1,NROUT)
     ,	,	NTOUT,(NAMET(J),J=1,NTOUT),(SCALET(J),J=1,NTOUT)
     ,	,	HRO,NB1,NSBR,NGR,NXOUT
C Append XLINE1 with spaces to avoid \0 at XLINE1(18:18) by usage from UPSTR
	   XLINE1=XLINE1(1:17)//"      "
	endif
	do	j=1,NRW
	   GRAL(j) = 0.
	   GRAP(j) = 1.E5
	enddo
C-------
C Read file ".exe/version" and define the current Astra version
	call	OPENRD(1,'.exe/version',0,IERR)
	if (IERR .gt. 1)	then
	   write(*,*)'>>> Warning: Current version is not determined'
	   goto	231
	endif
	do j=1,5
	   read(1,'(1A60)')STRI
	enddo
	j = index(STRI,'Version')
	RUNAST = STRI(j:j+30)//char(0)
	close(1)
C-------
	j0 = index(RUNAST,'.')
	if (j0 .eq. 0)	then
	   write(*,*)'>>> Warning: Current version is not determined'
	   goto	231
	endif
	read(RUNAST(j0-1:j0-1),*,ERR=231)AVERS 
	read(RUNAST(j0+1:j0+1),*,ERR=231)ARLEAS 
	j1 = INDEX(RUNAST(j0+1:),'.')
	if (j1 .eq. 0)	then
	   AEDIT = 0
	else
	   read(RUNAST(j0+j1+1:j0+j1+1),*,ERR=231)AEDIT 
	endif
	write(RUNAST,'(2(1I1,1A1),1I1)')AVERS,'.',ARLEAS,'.',AEDIT
C	write(*,'(A,2(1I1,1A1),1I1)')"Current  version: "
C     &,				AVERS,'.',ARLEAS,'.',AEDIT
	goto	232
 231	continue
	RUNAST = '5.3.?'//char(0)

 232	continue
	VVERS = '?.?.?'
	j0 = INDEX(VERSION,'.')
	if (j0 .eq. 0)	then
	   write(*,*)'>>> Warning: Unknown post-view file version'
	   goto	233
	endif
	read(VERSION(j0-1:j0-1),*,ERR=233)AVERS 
	read(VERSION(j0+1:j0+1),*,ERR=233)ARLEAS 
	j1 = INDEX(VERSION(j0+1:),'.')
	if (j1 .eq. 0)	then
	   AEDIT = 0
	else
	   read(VERSION(j0+j1+1:j0+j1+1),*,ERR=233)AEDIT 
	endif
	write(VVERS,'(2(1I1,1A1),1I1)')AVERS,'.',ARLEAS,'.',AEDIT
C	write(*,'(A,2(1I1,1A1),1I1)')"Run-time version: "
C     &,				AVERS,'.',ARLEAS,'.',AEDIT
 233	continue
	if (VVERS(1:3).ne.RUNAST(1:3) .and. RUNAST(1:3).ne.'6.2') then
	   write(*,*)char(7)		! Beep
	   if (VVERS(1:3) .eq. '?.?')	then
	      write(*,'(1A)')'  Astra version at run time UNKNOWN'
	   else
	      write(*,'(1A)')'  Astra version at run time: '//VVERS
	   endif
	   write(*,'(1A)')'  Current Astra version:     '//RUNAST(1:5)
	   write(*,*)">>> Warning >>> "
     >		,"The requested file might be not compatible with this "
	   write(*,*)"                "
     >  	,"Astra-viewer version. Some data can be corrupted." 
	   write(*,*)"    Please use appropriate viewer."
	   write(*,*)
C	   call	a_stop
	endif
C----------------------------------------------------------------------|
C !! NSBR is not used !!
	if (NXOUT .gt. 0 .and. NGR .gt. 0)	then
C Total length: 3*NGR*int+(3*NGR+GDEY(NGR)+NGRIDX(NGR)-1)*real+NARRX*int
	    read(3)
     .		(KTO(j),j=1,NGR),(NGRIDX(j),j=1,NGR),(NTYPEX(j),j=1,NGR)
     .	,	(TIMEX(j),j=1,NGR),(GDEX(j),j=1,NGR),(GDEY(j),j=1,NGR)
     .	,	(DATARR(j),j=1,GDEY(NGR)+NGRIDX(NGR)-1)
     .	,	(NAMEX(j),j=1,NARRX),(NWINDX(j),j=1,NARRX)
     .	,	(KOGDA(j),j=1,NARRX)
	endif
	DELOUT(13) = NB1
C	NRTIMS	=min(NRMAX,NCURVS/(NROUT+7))
	J4 = NCURVS
	NRTIMS	=min(NRMAX,J4*NRD/(NB1*(NROUT+7)))
	call APPRID(IDAY,IMONTH,IYEAR,IHOUR,IMINUT,J4)
C	write(*,*)IDAY,IMONTH,IYEAR,IHOUR,IMINUT
	TIME	=1.E32
	read(3,ERR=98,END=240)JJT
	read(3,ERR=98,END=240)TSTART
C----------------------------------------------------------------------|
	KROUT = 0
	KTOUT = 0
	goto	235
 234	continue
	read(3,ERR=98,END=240)JJT
	if (JJT .gt. 0)	read(3,ERR=98)YT		! TTOUT,TOUT   
 235	continue
	read(3,ERR=98,END=240)YT				! TIME
	KTOUT = KTOUT+JJT
	KROUT = KROUT+1
	if (KTOUT .le. NTMAX			.and.
     &	    KROUT .le. NRMAX			.and.
     &	    KROUT*(NROUT+7) .le. 5*NCURVS	.and.
     &	    KROUT*(NROUT+7)*NB1.le.NRD*NCURVS/4		! 4=sizeof(double)/2
     &	   )	then
	   j0 = KROUT
	   j4 = KTOUT
	endif
	read(3,ERR=98,END=240)YY				! CONSTF
	if (IRTYPE.ne.0) read(3,ERR=98,END=240) JJJ		! NA1,etc.
	do	J1=1,7+NROUT					! Equil. &
	   read(3,ERR=98,END=240)YY				! ROUT
	enddo
	goto	234
C----------------------------------------------------------------------|
 240	continue
	TEND = YT
C	write(*,*)
C	write(*,*)" Requested (r,t,p) ",KROUT,KTOUT,(NROUT+7)*KROUT*NB1
C     >			,KTOUT*NTOUT
C	write(*,*)" Allocated (r,t,p) ",NRMAX,NTMAX,NRD*NCURVS/4
C	write(*,*)"To be read (r,t,p) ",j0,j4,(NROUT+7)*j0*NB1
C	write(*,*)"Work space",NRD*(2*NRD+7)+NRD*2*NRD
C     >		,"    Xarray space",NRDX*NTARR,NRDX*NARRX
C	write(*,'(2(A,I12),1A,1F10.3)')" TOUT space",NTMAX*(NRW+1)
C     >		,"       ROUT space",NRD*NCURVS/4
C     >		,"     last time",TEND
	close(3)
	if (NCHL .eq. 0)	goto	250

	rewind(NCHL)
C Fill arrays PRNAME, DEVAR
	JINC = 0
	J = 0
 241	J = J+JINC
	read(NCHL,'(1A)',ERR=247)STRI
	if (STRI(1:10) .eq. ' Variables')	JINC = 1
	if (STRI(1:10) .eq. ' Constants')	goto 242
	if (J .gt. 0)	PRNAME(j) = STRI(1:6)
	if (J .gt. 0)	read(STRI(9:),*,ERR=241)DEVAR(j)
	goto	241
 242	J = J-1
	if (J+48 .ne. NPRNAM)	then
	 write(*,*)'>>> VIEW >>> Inconsistency in data file (variables)'
	endif
	if (J .gt. NCONST)	then
	   write(*,*)'>>> VIEW: More than NCONST variables'
	   call	a_stop
	endif
C Fill arrays CFNAME, CONSTF
	J = 0
 243	read(NCHL,'(1A)',ERR=247)STRI
	if (STRI(1:16) .eq. ' Control paramet')	goto 244
	J = J+1
	CFNAME(j) = STRI(1:6)
	read(STRI(9:),*,ERR=243)CONSTF(j)
	goto	243
 244	if (J .ne. NCFNAM)	then
	 write(*,*)'>>> VIEW >>> Inconsistency in data file (constants)'
	endif
	if (J .gt. NCONST)	then
	   write(*,*)'>>> VIEW: More than NCONST constants'
	   call	a_stop
	endif
C Fill arrays SRNAME, DELOUT (array DTNAME is set in BASTRA)
	j = index(STRI,':')
	read(STRI(j+1:),*)NSRNAM
	J = 0
 245	read(NCHL,'(1A)',ERR=247)STRI
	if (STRI(1:12) .eq. ' Color table')	goto 246
	J = J+1
	SRNAME(j) = STRI(1:6)
	read(STRI(9:),*,ERR=245)DELOUT(j)
	goto	245
 246	continue
	close(NCHL)
C Color table is not retrieved
	if (j .ne. NSRNAM)	then
	   write(*,*)'>>> VIEW >>> Inconsistency in data file (grids)'
	endif
	TIME = TSTART
C       write(*,*)NPRNAM,NCONST
C       write(*,*)(PRNAME(j),j=1,NPRNAM-48)
C       write(*,*)NCFNAM
C       write(*,*)(CFNAME(j),j=1,NCFNAM)
C       write(*,*)NSRNAM
C       write(*,*)(SRNAME(j),j=1,NSRNAM)
	goto	26
 247	write(*,*)'>>> Post-view file error after the string'
	write(*,*)STRI
	call	a_stop

 250	continue
C File model.log did not exist at the run time. Reading file const.inc
C to retrieve PRNAME, CFNAME, SRNAME
	FILENA='for/const.inc'
	call	OPENRD(1,FILENA,0,IERR)
	if(IERR.gt.0)	then
	   write(*,*)'>>> VIEW: File const.inc open error'
	   call	a_stop
	endif
C Read main variable list -> fill array PRNAME
		JINC=0
		J=0
 251		J=J+JINC
		if(J.gt.NCONST)	then
		   pause '>>> VIEW: More than NCONST variables'
		   goto 252
		endif
		read(1,'(2X,1A6)',ERR=254)VNAM
		if(VNAM.EQ.'Variab')	JINC=1
		if(VNAM.EQ.'End va')	goto 252
		if(J.gt.0)	PRNAME(J)=VNAM
		if(J.lt.NPRNAM)	goto 251
 252		JINC=0
C Read main constant list -> fill array CFNAME
		J=0
 253		J=J+JINC
		if(J.gt.NCONST)	then
		   pause '>>> VIEW: More than NCONST constants'
		   goto 255
		endif
		read(1,'(2X,1A6)',ERR=254)VNAM
		if(VNAM.EQ.'Consta')	JINC=1
		if(VNAM.EQ.'End co')	goto 256
		if(J.gt.0)		CFNAME(J)=VNAM
		if(J.lt.NCFNAM)		goto 253
		goto 255
 254		pause '>>> VIEW: File const.inc read error'

 255		JINC=0
C Read internal variable list -> fill array SRNAME
		NSRNAM=0
 256		NSRNAM=NSRNAM+JINC
		if(NSRNAM.gt.NCONST)	then
		   pause '>>> VIEW: More than NCONST internal variables'
		   goto 258
		endif
		read(1,'(2X,1A6)',ERR=257)VNAM
		if (VNAM .eq. 'Intern')	JINC=1
		if (VNAM .eq. 'End in')	then
		   NSRNAM = NSRNAM-1
		   goto 258
		endif
		if (NSRNAM .gt. 0)	SRNAME(NSRNAM)=VNAM
		goto 256
 257		pause '>>> VIEW: File const.inc read error'
 258		close(1)

 26	j1 = length(RSNAME)
	jj = 1
 27	j = nsymb(RSNAME(jj:j1),"/")
	if (j .ne. 0)	then
	   jj = jj+j
	   goto	27
	endif
	j4 = 0
	FILENA = RSNAME(jj:j1)//char(j4)
	j = j1-jj+2
	if (RSNAME(jj:j1) .eq. "profil.dat")	then
	   FILENA = "Last run"//char(j4)
	   j = 9
	endif
	FILENA = "Astra Viewer:  "//FILENA(1:j)
	j = j+16
C	write(*,*)'"',FILENA(1:j),'"',j

	jj = max(0,(15+NTOUT-64)/16)
	XWH = XWH+2*jj*(DYLET+2)
	call	initvm(XWX,XWY,XWW,XWH,COLTAB,FILENA,j)
	call	pcurso
	do	j=1,NB1
	   RHO(j) = (j-0.5)*HRO
	enddo
	jst = 0
	jsr = 0
C General (file start) information collected
	goto	40
C----------------------------------------------------------------------|
C____ Start reading after <PgDn> and <PgUp>
 30	continue
	jt = LOCTIM(TIME,TTOUT,NNTOUT)
	jr = LOCTIM(TIME,TROUT,NNROUT)
	if (NLOAD .lt. 0)	then
	   if (jt.lt.5.and.jst.eq.0 .or. jr.lt.2.and.jst.eq.0) then
	      STRI='>>> File at start. Loading data not possible'
	      write(*,'(/,1X,A)')STRI(1:45)
	      J0 = 0
	      J1 = XWH-124
	      call colovm(30)				! WarningColor
	      call textvm(j0,J1,STRI,45)
	      return
	   endif
C	   write(*,*)
C	   write(*,'(1A,F10.3)')" Load up.  Current time",TIME
	endif
	if (NLOAD .gt. 0)	then
	   if (KTOUT-jt.lt.5 .or. KROUT-jr.lt.1) then
	      STRI='>>> End of file. Loading data not possible'
	      write(*,'(/,1X,A)')STRI(1:42)
	      J0 = 0
	      J1 = XWH-124
	      call colovm(30)				! WarningColor
	      call textvm(j0,J1,STRI,42)
	      return
	   endif
C	   write(*,*)
C	   write(*,'(1A,F10.3)')" Load down.  Current time",TIME
	endif
	MTOUT	=0
	jst = 0
	jsr = 0
	jt = 0
	jr = 0
	j0 = 0
	call 	OPENRD(3,RSNAME,1,IERR)
	if(IERR.gt.1)	pause '>> VIEW: Profile file open error'
	j = 0
	IRTYPE = SKIPM(3,j)
	read(3)STR4
	if (NXOUT .gt. 0 .and. NGR .gt. 0) read(3)J	! KTO,...
 31	continue
	read(3,ERR=98,END=39)JJT
	if (jjt .gt. 0)	goto	32
	goto	33
 32	continue
	if (jt+jjt-jst.gt.NTMAX)	then
	   if (NLOAD.lt.0)	then
	      jst = jst+1
	      do	j=1,jt-jst
		 TTOUT(j) = TTOUT(j+1)
	      enddo
C	      write(*,'(1A,5I7,1F10.3)')"NTMAX overflow"
C     >		,jst,jjt,jt+jjt,NTMAX,jt+jjt,TTOUT(jt-jst+1)
	      goto	32
	   elseif (NLOAD.gt.0)	then
	      goto	39
	   else
	      write(*,*)">>> REVIEW >>> Logical error"
	      call	a_stop
	   endif
	endif
	read(3,ERR=98)(TTOUT(j-jst),(YY,j1=1,NTOUT),j=jt+1,jt+jjt)
	jt = jt+jjt					! Time point count
C	write(*,*)jst,jt,TTOUT(jt-jst)
	if (TIME .ge. TTOUT(jt-jst)) MTOUT = MTOUT+JJT
 33	continue
	jr = jr+1					! Radial point count
	j4 = (NROUT+7)*(jr-jsr)				! Profile count
	read(3,ERR=98,END=39)TROUT(jr-jsr)		! TIME
C	if (NLOAD .lt. 0)	goto	35
	if (NLOAD)35,34,34
 34	continue
C   NLOAD >= 0:
	if (TTOUT(jt-jst) .lt. TIME)	jst = jt
	if (TROUT(jr-jsr) .lt. TIME)	jsr = jr
	if (jr-jsr	  .ge. NRMAX)	goto	39
	if (j4     .gt. 5*NCURVS	.or.
     &      j4*NB1 .gt. NRD*NCURVS/4)	goto	39
	read(3,ERR=98,END=39)YY					! CONSTF
	if (IRTYPE.ne.0) read(3,ERR=98,END=39) JJJ		! NA1,NAB
	do	J1=1,7+NROUT					! Equil. &
	   read(3,ERR=98,END=39)YY				! ROUT
	enddo
	goto	31
C   NLOAD < 0:
 35	continue
	if (jr-jsr .ge. NRMAX)	then			! *(NRMAX) overflow
	   jsr = jsr+1
	   do	j=1,jr-jsr
	      TROUT(j) = TROUT(j+1)
	   enddo
	   write(*,*)"ARRAY(NRMAX) overflow",jsr,jr,TIME,TROUT(jr-jsr)
	endif
 36	continue
	j4 = (NROUT+7)*(jr-jsr)
	if (j4.gt.5*NCURVS .or. j4*NB1.gt.NRD*NCURVS/4)	then
	   jsr = jsr+1
	   do	j=1,jr-jsr
	      TROUT(j) = TROUT(j+1)
	   enddo
	   write(*,*)"RSCAL/RDOWN or int*2 RROUT array overflow",
     >	      j4,5*NCURVS,j4*NB1,NRD*NCURVS/4,jsr,jr,TIME,TROUT(jr-jsr)
	   goto	36
	endif
 38	continue
	read(3,ERR=98,END=39)YY					! CONSTF
	if (IRTYPE.ne.0) read(3,ERR=98,END=39) JJJ		! NA1,NAB
	do	J1=1,7+NROUT					! Equil. &
	   read(3,ERR=98,END=39)YY				! ROUT
	enddo
	if (TROUT(jr-jsr).gt.TIME .and.
     &      TTOUT(jt-jst).gt.TIME)	goto	39
	goto	31
 39	close(3)
CC	TEND = TROUT(jr-jsr)
CC	TINIT = max(TTOUT(max(1,jst)),TROUT(max(1,jsr)))
C	tst = max(TTOUT(1),TROUT(1))
C	ten = min(TTOUT(min(jt-jst,NTMAX)),TROUT(min(jr-jsr,NRMAX)))
C	write(*,'(1A,2F10.3)')"Time interval to plot",TST,TEN
C	write(*,'(1A,2F10.3)')" [TTOUT]"
C     >			,TTOUT(1),TTOUT(min(jt-jst,NTMAX))
C	write(*,'(3F10.3,1A,3F10.3)')(TTOUT(j),j=1,3),"  ...",
C     >		(TTOUT(j),j=min(jt-jst,NTMAX)-2,min(jt-jst,NTMAX))
CC	write(*,'(6F10.3,/,6F10.3)')(TTOUT(j),j=1,6),
CC     >		(TTOUT(j),j=min(jt-jst,NTMAX)-5,min(jt-jst,NTMAX))
C	write(*,*)"MTOUT,jst,jt,jsr,jr ",MTOUT,jst,jt,jsr,jr
C	write(*,*)"jt-jst, NTMAX ",jt-jst,NTMAX
C	write(*,'(1A,2F10.3)')" [TROUT]"
C     >			,TROUT(1),TROUT(min(jr-jsr,NRMAX))
C	write(*,'(3F10.3,1A,3F10.3)')(TROUT(j),j=1,3),"  ...",
C     >		(TROUT(j),j=min(jr-jsr,NRMAX)-2,min(jr-jsr,NRMAX))
C	write(*,*)"jr-jsr, NRMAX ",jr-jsr,NRMAX
C	j  = (jr-jsr)*(NROUT+7)
C	if (j .ge. 5*NCURVS) write(*,*)
C     >	   "Limiting factor:  (jr-jsr)*(NROUT+7) > 5*NCURV ",
C     >			j,5*NCURVS
C	if (j*NB1 .ge. NRD*NCURVS/4) write(*,*)
C     >	   "Limiting factor:  (jr-jsr)*(NROUT+7)*NB1, NRD*NCURVS/4",
C     >			j*NB1,NRD*NCURVS/4
	if(NLOAD.lt.0)	MTOUT = MTOUT+JJT		! max time # to read
	if(NLOAD.gt.0)	MTOUT = MTOUT-1
C----------------------- Main block start -----------------------------|
40	continue
C Start reading stored curves
C	if(NLOAD.eq.0) write(*,*)' VIEW: Loading data, please wait'
	NNROUT = 0				! Alias jr-jsr
	NNTOUT = 0				! Alias jt-jst=LTOUT
	IDENT  = 0			! Local use
	LTOUTO = 1			! Local use
	LTOUT1 = 1			! Local use
	call 	OPENRD(3,RSNAME,1,IERR)
	if(IERR.gt.1)	pause '>> VIEW: Profile file open error'
	j = 0
	IRTYPE = SKIPM(3,j)
	read(3)STR4
	if (NXOUT.gt.0 .and. NGR.gt.0)	read(3)JJ
 41	continue
	read(3,ERR=98,END=99)JTOUT
	if(NLOAD.lt.0 .and. LTOUT1.gt.MTOUT)	then
C	   write(*,'(1A,2I8)')"Exiting NLOAD < 0",LTOUT1,MTOUT
	   goto	99
	endif
	if(NLOAD.gt.0 .and. LTOUT1-1+JTOUT.lt.MTOUT)	then
	   if(JTOUT.gt.0)	then
	      read(3,ERR=98)YT
C	      read(3,ERR=98)(YT,(YY,JJ=1,NTOUT),J=1,JTOUT)
	      LTOUT1	=LTOUT1+JTOUT
	   endif
	   read(3,ERR=98,END=99)YT
	   read(3,ERR=98,END=99)YY		! (YY,j=1,NCFNAM+NPRNAM+3)
	   if (IRTYPE.ne.0)read(3,ERR=98,END=99)(JJ,j=1,12),(YY,j=1,10)
	   do	JJ=1,7+NROUT
	      read(3,ERR=98,END=99)YY
	   enddo
	   goto	41
	endif
C	write(*,*)NNROUT,NRTIMS,LTOUTO+JTOUT,NTMAX+1
	if(NNROUT.ge.NRTIMS.or.LTOUTO+JTOUT.ge.NTMAX+1)  goto	99
	if(JTOUT.eq.0)		goto 42
	LTOUT	=LTOUTO+JTOUT-1
	read(3,ERR=98)(TTOUT(J),(TOUT(J,JJ),JJ=1,NTOUT),J=LTOUTO,LTOUT)
	LTOUT1	=LTOUT1+JTOUT
	LTOUTO	=LTOUT+1
 42	read(3,ERR=98,END=99)YTIME
	NNROUT	=NNROUT+1
	if (IRTYPE.ne.0)	then
	   read(3,ERR=98,END=99)(CONS(J,NNROUT),J=1,NCFNAM)
     .		,(DDVAR(J,NNROUT),J=1,NPRNAM)
     .		,ABCT(NNROUT),ROCT(NNROUT),AVNE(NNROUT),QBND(NNROUT)
	   read(3,ERR=98,END=99)
     .	         INA1(NNROUT),INAB(NNROUT),(JJ,j=1,10),(YY,j=1,10)
	   J4	=(NROUT+7)*(NNROUT-1)
	   do	J1=1,7
	      J	=J1+J4
	      read(3,ERR=98,END=99)
     ,		RSCAL(J),RDOWN(J),(RROUT(JJ+IDENT),JJ=1,INAB(NNROUT))
	      IDENT = IDENT+INAB(NNROUT)
	   enddo
	   do	J1=8,NROUT+7
	      J	=J1+J4
	      read(3,ERR=98,END=99)
     ,		RSCAL(J),RDOWN(J),(RROUT(JJ+IDENT),JJ=1,INAB(NNROUT))
	      IDENT = IDENT+INAB(NNROUT)
	   enddo
	else
	   read(3,ERR=98,END=99)(CONS(J,NNROUT),J=1,NCFNAM)
     .		,(DDVAR(J,NNROUT),J=1,NPRNAM)
     .		,AVNE(NNROUT),ROCT(NNROUT),QBND(NNROUT)
	   J4	=(NROUT+7)*(NNROUT-1)
	   do	J1=1,7
	      J	=J1+J4
	      JJJ = NB1*(J1-1+J4)
	      read(3,ERR=98,END=99)
     ,			RSCAL(J),RDOWN(J),(RROUT(JJ+JJJ),JJ=1,NB1)
	   enddo
	   do	J1=8,NROUT+7
	      J	=J1+J4
	      JJJ = NB1*(J1-1+J4)
	      read(3,ERR=98,END=99)
     ,			RSCAL(J),RDOWN(J),(RROUT(JJ+JJJ),JJ=1,NB1)
	   enddo
	endif
	TROUT(NNROUT)	=YTIME
	if (NLOAD .eq. 0)	TIME = YTIME
	goto	41
C Main block end--------------------------------------------------------
98	close(3)
	write(*,*)'>>> View file read error: ',RSNAME(1:length(RSNAME))
	call	a_stop
C----------------------------------------------------------------------|
99	close(3)
	if(NNROUT.gt.NRTIMS) pause' profile # overflow '
	NNTOUT = LTOUT
	JJJ = (NROUT+7)*NNROUT*NB1
C	write(*,*)
C	write(*,*)' Total number of records (r,t): ',KROUT,KTOUT
C	write(*,*)' Number of records (r,t,p)  read: ',NNROUT,NNTOUT,JJJ
C	write(*,*)' Number of records (r,t,p) expected: ',jr-jsr,jt-jst
C	write(*,*)' nrtims, ltouto, JTOUT ',NRTIMS,LTOUTO,JTOUT,NLOAD
C	write(*,*)LTOUT,TTOUT(1),TTOUT(LTOUT),TSCALE,TTOUT(LTOUT)-TINIT
	if (NLOAD .eq. 0)		goto	5
	YY = TSCALE/5.
C	write(*,*)TTOUT(1)-YY,TINIT,TTOUT(1),TTOUT(LTOUT)
 3	continue
	if (TINIT .le. TTOUT(1))	goto	4
	TINIT = TINIT-YY
	goto	3
 4	continue
	if (TINIT .gt. TTOUT(1)-YY)	goto	5
	TINIT = TINIT+YY
	goto	4
 5	continue
C	write(*,*)TTOUT(1)-YY,TINIT,TTOUT(1)
	call	TIMOUT(TTOUT,TOUT)			! Define STOUT
C	do	JJ=1,NTOUT
C	   STOUT(JJ) = TOUT(LTOUT,JJ)
C	enddo
C	if (NLOAD.ne.0)		return
C	write(*,*)LTOUT,YY
	call	SETTSC(TTOUT(1),TTOUT(LTOUT))
	TSTART = TINIT
	do	j=1,NCONST
	   DEVAR(j)=DDVAR(j,1)
	enddo
	SCM = .1
	do	47	j=1,20
	if((RTOR+AWALL)/SCM8(j).gt.5..or.AWALL*ELONM/SCM8(j).gt.1.4)
     >						goto 47
	SCM = SCM8(j)
	return
 47	continue
	SCM=2.
	TIM7(3)	= abs(TSCALE)/8.
	TIM7(1)	= TINIT
C	write(*,*)TIM7(3),TINIT,TSCALE
	end
C======================================================================|
	subroutine	SETTSC(YTL,YTR)
C----------------------------------------------------------------------|
C This subroutine is incompatible with SETTSC(int) (file for/ifkey.f)
C----------------------------------------------------------------------|
C Input:  interval (YTL,YTR)
C Output: interval (TINIT,TSCALE)
C			TINIT < YTL < YTR < TINIT+TSCALE
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'		! TINIT,TSCALE,RTOR,AWALL
	integer	j,j0
	double precision YTL,YTR,Y,YS(10)
C----------------------------------------------------------------------|
	YS(1) = 1.0d-7
	YS(2) = 1.5d-7
	YS(3) = 2.0d-7
	YS(4) = 2.5d-7
	YS(5) = 5.0d-7
	YS(6) = 7.5d-7
	Y = YTR-YTL
 1	continue
	do	j=1,6
	   j0 = j
	   if (Y*0.97128378.le.YS(J)) goto 2
C		 (right_label_position)/(window_width)=575/592=0.97128378
	enddo
	do	j=1,6
	   YS(J)	=10.*YS(J)
	enddo
	goto 1
 2	continue
	TSCALE = YS(J0)
	Y = YS(j)/5.
	j = YTL/Y
	TINIT = j*Y
	if (j0.eq.2 .or. j0.eq.6)	then
	   Y = TSCALE/15
	else
	   Y = TSCALE/10
	endif
	if (YTL-TINIT.gt.Y)	TINIT = TINIT+Y
C	write(*,'(1A,1I4,6F10.3)')"SETTSC:",
C     >		j0,YTL,YTR,TINIT,TSCALE,TSCALE/5
	end
C======================================================================|
	subroutine	LININT(VARR,J1,IRTYPE,INAB,NROUT,NB1)
C VARR - radial array to retrieve
C J1   - No. of radial profiles to skip
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	integer	IRTYPE,INAB(*),J,J4,J40,J1,I1,I2,NROUT,NB1
	double precision	RSCAL,RDOWN,VARR(*)
	integer*2	RROUT(NRD*NCURVS)
	common	/AVIEW_RDAT/ RROUT,RSCAL(5*NCURVS),RDOWN(5*NCURVS)
C----------------------------------------------------------------------|
	if (IRTYPE.ne.0)	then
	   I1 = (J1-1)/(NROUT+7)
	   J4 = 0
	   if (I1.eq.0)	goto	1
	   do	J=1,I1
	      J4 = J4+(NROUT+7)*INAB(j)
	   enddo
 1	   continue
	   J40 = J4
	   I2 = J1-I1*(NROUT+7)-1
	   J4 = J4+I2*INAB(I1+1)
	   J40 = J40+(NROUT+7)*INAB(I1+1)
	   J40 = J40+I2*INAB(I1+2)
	   do 	J=1,INAB(I1+1)
	      VARR(J)=RDOWN(J1)+RSCAL(J1)*(RROUT(J+J4)+32768)/65535.
	   enddo
	else
	   J4 = J1-1
	   J4 = NB1*J4
	   do 	2	J=1,NB1
	      VARR(J)=RDOWN(J1)+RSCAL(J1)*(RROUT(J+J4)+32768)/65535.
 2	   continue
	endif
	end
C======================================================================|
	subroutine	RADOUT
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer J1,JJ,j,LOCTIM
	integer IROUT,ITOUT,NNROUT,NNTOUT,INA1,INAB,IRTYPE
	double precision AVNE,TROUT,CONS,DDVAR,ABCT,ROCT,QBND
	double precision YNE,QTIME,YCM,HA,YS,YS0,YS1,YDV
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     >			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
 	common	/AVIEW_CV/ CONS(NCONST,NRMAX),DDVAR(NCONST,NRMAX)
     .			,ABCT(NRMAX),ROCT(NRMAX),QBND(NRMAX)
C Get current record No.
	IROUT = LOCTIM(TIME,TROUT,NNROUT)
	J1 = (NROUT+7)*(IROUT-1)	! No. of preceeding radial records
	NA1 = INA1(IROUT)
	NAB = INAB(IROUT)

! No interpolation for constants and variables
	do	JJ=1,NCONST
	   CONSTF(JJ)= CONS(JJ,IROUT)
	   DEVAR(JJ) = DDVAR(JJ,IROUT)
	enddo

	if (IROUT.eq.1 .and. NNROUT.gt.1)	then
	   TAU = TROUT(IROUT+1)-TROUT(IROUT)
	elseif (NNROUT.eq.1)	then
	   TAU = TROUT(IROUT)
	else
	   TAU = TROUT(IROUT)-TROUT(IROUT-1)
	endif
	ABC = ABCT(IROUT)
	ROC = ROCT(IROUT)
	YNE = AVNE(IROUT)
	QTIME = QBND(IROUT)
	call	UPSTR(YNE,QTIME)
	call	TIMEDT(TIME,TAU)
C	write(*,*)"RADOUT:   ",ROC,QTIME,YNE
C	write(*,*)'IROUT=',IROUT,'    NROUT=',NROUT,'    NNROUT=',NNROUT
C     ,	,'   time:',time,"    J1=",J1
	call	LININT(AMETR,J1+1,IRTYPE,INAB,NROUT,NB1)
C	write(*,*)"AMETR:   ",(AMETR(j),j=1,NB1)
C Define NA,NA1
	do	j=1,NB1
	   RHO(j) = (j-0.5)*HRO
	enddo
	if (IRTYPE.eq.0)	then
	   do	jj=1,NB1
		if (RHO(jj) .lt. ROC)	NA=jj
	   enddo
	   NA1 = NA+1
	endif
	NA = NA1-1
	RHO(NA1) = ROC
	HROA	= ROC-RHO(NA)

C	write(*,*)
C	write(*,*)"RADOUT:   "
C	write(*,*)NA1,ABC,ROC,HRO
C	write(*,'(10f7.4)')(AMETR(j),j=NA1-8,NA1+1)
C	write(*,'(10f7.4)')(RHO(j),j=NA1-8,NA1)

	call	LININT(SHIF, J1+2,IRTYPE,INAB,NROUT,NB1)
	call	LININT(ELON, J1+3,IRTYPE,INAB,NROUT,NB1)
	call	LININT(TRIA, J1+4,IRTYPE,INAB,NROUT,NB1)
	call	LININT(FP,   J1+5,IRTYPE,INAB,NROUT,NB1)
	call	LININT(FP,   J1+6,IRTYPE,INAB,NROUT,NB1)
	call	LININT(FP,   J1+7,IRTYPE,INAB,NROUT,NB1)
C	write(*,'(4(10f7.3))')(FP(j),j=1,NB1)
C	write(*,'(4(10f7.3/))')(ROCT(j),j=1,NNROUT)
	J1	=J1+7
	do	J=1,NROUT
	   call	LININT(ROUT(1,J),J+J1,IRTYPE,INAB,NROUT,NB1)
	enddo
C	write(*,'(10f7.4)')(ROUT(j,1),j=NA1-8,NA1)
C	write(*,*)(ROUT(j,1),j=1,5)
C	write(*,*)NA1,NAB
C	write(*,*)ABC,ROC,AMETR(NA1),AMETR(NAB)
C Warning: accuracy in FP is not sufficient to define MU accurately
	YCM=1./(GP2*BTOR*HRO**2)
	do	JJ=1,NA
	   MU(JJ)	=YCM*(FP(JJ+1)-FP(JJ))/JJ
	enddo
	YCM = max(1.d-2,MU(NA-1)*(RHO(NA-1)/RHO(NA))**2)
	if (MU(NA) .lt. .8*YCM) MU(NA)=YCM
	MU(NA1)=MU(NA)*NA/(NA-.5+HROA/HRO)

	if (int(NEQUIL+0.1) .gt. 41)	goto	1
	YDV=0. 
	do J=2,NB1
	   VOLUM(J) = GP*GP2*AMETR(J)**2*ELON(J)*
     >		(RTOR+SHIF(J)-0.25*AMETR(J)*TRIA(J))
	   VRS(J-1) = (VOLUM(J)-YDV)/HRO
	   YDV	= VOLUM(J)
	   DRODA(J-1)=(RHO(J)-RHO(J-1))/(AMETR(J)-AMETR(J-1))
	enddo
	VOLUM(NA1)=GP*GP2*ABC**2*ELONG*(RTOR+SHIFT-0.25*ABC*TRIAN)
	VRS(NA)=(VOLUM(NA1)-VOLUM(NA))/(RHO(NA1)-RHO(NA1-1))
	if(NA1.lt.NB1)	
     >		VRS(NA1)=(VOLUM(NA1+1)-VOLUM(NA1))/(RHO(NA1+1)-RHO(NA1))
	DRODA(NA) = HROA/(AMETR(NA1)-AMETR(NA))

 1	continue
	if (IRTYPE.ne.0)	return
C       write(*,*)
C       write(*,*)NA1,NAB,ABC,AB,NB1,ROC,ROB
C       write(*,*)(RHO(j),j=NA-1,NB1)
C       write(*,'(10f7.3)')(AMETR(j),j=NA-1,NB1)
C piece of NEWGRD 
	if (AMETR(NA1) .ne. ABC)	AMETR(NA1) = ABC
	if (ABC .ge. AB)	then
	   NAB = NA1
	   ROB = RHO(NAB)
	   AMETR(NAB) = AB
	   return
	endif

	if (NA1 .lt. NB1)	then
	   NAB = NB1
	   HA = (AB-ABC)/(NB1-NA1)
	   if (HA .lt. ABC-AMETR(NA))	HA=ABC-AMETR(NA)
	   YS  = 2.*ELONG*RTOR/(RTOR+SHIFT)
	   YS0 = ABC**2/(1.+sqrt(1.-(ABC/(RTOR+SHIFT))**2))
	   do	j=NA1+1,NB1
C AMETR is taken from the profile.dat !
C             AMETR(j) = AMETR(j-1)+HA
	      YS1 = AMETR(j)**2/(1+sqrt(1-(AMETR(j)/(RTOR+SHIFT))**2))
	      RHO(j) = sqrt(ROC**2+YS*(YS1-YS0))
	      if (AMETR(j) .lt. AB)   NAB = j
	   enddo
	   if (NAB .lt. NB1)	NAB=NAB+1
	else
	   NAB = NA1
	endif
	ROB = RHO(NAB)
	AMETR(NAB) = AB

C       write(*,*)NA1,NAB,NB1,ABC,AB,ROC,ROB
C       write(*,'(10f7.4)')(RHO(j),j=NA-1,NB1)
C       write(*,'(10f7.4)')(AMETR(j),j=NA-1,NAB)
	end
C======================================================================|
	subroutine	TIMOUT(TTOUT,TOUT)
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	integer	icall,
     >		LOCTIM,NNTOUT,ITOUT,JJ,IROUT,NNROUT,INA1,INAB,IRTYPE
	double precision STOUT(NRMAX),TTOUT(NTMAX),TOUT(NTMAX,NRW)
	double precision TROUT,AVNE
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     >			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
	save	icall,STOUT
	data	icall/0/
	if (icall .eq. 0)	then
	   do	JJ=1,NTOUT
	      STOUT(JJ) = TOUT(LTOUT,JJ)
	   enddo
	   icall = 1
	   return
	endif

	ITOUT = LOCTIM(TIME,TTOUT,NNTOUT)
C	write(*,'(A,F10.3,3I5)')"TIMOUT: ",TIME,ITOUT,LTOUT,NNTOUT
C	if (ITOUT.gt.1)write(*,*)TOUT(ITOUT,1),TOUT(LTOUT,1),STOUT(1)
	if (NNTOUT .eq. 1)	then
	   TAU = 0.
	elseif (ITOUT .eq. 1)	then
	   TAU = TTOUT(ITOUT+1)-TTOUT(ITOUT)
	else
	   TAU = TTOUT(ITOUT)-TTOUT(ITOUT-1)
	endif

	if (MOD10 .eq. 6)	then
	   do	JJ=1,NTOUT
	      TOUT(LTOUT,JJ)	=STOUT(JJ)
	   enddo
	   return
	endif
	if (ITOUT .eq. 1 .or. NNTOUT .eq. 1)	then
	   do	JJ=1,NTOUT
	      TOUT(LTOUT,JJ) = TOUT(1,JJ)
	   enddo
	elseif (ITOUT .eq. LTOUT)	then
	   do	JJ=1,NTOUT
	      TOUT(LTOUT,JJ) = STOUT(JJ)
	   enddo
	else
	   do	JJ=1,NTOUT
	      TOUT(LTOUT,JJ) = TOUT(ITOUT-1,JJ)+(TIME-TTOUT(ITOUT-1))*
     .	      (TOUT(ITOUT,JJ)-TOUT(ITOUT-1,JJ))/TAU
	   enddo
	endif
C	write(*,'((4F10.3))')
C     .,	TTOUT(ITOUT),TIME,TTOUT(ITOUT-1),TOUT(ITOUT,1)-TOUT(ITOUT-1,1)
C     .,	TOUT(ITOUT,1),TOUT(ITOUT-1,1),TOUT(LTOUT,1)
	end
C======================================================================|
	subroutine	ARRPOS
C Draw arrow in a time mode
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
	integer	JXM,JYM,JYP,J,JX
	integer IROUT,ITOUT,NNROUT,NNTOUT,INA1,INAB,IRTYPE
	double precision AVNE,TROUT,XX,XSC0
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     >			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
	save	JXM
	data	JXM/-1/
	call	colovm(2)			! Red arrow
	JYM	=DYLET+2
	JYP	=JYM+3
C right_label_position=IDX*IDT=23*5*5=575 (see typdsp.f)
	XX	=575/abs(TSCALE)
	XSC0	=6*DXLET
C Mark radial time slices
	do	13	J=1,NNROUT
	JX	=XSC0+XX*(TROUT(J)-TINIT)
	if (JX.gt.XSC0)  call	drawvm(0,JX,JYM,JX,JYP)
 13	continue
C Erase & draw red arrow
	do	14	J=0,1
	if (j.eq.0)	call	colovm(0)
	if (j.eq.1)	call	colovm(2)
	call 	drawvm(0,JXM,JYM,JXM,JYM+8)
	call 	drawvm(0,JXM,JYM,JXM+2,JYM+3)
	call 	drawvm(0,JXM,JYM,JXM-2,JYM+3)
	JXM	=XSC0+XX*(TIME-TINIT)
 14	continue
	call	colovm(1)			! Return to black
	end
C======================================================================|
	integer	function LOCTIM(TIME,TIMES,NNOUT)
C The function returns
C   if NNOUT=1  then  LOCTIM=1,  otherwise
C   LOCTIM = an index of the array TIMES(1:NNOUT)
	implicit none
	double precision	TIME,TIMES(*)
	integer	NNOUT,j

	if (NNOUT .le. 1)	then
	   LOCTIM = 1
	   return
	endif
	if (TIME .le. TIMES(1))	then
	   LOCTIM = 1
	   TIME = TIMES(1)
	   return
	endif
	do	3	j=2,NNOUT
	   if (TIME .gt. TIMES(j))  goto	3
	   if (TIME .ge. 0.5*(TIMES(j)+TIMES(j-1)))	then 
	      LOCTIM = j
	   else
	      LOCTIM = j-1
	   endif
	   return
 3	continue
	LOCTIM = NNOUT
	TIME = TIMES(NNOUT)
	end
C======================================================================|
	subroutine	DRAWSPFLUX(jifnew,DYLET,YS0,YSC8)
C----------------------------------------------------------------------|
C Redraw magnetic surfaces:
	implicit	none
	integer		j,jj,jifnew,YS0,DYLET
	double precision YSC8
	character*132	WARNING
C----------------------------------------------------------------------|
	J = 48
	WARNING = "Sorry this output does not work in the view mode"
	call	textbf(0,550-104,Warning,J)
	end
C======================================================================|
	integer	function IFKEY(IFKL)
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer ICVMX
C NTMAX   max No. of stored time slices for output in mode 6
C TOUT   - Time variables output array
C TTOUT  - time-coordinate array for time output [s] TTOUT(1:LTOUT<=NTMAX)
	double precision TTOUT(NTMAX),TOUT(NTMAX,NRW)
	equivalence (WORK1D(1),TTOUT(1)),(WORK1D(NTMAX+1),TOUT(1,1))
C ICVMX    max No. of curves in one screen (<= 32)
	parameter(ICVMX=32)
C PT dimension: 4*NRD(Mode 5,8) 320(7) 2*NTMAX(Mode 6) 2*NRD(Modes 1-4)
	integer PT(4*NRD+2*NTMAX),ITO(NTMAX,ICVMX+2),OUTFIG(NRW)
	character*6	NAMEP(NTMAX)
	integer IFKL,IFKEY_,IROUT,ITOUT,NNROUT,NNTOUT,FLEN,LRDN,LWRN
	character*80 LHELP(30),CNSFIL*40,STR,STRB,PSNAME,STRI*132
     &		,UFLNAM*14,UNAMES(NRW)*10,DEFUNA*10,OUTNAME(NRW)*8
	integer MARK,J,NN0,NN1,NN80,JJ,J1,J2,NTP,LRJJ,NP1,JAB
	integer NSHOT,J4,INT4,IRET,LOCTIM,MODEX,length,ierr,ikey,KILLBL
	integer INA1,INAB,IRTYPE,WarningColor,NTVIEW,IM,XSC0,XSC,NX,NY
	logical MODADD
	double precision PRMARK(NTMAX),ARGT(NRMAX),
     .			 DELOLD(24),RADII(NRD),RARR(NRD),YARR(NRD)
	double precision AVNE,TROUT,CONS,DDVAR,ABCT,ROCT,QBND,ALFA,TIMEB
	double precision YNE,YSC,YQ,YDT,ABD
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     .			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
 	common	/AVIEW_CV/ CONS(NCONST,NRMAX),DDVAR(NCONST,NRMAX)
     .			,ABCT(NRMAX),ROCT(NRMAX),QBND(NRMAX)
C	data	UNAMES/NRW*'udb/tmp '/	MODADD/.FALSE./
	save	WarningColor,NN0,NN1,NN80,LRJJ,MARK,LHELP,ITO
	data	WarningColor/30/	PRMARK/NTMAX*0./
	data	DEFUNA/'      .tmp'/	UFLNAM/'udb/          '/
	data	NN0/0/  NN1/1/  NN80/80/ LRJJ/426/	MARK/0/
	data	(LHELP(j),j=1,10)/
     1	' ^C or / - STOP',
     2	'  1,2,3  - Radial profiles',
     3	'  4,5    - Time evolution of radial profiles',
     4	'  6      - Time evolution of local/global quantities',
     5	'  7      - Traces in a phase space',
     6	'  8      - Magnetic flux surfaces',
     7	'  A      - Adjust colors',
     8	'  B & N  - Backward & forward screen scan',
     9	'  C & V  - Constants & main variables',
     &	'  D      - Time step, output times'/
	data	(LHELP(j),j=11,22)/
     1	'  F & P  - Store displayed curves in a file',
     2	'  G & Q  - Save the screen as a PostScript file',
     3	'  L      - Model listing',
     4	'  M      - Change Markers or select time slice[s] in Mode 4',
     5	'  R      - Refresh screen',
     6	'  S      - New scales',
     7	'  T      - Type numerical values of current curves',
     8	'  U      - Write U-file; (1D in modes 1-3,6; 2D in modes 4,5)',
     9	'  W      - New windows',
     &	'  X      - X-axis info',
     1	'  Y      - Shift curve up/down',
     2	'  Z      - Restore original time scale'/
	data	(LHELP(j),j=23,30)/ 
     3	'  <-  /  ->   - Move up/down in time',
     4	'[Home]/[End]  - Move to minimal/maximal time',
     5	'[PgUp]/[PgDn] - Load data up/down',
     6	5*' '/
	LHELP(25)(1:) = LHELP(25)(1:33)//
     5			'with respect to the current time'
C ASCII codes: ^C  3  ^P 16  ^R 18  Esc 27    . 46  / 47 
C	  0 48  1 49  2 50  3 51  4 52  5 53  6 54  7 55  8 56  9 57
C	  A 65  B 66  C 67  D 68  E 69  F 70  G 71  H 72  I 73  J 74
C	  K 75	L 76  M 77  N 78  O 79  P 80  Q 81  R 82  S 83  T 84
C	  U 85  V 86  W 87  X 88  Y 89  Z 90
C----------------------------------------------------------------------|
	entry	IFKEY_(IFKL)
	NTVIEW = NTMAX
	IFKEY = 0
	IFKEY_= 0
C	write(*,*)'Entry IFKL ="',char(IFKL),'"',IFKL,KEY
C KEY =/= 0     =>    special action
	if(IFKL.ne.0)	then
		KEY=IFKL
		IFKL=0
		goto 9
			endif
C__________'KEY' is unrecognized	Condition: IFKL = 0
1	if (KEY .ne. 0)	then
Condition: (1) IFKL = 0, i.e. ASCII = 0 and (2) KEY =/= 256
	   KEY=0
	else
Condition: (1) IFKL = 0, i.e. ASCII = 0	and (2) KEY = 256
C	   call waitevent(KEY,ix,iy)
C	   if (NEVENT.eq.65000)   call   PUTXY(IX,IY,NTVIEW,TTOUT,TOUT)
C	   if (NEVENT.eq.65001)   call   PUTXY(IX,IY,NTVIEW,TTOUT,TOUT)
C	   call PUTXY(IX,IY)
	endif
 9	continue
C__________'P'
 10	if(KEY.ne.80)	goto 12
		KEY=0
C  Optional output format (G.W.Pacher)
		YNE = AVNE(IROUT)
		call	TIMOUT(TTOUT,TOUT)
		call	TYPDSP(1,YNE,NTVIEW,TTOUT,TOUT)
			return
C__________'A'
 12	if(KEY.ne.65)	goto 121
C		call SETCOL(COLTAB)
C		call colovm(1)
C			goto 20
C__________'N'
 121	if(KEY.ne.78)	goto 13
	   if (MOD10.ge.8 .or. MOD10.lt.0)	goto	20
		NSCR(MOD10)=NSCR(MOD10)+1
		J  = NROUT
		JJ = 8				! MOD10 = 2,3,6
		if (MOD10.eq.6 .or. MOD10.eq.7)	J = NTOUT
		if (MOD10 .eq. 1)	JJ=16
		if (MOD10 .eq. 4)	JJ=2
		if (MOD10 .eq. 5)	JJ=2
		if (MOD10 .eq. 7)	JJ=4
		if (J .gt. JJ*NSCR(MOD10))	goto	20
			NSCR(MOD10)=0
			goto 20
C__________'B'
 13	if(KEY.ne.66)	goto 14
	   if (MOD10.ge.8 .or. MOD10.lt.0)	goto	20
		NSCR(MOD10)=NSCR(MOD10)-1
		if (NSCR(MOD10)  .ge. 0)	goto	20
		J  = NROUT
		JJ = 8				! MOD10 = 2,3,6
		if (MOD10.eq.6 .or. MOD10.eq.7)	J = NTOUT
		if (MOD10 .eq. 1)	JJ=16
		if (MOD10 .eq. 4)	JJ=2
		if (MOD10 .eq. 5)	JJ=2
		if (MOD10 .eq. 7)	JJ=4
		NSCR(MOD10) = (J-1)/JJ
			goto 20
C__________ '1', '2', '3', '4', '5', '6', '7', '8'
 14	if(KEY.LT.49.OR.KEY.gt.56)	goto 16
		if (MOD10 .ge. 2 .and. MOD10 .lt. 6)	then
			if (MOD10 .eq. KEY-48)	then
				MODEY = -MODEY
			else
				MODEY = 1
			endif
		endif
		if (MOD10 .eq. 6)	then
			if (KEY.eq.54)	then
				MODEY = MODEY+1
				if (MODEY.eq.2)	MODEY = -1
			else
				MODEY = 1
			endif
		endif
		if (MOD10 .ne. KEY-48)	then
		   MOD10 = KEY-48
C		   NSCR=0
		endif
		goto 20
C__________'F'
 16	if(KEY.ne.70)	goto 18
		KEY=0
		YNE = AVNE(IROUT)
		call	TIMOUT(TTOUT,TOUT)
		call	TYPDSP(0,YNE,NTVIEW,TTOUT,TOUT)
			return
C__________'R'
 18	if(KEY.ne.82)	goto 21
 20		continue
C		write(*,*)KEY
		KEY=0
C		if(IFKL.eq.256)	return
		call	erasrw
		IM = 1
		NST = 0
		if (MOD10 .eq. 7)	call NEGA(IM,NTVIEW,TOUT)
		call	SETFRM(IM,XSC0,XSC,NX,NY)
		call	AFRAME(IM,XSC0,XSC,NX,NY)
		call	TIMOUT(TTOUT,TOUT)
		call	RADOUT
		if (NXOUT .gt. 0 .and. NGR .gt. 0) call	SETARX(0)
		J4 = XOUT+.49
		call	VIEWMN(J4)
		J4 = 80
		call	textbf(0,XWH-104,RUNID,J4)	! Task ID
		J4 = 0
		if(MOD10.le.3)	then
			TIME = TROUT(IROUT)
			call	OUTDSP(MARK,11,PT,ITO,NTVIEW,TTOUT,TOUT)
			call	DNSTR(J4,NTVIEW,TOUT)
		endif
		if(MOD10.eq.4 .or. MOD10.eq.5)	then
			YQ  = QBND(IROUT)
			YNE = AVNE(IROUT)
			call	UPSTR(YNE,YQ)
			call	TIMEDT(TIME,TAU)
			call	TIMOUT(TTOUT,TOUT)
			call	SMODE4(MARK,PT,PRMARK,NAMEP)
			call	DNSTR(J4,NTVIEW,TOUT)
		endif
		if(MOD10.eq.6)	then
			YQ  = QBND(IROUT)
			YNE = AVNE(IROUT)
			call	UPSTR(YNE,YQ)
			call	OUTDSP(MARK,11,PT,ITO,NTVIEW,TTOUT,TOUT)
			j4    = LOCTIM(TIME,TTOUT,LTOUT)
			YDT   = 0.
			if (j4 .lt. LTOUT)  YDT = TTOUT(j4+1)-TTOUT(j4)
			call	TIMEDT(TIME,YDT)
			call	DNSTR(j4,NTVIEW,TOUT)
C			call	colovm(2)	! Time in red
C			STRI(1:6)   = 'Time ='
C			STRI(12:14) = ' s'
C			call	FMTF50(STRI(7:11),TIME)
C			call	textvm(DXLET,426-4*DYLET+DYLET/3,STRI,14)
			if (KPRI .eq. 0)	call	ARRPOS
		endif
		if(MOD10.eq.7)	then
			IROUT = LOCTIM(TIME,TROUT,NNROUT)
			YQ    = QBND(IROUT)
			YNE   = AVNE(IROUT)
			call	UPSTR(YNE,YQ)
			call	TIMEDT(TIME,TAU)
			call	OUTDSP(MARK,11,PT,ITO,NTVIEW,TTOUT,TOUT)
			call	MARKTI(TIME,NTVIEW,TTOUT,TOUT)
			call	DNSTR(J4,NTVIEW,TOUT)
		endif
		if(MOD10.eq.8)	then
			TIME = TROUT(IROUT)
			call	OUTDSP(MARK,11,PT,ITO,NTVIEW,TTOUT,TOUT)
		endif
		if(KPRI.ge.1 .and. KPRI.le.2)	then
			STRI='>>>  The figure is stored in the file: '
     /					//PSNAME(1:FLEN)
			write(*,*)STRI(1:39+FLEN)
			call	colovm(1)
			if (KPRI .eq. 1)	call	PCOVA
			STRI='The figure is stored in the file: '
     /					//PSNAME(1:FLEN)
			call	textvm(40,980,STRI,34+FLEN)
			call	PSCLOS
			KPRI = 0
		endif
	return
C__________'U'
 21	if (KEY.ne.85)	goto	220
	if (MOD10.ge.7)	goto	220
	if (NUF.gt.NRD) goto	91
	call	TIMOUT(TTOUT,TOUT)
C Triangularity corrected MHD q (accoding to ITER guidelines)
C	YQ	=ELON(NA)**2
C	YD	=TRIA(NA)
C	YQ=(1.+YQ*(1.+YD**2*(2.-1.2*YD)))/(MU(NA)*(1.+YQ))
C----------------------------------------------------------------------
      if(MOD10.lt.1 .or. MOD10.gt.3)   goto   211
	INT4 = NROUT
	call	ASKUNA(INT4,UNAMES,NAMER,DEFUNA)
	call	RADOUT
	YNE = AVNE(IROUT)
	MODEX = XOUT+.49
	if	(MODEX .eq. 1)	then
C			 Write up to ABC against "a"
	   NP1 = NA1
	   ABD = ABC
	   do	j = 1,NP1
		RARR(j) = AMETR(j)/ABC
	   enddo
	elseif	(MODEX .eq. 2)	then
C			 Write up to ROC against "rho"
	   NP1 = NA1
	   ABD = 1.
	   do	j = 1,NP1
		RARR(j) = RHO(j)/ROC
	   enddo
	else
C			 Write up to AB against "a"	(default)
	   if (MODEX .ne. 0)	write(*,*)
     ,	   ">>>  Warning: Unknown radial mode.  Writing anyway [0,AB]"
	   NP1 = NAB
	   ABD = AB
	   do	j = 1,NP1
		RARR(j) = AMETR(j)/AB
	   enddo
	endif
C Radial coordinate in a U-file in [m] (presently)
C	if (MOD10.eq.1) MODADD = .FALSE.
C	if (MOD10.ne.1) MODADD = .TRUE.
	do	210	jj=1,NROUT
C	    if(UNAMES(jj)(1:8).eq.DEFUNA(1:7)//' ')	goto 210
		if(UNAMES(jj).eq.DEFUNA)	goto 210
	        if (IRTYPE.eq.0)	then
		   do	j=1,NB1
		      if (RARR(j) .ge. 1.)	goto	232
		      JAB = j
		   enddo
 232		   continue
		   JAB = JAB+1
	        else
	           JAB = INAB(IROUT)
	        endif
		RARR(JAB) = 1.

		do	j=1,NUF
		   RADII(j)	=(j-1.)/(NUF-1.)
		enddo
		ALFA	=.0001
C Transfer to an equidistant radial grid
	CALL	SMOOTH(ALFA,JAB,ROUT(1,jj),RARR,NUF,YARR,RADII)
C Radial coordinate in a U-file in [m] (presently)
C or normalized (if the following cycle is suppressed)
		do	j=1,NUF
		   RADII(j) = RADII(j)*ABD
		enddo
	call	UF1DWA('AUGD',RUNID,UNAMES(jj),TIME,NAMER(jj),NUF,1,4,
     &		  RADII,YARR,RTOR,AB,BTOR,IPL,YNE,MODADD,MODEX)
 210	continue
	goto   219
 211	continue
C----------------------------------------------------------------------
      if (MOD10 .ne. 4)   goto   216
	IROUT =LOCTIM(TIME,TROUT,NNROUT)
	INT4 = NROUT
	call ASKUNA(INT4,UNAMES,NAMER,DEFUNA)
C NTP - total number of time slices: 
	NTP = 0
	do	212	j = 1,NNROUT
	if (PRMARK(j) .lt. 0.)	goto	212
		NTP = NTP+1
		ARGT(NTP) = TROUT(j)
 212	continue
	do	j = 1,NUF
	   RADII(j)	=(j-1.)/(NUF-1.)
	enddo

	do	2120	jj=1,NROUT
C		if(UNAMES(jj)(1:8).eq.DEFUNA(1:7)//' ')	goto 2120
C	NSHOT = 0
C	read(UNAMES(jj)(2:6),*)NSHOT
C	open(3,file='udb/'//UNAMES(jj),status='UNKNOWN')
		if(UNAMES(jj).eq.DEFUNA)	goto 2120
	read(UNAMES(jj)(2:6),*)NSHOT
	UFLNAM = 'udb/'//UNAMES(jj)
	open(3,file=UFLNAM,status='UNKNOWN')
	write(3,1001)NSHOT,'AUGD',2,0
	write(3,1002)
	write(3,1003)0
	if (NUF.lt.0)	goto	2121
C NOTE! a radially contiguous u-file is written
	write(3,1007)
	write(3,1005)
	write(3,1008)NAMER(jj)
	write(3,1009)4
	write(3,1012)NUF
	write(3,1011)NTP
C Radial coordinate in a U-file in [m]
	write(3,'(1P,6E13.5)')(RADII(J)*AB,J=1,NUF)
C or normalized
C	write(3,'(1P,6E13.5)')(RADII(J),J=1,NUF)
	write(3,'(1P,6E13.5)')(ARGT(J),J=1,NTP)

	do	J2 = 1,NNROUT
	    if (PRMARK(j2) .ge. 0.)   then
		J1	=(NROUT+7)*(j2-1)
		call	LININT(RARR,J1+1,IRTYPE,INAB,NROUT,NB1)
	        if (IRTYPE.eq.0)	then
		   do	j=1,NB1
		      RARR(J) = RARR(j)/AB
		      if (RARR(j) .ge. 1.)	goto	233
		      JAB = j
		   enddo
 233		   continue
		   JAB = JAB+1
	        else
	           JAB = INAB(j2)
	        endif
		RARR(JAB) = 1.
		J1  =(NROUT+7)*(J2-1)+7
		call  LININT(ROUT(1,jj),jj+J1,IRTYPE,INAB,NROUT,NB1)
C Transfer to an equidistant radial grid
		ALFA = 0.0001
		CALL  SMOOTH(ALFA,JAB,ROUT(1,jj),RARR,NUF,YARR,RADII)
		write(3,'(1P,6E13.5)')(YARR(J),J=1,NUF)
	    endif
	enddo
	write(3,*)' ;----END-OF-DATA---------------COMMENTS:-----------'
	write(3,'(1A80)')RUNID
	write(3,'(1A4,1F6.2,3(1A8,1F6.2),1A13,1F6.2,1A7)')
     +	    ' R =',RTOR,'m,   a =',AB,'m,   B =',BTOR,
     +	    'T,   I =',IPL,'MA,   <n_e> =',.1*YNE,'E20m^-3'
	goto	215

 2121	continue
C NOTE! a time contiguous u-file is written
	write(3,1005)
	write(3,1007)
	write(3,1008)NAMER(jj)
	write(3,1009)4
	write(3,1011)NTP
	write(3,1012)NUF
	write(3,'(1P,6E13.5)')(ARGT(J),J=1,NTP)
C Radial coordinate in a U-file in [m]
	write(3,'(1P,6E13.5)')(RADII(J)*AB,J=1,NUF)
C or normalized
C	write(3,'(1P,6E13.5)')(RADII(J),J=1,NUF)
	do  213  J2 = 1,NNROUT
	    if (PRMARK(j2) .lt. 0.)   goto   213
	    J1	=(NROUT+7)*(j2-1)
	    call   LININT(RARR,J1+1,IRTYPE,INAB,NROUT,NB1)
	    do	j=1,NB1
		RARR(J) = AMETR(j)/AB
		if (RARR(j) .lt. 1.)	then
		   JAB = j
		else
		   goto	234
		endif
	    enddo
 234	    continue
	    JAB = JAB+1
	    RARR(JAB) = 1.
	    J1  =(NROUT+7)*(J2-1)+7
	    call   LININT(ROUT(1,jj),jj+J1,IRTYPE,INAB,NROUT,NB1)
C Transfer to an equidistant radial grid
	    ALFA = 0.0001
	    CALL   SMOOTH(ALFA,JAB,ROUT(1,jj),RARR,NUF,YARR,RADII)
	    write(3,'(1P,6E13.5)')(YARR(J),J=1,NUF)
 213    continue
	write(3,*)' ;----END-OF-DATA---------------COMMENTS:-----------'
	write(3,'(1A80)')RUNID
	write(3,'(1A4,1F6.2,3(1A8,1F6.2),1A13,1F6.2,1A7)')
     +	    ' R =',RTOR,'m,   a =',AB,'m,   B =',BTOR,
     +	    'T,   I =',IPL,'MA,   <n_e> =',.1*YNE,'E20m^-3'
	if(.TRUE.)	goto 215
C The following part has to be replaced so that the model is read from
C 						the current view file
C The model listing can be optionally appended to the U-file:
	LWRN	= length(RUNID(69:77))
	write(3,*)'>>>  Model  "',RUNID(69:68+LWRN),'"  listing:'
	call	OPENRD(1,'equ/'//RUNID(69:68+LWRN),0,IERR)
	if(IERR.gt.0)   then
		write(*,*) '>>> U-file write error:  Model "equ/',
     +		RUNID(69:68+LWRN),'" not found'
		goto   214
	endif
	do  J1=1,10000
	   read (1,'(1A132)',END=214)STRI
	   write(3,*)STRI(1:length(STRI))
	enddo
 214	close(1)
 215	write(*,*)'>>>  Data "',NAMER(jj),
     &		'" are written in the 2D U-file "udb/',UNAMES(jj),'"'
	close(3)
 2120	continue
	goto   219
 1001	format(1I7.5,1A4,2I2,16X,'; Shot #  ')
 1002	format(1X,9HDD-MMM-YY,21X,
     +	'; Shot date -  Ufiles ASCII file system')
 1003	format(1I4,27X,'; # of associated scalars')
 1004	format(1P,1E13.6,18X,'; Scalar value, LABEL FOLLOWS:')
 1005	format(' Time:',15X,'sec',7X,'; Independent variable: time')
 1006	format(' Time:',15X,'sec',7X,';')
 1007	format(' MINOR RAD',11X,'m',9X,'; Independent variable: Rho')
 1008	format( ' Function            ??        ;',
     +		' Function name:   "',1A4,'"')
 1009	format(1I2,29X,'; Processing code: ..., 4 - Astra output')
 1011	format(1I11,20X,';-# of  time  pts  t, F(t) (data follow)')
 1012	format(1I11,20X,';-# of radial pts  X, F(X) (data follow)')
 216	continue
C----------------------------------------------------------------------
	if (MOD10 .ne. 5)   goto   217
	write(*,*)'>>>  WARNING: No U-file can be written in this mode'
C-----------------------------------------------------------------------
 217    if(MOD10.ne.6)   goto   219
C	call	TIMOUT(TTOUT,TOUT)
C	INT4 = NTOUT
	call ASKUNA(NTOUT,UNAMES,NAMET,DEFUNA)
	do	218	jj=1,NTOUT
C	if(UNAMES(jj)(1:8).eq.DEFUNA(1:7)//' ')	goto 218
	if(UNAMES(jj).eq.DEFUNA)	goto 218
C	MODADD = .TRUE.
	call	UF1DWA('AUGD',RUNID,UNAMES(jj),TIME,NAMET(jj),LTOUT,0,
     &		4,TTOUT(1),TOUT(1,jj),RTOR,AB,BTOR,IPL,YNE,MODADD,MODEX)
 218	continue
 219	continue
                  goto 20
C__________'W'
 220	if(KEY.ne.87)	goto 22
	if(MOD10.EQ.1)			call ASKINT(NROUT,NWIND1,NAMER)
	if(MOD10.EQ.2 .or. MOD10.EQ.3)	call ASKINT(NROUT,NWIND2,NAMER)
	if(MOD10.EQ.4 .or. MOD10.EQ.5)	call ASKINT(NROUT,NWIND4,NAMER)
	if(MOD10.EQ.6)			call ASKINT(NTOUT,NWIND3,NAMET)
	if(MOD10.EQ.7)			call ASKINT(NTOUT,NWIND7,NAMET)
			goto 20
C__________'L'
 22	if(KEY.ne.76)	goto 23
		YSC=1.333*YSCMAX
		call colovm(1)
		LWRN	= length(EQNAME)
		CNSFIL='equ/txt/'//EQNAME(1:LWRN)
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
		if(IERR.gt.0)	pause ' >>> IFKEY: Constants file error'
		if(IERR.gt.0)	goto 7
		IKEY=0
 224		IF(IKEY.LT.0)	goto 226
		do 221 J1=1,int(YSC/DYLET)-1
		IF(IKEY.LT.0)	goto 225
		read(1,'(1A80)',END=222)STR
			goto 225
 222		IKEY=-1
 225		if(IKEY.LT.0)		then
			STRB=' '
			call prntxt(STRB)
					else
			call prntxt(STR)
					endif
 221		continue
		KEY=0
		goto 224
 226		close (1)
			goto 20
C__________'S'
 23	if(KEY.ne.83)	goto 230
		if (MOD10.le.5)	then
			call ASKLIS(NROUT,SCALER,NAMER,5)
		else
			call ASKLIS(NTOUT,SCALET,NAMET,5)
		endif
			goto 20
C__________'Y'
 230	if(KEY.ne.89)	goto 231
		if (MOD10.le.5)	then
			call ASKLIS(NROUT,OSHIFR,NAMER,6)
		else
			call ASKLIS(NTOUT,OSHIFT,NAMET,6)
		endif
			goto 20
C__________'X'
 231	if (KEY.ne.88)	goto 24
	MODEX = XOUT+.49
	if	(MOD10 .eq. 0)	then
	   write(*,*)'X-axis:   none'
	elseif  (MOD10 .eq. 3)	then
	   write(*,101)'poloidal flux,  1 < j < NA1 =',NA1
	elseif  (MOD10 .eq. 4)	then
	   write(*,102)'0 < a < AB =',AB,'m     1 < j < NAB =',NAB
	elseif  (MOD10 .eq. 5)	then
	   write(*,101)'major radius in the mid-plane'
	elseif  (MOD10 .eq. 6)	then
	   write(*,*)'X-axis:   time [s]'
	elseif  (MOD10 .eq. 7)	then
	   write(*,*)'X-axis:   phase space'
	elseif  (MOD10 .eq. 8)	then
	   write(*,*)'X-axis:   major radius [m]'
	elseif  (MOD10 .eq. 9)	then
	   write(*,*)"User's plot"
	elseif  (MODEX .eq. 0)	then
	   write(*,102)'0 < a < AB =',AB,'m,     1 < j < NAB =',NAB
	elseif  (MODEX .eq. 1)	then
	   write(*,102)'0 < a < ABC =',ABC,'m,    1 < j < NA1 =',NA1
	elseif  (MODEX .eq. 2)	then
	   write(*,102)'0 < rho < ROC =',ROC,'m,    1 < j < NA1 =',NA1
	elseif  (MODEX .eq. 3)	then
	   write(*,101)'0 < Psi < FP(NA1),  1 < j < NA1 =',NA1
	else
	   write(*,*)'X-axis:   Unknown option'
	endif
 101	format('X-axis:   ',1A29,1I3)
 102	format('X-axis:   ',1A,F5.2,1A,1I3)
			goto 20
C__________'D'
 24	if(KEY.ne.68)	goto 25
		NDTNAM=1024
		do	j=1,24
		    DELOLD(j)=DELOUT(j)
		enddo
		TIMEB=TIME
		MODEX = DELOUT(15)+.49
		call ASKLIS(NDTNAM,DELOUT,DTNAME,3)
		do	j=1,3
		    DELOUT(j)=DELOLD(j)
		enddo
C Skip time
		do	j=5,10
		    DELOUT(j)=DELOLD(j)
		enddo
		DELOUT(13) = NB1
		NUF = DELOUT(14)
		do	j=16,24
		    DELOUT(j)=DELOLD(j)
		enddo
		j = DELOUT(15)+.49
		if (j.lt.0 .or. j.gt.3)	then
		   write(*,*)">>> Unknown X-axis. Re-definition ignored"
		   j = MODEX
		   DELOUT(15) = MODEX
		endif
		IF(TIME.GE.TIMEB)	goto 20
		TTOUT(LTOUT-1)=TIME
			goto 20
C__________'C'
 25	if(KEY.ne.67)	goto 26
	INT4 = NCFNAM+1000
	call ASKLIS(INT4,CONSTF,CFNAME,2)
			goto 20
C__________'V'
 26	if(KEY.ne.86)	goto 27
	INT4 = NPRNAM-48+1000
	call ASKLIS(INT4,DEVAR,PRNAME,1)
			goto 20
C__________'M'
 27	if(KEY.ne.77)	goto 28
	if (MOD10.eq.1)			call
     &	  ASXWIN(NROUT,NWIND1,NAMER,SCALER,OSHIFR,GRAL,GRAP,MOD10,MODEY)
	if (MOD10.eq.2.or.MOD10.eq.3)	call
     &	  ASXWIN(NROUT,NWIND2,NAMER,SCALER,OSHIFR,GRAL,GRAP,MOD10,MODEY)
	if (MOD10.eq.6)			call
     &	  ASTWIN(NTOUT,NWIND3,NAMET,SCALET,OSHIFT,MOD10,MODEY)
	if (MOD10.le.3 .or. MOD10.eq.6)		goto 20
	IF(MOD10.ne.7)	goto 271
	INT4 = 4
	call ASKLIS(INT4,TIM7,NAM7,7)
			goto 20
271	IF(MOD10.ne.4 .and. MOD10.ne.5)	goto 28
	INT4 = -MAX(4,NNROUT)
	call ASKLIS(INT4,PRMARK,NAMEP,8)
		goto 20
C__________'.'
 28	if(KEY.ne.46)	goto 30
		MARK = MARK+1
		if (MARK .eq. 2) MARK = -1
			goto 20
C__________'H'
 30	if(KEY.ne.72)	goto	31
		write(*,*)
		write(*,*)"The following keys are operable in ",
     >				"this mode:"
		do	J=1,22
		call prntxt(LHELP(J))
		enddo
		write(*,*)
		do	J=23,26
		call prntxt(LHELP(J))
		enddo
		KEY=0
		write(*,*)
		goto 20
C__________'/' & '^C' & 'Esc'
C 31	if(KEY.ne.47.AND.KEY.ne.3.AND.KEY.ne.27)	goto 32
 31	if(KEY.ne.47.AND.KEY.ne.27)	goto 32
		call endvm
		CALL	A_STOP
C__________'T'
 32	if(KEY.ne.84)	goto	42
		YNE = AVNE(IROUT)
		call	TIMOUT(TTOUT,TOUT)
		call	TYPDSP(5,YNE,NTVIEW,TTOUT,TOUT)
		goto 20
C__________<space>
 42	if(KEY.ne.32)	goto	43
		KEY=0
		call	TIMOUT(TTOUT,TOUT)
		call	RADOUT
		if (NXOUT .gt. 0 .and. NGR .gt. 0) call	SETARX(0)
		J4 = XOUT+.49
		call	VIEWMN(J4)
		J4 = 80
		J4 = 0
		call	textbf(0,XWH-104,RUNID,J4)	! Task ID
		if(MOD10.le.3)	then
			call	OUTDSP(MARK,10,PT,ITO,NTVIEW,TTOUT,TOUT)
			call	DNSTR(J4,NTVIEW,TOUT)
				endif
		if(MOD10.eq.4 .or. MOD10.eq.5)	then
			YQ    = QBND(IROUT)
			YNE   = AVNE(IROUT)
			call	UPSTR(YNE,YQ)
			call	TIMEDT(TIME,TAU)
			call	DNSTR(J4,NTVIEW,TOUT)
				endif
		if(MOD10.eq.6)	then
			IROUT = LOCTIM(TIME,TROUT,NNROUT)
			YQ    = QBND(IROUT)
			YNE   = AVNE(IROUT)
			call	UPSTR(YNE,YQ)
			call	ARRPOS
			j4    = LOCTIM(TIME,TTOUT,LTOUT)
			YDT   = 0.
			if (j4 .lt. LTOUT)  YDT = TTOUT(j4+1)-TTOUT(j4)
			call	TIMEDT(TIME,YDT)
			call	DNSTR(j4,NTVIEW,TOUT)
			call	colovm(2)
			STRI(1:5)   = 'Time='
			STRI(11:12) = 's'
			call	FMTF50(STRI(6:10),TIME)
			j2	= YSCMAX+3*DYLET+15
			call	textvm(DXLET,j2-3*DYLET+DYLET/2,STRI,12)
				endif
		if(MOD10.eq.7)	then
			IROUT = LOCTIM(TIME,TROUT,NNROUT)
			YQ    = QBND(IROUT)
			YNE   = AVNE(IROUT)
			call	UPSTR(YNE,YQ)
			call	TIMEDT(TIME,TAU)
	 		call	MARKTI(TIME,NTVIEW,TTOUT,TOUT)
			call	DNSTR(J4,NTVIEW,TOUT)
				endif
		if(MOD10.eq.8)	then
			call	OUTDSP(MARK,10,PT,ITO,NTVIEW,TTOUT,TOUT)
				endif
		goto	99
C__________'G' & 'Q'
 43	if(KEY.ne.71 .and. KEY.ne.81)	goto	44
	   J4 = 0
	   J1 = length(RDNAME)
	   J2 = length(EQNAME)
	   PSNAME='dat/'//RDNAME(1:J1)//'-'//EQNAME(1:J2)//char(J4)
	   FLEN = KILLBL(PSNAME,NN80)
	   if (KEY.eq.71)	J4 = NN0
	   if (KEY.eq.81)	J4 = NN1
 	   call PSOPEN (PSNAME, J4, IRET)
	   FLEN = KILLBL(PSNAME,NN80)
	   if(IRET.eq.0)	then
	      if (KEY .eq. 71)	KPRI = 1
	      if (KEY .eq. 81)	KPRI = 2
	      goto	20
	   endif
	   if(IRET.eq.1)	then
	      STRI='>>>  Can not open file: '//PSNAME(1:FLEN)//' '
	      call colovm(WarningColor)
	      call textvm(NN0,LRJJ,STRI,24+FLEN)
	   endif
	   KEY = 0
	   goto	7
C__________'O'
 44	if (KEY.ne.79)	goto	45
	if (MOD10 .eq. 6)  then
	   call	ASKXGR(NTOUT,NWIND3,NAMET,MODEY,jj,j1,j2,OUTFIG,OUTNAME)
	   call	TIMOUT(TTOUT,TOUT)
	   call	WRFIGS(jj,j1,j2,OUTFIG,OUTNAME,NWIND3,NTVIEW,TTOUT,TOUT)
	elseif (MOD10 .le. 3)	then
	   call	ASKXGR(NROUT,NWIND1,NAMER,MODEY,jj,j1,j2,OUTFIG,OUTNAME)
C	   write(*,'(1A8,1I6)')(outname(j),outfig(j),j=1,NROUT)
	   call	WRFIGS(jj,j1,j2,OUTFIG,OUTNAME,NWIND1,NTVIEW,TTOUT,TOUT)
	endif
	goto	7
C__________'Z'
 45	if (KEY.ne.90)	goto	7
	call	SETTSC(TTOUT(1),TTOUT(LTOUT))
C	write(*,*)"times",TINIT,TSCALE
	if (MOD10 .eq. 6)	goto	20
	goto	7
C__________
 7	continue
 99	return
 91	write(*,*)'>>> No. of radial points in a U-file is too large:'
	write(*,*)'           Decrease NUF or increase NB1'
	end
C======================================================================|
 	subroutine	SMODE4(MARK,PT,PRMARK,NAMEP)
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'		! NAB,NB1,AB
	include	'for/status.inc'	! SHIF,AMETR
	include	'for/outcmn.inc'
	double precision 	SC(NRW),TROUT,YS,YROUT,AVNE,PRMARK(*)
	integer MARK,PT(*),IYMN,IYMX,JDSP,STYL,JPOS,JX,IDSP(2*NRD),JYO,
     1		IX(2*NRD),JXO,IK(5),FSHIFT,j,IROUT,ITOUT,NP,IS,j1,jj,
     2		JAB,NP1,NNROUT,NNTOUT,JP,INA1,INAB,IRTYPE,lonlen,JTIM
	character NAMEP(*)*6,ST*9,CHAR4*4,STRI*80
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     >			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
	data 	FSHIFT/10/
C----------------------------------------------------------------------|
C Mode 4:
	IYMN=YSCMAX-IYM
	IYMX=YSCMAX-IY0
C Scale is determined by the current radial distribution
c	call	SCAL(NROUT,SC,SCALER,ROUT(4,1),NAB-3,NRD)
	JP = 0
C	write(*,*)(PRMARK(JTIM),JTIM=1,NNROUT)
	do   1	JTIM=1,NNROUT
	   if(PRMARK(JTIM).le.0.)	goto 1
	   JP = JP+1
 1	continue
C	write(*,*)"JP =",JP
	call	SCAL(NROUT,SC,SCALER,ROUT(4,1),NAB-3,NRD)
	do	10	JTIM=1,NNROUT
	CALL	FMTF4(NAMEP(JTIM),TROUT(JTIM))
	if(PRMARK(JTIM).lt.0.)	goto 10
	do	9	J=1,NROUT
	   NP	=NWIND4(J)-2*NSCR(MOD10)
	   JX	=(NP-1)*DXM
	   if(NP.ne.1.and.NP.ne.2)	goto 9
	   do	IS=1,5
	      IK(IS)	=0
	   enddo
	   YS = OSHIFR(j)
	   if (YS)	2,5,3
 2	   YS = -YS
	   call FMTF4(CHAR4,YS)
	   ST = '-'//CHAR4
	   goto	4
 3	   call FMTF4(CHAR4,YS)
	   ST = '+'//CHAR4
 4	   call textvm((NP-1)*DXM+4*DXLET,FSHIFT+2*DYLET+2,ST,5)
 5	   continue
C	write(*,*)"From smode4"

	   J1 = (NROUT+7)*(JTIM-1)
	   call	LININT(AMETR,J1+1,IRTYPE,INAB,NROUT,NB1)
	   call	LININT(SHIF, J1+2,IRTYPE,INAB,NROUT,NB1)
C	   call	LININT(ELON, J1+3,IRTYPE,INAB,NROUT,NB1)
C	   call	LININT(TRIA, J1+4,IRTYPE,INAB,NROUT,NB1)
C	   call	LININT(FP,   J1+7,IRTYPE,INAB,NROUT,NB1)
	   J1 = J1+7
	   call	LININT(ROUT(1,J),J+J1,IRTYPE,INAB,NROUT,NB1)
	   if (IRTYPE.eq.0)	then
	      do	jj=1,NB1
		 if (AMETR(jj) .ge. AB)	goto	6
		 JAB = jj
	      enddo
 6	      continue
	      if (JAB.lt.NB1
     &		 .and. 2.*AB.gt.AMETR(JAB+1)+AMETR(JAB)) JAB = JAB+1
	   else
	      JAB = INAB(JTIM)
	   endif

C	   if (JP.eq.0 .or. PRMARK(JTIM).ne.0)
C     &
C	       call	SCAL(NROUT,SC,SCALER,ROUT(4,1),NAB-3,NRD)
	   call	FMTF4(ST(1:4),SC(J))
	   ST(5:9)	=NAMER(J)
	   call	colovm(1)
	   call	textvm((NP-1)*DXM,DYLET+FSHIFT,ST,9)
	   do	JJ=1,JAB
	      YROUT=min(max((ROUT(JJ,J)+OSHIFR(j))/SC(J),-7.d0),7.d0)
	      JDSP	=10*(DYM*YROUT+IYMN)
	      if (MODEY .eq. -1)	JDSP = JDSP+10*DYM
	      IDSP(JJ) = 3500-min(max(JDSP,10*IYMN),10*IYMX)
	   enddo

	   if (MOD10 .eq. 4)	then
	      NP1 = JAB
	      do	JJ=1,JAB
		 IX(jj) = JX+320*AMETR(jj)/AMETR(JAB)
	      enddo
	   else
C Mode 5
C Take AMETR(x,t) and SHIF(x,t) from "profile.dat"
	      NP1 = 2*JAB
	      do jj = 1,JAB
		 IDSP(JAB+jj) = IDSP(jj)
	      enddo
	      do jj = 1,JAB
	      IX(JAB+jj) = JX+160*(1.+min(1.d0,(SHIF(jj)+AMETR(jj))/AB))
	     IX(JAB+1-jj)=JX+160*(1.+max(-1.d0,(SHIF(jj)-AMETR(jj))/AB))
		IDSP(JAB+1-jj) = IDSP(JAB+jj)
	      enddo
	      IX(1) = min(IX(1),JX)
	      IX(NP1) = max(IX(NP1),JX+320)
	   endif

	   j1 = 31
	   do	IS=1,5
	      if(PRMARK(JTIM) .eq. IS)	then
		 if (IK(IS) .ne. 0)	goto 8
		 IK(IS) = 1
		 j1 = IS+1
	      endif
	   enddo
 8	   continue
	   STYL = (j1-1)*MARK
	   if (PRMARK(JTIM) .eq. 0)	STYL=0
	   if (KPRI.ge.1 .and. KPRI.le.2)	then
	      write(STRI,'(1A6,1A6,1A1)')'Plot "',NAMEP(JTIM),'"'
	      call	pscom(STRI,lonlen(STRI))
	   endif
	   call	PLOTCR(NP1,-1,IX,JXO,IDSP,JYO,j1,STYL,PT)
	   if(STYL.eq.0)	goto 9
	   JPOS	=3*DXLET+PRMARK(JTIM)*DXM/6.5+DXM*(NP-1)
	   call	CMARKP(DYLET+FSHIFT,JPOS,NAMEP(JTIM),STYL)
 9	continue
 10	continue
	end
C======================================================================|
 	subroutine	MARKTI(TIME,ITIMES,TTOUT,TOUT)
	implicit none
	include	'for/parameter.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer ITIMES,NUM(4),JMIN,JMAX,NKL1,NKL2,MODK,IYMN,IYMX,JDSP
	integer	I,J,JJ,J1,JX,JY,JPOSX,JPOSY,PT(2),WarningColor
	common /AC_NEG1/ NUM,JMIN,JMAX,NKL1,NKL2,MODK(2)
	double precision
     1		TIME,SC(NRW),XROUT,YROUT,TTOUT(ITIMES),TOUT(ITIMES,NRW)
	save	WarningColor
	data	WarningColor/30/
	if(LTOUT.lt.6)		return
	IYMN=YSCMAX-IYM
	IYMX=YSCMAX-IY0
	I=0
	do 58 JJ=1,4
	   NUM(JJ)=0
	   do 59 J=1,NTOUT
	      if(NWIND7(J).EQ.JJ)	then
		 NUM(JJ)=J
		 goto 58
	      endif
 59	   continue
 58	continue
	NKL1=2
	NKL2=1
	IF(NUM(1).ne.0.AND.NUM(2).ne.0)NKL1=1
	IF(NUM(3).ne.0.AND.NUM(4).ne.0)NKL2=2
	IF(NKL1.EQ.2.AND.NKL2.EQ.1)	return
c* Time interval : Tmin7 < t < Tmax7, Mark interval - Tmet
	JMIN	=1
	JMAX	=0
	do	J=1,LTOUT-1
	   if(TTOUT(J) .le. TIM7(1))	JMIN=J
	   IF(TTOUT(J) .le. TIM7(2))	JMAX=J
	enddo
	if(JMIN.gt.JMAX)	return
	call SCAL(NTOUT,SC,SCALET,TOUT,ITIMES,ITIMES)
	call	colovm(WarningColor)
	do	JJ=nkl1*2,nkl2*2,2
	   JPOSY=MODK(JJ/2)*(IYM-IY0)/2	
	   JPOSX=JJ/4*XSCMAX/2+MODK(JJ/2)*XSCMAX/4
	   do	J=JMIN,JMAX
	      YROUT = max(TOUT(J,NUM(JJ-1))/SC(NUM(JJ-1)),-1.d0)
	      YROUT = min(YROUT,1.d0)
	      JDSP = (2-MODK(JJ/2))*(IYM-IY0)/2*YROUT+IYMN+JPOSY
	      JDSP = max(JDSP,-MODK(JJ/2)*(IYM-IY0)/2+IYMN)
	      J1 = J-JMIN+1
	      JY = min(JDSP,IYMX)
	      XROUT = min(max(TOUT(J,NUM(JJ))/SC(NUM(JJ)),-1.d0),1.d0)
	      JDSP = (2-MODK(JJ/2))*XSCMAX/4*XROUT
	      JDSP = max(JDSP,XSCMAX/4*(-MODK(JJ/2)))
	      JX = min(JDSP,XSCMAX/2)
	      if (TTOUT(J).le.TIME .and. TIME.lt.TTOUT(J+1))	then
		 PT(1) = JPOSX+JX
		 PT(2) =  350 -JY
		 call NMARK(PT,2)
		 return
	      endif
	   enddo
	enddo
	end
C======================================================================|
	blockdata BASTRA
	implicit none
	include	'for/parameter.inc'
	include	'for/pareview.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer	IROUT,ITOUT,NNROUT,NNTOUT,INA1,INAB,IRTYPE,NCONST48,j
	double precision AVNE,TROUT
	common	/AVIEW_MAIN/ AVNE(NRMAX),TROUT(NRMAX),IROUT,ITOUT,
     >			 NNROUT,NNTOUT,INA1(NRMAX),INAB(NRMAX),IRTYPE
	parameter(NCONST48=NCONST-48)
	data	IROUT/1/ ITOUT/1/ NNROUT/0/ NNTOUT/0/
C Default constants
	data	DTNAME/
C Time interval control names
     1		'dRout ','dTout ','dPout ','Time  ',
     2		'TAUmin','TAUmax','TAUinc','DELvar',
     3		'Iterex','N/used','Tinit ','Tscale',
     4		'NB1   ','NUF   ','Xaxis ','Xdeflt',
     5		'NB2EQL','NEQUIL','NBND  ','Xflag ',
     6		'DTeql ','MEQUIL','Tpause','Tend  ',
     7		80*'      '/
	data DELOUT/
C Output and time step control values
     1		.01,     .01,     .01,     0.,
     2		.00003,  .05,     1.1,     .1,
     3		1.,      41.,      0.,     1.,
     4		41.,     41.,      1.,     1.,
     5		1.,       1.,      1.,     1.,
     6		1.,       1.,      1.,     1.,
     7		80*.0/
C Initial times
	data	TAU/.00003/  TEQ/20*-99999./  LTOUT/1/  IPOUT/1/
C Screen parameters
	data	XSCMAX/640/ YSCMAX/350/ DXLET/8/ DYLET/13/ IDT/5/
	data	XWX/470/    XWY/10/     XWW/660/ XWH/550/
C Astra colors: #0 - background, ##1-6 - plots 1-6
C		#7, #8 - not used
C		#9 - delete curve, #10 - standard text color (black)
C		#11 - Black
C		#12 - warning messages, axis upper marks in View
C		#13 -> #15 -reserved (black)
C	data	COLTAB/
C     1		0,0,	50,0,	3, 0,	62,0,	37,0,
C     2		5,0,	26,0,	30,0,	50,0,	0, 0,
C     3		1,0,	1,0,	50,0,	1,0,	1,0,
C     4		1,0,	1,0,	1,0,	1,0,	1,0,
C     5		1,0,	1,0,	1,0,	1,0,	1,0,
C     6		1,0,	1,0,	1,0,	1,0,	1,0,
C     7		1,0,	1,0/
	data	TIM7/0,9999,9999,1/
     5		NAM7/'Tmin','Tmax','Tmark','Style'/  NSCR/9*0/	MOD10/1/
     6		IFDFVX /NCONST*-1./ MODEY/1/
	data	NSBR /0/	NXOUT /0/
	data    EXARNM /
     1		'CAR1X ','CAR2X ','CAR3X ','CAR4X ',
     2		'CAR5X ','CAR6X ','CAR7X ','CAR8X ',
     3		'CAR9X ','CAR10X','CAR11X','CAR12X',
     4		'CAR13X','CAR14X','CAR15X','CAR16X',
     5		'CAR17X','CAR18X','CAR19X','CAR20X',
     6		'CAR21X','CAR22X','CAR23X','CAR24X',
     7		'CAR25X','CAR26X','CAR27X','CAR28X',
     8		'CAR29X','CAR30X','CAR31X','CAR32X',
     9		'F0X   ','F1X   ','F2X   ','F3X   ','F4X   ',
     &		'F5X   ','F6X   ','F7X   ','F8X   ','F9X   ',
     1		'MUX   ','MVX   ','GNX   ','SNX   ','PEX   ',
     2		'PIX   ','PRADX ','TEX   ','TIX   ','NEX   ',
     3		'CUX   ','ZEFX  ','VRX   ','SHX   ','ELX   ',
     4		'TRX   ','G11X  ','G22X  ','G33X  ','DRODAX',
     5		'IPOLX ',',NIX  ','VPOLX ','VTORX ','SLATX'/
C Metric array data
C	data	VRX   /41*1./     G11X  /41*1./      G22X  /41*0./
C     1		SHX   /41*1./     ELX   /41*1./      TRX   /41*1./
C     2		DRODAX/41*1./     IPOLX /41*1./      G33X  /41*1./
C     3         RESX  /41*0./
C General array data
C	data	ELON  /41*1./     TRIA  /41*0./      SHIF  /41*0./
C     2		TE    /41*.01/    TI    /41*.01/     NE    /41*.1/
C     3		CU    /41*.1/
C     4		MU    /41*.3/     ZEF   /41*1./      AMAIN /41*1./
C     5		ZMAIN /41*1./     GP    /3.1415926/  GP2   /6.283185/
C     7		G33 /41*1./	  IPOL  /41*1./      SQEPS /41*0.2/
C     8		VOLUM /41*1./     NI    /41*.1/      ZEFX  /41*1./
C	data	FP    /0.,1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,
C     &			13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,
C     &			23.,24.,25.,26.,27.,28.,29.,30.,31.,32.,
C     &			33.,34.,35.,36.,37.,38.,39.,40./
C Main variables data
	data	AB    	/.3/	ABC	/.3/	AMJ	/2./
     1		AWALL	/.4/	BTOR	/3./	ELONG	/1./
     2		ELONM	/1./	ENCL	/.002/	ENWM	/.02/
     3		IPL	/.3/	RTOR	/1.5/
     4		SHIFT	/0./	NNCL	/.001/	NNWM	/.0001/
     5		TRIAN	/.0/	TRICH	/.0/	ZMJ	/1./
     6		NB1	/41/	NA1	/41/	NUF	/41/
     7		NAB	/41/	GP  /3.1415926/	GP2  /6.283185/
C Constants data
	data	(CONSTF(j),j=1,16)	/16*1./
     1		(CONSTF(j),j=17,32)	/16*0./
     2		(CONSTF(j),j=33,40)	/8*1./
     3		(CONSTF(j),j=41,48)	/8*0./
     4		(CONSTF(j),j=49,NCONST)	/NCONST48*1./
	data	CNFILE/'***'/	NBFILE/'***'/	ECFILE/'***'/
	data	ICFILE/'***'/	MSFILE/'***'/	TASK/'VIEW'/
	end
C======================================================================|
