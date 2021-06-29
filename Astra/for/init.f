C Modules:	>>>>>> READAT, INTVAR, BASTRA <<<<<<<
C=======================================================================
	subroutine	READAT
C-----------------------------------------------------------------------
C  NCONST   amount of simple variables readable (initiated)
C	    from a data file (List in the file "const.inc") 
C  NTVAR    maximal number of time slices for all variables
C  IVAR	    number of actually defined variables
C  NTARR    maximal number of time slices for all arrays 
C  NGR	    number of actually defined groups
C  NARRX    maximal number of arrays recognizable from a data file 
C-----------------------------------------------------------------------
C The subroutine is called once at the start-up, it reads the "exp" file
C and stores the time evolution of all input data in the arrays 
C VARDAT (simple variables) and 
C DATARR (sequential data in groups [grid,quantity]), TIMEX (group time)
C NGRIDX (grid size), NTYPEX(grid type), KTO (internal group name)
C FILTER (smoothing/transfer parameter)
C GDEX   (relative position of grid in DATARR)
C GDEY   (relative position of profile in DATARR)
C KOGDA  (relative position of time in TIMEX)
C 	      
C	VARDAT(1,NTVAR)	- time
C	VARDAT(2,NTVAR)	- value
C	VARDAT(3,NTVAR)	- error (not used)
C	DATARR(NRDX*NTARR) - data array
C  	    Let   1 <= j <= NTARR is an ordinal number of array in DATARR
C	TIMEX(j)  - time for this array
C	NGRIDX(j) - number of grid points
C	NTYPEX(j) - type of grid
C	KTO(j)	  - pointer to a position in the array EXARNM
C		    so that EXARNM(KTO(j)) gives the name of the quantity
C	GDEX(j)   - pointer to grid in DATARR
C	GDEY(j)   - pointer to data in DATARR
C	XAXES(NRDX,j) - not used here
C  	    Let   1 <= jx <= NARRX is an ordinal number of q-ty EXARNM(jx)
C	KOGDA(jx)  - pointer to a position in the array TIMEX
C		     KOGDA(KTO(j)) gives pointer to 1st time for KTO(j)
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer     length,nsymb,noblan,killbl,nextwd,npnt,IFKEY,getpid
	integer     IFDEFX
	integer     INDVAR,IVAR,jarr,INDEX,INTYPE,lname,ishot,jtype,IM
	integer	    NCH,ERCODE,ITYPE,SYSTEM,lonlen,jbdry,jdim,jtim,jtim1
	integer     jj,jy,jt,j,j0,j1,j2,j3,JINC,IERR,JRUNN,jrad,jexar,NY
	integer	    KAB,KABC,KAWALL,KRTOR,KELONM,KTRICH,iarg,XSC0,XSC,NX
	logical     EXI,EXILOG,TCONTG,IFREE,EXIRAD,EXITIM,XFILE
	character   RUNAST*79,STRI*132,STRJ*132,STRAD*132,STRAD2*32,CH*1
	character*6 VNAM,VNAMO,VNAMU,VNAMX,VTIM,VDAT,VERR,VARNAM,COMNAM
	character*6 NAMECL(7),ARRNAM,KEYWRD,IDEV*4,FILENA*40,STR40*40
	double precision CHORDN,LINEAV
	double precision 
     1		YS(10),XBDRY,SCM8(20),YB,YB1,YXB,YXB1,YY,ALFA,YSMAX,
     2		VRDATA,VARDAT,FACTOR,TIMEVR,VRERR,ROC3A,YTP
	common	/EXPDAT/ VARDAT(3,NTVAR),INDVAR(NTVAR),IVAR
C----------------------------------------------------------------------|
C Ex/ext-file requisites:
	integer		MPEX,MSIGEX,MTEX,JSHOT,JPEX,MSIG,JTEX,MEXT
	parameter	(MPEX=101, MSIGEX=1, MTEX=50, MSIG=1)
	parameter	(MEXT=MPEX*MTEX)
	integer		IERSIG(MSIGEX)
	double precision TIMJEX(MTEX),
     &			 JEXSIG(MPEX,MSIGEX,MTEX),JEXTSG(MSIGEX,MEXT),
     &			 JEXRAD(MPEX,MSIGEX,MTEX),JEXTIM(MSIGEX,MEXT)
	equivalence	 (JEXSIG,JEXTSG),(JEXRAD,JEXTIM)
	character	NAMREQ(MSIG)*8,VERJEX*20
C----------------------------------------------------------------------|
	save	SCM8,NAMECL,IFREE
	data	SCM8/.1,.15,.2,.25,.3,.4,.5,.6,.8,1.,1.25,1.5,
     &			2.,2.5,3.,4.,5.,6.,8.,10./
	data	NAMECL/ 'POINTS','NAMEXP','GRIDTY','NTIMES',
     &			'FILTER','FACTOR','END   '/
	data	IFREE/.TRUE./	! Default input mode: read(u-file,*)
C----------------------------------------------------------------------|
C	call	STARTCLOCK
	call	ADDTIME(CPT)		! Initialize timer
C	do j=0,iargc()
C	   call	getarg(j,STRI)
C	   write(*,*)'"',STRI(1:length(STRI)),'"',j,iargc()
C	enddo
	AS_PID = getpid()
C	write(*,*)'Astra PID =',AS_PID
	call	getarg(0,STRI)
C	write(*,*)STRI
	j = min(len(TASKID),length(STRI))
	TASKID = STRI(1:j)
C	write(*,*)TASKID
C	if (iargc() .gt. 0)	then
C          call	getarg(1,STRI)
C	   j = length(STRI)
C	   if (STRI(1:j) .eq. 'background')	TASK = 'BGD '
C	endif
C----------------------------------------------------------------------|
	call	markloc("READAT"//char(0))
	do	j=1,NRW
	   NWIND1(j) = j
	   NWIND2(j) = j
	   NWIND3(j) = j
	   NWIND7(j) = j
	enddo
C	data	NWIND4/
C     1		1, 5, 9, 13,3, 7, 11,15,2, 6, 10,14,4, 8, 12,16,
C     2		17,21,25,29,19,23,27,31,18,22,26,30,20,24,28,32,
C     3		33,37,41,45,35,38,43,47,34,38,42,46,36,40,44,48,
C     4		49,53,57,61,51,55,59,63,50,54,58,62,52,56,60,64,
C     5		65,69,73,77,67,71,75,79,66,70,74,78,68,72,76,80,
C     6		81,85,89,93,83,87,91,95,82,86,90,94,84,88,92,96/
C     2		17,0,0,0,19,0,0,0,18,0,0,0,20,67*0./
C     &		1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,44*0./
	do	j=1,4
	   NWIND4(j) = 4*j-3
	   NWIND4(j+4) = NWIND4(j)+2
	   NWIND4(j+8) = NWIND4(j)+1
	   NWIND4(j+12) = NWIND4(j)+3
	enddo
	do	j=17,NRW
	   NWIND4(j) = NWIND4(j-16)+16
	enddo
	JRUNN = -1				! JAMS run ID (-1 none)
	do	j=1,NRD
	    FP(j) = 0.1*j
	    QE(j) = FP(j)
	    QI(j) = FP(j)
	    QN(j) = FP(j)
	    GN(j) = FP(j)
	    do	j1=1,2*NRD
	       WORK(j,j1) = 0.d0
	    enddo
	    do	j1=1,2*NRD+7
	       WORK1(j,j1) = 0.d0
	    enddo
	    do	j1=1,32
	       CAR(j,j1) = 0.d0
	    enddo
	    do	j1=1,NARRX
	       EXT(j,j1) = 0.d0
	    enddo
	enddo
	do	j=1,NTARR
	    FILTER(j) = 0.001
	    TIMEX(j) = 0.
	enddo
	do	j=1,NARRX
	    IFDFAX(j) = -1.
	enddo
C----------------------------------------------------------------------|
C Read file ".exe/version"
	call	OPENRD(1,'.exe/version',0,IERR)
	if (IERR .gt. 1)	then
	   write(*,*)'>>> Warning: Unknown version'
	   goto	432
	endif
	do j=1,5
	   read(1,'(1A60)')STRI
	enddo
	j = index(STRI,'Version')
	VERSION = STRI(j:j+30)//char(0)
	close(1)
	j0 = index(VERSION,'.')
	if (j0 .eq. 0)	then
	   write(*,*)'>>> Warning: Unknown version'
	   goto	432
	endif
	read(VERSION(j0-1:j0-1),*,ERR=432)AVERS 
	read(VERSION(j0+1:j0+1),*,ERR=432)ARLEAS 
	j1 = INDEX(VERSION(j0+1:),'.')
	if (j1 .eq. 0)	then
	   AEDIT = 0
	else
	   read(VERSION(j0+j1+1:j0+j1+1),*,ERR=432)AEDIT 
	endif
 432	continue
C	write(*,'(1I1,1A1,1I1,1A1,1I1)')AVERS,'.',ARLEAS,'.',AEDIT
C----------------------------------------------------------------------|
C Read file const.inc
	FILENA='for/const.inc'
	call	OPENRD(1,FILENA,0,IERR)
	if (IERR .gt. 0)	goto 900
C Read main variables list: Character*6 array PRNAME(1:NPRNAM) is filled
	JINC=0
	J=0
 1	J=J+JINC
	if(J .gt. NCONST)	then
	   write(*,*)'>>> READAT: File const.inc error: ',
     >			'            more than',NCONST,' variables'
	   goto 2
	endif
	read(1,100,ERR=924)VNAM
 100	format(2X,1A6)
	if(VNAM .eq. 'Variab') JINC=1
	if(J .le. 0)		goto 1
	if(VNAM .eq. 'End va')	goto 2
C Find ordinal numbers of some variables on the "for/const.inc" list
	if(VNAM .eq. 'AB    ')	KAB=J
	if(VNAM .eq. 'ABC   ')	KABC=J
	if(VNAM .eq. 'AWALL ')	KAWALL=J
	if(VNAM .eq. 'RTOR  ')	KRTOR=J
	if(VNAM .eq. 'ELONM ')	KELONM=J
	if(VNAM .eq. 'TRICH ')	KTRICH=J
	PRNAME(J)= VNAM
				goto 1
 2	NPRNAM=J-1
C Read main constant list: Character*6 array CFNAME(1:NCFNAM) is filled
	JINC=0
	J=0
 3	J=J+JINC
	if(J .gt. NCONST)	then
	   write(*,*)'>>> READAT: File const.inc error: ',
     >			'            more than',NCONST,' constants'
	   goto 4
	endif
	read(1,100,ERR=924)VNAM
	if(VNAM .eq. 'Consta') JINC=1
	if(J .le. 0)		goto 3
	if(VNAM .eq. 'End co')	goto 4
	CFNAME(J)=VNAM
				goto 3
 4	NCFNAM=J-1
C Read internal variable list: Character*6 array SRNAME(1:NSRNAM) is filled
	JINC=0
	J=0
	NINVAR=0
 13	J=J+JINC
	if (J .gt. NCONST)	then
	   write(*,*)'>>> READAT: File const.inc error: ',
     >			'            more than',NCONST,' parameters'
	   goto 14
	endif
	read(1,'(A)',ERR=924)STRI
	VNAM = STRI(3:8)
	if(VNAM .eq. 'Intern') JINC=1
	if(J .le. 0)		goto 13
	j1 = index(STRI,'- ASTRA function:')
	if (j1.ne.0 .and. NINVAR.eq.0) NINVAR = J-1
C	write(*,*)j1,VNAM,j
	if(VNAM .eq. 'End in')	goto 14
	SRNAME(J)=VNAM
				goto 13
 14	NSRNAM=J-1
	close(1)
C Next 9 lines provide an output which forms a body of 
C the subroutine APPST, file for/defarr.f
!	do	j=1,NPRNAM
!	   write(*,*)'	if (JFALT .eq.',j,')  YARR = ',PRNAME(j)
!	enddo
C Exception for DTEQ(4,20)
!	do	j=1,NINVAR
!	   write(*,*)'	if (JFALT .eq.',j,'+NCONST) YARR = ',SRNAME(j)
!	enddo
!	do	j=1,NCFNAM
!	   write(*,*)'	if (JFALT .eq.',j,'+2*NCONST) YARR = ',CFNAME(j)
!	enddo
C	write(*,*)NSRNAM,NINVAR
C	write(*,'(10(A7))')(SRNAME(j),j=1,NSRNAM)
	TIME = TSTART		! Here TSTART=0
C----------------------------------------------------------------------|
C	write(*,*)'File "for/const.inc" OK'
C----------------------------------------------------------------------|
C Read file status.inc
	FILENA='for/status.inc'
	call	OPENRD(1,FILENA,0,IERR)
	if (IERR .gt. 0)	goto 901
C Read global array list: Character*6 array ARNAME(1:NARNAM) is filled 
	JINC=0
	J=0
 11	J=J+JINC
	if(J .gt. 512)	then
	   write(*,*)'>>> READAT: File for/status.inc error: ',
     >			'            more than',512,' arrays'
	   goto 12
	endif
	read(1,100,ERR=924)VNAM
	if(VNAM .eq. 'Array ') JINC=1
	if(J .le. 0)		goto 11
	if(VNAM .eq. 'End ar')	goto 12
	ARNAME(J)= VNAM
				goto 11
 12	NARNAM=J-1
	close(1)
C	write(*,*)NARNAM
C	write(*,'(10(A7))')(ARNAME(j),j=1,NARNAM)
C	write(*,*)'File "for/status.inc" OK'
C Next 4 lines provide an output which forms a body of 
C the subroutine APPSR, file for/defarr.f
!	do	j=1,NARNAM
!	   write(*,*)'     if (NUMALR .eq.',j,
!     >		') write(JCH,*)(',ARNAME(j),"(j),j=1,JA1)"	
!	enddo
C----------------------------------------------------------------------|
C	write(*,*)'File "for/status.inc" OK'
C----------------------------------------------------------------------|
C Read file 'tmp/astra.log'
	FILENA='tmp/astra.log'
	XFILE=.FALSE.
	call	OPENRD(1,FILENA,0,IERR)
	if (IERR .gt. 0)	goto 902
 103	read(1,fmt='(1A79)',ERR=902,END=104)STRJ
	COMNAM=STRJ(1:6)
C TASK, TSTART=TIME, TPAUSE, TEND are defined by the 2nd reading
C User's home directory: $HOME
	if (COMNAM .eq. 'WHOME:')	then
	   j1 = length(STRJ(7:))
	   if (j1 .eq. 0)	goto	942
	   WHOME=STRJ(7:)
	endif
C AROOT:
	if (COMNAM .eq. 'AROOT:')	then
	   j1 = length(STRJ(7:))
	   if (j1 .eq. 0)	goto	942
	   AROOT=STRJ(7:)
	endif
C AWD:
	if (COMNAM .eq. 'AWD:  ')	then
	   j1 = length(STRJ(7:))
	   if (j1 .eq. 0)	goto	942
	   AWD=STRJ(7:)
	endif
C VAR:   (Variant name)
	if (COMNAM .eq. 'VAR:  ')	then
	   j1 = length(STRJ(7:))
	   RDNAME=STRJ(7:)
	   if (length(RDNAME) .eq. 0)	goto	942
	endif
C MOD:   (Model name)
	if (COMNAM .eq. 'MOD:  ')	then
	   EQNAME=STRJ(7:)
	   if (length(EQNAME) .eq. 0)	goto	942
	endif
C MACH:  (Machine name - presently not used)
	if (COMNAM .eq. 'MACH: ')	then
	   j1 = length(STRJ(7:))
	   if (j1 .ne. 0)	XLINE1=STRJ(7:)
	endif
C SHOT:  (Shot No. - presently not used)
	if (COMNAM .eq. 'SHOT: ')	then
	   j1 = length(STRJ(7:))
C	   if (j1 .ne. 0)	XLINE2=STRJ(7:)
	endif
C XFILE: (Ex-file name)
	if (COMNAM .eq. 'XFILE:' .and. length(STRJ(7:)).gt.0)	then
	   j1 = length(STRJ(7:))
C	   write(*,*)j1,'  "',STRJ(7:7),'"'
	   if (j1 .ne. 0)	XFILE=.TRUE.
C	   if (XFILE)		write(*,*)'Ex-file requested'
	endif
C RUNN:  (JAMS run No.)
	if (COMNAM .eq. 'RUNN: ' .and. length(STRJ(7:)).gt.0)	then
	   j1 = length(STRJ(7:))
	   j2 = length(WHOME)
	   STRI = WHOME(1:j2)//'/astra/runs/run'//STRJ(7:6+j1)//'/'
	   read(STRJ(7:6+j1),*)JRUNN
	   FILEX = STRI(1:j1+16+j2)
C	   write(*,*)'Exchange directory: "',FILEX(1:15+j1+j2),'"'
	endif
C PROFT: Profile file name
	if (COMNAM.eq.'RSNAM:')	FILENA=STRJ(7:)
	if (COMNAM.eq.'PROFT:')	RSNAME='.res/'//STRJ(7:)
C RTYPE: Task type & Pause mode setting time
	if (COMNAM(1:6) .eq. 'RTYPE:') then
	   TASK=STRJ(7:9)//' '
	   call UPCASE(3,TASK)
	   STR40 = STRJ(1:40)
	   j = nextwd(STRJ(7:))	! 2nd word, TPAUSE, is present
C	   write(*,*)j,STRJ(7:)
	   YTP = -1.d20
	   if (j.ne.0) read(STRJ(6+j:),*)YTP
C	   write(*,*)"Wtime =",YTP
	endif
			goto 103
 104	close(1)
C Uncomment the next line to set the "DSP" mode at the start
C	TASK(1:3)='DSP' 	! STOP, STEP, WAIT
	j = length(RSNAME)
C	write(*,*)"RSNAME =",RSNAME(1:j)
C	write(*,*)"FILENA =",FILENA(1:)
	if (RSNAME(j:j).ne."/")	goto	105		! No name specified
	if (FILENA(1:5) .eq. "/tmp/")	then
	   j = length(FILENA)
	   RSNAME = FILENA(1:j)//char(0)
	else
	   j = length(RSNAME)
	   if(RSNAME(j:j).eq."/")RSNAME='.tsk'//RSNAME(5:j)//'profil.dat'
	endif
 105	continue
C	write(*,*)'Review file: "',RSNAME(1:length(RSNAME)),'"'
C----------------------------------------------------------------------|
C	goto	107
C  A subroutine APPRID forms RUNID, IDAY,IMONTH,IYEAR,IHOUR,IMINUT,ISEC
	IDAY = 0
	call	APPRID(IDAY,IMONTH,IYEAR,IHOUR,IMINUT,ISEC)
C	write(*,*)IDAY,IMONTH,IYEAR,IHOUR,IMINUT
C	write(*,*)'"',RUNID(1:lonlen(RUNID)),'"'
	j2 = length(RSNAME)
	j = 1
	jj = 1
 106	j1 = j
	j = jj
	jj = j+index(RSNAME(j+1:j2),'/')
	if (jj .ne. j)	goto	106
	STRI = "ASTRA @ "//RSNAME(j1+1:jj-1)//char(0)
	j = jj-j1+8
C       write(*,*)j,jj,j1,j2,'"',RSNAME(j1+1:jj-1),'"',STRI(1:j),'"'
C       write(*,*)'"',TASK,'"'
	if (TASK(1:3) .eq. 'BGD')	then
	   STRI = 'BGD'//char(0)
	   call	initvm(XWX,XWY,XWW,XWH,COLTAB,STRI(1:3),3)
	   goto	107
	endif
	jj = max(0,(15+NTOUT-64)/16)
	XWH = XWH+2*jj*(DYLET+2)
	call	initvm(XWX,XWY,XWW,XWH,COLTAB,STRI(1:j),j)
	IM  = 1
	NST = 0
	MOD10 = 1
	call	SETFRM(IM,XSC0,XSC,NX,NY)
	call	AFRAME(IM,XSC0,XSC,NX,NY)
C	call	AFRAME
	j = XOUT+.49
	call	ASRUMN(j)	! Task menu
	call	textbf(0,XWH-104,RUNID,80) ! Task ID
C       CHORDN = LINEAV(1)
C       call	UPSTR(CHORDN,1./MU(NA))
Configfile exists?
C       if (CNFILE(1:1) .ne. '*')	call	createpixmap(1)
C----------------------------------------------------------------------|
C Read file equ/MODEL.log
 107	jj = length(EQNAME)
	if (jj .eq. 0)	goto	904
	j1 = 0
	do	j3=1,jj
	   if (EQNAME(j3:j3) .eq. '/')	j1 = j3
	enddo
	if (j1 .eq. 1 .or. j1 .ge. jj)	goto	904
	if (j1 .eq. 0)	then
	   inquire(file='equ/'//EQNAME(1:jj)//'.log',exist=EXILOG)
	   if ( EXILOG ) 	then
	      write(*,*)'>>> Warning: File "',
     &		'equ/'//EQNAME(1:jj)//'.log" has been found'
	      write(*,*)'             Most probably it should be ',
     &		'moved to "','equ/log/'//EQNAME(1:jj),'"'
	   endif
	   FILENA = 'equ/log/'//EQNAME(1:jj)//char(0)
	else
	   FILENA = 'equ/'//EQNAME(1:j1)//'log/'//EQNAME(j1+1:jj)
	endif
C	write(*,*)'"',FILENA(1:length(FILENA)),'"'
C	write(*,*)'"',FILENA(1:length(FILENA))//'.log','"'
	inquire(file=FILENA,exist=EXILOG)
	if (.not. EXILOG) 	then
	   inquire(file=FILENA(1:length(FILENA))//'.log',exist=EXILOG)
	   if ( EXILOG ) 	then
	      write(*,*)'>>> Warning: File "',
     &		FILENA(1:length(FILENA))//'.log," has been found'
	      write(*,*)'             Most probably it should be ',
     &		'moved to "',FILENA(1:length(FILENA)),'"'
	   endif
	   goto	33
	endif
	call	OPENRD(1,FILENA,0,IERR)
	if (IERR .gt. 0)	goto 914
	read(1,132,ERR=913,END=913)STRI
C What version?
	j1 = index(STRI,'version')
	if (j1 .ne. 0)	goto	109
	read(1,132,ERR=913,END=913)STRI
	j1 = nextwd(STRI(1:))
	j3 = 92
	if(j1.lt.132 .and. j1.gt.0) read(STRI(j1:),*,ERR=108)j3
 108	if(j3.eq.92.or.j3.eq.60)	then
C version 5.0 and older
	   write(*,102)char(7)				! Beep
 102	   format(1X,A,'>>> Warning! Backward compatibility with ',
     &		'old type "equ/*.log" files'/'             is not ',
     &		'possible. Please readjust tuning parameters.')
	   read(1,*,ERR=913,END=913)(DEVAR(jy),jy=1,22),
     &		(DEVAR(jy+2),jy=23,40),(DEVAR(jy+4),jy=41,j3-4),
     &		(YY,jy=1,4)
	else
C version 5.1: NPRNAM=80; NCFNAM=88
	   write(*,102)char(7)				! Beep
	   read(1,*,ERR=913,END=913)		!(DEVAR(jy),jy=1,j3)
     1		AB,    ABC,   YY,    AIM1,  AIM2,  AIM3,  AMJ,   AWALL,
     2		BTOR,  YY,    YY,    YY,    YY,    YY,    ELONG, ELONM,
     3		ENCL,  ENWM,  FECR,  FFW,   FICR,  FLH,   GN2E,  GN2I,
     4		YY,    IPL,   LEXT,  NNCL,  NNWM,  QNBI,  QECR,  QFW,
     5		QICR,  QLH,   YY,    YY,    RTOR,  SHIFT, TRIAN, TRICH,
     6		UEXT,  UPDWN, YY,    YY,    WNE,   WTE,   WTI,   ZMJ,
     7		(DEVAR(jy),jy=49,j3)
	endif
	read(1,132,ERR=913,END=913)STRI
	j1 = nextwd(STRI(1:))
	if (j1 .lt. 132)	read(STRI(j1:),*)j3
	read(1,*,ERR=913,END=913)			!(CONSTF(jy),jy=1,j3)
     1		CF1,   CF2,   CF3,   CF4,   CF5,   CF6,   CF7,   CF8,
     2		CF9,   CF10,  CF11,  CF12,  CF13,  CF14,  CF15,  CF16,
     3		CV1,   CV2,   CV3,   CV4,   CV5,   CV6,   CV7,   CV8,
     4		CV9,   CV10,  CV11,  CV12,  CV13,  CV14,  CV15,  CV16,
     5		CBND1, CBND2, CBND3, CBND4, CNB1,  CNB2,  CNBI3, CNBI1,
     6		(YY,jy=1,12),               CNB4,  CNBI2, CNB3,  CNBI4,
     8		CNEUT1,CNEUT2,CNEUT3,CNEUT4,CCD1,  CCD2,  CCD3,  CCD4,
     9		CHE1,  CHE2,  CHE3,  CHE4,  CIMP1, CIMP2, CIMP3, CIMP4,
     1		CRAD1, CRAD2, CRAD3, CRAD4, CMHD1, CMHD2, CMHD3, CMHD4,
     2		CPEL1, CPEL2, CPEL3, CPEL4, CFUS1, CFUS2, CFUS3, CFUS4
	read(1,132,ERR=913,END=913)STRI
	j1 = nsymb(STRI,":")
	if ( j1 .eq. 0 )	then
	   j1 = 16
	else
	   read(STRI(j1+1:),*)j1
	endif
	read(1,*,ERR=913,END=229)(DELOUT(jy),jy=1,j1)
	read(1,132,ERR=913,END=913)STRI
	goto	228
 109	continue
C version 5.3 (namelist-like format)
C	write(*,*)"New format detected"
	read(1,132,ERR=913,END=913)STRI
	do	222	j1=1,NCONST+1
	   read(1,132,ERR=913,END=913)STRI
	   if (STRI(1:10) .eq. ' Constants')	goto	223
	   do   221	jy=1,NPRNAM
	      if (PRNAME(jy) .eq. 'ZRD1  ')	goto	222
	      if (PRNAME(jy) .ne. STRI(1:6))	goto	221
	      read(STRI(9:),*)DEVAR(jy)
C	      write(*,*)PRNAME(jy),'=',DEVAR(jy)
 221	   continue
 222	continue
 223	continue
	do	j1=1,NCONST+1
	   read(1,132,ERR=913,END=913)STRI
C	   if (STRI(1:10) .eq. ' Time cont')	goto	227
	   if (STRI(1:10) .eq. ' Control p')	goto	225
	   do   224	jy=1,NCFNAM
	      if (CFNAME(jy) .ne. STRI(1:6))	goto	224
	      read(STRI(9:),*)CONSTF(jy)
C             write(*,*)CFNAME(jy),'=',CONSTF(jy)
 224	   continue
	enddo
 225	continue
	j2 = nsymb(STRI,":")
	if ( j2 .eq. 0 )	goto	913
	read(STRI(j2+1:),*,err=913)j2
	do	j1=1,NCONST+1
	   read(1,132,ERR=913,END=913)STRI
	   if (STRI(1:10) .eq. ' Color tab')	goto	228
	   do   226	jy=1,NSRNAM
	      if (SRNAME(jy) .ne. STRI(1:6))	goto	226
C	      write(*,*)SRNAME(jy),'=',DELOUT(jy),jy		! Before
	      read(STRI(9:),*)DELOUT(jy)
C	      write(*,*)SRNAME(jy),'=',DELOUT(jy),jy,nsrnam	! After
 226	   continue
	enddo
C 227	continue
C	j1 = nsymb(STRI,":")
C	if ( j1 .eq. 0 )	then
C	   j1 = 16
C	else
C	   read(STRI(j1+1:),*)j1
C	endif
C	read(1,*,ERR=913,END=229)(DELOUT(jy),jy=1,j1)
C	read(1,132,ERR=913,END=229)STRI

 228	continue
	NB1 = DELOUT(13)
	NUF = DELOUT(14)
C	NBND = DELOUT(19)
	XFLAG = DELOUT(20)
	j1 = nsymb(STRI,":")
	if ( j1 .eq. 0 )	then
	   j1 = 24
	else
	   read(STRI(j1+1:),*)j1
	endif
	read(1,*,ERR=913,END=229)(COLTAB(jy),jy=1,j1)
 229	close(1)
C	write(*,*)(DEVAR(jy),jy=1,80)
C	write(*,*)(CONSTF(jy),jy=1,88)
C	write(*,*)(COLTAB(jy),jy=1,j1,2)
C	write(*,*)'File "',FILENA(1:length(FILENA)),'" OK'
	goto 36
 33	continue
C 14-02-97
C Read file equ/MODEL.var
	FILENA='equ/'//EQNAME(1:jj)//'.var'
	inquire(FILE=FILENA,EXIST=EXILOG)
	if(EXILOG) 	then
	   write(*,*)char(7)				! Beep
	   write(*,*)">>> Warning! Backward compatibility with ",
     &	   	'old type "equ/*.var" files'
	   write(*,*)"             is not possible. Please readjust",
     &			' tuning parameters.'
C Old standard
C		call	OPENRD(1,FILENA,1,IERR)
C		if(IERR.gt.0)	goto 914
C	read(1,ERR=35)DEVAR,CONSTF,(DELOUT(jy),jy=1,12),(YY,jy=1,80)
C 35		close(1)
C		goto 36
			endif
C Read file equ/MODEL.cns
C	FILENA='equ/'//EQNAME(1:jj)//'.cns'
C	inquire(FILE=FILENA,EXIST=EXI)
C	if(EXI) 	then
C Oldest standard
C	jj	= length(FILENA)
C	INT4 = 0
C	FILENA(jj+1:)=char(INT4)
C		call	OPENRD(1,FILENA,0,IERR)
C		if(IERR.gt.0)	goto 914
C		do	J1=1,90,6
C		read(1,'(6F13.6)',ERR=35)(CONSTF(J),J=J1,J1+5)
C		enddo
C 35		close(1)
C			endif

 36	continue
C	write(*,*)'File "for/MODEL.log" OK'
C----------------------------------------------------------------------|
C Start time signal reading
C----------------------------------------------------------------------|
C	write(*,*)XFILE,JRUNN
	IVAR=0
	if (JRUNN .le. 0)	goto	368
	if ( .not. XFILE)	goto	368
C Read Ex-file according to the convertion table from '~/astra/runs/ext2a'
C	write(*,*)'JAMS run ID: ',JRUNN
C	write(*,*)FILEX(1:index(FILEX,'runs')+4)//'ext2a'
	call	
     & OPENRD(1,FILEX(1:index(FILEX,'runs')+4)//'ext2a'//char(0),0,IERR)
	if (IERR .ne. 0)	goto	9411
	jj = length(FILEX)
	call	OPENRD(2,FILEX(1:jj)//'astra.ext'//char(0),0,IERR)
	if (IERR .ne. 0)	then
	   write(*,*)'>>> READAT: Warning. File "',FILEX(1:jj)//
     &	       'astra.ext"'
	   write(*,*)'    (time traces) is not found. Continue anyway.'
	   goto	368
C	   call	a_stop
	endif
C	write(*,*)'Reading Ext-file: "',FILEX(1:jj),'astra.ext"'
 361	read(1,fmt='(1A79)',ERR=940,END=367)RUNAST
	jj = ichar(RUNAST(1:1))
	if (jj .lt. 65 .or.  jj .gt. 122)	goto	361	! If non-
	if (jj .gt. 91 .and. jj .lt. 97 )	goto	361	! letter
	j3 = killbl(RUNAST(1:79),79)
	call	UPCASE(j3,RUNAST(1:j3))
	j0 = index(RUNAST(1:j3),'=')
	j1 = index(RUNAST(1:j3),'<')
	j2 = max(j0,j1)
	if (j2 .eq. 0)			goto	361
	jj  = min(j0,j1)
	if (jj .eq. 0)	jj = j2
	VNAM = RUNAST(1:jj-1)
	if (VNAM(1:3) .eq. 'END')	goto	367
	do	j2=1,NPRNAM
	   if (PRNAME(j2) .ne. VNAM)	goto	362
	   jtim1 = j2
	   goto	363
 362	   continue
	   if (j2 .eq. NPRNAM)	then
	    write(*,*)'>>> Warning >>> Error in file ~/astra/runs/ext2a'
	    write(*,*)'    "',VNAM,'" is not recognized'
	   endif
	enddo
 363	continue
C	write(*,*)'VNAM="',VNAM,'",   PRNAME="',PRNAME(jtim1),'",   No.',jtim1
	j1 = index(RUNAST(1:j3),':')
	j2 = 0
	if (j1 .eq. 0)	then
	   j1 = j3
	   factor = 1.
	   goto	364
	else
	   j2 = index(RUNAST(j1+1:j3),':')
	endif
	if (j2 .gt. 1)	then
	   read(RUNAST(j1+1:j1+j2-1),*,ERR=932)factor
	endif
 364	continue
	NAMREQ(1) = RUNAST(jj+1:j1-1)//'       '
	write(*,'(5A,1p,1E12.3)')'"',VNAM,'" is taken from "'
     >			,RUNAST(jj+1:j1-1),'"',factor
C	   call	jexbase
C	   call	a_stop
C Input:
C    2      logical unit for file .../astra.ext
C    6      logical unit for message
C    MSIGEX,MEXT dimensions of JEXTSG(MSIGEX,MEXT)
C                              JEXSIG(MSIGEX,MEXT)
C    MSIG   number of signals required (currently = 1)
C    NAMREQ names  of signals required (character*8)
C Output:
C    JTEX   number of time points for each trace
C    IERSIG error (MSIGEX) code for each signal
C	    = 0	Signal found OK
C	    = 1	Signal not in Ex-file
C    JEXTSG signal (MSIGEX,MEXT)
C    JEXTIM signal times (MSIGEX,MEXT)
C    VERJEX the Ex-file version
C    JSHOT  the shot number
C    IERR   error code  0 OK
C			1 - error reading version, shot
C			2 - error reading data
C			3 - MTEX too small
C			4 - MPEX too small
C----------------------------------------------------------------------|
	JTEX = -1
	call	JEXFILET(2,6,MSIGEX,MEXT,MSIG,NAMREQ,
     >		JTEX,IERSIG,JEXTSG,JEXTIM,VERJEX,JSHOT,IERR)
	if (IERR      .ne. 0)	goto	9412
	if (IERSIG(1) .eq. 0)	goto	365
	if (RUNAST(jj:jj) .eq. '=')	then
	   write(*,*)'>>> READAT: Ext-file reading error'
	   write(*,'(A,$)')'             Signal "'
	   write(*,*)NAMREQ,'" not found'
	   call	a_stop
	endif
	if (RUNAST(jj:jj) .eq. '<')	then
	   write(*,'(A)')' >>> READAT: Ext-file reading error'
	   write(*,'(A,$)')'             Signal '
	   write(*,*)'"',NAMREQ,'" not found, continue anyway'
	   goto	361
	endif
 365	continue				! Reading OK
	jtim = JTEX
	if (jtim .le. 0)	then
	   IFDFVX(jtim1)=-1
	   write(*,*)'>>> READAT: Error reading Ex-file. Signal "'
     &				,NAMREQ(1),'" ignored.'
	   goto	361
	elseif (jtim .eq. 1)	then
	   IFDFVX(jtim1)=0
	else
	   IFDFVX(jtim1)=1
	endif
	if (IVAR+jtim .gt. NTVAR)	then
	   write(*,*)'>>> READAT: Too many time points >',NTVAR
	   call	a_stop
	endif
	DEVAR(jtim1) = factor*JEXTSG(1,1)
	do	j = 1,jtim
	   IVAR = IVAR+1
	   INDVAR(IVAR) = jtim1
	   VARDAT(1,IVAR) = JEXTIM(1,j)
	   VARDAT(2,IVAR) = factor*JEXTSG(1,j)
	   VARDAT(3,IVAR)=0.
	enddo
C	write(*,*)(JEXTIM(1,j),j=1,jtim)
C	write(*,*)(JEXTSG(1,j),j=1,jtim)
C	write(*,*)'VNAM="',VNAM,'",   PRNAME="',PRNAME(jtim1),
C     &'",   No.',jtim1
C	write(*,*)NAMREQ(1),MEXT,jtim
C	write(*,*)(VARDAT(1,IVAR-jtim+j),j=1,jtim)
C	write(*,*)(VARDAT(2,IVAR-jtim+j),j=1,jtim)
	goto	361
 367	close(1)
	close(2)
C	goto	20				! Skip variable reading
C End reading Ext-file
C----------------------------------------------------------------------|
 368	continue
C Read experimental data:
	jj	=length(RDNAME)
	FILENA='exp/'//RDNAME(1:jj)//char(0)
	call	OPENRD(2,FILENA,0,IERR)
	if (IERR .gt. 0)	goto 903
C Read XLINE1 and append with spaces if shorter than 132
	read(2,132,ERR=905)XLINE1
	j1=1
 37	j2=noblan(XLINE1(j1:))
	if (j2 .eq. 0)		goto	38
	j1=j1+j2-1
	j1=j1+length(XLINE1(j1:))
C j1 is position_of_the_word_end + 1
	if (j1 .lt. 132)	goto	37
 38	continue
	if (j1 .lt. 132) write(XLINE1(j1:),'(132A1)')(" ",j=j1,132)
	read(2,132,ERR=905)XLINE2
 132	format(1A132)

C----------------------------------------------------------------------|
C Read simple variable loop (between the labels "5" and "10"):
C
	VNAMO=' '
	VNAMU=' '
 5	continue
	read(2,132,END=39)STRI
C	write(*,*)STRI
	STRJ = STRI
	call	UPCASE(132,STRJ)
	if (STRJ(1:6) .eq. 'PROFIL')	goto	20
	if (STRJ(1:3) .eq. 'END')	goto	20
C NAMECL:/'POINTS','NAMEXP','GRIDTY','NTIMES','FILTER','FACTOR','END   '/
	if (STRJ(1:6) .eq. NAMECL(1) .or.
     +	    STRJ(1:6) .eq. NAMECL(3) .or.
     +	    STRJ(1:6) .eq. NAMECL(5) .or.
     +	    STRJ(1:6) .eq. NAMECL(7) )	goto	20
C Variable name, time, value, error
	VNAM = VARNAM(STRI(1:6),jt)
	call	UPCASE(6,VNAM)
	VTIM = STRI(9:14)
	VDAT = STRI(17:22)
	VERR = STRI(25:30)
C A special tretment required for NB1,TSTART,TEND
	if (VNAM .eq. 'NB1   ')	then
	   if (jt .ne. 0 )		goto	927
	   if (VNAMO.eq.'NB1   ')	goto	928
	   call	READF6(VDAT,VRDATA,IERR)
	   NB1 = VRDATA
	   VNAMO = VNAM
	   goto	5
	endif
	if (VNAM .eq. 'TSTART')	then
	   if (jt .ne. 0 )		goto	927
	   if (VNAMO.eq.'TSTART')	goto	928
	   call	READF6(VDAT,TSTART,IERR)
	   TIME  = TSTART
	   VNAMO = VNAM
	   goto	5
	endif
	if (VNAM .eq. 'TEND  ')	then
	   if (jt .ne. 0 )		goto	927
	   if (VNAMO.eq.'TEND  ')	goto	928
	   call	READF6(VDAT,VRDATA,IERR)
	   TEND = VRDATA
	   VNAMO = VNAM
	   goto	5
	endif
	jtim = 0
	jtim1 = 0
	factor = 1.
	j1 = 1
	j2 = 1
	do	211	j3=1,10
	   if (j1 .gt. 1)	j2 = nextwd(STRI(j1:))
	   j1 = j1-1+j2
	   if (j2 .eq. 0)	goto	212
	   if (j1 .ge. 132)	goto	212
	   do	j=1,6
	      KEYWRD = STRI(j1:j1+5)
	      call	UPCASE(6,KEYWRD)
	      if (KEYWRD.eq.NAMECL(j))	then
		 j1 = j1+6
		 goto	(20,213,20,214,20,216)	j
	      endif
	   enddo
	   goto		211
 213	   j1 = j1-1+noblan(STRI(j1:))
	   VNAM = VARNAM(STRI(j1:j1+5),j)
	   call		UPCASE(6,VNAM)
	   jtim1 = 1
	   goto		211
 214	   read(STRI(j1:),*,ERR=932)jtim
	   jtim1 = 1
	   goto		211
 216	   read(STRI(j1:),*,ERR=932)factor
	   jtim1 = 1
 211	continue
 212	continue
C	write(*,*)VNAM,VNAMU,VNAMO,VTIM,VDAT,VERR,jtim,factor

C Special treatment required for BND (unimplemented)
	if (VNAM .eq. 'BND   ')	then
	   j = INDEX(STRJ,NAMECL(1))
	   if (j .ne. 0)	goto	20
	   write(*,*)'>>> Control line "',STRI(1:20),'... " error'
	   write(*,*)'    "BND" definition ignored,',
     &		     '  Please, move the group to array input'
	   VNAM = ' '
	   goto	5
	endif

	do	10 J=1,NPRNAM
	if (VNAM .ne. PRNAME(J))	goto 10
C	write(*,*)'"',VNAM,'"',j,'  "',PRNAME(J),'"',VTIM,VDAT,VERR
	if (IFDFVX(j) .ge. 0 .and. XFILE)	then
	   if (VNAM .eq. VNAMO)	goto	9
	   write(*,*)'>>> READAT: Warning. ',
     &		'Signal "',PRNAME(J),'" has been defined by Ex-file.'
	   write(*,*)'            ',
     &		'Repeated definition by the Astra file is ignored.'
	   goto	9
	endif

C U-file name duplicated:
	if(VNAM .eq. VNAMU)	goto 923

C Exp_file (start reading)
C IFDFVX: = 0 - determined by the data file (time independent),
C	  = 1 - determined by the data file (time dependent),
C	  = 2 - determined by the MODEL,
C	  = 3 - keyboard
C	  = 4 - any changes are forbidden (AB, RTOR, ELONM, TRICH)

	j1 = index(STRI,':')
	if (j1 .ne. 0)	goto	8		! Pointer to Ex/U-file
	if (jtim1 .gt. 0)	then
	    if (IVAR+jtim .gt. NTVAR)	goto	909
	    if (jtim .le. 1)	then
		jtim = 1
		IFDFVX(J)=0
	    else
		IFDFVX(J)=1
	    endif
C			 Read "time" array & function array
	    read(2,*,ERR=906)(VARDAT(1,IVAR+jj),jj=1,jtim)
	    read(2,*,ERR=906)(VARDAT(2,IVAR+jj),jj=1,jtim)
	    DEVAR(J) = factor*VARDAT(2,IVAR+1)
	    do	jj = 1,jtim
		IVAR = IVAR+1
		INDVAR(IVAR) = J
		VARDAT(2,IVAR) = factor*VARDAT(2,IVAR)
		VARDAT(3,IVAR)=0.
	    enddo
	    goto	9
	endif

C Tabulation encountered in the old-standard line:
	if (jt .ne. 0)	goto	927

	call READF6(VTIM,TIMEVR,IERR)
	if (IERR .ne. 0)	goto	913
	call READF6(VDAT,VRDATA,IERR)
	if (IERR .ne. 0)	goto	913
	call READF6(VERR,VRERR,IERR)
	if (IERR .ne. 0)	goto	913
	IFDFVX(J) = 0
	DEVAR(J) = factor*VRDATA
	if(VNAM .eq. VNAMO)	IFDFVX(J) = 1
C Repeated name
	IVAR = IVAR+1
	if (IVAR .gt. NTVAR)		goto 909
	INDVAR(IVAR) = J
	VARDAT(1,IVAR) = TIMEVR
	VARDAT(2,IVAR) = factor*VRDATA
	VARDAT(3,IVAR) = VRERR
	goto	9
C Old-standard exp_file (end reading)

 8	continue
C ":" found in the input string "STRI"
	if(VNAM .eq. VNAMO)	goto	923
	call	ANLSTR(STRI,FILENA,lname,STR40,j0,factor,ierr)
C	if (j0 .ne. 0)write(*,*)'"',STR40(1:j0),'"',j0
C	write(*,*)'"',FILENA(1:lname),'"',factor,ierr
	if (ierr .eq. 1)	goto	930
	if (ierr .eq. 2)	goto	931
	if (ierr .eq. 3)	goto	908
	if (ierr .eq. -1)	then
	   write(*,*)STRI
	   write(*,*)FILENA
	   goto	20		! Ex-file pointer
	endif
	inquire(FILE=FILENA(1:lname),EXIST=EXI)
	if(.not.EXI)	goto	915
C Start reading 1D U-file
	VNAMX = VNAM
	do	j3=1,6
	    if (VNAMX(j3:j3) .eq. ' ')	then
		VNAMX(j3:j3)='X'
		goto	84
	    endif
	enddo
 84	continue
C	write(*,*)VNAMX,' is taken from the U-file "',
C     +		FILENA(1:lname),'".   Factor:',factor
	VNAMU = VNAM
	call	OPENRD(1,FILENA,0,IERR)
C unit 1 will be closed in UF1RD()
	if (IERR .gt. 0)	goto 916
	if (j0 .eq. 0)	STR40 = ' '
C	if (j0 .ne. 0)	write(*,*)'Internal name: "',STR40(1:j0),'"',j0
	ERCODE = 0
	NCH = 1
	j0 = min(6,j0)
	if (j0 .eq. 0)	j0 = 1
	call	UF1RD(NCH,J,DEVAR,FACTOR,ISHOT,IDEV,STR40(1:j0),ERCODE)
	goto (909,921,925,929,926),ERCODE
 9	VNAMO=VNAM
C	write(*,*)'Var: ',VNAM,'    IVAR =',IVAR,J,'  ',PRNAME(J)
 10	continue
	goto	5
C End reading variable loop (labels 5 <-> 10):
C----------------------------------------------------------------------|
C Start Astra standard assignments
 20	continue
	if (NB1.gt.NRD)	then
	write(*,*)'>>> FATAL ERROR: The radial grid size out of range.'
	write(*,*)'                 Parameter "NB1" cannot exceed',NRD
	call	a_stop
	endif
C 2nd reading file 'tmp/astra.log': override TPAUSE & TEND if given
C				    in a command line
C	FILENA='tmp/astra.log'
	call	OPENRD(1,'tmp/astra.log',0,IERR)
 434	read(1,fmt='(1A79)',ERR=902,END=433)RUNAST
	COMNAM=RUNAST(1:6)
C Start time
	if (COMNAM.eq.'TIME: ')	then
	   VTIM = RUNAST(7:)
	   if (length(VTIM) .ne. 0)	then
	      call READF6(VTIM,TSTART,IERR)
	      if (IERR .ne. 0)write(*,*)
     >		'>>> Error in file "tmp/astra.log", line TIME:'
	   endif
	endif
C Stop time
	if (COMNAM.eq.'TEND: '.and. length(RUNAST(7:)).gt.0)
     >	      read(RUNAST(7:),*)TEND
C Task type & Pause mode setting time
	if (COMNAM(1:6) .eq. 'RTYPE:') then
	   TASK=RUNAST(7:9)//' '
	   call UPCASE(3,TASK)
	   STR40 = RUNAST(1:40)
	   j = nextwd(RUNAST(7:))	! 2nd word is present, Pause time is set
C	   write(*,*)j,RUNAST(7:)
	   YTP = -1.d20
	   if (j.ne.0) read(RUNAST(6+j:),*)YTP
C	   write(*,*)"Wtime =",YTP
	endif
			goto 434
 433	close(1)
C----------------------------------------------------------------------|
C Uncomment the next line to set the "WAIT" mode at the start
C	TASK='WAIT'
	TIME = TSTART
C	j = KILLBL(STR40,40)
C	if (j.ge.10) read(STR40(10:j),*,ERR=434)TPAUSE
	if (YTP .gt. -9.d19)	TPAUSE = YTP
C	write(*,*)"Tpause =",TPAUSE

C	call	markloc("INTVAR"//char(0))
C	call	add2loc("Calling INTVAR from READAT"//char(0))
	call	INTVAR
C	call	markloc("READAT"//char(0))
	if (.not. EXILOG)				then
C Determine ABC & AB if not defined by "exp" file
	  if (IFDFVX(KABC).lt.0 .and. IFDFVX(KAB).lt.0)	then
	    write(*,*)
     +		'>>> The minor radius AB is not defined in "exp/'//
     +		RDNAME(1:length(RDNAME))//'"'
	    write(*,*)'  This can cause equilibrium convergence problem'
	  endif
	  if (IFDFVX(KAB).gt.0)
     +		write(*,*)'>>> Warning: AB cannot depend on time'
	  if (IFDFVX(KABC).ge.0 .and. IFDFVX(KAB).lt.0)   AB = ABC
	  if (IFDFVX(KABC).lt.0 .and. IFDFVX(KAB).ge.0)	then
	    ABC = AB
	    IFDFVX(KABC) = 0
	  endif
C Determine AWALL if not defined by an "exp" file
	  if(IFDFVX(KAWALL) .lt. 0)	AWALL=AB
	endif
	if (AWALL .lt. AB)	then
	    if (IFDFVX(KAWALL) .ge. 0)	write(*,*)
     >	    '>>> Warning: AWALL < AB.  Setting AWALL = AB'
	    AWALL=AB
	endif
	if (IFDFVX(KAWALL).lt.0 .and.
     +	    (AWALL.lt.AB .or. AWALL.gt.1.2*AB))	AWALL = AB

	if (ABC .gt. AB)	then
	    write(*,*)
     +		'>>> Warning: ABC cannot exceed AB. Setting ABC = AB'
	    if (IFDFVX(KABC) .lt. 1)	ABC = AB
	endif

	if (AWALL .gt. 1.2*AB .or. AWALL .gt. RTOR)	then
	    write(*,*)
     +		'>>> Warning: AWALL is set unreasonably large'
	    write(*,*)
     +		'    Check settings in data and log files'
	endif
 	IFDFVX(KAB)=4
 	IFDFVX(KELONM)=4
 	IFDFVX(KRTOR)=4
 	IFDFVX(KTRICH)=4
	IFDFVX(KAWALL)=4
C	VOLUME	=GP2*GP*(RTOR+SHIFT-.25*ABC*TRIAN)*AB*AB*ELONG
	VOLUME	=GP2*GP*RTOR*AB*AB*ELONG
C	ROC = ABC/(RTOR+SHIFT)*
C     *		SQRT(ELONG*RTOR*(RTOR+SHIFT-.25*ABC*TRIAN))
	ROC	= ROC3A(RTOR,SHIFT,ABC,ELONG,TRIAN)
C	write(*,*)RTOR,SHIFT,ELONG,TRIAN
C	write(*,*)ABC,ROC3A(1.d0,SHIFT,ABC,ELONG,TRIAN)
	if (AWALL .gt. RTOR+SHIFT)	then
	   ROWALL = AWALL*SQRT(max(ELONM,ELONG))
	elseif (ELONM .gt. ELONG)	then
	   ROWALL = ROC3A(RTOR,SHIFT,AWALL,ELONM,TRICH)
	else
	   ROWALL = ROC3A(RTOR,0.d0,AWALL,ELONG,TRIAN)
	endif
C	write(*,*)ROWALL,ROC3A(RTOR,SHIFT,AWALL,ELONG,TRIAN)
	HRO	= 1.1*ROWALL/(NB1-0.5)
	HROA    =HRO
C	HRO	=1.5*ROWALL/(NB1-0.5)
	do	j=1,NB1
	    RHO(J) = (J-0.5)*HRO
	    G22(J) = RHO(J)
	    VR(J)  = (GP2*(RTOR+SHIFT))**2*RHO(J)/RTOR
	    VRS(J) = (GP2*(RTOR+SHIFT))**2*J*HRO/RTOR
	    G11(J) = VRS(J)
	    if (RHO(j) .lt. ROC)	NA = j
	enddo
	call	SETGEO(0)
	VOLUM(1) = HRO*VR(1)
	do	J=2,NA
	   VOLUM(J) = VOLUM(J-1)+HRO*VR(J)
	enddo
	NA1 = NA+1
	NA1N = NA1
	NA1E = NA1
	NA1I = NA1
	NA10 = NA1
	NA11 = NA1
	NA12 = NA1
	NA13 = NA1
	NA14 = NA1
	NA15 = NA1
	NA16 = NA1
	NA17 = NA1
	NA18 = NA1
	NA19 = NA1
	VOLUM(NA1) = VOLUM(NA)+(HROA-0.5*HRO)*VR(NA)
	VOLUME = VOLUM(NA1)

	NAB = NA1
	if (NA1 .lt. NB1 .and. AB .gt. ABC)	then
	   do	j=NA1+1,NB1
		if (AMETR(j) .lt. AB)   NAB = j
	   enddo
	   if (NAB .lt. NB1)	NAB=NAB+1
	endif
	ROB = RHO(NAB)
	RHO(NA1) = ROC
	AMETR(NAB) = AB
	AMETR(NA1) = ABC
	call	SETGEO(NA1)

C	write(*,*)NB1,AB,ABC,AWALL,NA1,JRUNN
C	call	a_stop
C	write(*,'(1P,6E13.5)')(AMETR(j3),j3=NA-4,NB1)
C	write(*,*)ABC,NA1,AMETR(NA-1),AMETR(NA),AMETR(NA1),ROC
C	write(*,'(1P,6E13.5)')(AMETR(j3),j3=1,NA1)
C	write(*,'(1P,6E13.5)')(AMETR(j3),j3=NA1,NB1)
C----------------------------------------------------------------------|
C Input information is taken from JAMES through Ex-files
C	write(*,*)XFILE,JRUNN
	NGR = 0
	jarr = 0
	if (JRUNN .le. 0)	goto	110
	if ( .not. XFILE)	goto	110
C Read Ex-file according to the convertion table from '~/astra/runs/ex2a'
C       RDNAME = 'JAMS'			! Indication of mode
C Warning: unit 2 is opened! exp/file is connected
C This will be changed when reading .ext-file is enabled
	j0 = 0
	j1 = 0
	j2 = 0
	j3 = 0
	call
     >	OPENRD(1,FILEX(1:index(FILEX,'runs')+4)//'ex2a'//char(0),0,IERR)
	if (IERR .ne. 0)	goto	941
	jj = length(FILEX)
	call	OPENRD(3,FILEX(1:jj)//'astra.ex'//char(0),0,IERR)
	if (IERR .ne. 0)	call	a_stop
C	write(*,*)'Reading Ex-file: "',FILEX(1:jj),'astra.ex"'
 351	read(1,fmt='(1A79)',ERR=940,END=357)RUNAST
	jj = ichar(RUNAST(1:1))
	if (jj .lt. 65 .or.  jj .gt. 122)	goto	351	! If non-
	if (jj .gt. 91 .and. jj .lt. 97 )	goto	351	! letter
	j3 = killbl(RUNAST(1:79),79)
	call	UPCASE(j3,RUNAST(1:j3))
	j0 = index(RUNAST(1:j3),'=')
	j1 = index(RUNAST(1:j3),'<')
	j2 = max(j0,j1)
	if (j2 .eq. 0)			goto	351
	jj  = min(j0,j1)
	if (jj .eq. 0)	jj = j2
	VNAM = RUNAST(1:jj-1)				! Signal to read
	if (VNAM(1:3) .eq. 'END')	goto	357
	do	j2=1,NARRX
	   if (EXARNM(j2) .ne. VNAM)	goto	352
	   jexar = j2
	   goto	353
 352	   continue
	   if (j2 .eq. NARRX)	then
	     write(*,*)'>>> Warning >>> Error in file ~/astra/runs/ex2a'
	     write(*,*)'    "',VNAM,'" is not recognized as X-vector'
	   endif
	enddo
 353	continue
C	write(*,*)'VNAM="',VNAM,'",   EXARNM="',EXARNM(j2),'",   No.',j2
	j1 = index(RUNAST(1:j3),':')
	j2 = 0
	if (j1 .eq. 0)	then
	   j1 = j3
	   factor = 1.
	   goto	354
	else
	   j2 = index(RUNAST(j1+1:j3),':')
	endif
	if (j2 .gt. 1)	then
	   read(RUNAST(j1+1:j1+j2-1),*,ERR=932)factor
	else
	   factor = 1.
	endif
 354	continue
	NAMREQ(1) = RUNAST(jj+1:j1-1)//'       '
C	write(*,'(5A,1p,1E12.3)')'"',VNAM,'" is taken from "'
C     >			,RUNAST(jj+1:j1-1),'"',factor
	if (NGR .ne. 0)	goto	355
C First, read the radial grid
	INTYPE = 12
	call	JEXFILE(3,6,MPEX,MSIGEX,MTEX,MSIG,'XRHO    ',
     >		JPEX,JTEX,IERSIG,JEXRAD,TIMJEX,VERJEX,JSHOT,IERR)
C	write(*,*)"XRHO"
C	write(*,*)(JEXRAD(j1,1,1),j1=1,JPEX)
	if (IERR      .ne. 0)	goto	9413
	if (IERSIG(1) .eq. 0)	goto	355
C Presently is supposed that the signal "XRHO" (normalized rho toroidal)
C   is always present in Exfile.
C A possibility is forseen to treat the radial grid as other signals.
C Then setting "<" as the first delimiter will result in equidistant grid
C while "=" at the same place will require definition through ex-file.
C Remove lines with "C<" to enable this option.
	write(*,*)'>>> READAT: Ex-file reading error:'
	write(*,'(A,$)')'             Radial grid not found.'
	write(*,*)' Equidistant grid is used.'
	INTYPE = 2
 355	continue
C	if (NGR .eq. 0)	then
C	   write(*,*)"XRHO",JPEX
C	   do	j=1,JTEX
C	      write(*,*)"t = ",TIMJEX(j)
C	      write(*,'(1p,5E12.3)')(JEXRAD(j1,1,j),j1=1,JPEX)
C	   enddo
C	endif
C Input:
C    3      logical unit for file FILENA
C    6      logical unit for message
C    MPEX,MSIGEX,MTEX dimensions of JEXSIG(MPEX,MSIGEX,MTEX)
C    MSIG   number of data signals required (currently = 1)
C    NAMREQ names of data signals required (character*8)
C Output:
C    JPEX   actual length of profile <= MPEX
C    JTEX   number of times returned
C    IERSIG(MSIGEX) error code for each signal
C	0   normal exit
C	1   error: signal not in Ex-file
C    JEXSIG profile signal
C    TIMJEX profile times
C    JSHOT  the shot number
C    IERR   error code  0 OK
C			1 - error reading version, shot
C			2 - error reading data
C			3 - MTEX too small
C			4 - MPEX too small
C    VERJEX the Ex-file version
C----------------------------------------------------------------------|
	call	JEXFILE(3,6,MPEX,MSIGEX,MTEX,MSIG,NAMREQ,
     >		JPEX,JTEX,IERSIG,JEXSIG,TIMJEX,VERJEX,JSHOT,IERR)
	if (IERR      .ne. 0)	goto	9413
	if (IERSIG(1) .eq. 0)	goto	356
	if (RUNAST(jj:jj) .eq. '=')	then
	   write(*,*)'>>> READAT: Ex-file reading error'
	   write(*,'(A,$)')'             Signal "'
	   write(*,*)NAMREQ,'" not found'
	   call	a_stop
	endif
	if (RUNAST(jj:jj) .eq. '<')	then
	   write(*,'(A)')' >>> READAT: Ex-file reading error'
	   write(*,'(A,$)')'             Signal '
	   write(*,*)'"',NAMREQ,'" not found, continue anyway'
	   goto	351
	endif
 356	continue
C	if (VNAM .eq. 'PIX   ')	then
C	write(*,'(6A,1I4,1A,1I3,1A,1p,1E12.3)')'"',VNAM,'" has been '
C     >			,'taken from "',NAMREQ(1),'"',JPEX,' x',JTEX
C     >			,'      Factor =',factor
C	do	j=1,JTEX
C	   write(*,*)"t = ",TIMJEX(j)
C	   write(*,'(1p,5E12.3)')(JEXSIG(j1,1,j),j1=1,JPEX)
C	   write(*,'(1p,5E12.3)')(factor*JEXSIG(j1,1,j),j1=1,JPEX)
C	   write(*,*)(factor*JEXSIG(j1,1,j),j1=1,JPEX)
C	enddo
C	endif
	jbdry = JPEX
	jtim  = JTEX
	jtim1 = max(JTEX,1)
	if(NGR+jtim1 .gt. NTARR)		goto 910
	if(jarr+jbdry+1 .gt. NRDX*NTARR)	goto 939
	KOGDA(jexar) = NGR+1
	IFDFAX(jexar) = 0

	do	j=1,jtim1
	   NGR = NGR+1
	   KTO(NGR) = jexar
	   TIMEX(NGR) = TIMJEX(j)
	   NGRIDX(NGR) = jbdry
	   FILTER(NGR) = 0.0001		! Adjust !
C The group of commands marked as "C2" below serves for setting
C     an equidistant grid
C2	   NTYPEX(NGR) = 2
C2	   GDEX(NGR) = jarr+1
C2	   GDEY(NGR) = jarr+1
C2	   do	jj=1,jbdry
C2	      DATARR(jarr+jj)=factor*JEXSIG(jj,1,j)
C2	   enddo
C2	   jarr = jarr+jbdry
	   NTYPEX(NGR) = INTYPE
	   if (j .eq. 1)	then
C These GDEX and GDEY definitions work for 10 <= INTYPE <= 16
C Define radial grid 
	      GDEX(NGR) = jarr+1
	      do	jj=1,jbdry
		 DATARR(jarr+jj)=JEXRAD(jj,1,j)
	      enddo
	      jarr = jarr+jbdry
	      XBDRY = DATARR(jarr)
C Define data values
	      GDEY(NGR) = jarr+1
	      do	jj=1,jbdry
		 DATARR(jarr+jj)=factor*JEXSIG(jj,1,j)
	      enddo
	      jarr = jarr+jbdry
	   else
	      GDEX(NGR) = GDEX(NGR-1)
	      GDEY(NGR) = jarr+1
	      do	jj=1,jbdry
		 DATARR(jarr+jj)=factor*JEXSIG(jj,1,j)
	      enddo
	      jarr = jarr+jbdry
	   endif
	enddo
C       write(*,*)(DATARR(jarr+jj),jj=1,jbdry)
	goto	351
 357	close(1)
	close(3)
C	goto	40				! Skip profile reading
C End reading Ex-file
C----------------------------------------------------------------------|
C Start reading radial profiles:		Radial grid: jbdry 
C 			Option for different input grids can be added
 110	continue
C	NGR=0
	if (STRI(1:3) .eq. 'END')	goto	39
C	jarr=0
	VNAM=' '
	VNAMO=' '
 111	continue
C	backspace(2)
C	read(2,132)STRI
	INTYPE = -1
	TIMEVR = .0
	jtim = 0
	jbdry = 0
	j1 = 1
	j2 = 1
	factor = 1.
	ALFA = 0.001
	if (STRJ(1:6).eq.'PROFIL') goto 119
	if (STRI(1:6).eq.NAMECL(1) .and. nextwd(STRI(7:)).eq.0)	then
C          the string consists of "POINTS Number" only:
	   read(STRI(7:),*,ERR=932)NPNT
	   jbdry = NPNT
	   goto		119
	endif

C parse floating format control string:	! 12=6x2
	do	118	jj=1,12
	   if (j1 .gt. 1)	j2 = nextwd(STRI(j1:))
	   j1 = j1-1+j2
	   if (j2 .eq. 0)	goto	119
	   if (j1 .ge. 132)	goto	119
	   do	j=1,6
	      if (STRI(j1:j1+5) .eq. NAMECL(j))	then
		 j3 = nextwd(STRI(j1:))
		 if (j3 .eq. 0)  goto	937
		 j1 = j1+j3-1
		 goto (112,113,114,115,116,117)	j
	      endif
	   enddo
	   write(*,*)'>>> READAT error unknown key word in string:'
	   write(*,*)STRI(1:length(STRI))
 112	   read(STRI(j1:),*,ERR=932)NPNT
	   jbdry = NPNT
	   goto		118
 113	   j1 = j1-1+noblan(STRI(j1:))
	   VNAM = ARRNAM(STRI(j1:j1+5))
	   call		UPCASE(6,VNAM)
	   goto		118
 114	   read(STRI(j1:),*,ERR=932)INTYPE
	   goto		118
 115	   read(STRI(j1:),*,ERR=932)jtim
	   goto		118
 116	   read(STRI(j1:),*,ERR=932)ALFA
	   goto		118
 117	   read(STRI(j1:),*,ERR=932)factor
 118	continue
 119	continue

	if (VNAM .eq. 'BNDX  ')	then
	   if (NBNT .ne. 0)			goto	934
	   j = INDEX(STRJ,NAMECL(1))
	   if (j .ne. 0) read(STRI(j+6:),*) NBND
	   if (j .eq. 0) 			goto	935
	   NBNT = max(jtim,1)
	   if ((2*NBND+1)*NBNT .gt. 25000) 	goto	936
	   read(2,*,ERR=938)(BNDARR(j),j=1,(2*NBND+1)*NBNT)
	   VNAMO = VNAM
	   VNAM = ' '
	   goto		24
	endif

C	write(*,*)'Quantity: "',VNAM,'",   Grid =',jbdry, 
C     >	',   Type =',INTYPE,',   Slices =',jtim
C	write(*,'(2(a,f10.3)')
C     >	'          Factor =',factor,',     Filter =',ALFA
 24	continue
	if (VNAM .eq. ' ')	then
	   STRI=' '
	   read(2,132,ERR=906,END=39)STRI
	   STRJ = STRI
	   call	  UPCASE(132,STRJ)
C	   write(*,*)"String read:",STRI(1:20),'"'
	   if (STRI(1:6) .eq. NAMECL(1) .or.
     +	       STRI(1:6) .eq. NAMECL(2) .or.
     +	       STRI(1:6) .eq. NAMECL(3) .or.
     +	       STRI(1:6) .eq. NAMECL(4) .or.
     +	       STRI(1:6) .eq. NAMECL(5) .or.
     +	       STRI(1:6) .eq. NAMECL(6) )	goto 111 ! New input group
	   VNAM = ARRNAM(STRI(1:6))
	   call	UPCASE(6,VNAM)
C	   if(VNAM .ne. VNAMO .and. VNAMO(1:1) .ne. ' ') goto 111
	   if(VNAM .eq. 'ENDX  ')		goto 39
	   INTYPE = -1
	   jbdry = NPNT
	endif

C	write(*,*)'"',STRI(1:6),'",      "',VNAM,'"',INTYPE
Cycle over all arrays and all types of input data
      do   30   jexar=1,NARRX
	 if (EXARNM(jexar) .ne. VNAM)	goto 30
	 if (IFDFAX(jexar) .ge. 0 .and. XFILE)	then
C	    write(*,*)'"',VNAMO,'"',VNAM,'"',EXARNM(jexar),'"',XFILE
	    if (VNAM .eq. VNAMO)	goto	29
C	    call system('tput smso')		! enter standout mode
C	    call system('tput rmso')		! exit standout mode
C	    call system('tput smul')		! enter underline mode
C	    call system('tput rmul')		! exit underline mode
C	    call system('tput bold')
	    write(*,'(1a,$)')' >>> READAT: Warning. Signal "'
	    write(*,'(1a,$)')EXARNM(jexar)
	    write(*,'(1a)')'" has been defined by Ex-file.'
C	    call system('tput rmso')
	    write(*,*)'            ',
     &		'Repeated definition by the Astra file is ignored.'
	    goto	29
	endif
C 							New quantity
C	 write(*,*)'"',VNAMO,'"',VNAM,'"',EXARNM(jexar),'"',XFILE
C	 write(*,*)jexar,IFDFAX(jexar),INTYPE,jbdry,factor
	 if(VNAM .ne. VNAMO)	then
	    if (IFDFAX(jexar) .ge. 0)	goto	933
	    KOGDA(jexar) = NGR+1
	 endif
C --------------------- Check for U-file command string
	 do   j1=1,132
	    if (STRI(j1:j1) .eq. ':')	then
	       STRJ = STRI
	       j = KILLBL(STRJ,j1)
	       if (j .gt. 6)	then
		  call UPCASE(7,STRJ(j-6:j))
		  if (STRJ(j-6:j) .eq. 'U-FILE:')	goto	25
		  if (STRJ(j-6:j) .eq. 'J-FILE:')	goto	25
	       endif
	       if (j .gt. 7)	then
		  call UPCASE(8,STRJ(j-7:j))
		  if (STRJ(j-7:j) .eq. 'EX-FILE:')	goto	25
	       endif
	    endif
	 enddo
	 if (jbdry .le. 1)	goto	912
	 if (jbdry .gt. NRDX)	goto	911

	 jtype = 1+INTYPE/10
	 jtim1 = max(jtim,1)
	 if(NGR+jtim1 .gt. NTARR)		goto 910
	 if(jarr+jtype*jbdry+1 .gt. NRDX*NTARR)	goto 939
C----------------------------------------------------------- INTYPE = -1
C 					Reading old-standard data
	if (INTYPE .eq. -1)	then
	   NGR = NGR+1
	   read(STRI(7:11),*)TIMEX(NGR)
	   GDEX(NGR) = jarr+1
	   GDEY(NGR) = jarr+1
	   KTO(NGR) = jexar
	   NGRIDX(NGR) = jbdry
	   FILTER(NGR) = ALFA
	   NTYPEX(NGR) = XINPUT+.49
	   do	jj=1,jbdry
	      j3 = 12+5*(jj-1)
	      read(STRI(j3:j3+4),*,ERR=913)DATARR(jarr+jj)
	   enddo
C	   write(*,*)(DATARR(jarr+jj),jj=1,jbdry)
	   jarr = jarr+jbdry
	elseif
C------------------------------------------------------ INTYPE = 0 -> 20
     +	 (  0 .le. INTYPE .and. INTYPE .le. 20)	then
	   if (jtim .gt. 0)	read(2,*,ERR=932)(TIMEX(NGR+j),j=1,jtim)
	   do	j=1,jtim1
	      NGR = NGR+1
	      KTO(NGR) = jexar
	      NGRIDX(NGR) = jbdry
	      NTYPEX(NGR) = INTYPE
	      FILTER(NGR) = ALFA
	      if (j .eq. 1)	then
		 GDEX(NGR) = jarr+1
		 if (INTYPE .eq. 18 .or. INTYPE .eq. 19)	then
		    jarr = jarr+1
		    read(2,*,ERR=932)DATARR(jarr)
		 endif
		 do	j1=1,jtype
		    read(2,*,ERR=907)(DATARR(jarr+jj),jj=1,jbdry)
		    jarr = jarr+jbdry
		    if (jtype.eq.2 .and. INTYPE.le.17 .and. j1.eq.1)
     +					XBDRY = DATARR(jarr)
		 enddo
		 GDEY(NGR) = jarr-jbdry+1
		 do	jj=1,jbdry
		    DATARR(jarr-jbdry+jj)=factor*DATARR(jarr-jbdry+jj)
		 enddo
	      else
		 GDEX(NGR) = GDEX(NGR-1)
		 GDEY(NGR) = jarr+1
		 read(2,*,ERR=907)(DATARR(jarr+jj),jj=1,jbdry)
		 do	jj=1,jbdry
		    DATARR(jarr+jj)=factor*DATARR(jarr+jj)
		 enddo
		 jarr = jarr+jbdry
	      endif
	   enddo
C          if ((INTYPE.ge.12 .and. INTYPE.le.15)
C     +		 .and. abs(XBDRY-1.).gt.1.e-6)	then
C		write(*,*)">>> ERROR: Inconsistency in the input data ",
C     +			" GRIDTYPE =",INTYPE
C		write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'": ',
C     >	   	'    Input radial grid is expected to be normalized'
C		call	a_stop
C	    endif
	   if (INTYPE.eq.17) write(*,*)
     +	     'X-axis analysis is not implemented GRIDTYPE =',INTYPE
	else
C-------------------------------------------------------- INTYPE unknown
	   write(*,*)'>>> ERROR: Unknown input type =',
     >		INTYPE,',  ignored'
	   call	a_stop
	endif
C----------------------------------------------------------------------|
	goto	29

C^^^^^^^^^^^^^^^^^^ Old-standard exp_file (end reading) ^^^^^^^^^^^^^^^
C----------------------------------------------------------------------|
 25	continue
C	if (IFDFAX(jexar) .ge. 0)	goto	923
	call	ANLSTR(STRI,FILENA,lname,STR40,j0,factor,ierr)
	if (ierr .eq. 3)	goto	908
	if (ierr .eq. 1)	goto	930
	if (ierr .eq. 2)	goto	931
	TIMEVR = 0.
	inquire(FILE=FILENA(1:lname),EXIST=EXI)
	if (.not.EXI)		goto	915
	if (ierr .eq. 0)	goto	251		! -> U-file
C Start reading Ex-file
C	write(*,*)STRI
C	write(*,*)'"',FILENA(1:lname+1),'"',factor,lname,ierr
	j = nsymb(STRI,"<")
	if (j .eq. 0)		goto	919
	j1 = noblan(STRI(j+1:))
C	write(*,*)j,j1,j+j1,'"',STRI(j+j1:j+j1+1),'"'

C	j2 = nextwd(STRI(j+j1:))
C	if (j2 .eq. 0)	NAMREQ(1) = STRI(j+j1:)//'       '
C	if (j2 .ne. 0)	NAMREQ(1) = STRI(j+j1:j+j1+j2-2)//'       '
	j3 = length(STRI(j+j1:))
	NAMREQ(1) = STRI(j+j1:j+j1+j3-1)//'       '
	call	UPCASE(8,NAMREQ(1))
C	write(*,*)'"',NAMREQ,'"'
	call	OPENRD(1,FILENA,0,IERR)
C	open(UNIT=1,FILE=FILENA,STATUS='OLD')
	if (ierr .ne. 0)	goto	918
C Ex-file requisites:
C	integer		MPEX,MSIGEX,MTEX,JSHOT,JPEX,MSIG,JTEX
C	parameter	(MPEX=51, MSIGEX=1, MTEX=50, MSIG=1)
C	integer		IERSIG(MSIGEX)
C	double precision JEXSIG(MPEX,MSIGEX,MTEX),TIMJEX(MTEX)
C	character	NAMREQ(MSIG)*8,VERJEX*20
C----------------------------------------------------------------------|
C emulation JEXFILE:
C	JPEX = 12
C	JTEX = 10
C	do	j=1,JTEX
C	   TIMJEX(j) = 0.1*j
C	   do	jj=1,JPEX
C	      JEXSIG(jj,1,j) = 1.-10*j*(jj-1.)/(JPEX-1.)**2
C	   enddo
C	enddo
	call	JEXFILE(1,6,MPEX,MSIGEX,MTEX,MSIG,NAMREQ,
     >		JPEX,JTEX,IERSIG,JEXSIG,TIMJEX,VERJEX,JSHOT,IERR)
C Input:
C    1      logical unit for file FILENA
C    6      logical unit for message
C    MPEX,MSIGEX,MTEX dimensions of JEXSIG(MPEX,MSIGEX,MTEX)
C    MSIG   number of data signals required (currently = 1)
C    NAMREQ names of data signals required (character*8)
C Output:
C    JPEX   actual length of profile <= MPEX
C    JTEX   number of time returned
C    IERSIG(MSIGEX) error code for each signal
C	0   normal exit
C	1   error: signal not in Ex-file
C    JEXSIG profile signal
C    TIMJEX profile times
C    JSHOT  the shot number
C    IERR   error code  0 OK
C			1 - error reading version, shot
C			2 - error reading data
C			3 - MTEX too small
C			4 - MPEX too small
C    VERJEX the Ex-file version
	close(1)
C----------------------------------------------------------------------|
C      write(*,*)VNAM,' is taken from the Ex-file "',FILENA(1:lname),'"'
C      write(*,*)'       Factor:',factor,',    JETTO_NAME  "',NAMREQ,'"'

	if (IERR .ne. 0)			goto 917
	if (IERSIG(MSIG) .ne. 0)	then
	   write(*,*)'>>> READAT: Ex-file "',FILENA(1:lname),
     >		'" reading error'
	   write(*,'(A,$)')'            '
	   write(*,*)' Signal "',NAMREQ,'" not found, continue anyway'
	endif
	jbdry = JPEX
	jtim  = JTEX
	jtim1 = max(JTEX,1)
	if(NGR+jtim1 .gt. NTARR)		goto 910
	if(jarr+jbdry+1 .gt. NRDX*NTARR)	goto 939

	do	j=1,jtim1
	   NGR = NGR+1
	   TIMEX(NGR) = TIMJEX(j)
	   KTO(NGR) = jexar
	   NGRIDX(NGR) = jbdry
	   NTYPEX(NGR) = 2		! RHO grid: \rho_j=(j-1)*h
	   FILTER(NGR) = 0.0001		! Adjust !
	   GDEX(NGR) = jarr+1
	   GDEY(NGR) = jarr+1
	   do	jj=1,jbdry
	      DATARR(jarr+jj)=factor*JEXSIG(jj,MSIG,j)
	   enddo
	   jarr = jarr+jbdry
	enddo
C       write(*,*)(DATARR(jarr+jj),jj=1,jbdry)
	goto	29
C End reading Ex-file
C----------------------------------------------------------------------|
 251	continue
C Start reading U-file
C	write(*,*)VNAM,' is taken from the U-file "',
C     +		FILENA(1:lname),'".   Factor:',factor
	call	OPENRD(1,FILENA,0,IERR)
	if (IERR .gt. 0)	goto 916
C Shot #, Device ID, U-file dimension,
	read(1,132,ERR=925)STRI
	read(STRI(2:6),  fmt='(1I5)',ERR=925)ISHOT
	read(STRI(8:11), fmt='(1A4)',ERR=925)IDEV
	read(STRI(13:13),fmt='(1I1)',ERR=925)jdim
C	read(STRI(15:15),fmt='(1I1)',ERR=925)
	if(jdim.le.0 .or. jdim.gt.2)	goto	921
C Shot date
	read(1,132,ERR=925)STRI
C Number of associated scalar quantities
	TIMEVR = 0
	read(1,*,ERR=925)j3
	if (j3 .gt. 0)	then
	   do	j=1,j3
	      read(1,*,ERR=925)YY
	      read(1,132,ERR=925)STRI
	      jj = KILLBL(STRI,80)
	      call UPCASE(8,STRI)
	      if (STRI(1:4).eq.'TIME')	jj = 5
	      if (STRI(5:5).eq.':')	jj = jj+1
	      if (STRI(jj:jj) .eq. 'S')	TIMEVR = YY
	   enddo
	endif
C	write(*,*)jdim,NGR,jarr
	if (jdim .eq. 2)	goto	256
C Continue reading 1D U-file
C Independent variable label
	    read(1,132,ERR=925)STRJ
C Dependent variable label
	    read(1,132,ERR=925)STRI
C Processing code
	    read(1,132,ERR=925)STRI
	    read(1,*,ERR=925)jrad
	    if (jrad .gt. NRDX)			goto	922
C Read "radial" coordinate array
	    if (jarr+jrad .gt. NRDX*NTARR)	goto	939
	    if (NGR+1 .gt. NTARR)		goto	910
	    TIMEX(NGR+1) = TIMEVR
	    if (IFREE)	
     &		read(1,*,ERR=253)(DATARR(jarr+j),j=1,jrad)
	    if (.not.IFREE)
     &		read(1,'(1X,6F13.1)',ERR=926)(DATARR(jarr+j),j=1,jrad)
	    call    CHECKU(INTYPE,ABC,AB,XBDRY,
     &			DATARR(jarr+1),jrad,jbdry,STRJ,FILENA)
	    if (INTYPE .eq. 18 .or. INTYPE .eq. 19)	then
		jarr = jarr+1
		do	j=jarr+jbdry,jarr,-1
		    DATARR(j+1) = DATARR(j)
		enddo
		if (INTYPE .eq. 18)	DATARR(jarr) = RTOR
		if (INTYPE .eq. 19)	DATARR(jarr) = 0.
	    endif
C	    write(*,*)jarr,jbdry
C	    write(*,*)(DATARR(jarr+j),j=1,jbdry)
C	    call	a_stop
	    jarr = jarr+jbdry
C Read function, first in a free format, then in a fixed.
	    if(jarr+jbdry .gt. NRDX*NTARR)	goto	939
	    if (.not.IFREE)	goto	254
	    read(1,*,ERR=252)(DATARR(jarr+j),j=1,jbdry)
	    goto	255
 252	    if (INTYPE .eq. 18 .or. INTYPE .eq. 19)	jarr = jarr-1
	    jarr = jarr-jbdry
 253	    close(1)
	    IFREE = .FALSE.
	    goto	251
 254	    read(1,'(1X,6F13.1)',ERR=926)(DATARR(jarr+j),j=1,jbdry)
 255	    IFREE = .TRUE.
C	    write(*,*)jarr,jbdry
C	    write(*,*)(DATARR(jarr+j),j=1,jbdry)
C End reading 1D U-file
	    jtim = 1
	    goto	266
 256	continue
C Continue reading 2D U-file
C 1st independent variable label: X-
	   read(1,132,ERR=925)STRI
	   STRAD2 = STRI(1:32)
	   jj = KILLBL(STRI,80)
	   call UPCASE(8,STRI)
C 2nd independent variable label: Y-
	   read(1,132,ERR=925)STRJ
	   STRAD = STRJ(1:32)
	   jj = KILLBL(STRJ,80)
	   call UPCASE(8,STRJ)
	   if (STRI(1:4).ne.'TIME' .and. STRJ(1:4).ne.'TIME')	then
		write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'" ',
     >		'error:  unrecognized 2D U-file type. Input ignored'
		close(1)
		goto 29
	   endif
	   TCONTG = .FALSE.
	   j = 5
	   if (STRI(1:4).eq.'TIME')	TCONTG = .TRUE.
	   if (TCONTG .and. STRI(1:5).eq.'TIME:')	j = 6
	   if (.not.TCONTG .and. STRJ(1:5).eq.'TIME:')	j = 6
C	   if(TCONTG)	  write(*,*)'                      2D U-file "'
C     +		,FILENA(1:length(FILENA)),'" is time contiguous'
C	   if(.not.TCONTG)write(*,*)'                      2D U-file "'
C     +		,FILENA(1:length(FILENA)),'" is radially contiguous'
	   if (TCONTG .and. (STRI(j:j)   .ne. 'S')
     +		      .and. (STRI(j:j+2) .ne. 'SEC')
     +		      .and. (STRI(j:j+5) .ne. 'SECOND')
     +		      .and. (STRI(j:j+6) .ne. 'SECONDS'))
     +	   write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'"',
     +				' time unit should be second'
	   if (.not.TCONTG .and. (STRJ(j:j)   .ne. 'S')
     +			   .and. (STRJ(j:j+2) .ne. 'SEC')
     +			   .and. (STRJ(j:j+5) .ne. 'SECOND')
     +			   .and. (STRJ(j:j+6) .ne. 'SECONDS'))
     +	   write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),'"',
     +				' time unit should be second'
	   if (.not.TCONTG) STRAD = STRAD2

C Dependent variable label
	   read(1,132,ERR=925)STRI
C Processing code
	   read(1,132,ERR=925)STRI
	   read(1,*,ERR=925)j
	   read(1,*,ERR=925)jj
	   if (TCONTG)	then
		jrad = jj
		jtim = j
			else
		jrad = j
		jtim = jj
			endif

	   if (jrad .gt. NRDX)			goto	922
	   if (NGR+jtim+1 .gt. NTARR)		goto	910
	   if (jarr+jrad  .gt. NRDX*NTARR)	goto	939

	   if (.not.IFREE)	goto	262
	   if (TCONTG)	then
C Read "time" coordinate array then read "radial" coordinate array 
		read(1,*,ERR=261)(TIMEX(NGR+j),j=1,jtim)
		read(1,*,ERR=261)(DATARR(jarr+j),j=1,jrad)
	   else
C Read "radial" coordinate array then read "time" coordinate array
		read(1,*,ERR=261)(DATARR(jarr+j),j=1,jrad)
		read(1,*,ERR=261)(TIMEX(NGR+j),j=1,jtim)
	   endif
	   goto		263
 261	   close(1)
	   IFREE = .FALSE.
	   goto		251
 262	   if (TCONTG)	then
C Read "time" coordinate array then read "radial" coordinate array 
		read(1,'(1X,6F13.1)',ERR=926)(TIMEX(NGR+j),j=1,jtim)
		read(1,'(1X,6F13.1)',ERR=926)(DATARR(jarr+j),j=1,jrad)
	   else
C Read "radial" coordinate array then read "time" coordinate array
		read(1,'(1X,6F13.1)',ERR=926)(DATARR(jarr+j),j=1,jrad)
		read(1,'(1X,6F13.1)',ERR=926)(TIMEX(NGR+j),j=1,jtim)
	   endif
 263	   continue

C	write(*,*)jarr,jtim
C	write(*,'(1P,6E13.5)')(TIMEX(NGR+j),j=1,jtim)
C	write(*,*)jrad,INTYPE
C	write(*,'(1P,6E13.5)')(DATARR(jarr+j),j=1,jrad)
	   call	   CHECKU(INTYPE,ABC,AB,XBDRY,
     >			DATARR(jarr+1),jrad,jbdry,STRAD,FILENA)
	   if (INTYPE .eq. 18 .or. INTYPE .eq. 19)	then
	      jarr = jarr+1
	      do	j=jarr+jbdry,jarr,-1
		 DATARR(j+1) = DATARR(j)
	      enddo
	      if (INTYPE .eq. 18)	DATARR(jarr) = RTOR
	      if (INTYPE .eq. 19)	DATARR(jarr) = 0.
	   endif
C	   write(*,*)jarr,jbdry,INTYPE
C	   write(*,'(1P,6E13.5)')(DATARR(jarr+j),j=0,jbdry)
	   jarr = jarr+jbdry
C Read function
	   if (jarr+jrad+(jtim-1)*jbdry .gt. NRDX*NTARR)  goto	939

	   if (.not.IFREE)	goto	265
	   if   (TCONTG)	then
C	      write(*,*)' *:	Reading time contiguous function'
	      read(1,*,ERR=264)
     &            ((DATARR(jarr+j+jbdry*j0),j0=0,jtim-1),j=1,jbdry)
	   else
C 	      write(*,*)' *:	Reading radially contiguous function'
	      read(1,*,ERR=264)
     &            ((DATARR(jarr+j0+jbdry*j),j0=1,jrad),j=0,jtim-1)
	   endif
	   goto		266
 264	   continue
	   if (INTYPE .eq. 18 .or. INTYPE .eq. 19)	jarr = jarr-1
	   jarr = jarr-jbdry
	   close(1)
	   IFREE = .FALSE.
	   goto	251
 265	   if   (TCONTG)	then
C	      write(*,*)'F13.1:	Reading time contiguous function'
	      read(1,'(1X,6F13.1)',ERR=926)
     &            ((DATARR(jarr+j+jbdry*j0),j0=0,jtim-1),j=1,jbdry)
	   else
C 	      write(*,*)'F13.1:	Reading radially contiguous function'
	      read(1,'(1X,6F13.1)',ERR=926)
     &            ((DATARR(jarr+j0+jbdry*j),j0=1,jrad),j=0,jtim-1)
	   endif
C	   write(*,*)jbdry
C	   do	j0=0,jtim-1
C	       write(*,'(1P,6E13.5)')(DATARR(jarr+j+jbdry*j0),j=1,jbdry)
C	   enddo
C End reading 2D U-file
C----------------------------------------------------------------------|
Close U-file:
 266	close(1)
	IFREE = .TRUE.
	if (INTYPE .gt. 13 .and. INTYPE .ne. 19)
     >		write(*,*)"Unknown U-file type"
	YXB = DATARR(jarr)
	YXB1 = DATARR(jarr-1)
	do	j=1,jtim
	   NGR = NGR+1
	   if (j .eq. 1)	then
	      GDEX(NGR) = jarr-jbdry+1
	      GDEY(NGR) = jarr+1
	   else
	      GDEX(NGR) = GDEX(NGR-1)
	      GDEY(NGR) = jarr+1
	   endif
	   do	j0=1,jbdry
	      DATARR(jarr+j0) =DATARR(jarr+j0)*factor
	   enddo
	   jarr = jarr+jbdry
	   if (jbdry .ne. jrad)	then
	      YB = DATARR(jarr)
	      YB1 = DATARR(jarr-1)
	      YB = (YB*(XBDRY-YXB1)-YB1*(XBDRY-YXB))/(YXB-YXB1)
	      DATARR(jarr) = YB
	   endif
	   KTO(NGR) = jexar
	   NGRIDX(NGR) = jbdry
	   NTYPEX(NGR) = INTYPE
	   FILTER(NGR) = ALFA
C	    write(*,*)'End 2D:  NGR, jarr',NGR,jarr,KTO(NGR)
C	    write(*,'(1P,6E13.5)')(DATARR(jarr-jbdry+j0),j0=1,jbdry)
C 28	continue
	enddo
 29	VNAMO = VNAM
	VNAM = ' '
	IFDFAX(jexar) = 0
 30   continue
	VNAM = ' '
	goto 24
 39	continue
	close(2)
C End reading "exp" file
 40	continue
C	write(*,*)' X-array[s] in use'
	j3 = -1
	do  42	j=1,NARRX
	   if (ARXUSE(j) .eq. 0)	goto	42
	   j2 = -1
	   do   41  j1=1,NARRX		! Find X-array place in EXARNM
	      if (IFDFAX(j1) .lt. 0)	goto	41
	      if (ARNAME(ARXUSE(j)) .ne. EXARNM(j1))	goto	41
	      j2 = 0			! X-array No. "j1" has been defined
C	      write(*,*)"Defined",j1
 41	   continue
	   if (j2 .eq. 0)	goto	42
	   j3 = j3+1
	   j0 = min(length(ARNAME(ARXUSE(j))),6)
	   if (j3 .eq. 0)	write(*,*)
	   write(*,*)'>>> Warning >>> X-array used but not defined: "',
     >		ARNAME(ARXUSE(j))(1:j0),'"'
 42	continue
	if (IFDEFX("TEX").eq.0)	then
	   j = system("grep HEXP= ./tmp/*.tmp | grep TEX > /dev/null")
	   if (j .eq. 0) write(*,*)
     >		'>>> Warning >>> X-array "TEX" is used but not defined'
	endif
	if (IFDEFX("TIX").eq.0)	then
	   j = system("grep XEXP= ./tmp/*.tmp | grep TIX > /dev/null")
	   if (j .eq. 0) write(*,*)
     >		'>>> Warning >>> X-array "TIX" is used but not defined'
	   j = system("grep SVCXX ./tmp/*.tmp | grep TIX > /dev/null")
	   if (j .eq. 0) write(*,*)
     >		'>>> Warning >>> X-array "TIX" is used but not defined'
	endif

	DELOUT(13) = NB1
	DELOUT(14) = NUF
	DELOUT(19) = NBND
	DELOUT(20) = XFLAG
	TIMEQL = TIME-DTEQL-1.d-7

C Define TAU,TAUMIN,TAUMAX,TSCALE,DROUT,DTOUT,DPOUT
	j = IFKEY(258)

	do	j=1,NRW
	   GRAP(j) = AB
	enddo
	TIM7(1)	= TINIT
C	if (TIME .lt. TINIT)	TINIT = TSTART
	if (TIME .gt. TINIT+1.025*abs(TSCALE))	TINIT = TSTART
C	TSCALE	= -abs(TSCALE)
	TIM7(3)	= abs(TSCALE)/8.

	STRI = XLINE1
	call	UPCASE(len(STRI),STRI)

	if (TASK(1:3) .ne. 'BGD')	then
	   CHORDN = LINEAV(1)
	   call	UPSTR(CHORDN,1./MU(NA))
	endif
C Define configuration file name
	j1 = length(RDNAME)
	CNFILE = 'exp/cnf/'//RDNAME(1:j1)	! Check default name
	j1 = j1+8
	if (j1 .gt. 32)	then
	   write(*,*)'>>> Configuration file name is too long'
	   call	a_stop
	endif
	inquire(FILE=CNFILE(1:j1),EXIST=EXI)
	if ( EXI )	goto	47
C Check generic names: AUG, FEAT, JET
	j  = index(STRI,'AUG')

	if (j .ne. 0)	then
	   CNFILE = 'exp/cnf/aug'
	   j1  = 11
	endif
	j  = index(STRI,'FEAT')
	if (j .ne. 0)	then
	   CNFILE = 'exp/cnf/feat'
	   j1  = 12
	endif
	j  = index(STRI,'JET')
	if (j .ne. 0)	then
	   CNFILE = 'exp/cnf/jet'
	   j1  = 11
	endif
	inquire(FILE=CNFILE(1:j1),EXIST=EXI)
	if ( .not. EXI )  CNFILE = '***'
 47	continue

C Define NBI configuration file name
	j1 = length(RDNAME)
	NBFILE = 'exp/nbi/'//RDNAME(1:j1)		! Check default name
	j1 = j1+8
	if (j1 .gt. 32)	then
	   write(*,*)'>>> NB configuration file name is too long'
	   call	a_stop
	endif
	inquire(FILE=NBFILE(1:j1),EXIST=EXI)
	if ( EXI )	goto	48
C Check generic names: AUG, FEAT, JET
	j  = index(STRI,'AUG')+index(STRI,'ASD')
	if (j .ne. 0)	then
	   NBFILE = 'exp/nbi/aug'
	   j1  = 11
	   CNB1 = 8
	endif
	j  = index(STRI,'FEAT')
	if (j .ne. 0)	then
	   NBFILE = 'exp/nbi/feat'
	   j1  = 12
	endif
	j  = index(STRI,'JET')
	if (j .ne. 0)	then
	   NBFILE = 'exp/nbi/jet'
	   j1  = 11
	endif
	inquire(FILE=NBFILE(1:j1),EXIST=EXI)
	if ( .not. EXI )  NBFILE = '***'
 48	continue
C	write(*,*)'Input data:   "',RDNAME(1:length(RDNAME)),'"'
C	write(*,*)'Geometry file:  "',CNFILE(1:length(CNFILE)),'"'
C	write(*,*)'NBI parameters: "',NBFILE(1:length(NBFILE)),'"'
Configfile exists?
	if (TASK(1:3) .ne. 'BGD' .and. CNFILE(1:1) .ne. '*')
     >				call	createpixmap(1)

C Define MSE configuration file name
	j1 = length(RDNAME)
	MSFILE = 'exp/mse/'//RDNAME(1:j1)	! Check default name
	j1 = j1+8
	if (j1 .gt. 32)	then
	   write(*,*)'>>> MSE configuration file name is too long'
	   call	a_stop
	endif
	inquire(FILE=MSFILE(1:j1),EXIST=EXI)
	if ( EXI )	goto	481
C Check generic names: AUG, FEAT, JET
	j  = index(STRI,'AUG')
	if (j .ne. 0)	then
	   MSFILE = 'exp/mse/aug'
	   j1  = 11
	endif
	j  = index(STRI,'FEAT')
	if (j .ne. 0)	then
	   MSFILE = 'exp/mse/feat'
	   j1  = 12
	endif
	j  = index(STRI,'JET')
	if (j .ne. 0)	then
	   MSFILE = 'exp/mse/jet'
	   j1  = 11
	endif
	inquire(FILE=MSFILE(1:j1),EXIST=EXI)
	if ( .not. EXI )  MSFILE = '***'
 481	continue
C	write(*,*)'MSE geometry: "',MSFILE(1:length(MSFILE)),'"'

C Define ECR configuration file name
	j1 = length(RDNAME)
	ECFILE = 'exp/ecr/'//RDNAME(1:j1)	! Check default name
	j1 = j1+8
	if (j1 .gt. 32)	then
	   write(*,*)'>>> ECR configuration file name is too long'
	   call	a_stop
	endif
	inquire(FILE=ECFILE(1:j1),EXIST=EXI)
	if ( EXI )	goto	482
C Check generic names: AUG, FEAT, JET
	j  = index(STRI,'AUG')
	if (j .ne. 0)	then
	   ECFILE = 'exp/ecr/aug'
	   j1  = 11
	endif
	j  = index(STRI,'FEAT')
	if (j .ne. 0)	then
	   ECFILE = 'exp/ecr/feat'
	   j1  = 12
	endif
	j  = index(STRI,'JET')
	if (j .ne. 0)	then
	   ECFILE = 'exp/ecr/jet'
	   j1  = 11
	endif
	inquire(FILE=ECFILE(1:j1),EXIST=EXI)
	if ( .not. EXI )  ECFILE = '***'
 482	continue

	SCM = .1
	call	markloc("READAT exit"//char(0))
	do	49	J=1,20
	if((RTOR+AWALL)/SCM8(J).gt.5..or.AWALL*ELONM/SCM8(J).gt.1.4)
     >						goto 49
	SCM = SCM8(J)
C	write(*,*)CNFILE(1:1),SCM
C	if (CNFILE(1:1).ne.'*' .and. J.lt.20)	SCM = SCM8(J+1)
	return
 49	continue
	SCM=2.
	return

 900	write(*,*)'>>> READAT: File for/const.inc error'
	goto	999
 901	write(*,*)'>>> READAT: File for/status.inc error'
	goto	999
 902	write(*,*)'>>> READAT: File tmp/astra.log error'
	goto	999
 903	write(*,*)'>>> READAT: No such experimental variant "',
     .		RDNAME(1:length(RDNAME)),'"'
	goto	999
 904	write(*,*)'>>> READAT: Error in model file name: "',
     .		EQNAME(1:length(EQNAME)),'"'
	goto	999
 905	write(*,*)'>>> READAT: Data file "',
     .		RDNAME(1:length(RDNAME)),'" error in header'
	goto	999
 906	write(*,*)'>>> READAT: Data file "',
     .		RDNAME(1:length(RDNAME)),'".  Format error in group '
	write(*,*)STRI(1:132)
Clength(STRI))
	goto	999
 907	write(*,*)'>>> READAT: Data file "',
     .		RDNAME(1:length(RDNAME)),'".  Format error in group '
	write(*,*)STRI(1:132)
	write(*,*)'Last data read:'
	write(*,'(1p,6e12.4)')(DATARR(jarr+jj),jj=1,jbdry)
Clength(STRI))
	goto	999
 908	write(*,*)'>>> READAT: Data file "',
     .		RDNAME(1:length(RDNAME)),'".   Format error in line '
	write(*,*)STRI(1:length(STRI))
	goto	999
 909	write(*,*)'>>> READAT: Time dependent variable strings >',NTVAR
	goto	999
 910	write(*,*)'>>> READAT: Number of time dependent arrays cannot'//
     .			' exceed',NTARR
	goto	999
 911	write(*,*)'>>> READAT: Data file "',
     .		RDNAME(1:length(RDNAME)),'" error'
	write(*,'(10X,1A26,1I2)')'Number of radial points > ',NRDX
	goto	999
 912	write(*,*)'>>> READAT: Data file "',
     .		RDNAME(1:length(RDNAME)),'" error'
	write(*,'(8X,1A16,1A6,$)')'Input quantity: ',VNAM
	write(*,*)'Number of grid points must be > 1'
	goto	999
 913	write(*,*)'>>> READAT: File "',FILENA(1:length(FILENA)),
     .	'" reading error'
	goto	999
 914	write(*,*)'>>> READAT: File "',FILENA(1:length(FILENA)),
     .	'" open error'
	goto	999
 915	if (ierr .eq. 0)
     .	write(*,*)'>>> READAT: U-file "',FILENA(1:lname),
     .	'" is not found'
	if (ierr .eq. -1)
     .	write(*,*)'>>> READAT: Ex-file "',FILENA(1:lname),
     .	'" is not found'
	goto	999
 916	write(*,*)'>>> READAT: U-file "',FILENA(1:lname),
     .	'" reading error'
	goto	999
 917	write(*,*)'>>> READAT: Ex-file "',FILENA(1:lname),
     .	'" reading error'
	write(*,'(A,$)')'            '
	if (ierr .eq. 1) write(*,*)' Version, shot reading error'
	if (ierr .eq. 2) write(*,*)' Data reading error'
	if (ierr .eq. 3) write(*,*)' Time dimension is too small'
	if (ierr .eq. 4) write(*,*)' Radial dimension is too small'
	goto	999
 918	write(*,*)'>>> READAT: Ex-file "',FILENA(1:lname),
     .	'" open error'
	goto	999
 919	write(*,*)'>>> READAT: Ex-file "',FILENA(1:lname),
     .	'" reading error'
	write(*,'(A,$)')'            '
	write(*,*)' Required parameter "JETTONAME" is missed'
	goto	999
C 920	write(*,*)'>>> READAT: Ex-file "',FILENA(1:lname),
C     .	'" reading error'
C	write(*,'(A,$)')'            '
C	write(*,*)' Error reading data'
C	goto	999
 921	write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),
     .	'" error: wrong dimensionality'
	goto	999
 922	write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),
     .	'" error: radial grid size',jrad,' is larger than',NRDX
C	write(*,*)'>>> Re-compilation of all libraries is required'
	goto	999
 923	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),
     .	'" error: ambiguous ',VNAM(1:length(VNAM)),' definition'
	goto	999
 924	write(*,*)'>>> READAT: File const.inc read error'
	goto	999
 925	write(*,*)'>>> U-file "'
     .	,FILENA(1:length(FILENA)),'" read error'
	goto	999
 926	write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),
     .  '" read error. Possible reason: unrecognized format'
	goto	999
 927	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),
     &	     '",    line "',VNAM(1:length(VNAM)),'" error:'
	write(*,*)'>>> Tabulation is not is not permitted',
     &			  ' in this type of input'
	goto	999
 928	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),'" error:'
	write(*,*)'    "',VNAM,'" cannot vary in time'
	goto	999
 929	write(*,*)'>>> U-file "',FILENA(1:length(FILENA)),
     +	   '" unknown U-file type'
	goto	999
 930	write(*,*)'>>> File "',RDNAME(1:length(RDNAME)),'" error:'
     +			,' Wrong pointer to the U-file name'
	goto	999
 931	write(*,*)'>>> Warning <<< Wrong U-file name in  '
     +		,'"exp/'//RDNAME(1:length(RDNAME))//'"'
	goto	999
 932	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),
     &	    '",    format error in the group started with: '
	write(*,*)STRI(1:132)
	goto	999
 933	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),
     .	'" error: ',VNAM(1:length(VNAM)),' out of order'
	goto	999
 934	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),'" error:'
	write(*,*)'    Boundary must be defined in a single group'
	goto	999
 935	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),'" error:'
	write(*,*)'    Number of boundary points must be defined'
	goto	999
 936	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),'" error:'
	write(*,*)'    Boundary data array cannot exceed 25000'
	goto	999
 937	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),'" error:'
	write(*,*)'    No data input present for ',STRI(j1:j1+5)
	goto	999
 938	write(*,*)'>>> Data file "',RDNAME(1:length(RDNAME)),'" error:'
	write(*,*)'    More data items than data values for BND group'
	goto	999
 939	write(*,*)'>>> READAT: Buffer size exceeded'
	goto	999
 940	write(*,*)'>>> READAT: File "~/astra/runs/ex2a" read error'
	goto	999
 941	write(*,*)'>>> READAT: File "~/astra/runs/ex2a" open error'
	goto	999
 9411	write(*,*)'>>> READAT: File "~/astra/runs/ext2a" open error'
	goto	999
 9412	write(*,*)'>>> READAT: Ext-file reading error'
	goto	9414
 9413	write(*,*)'>>> READAT: Ex-file reading error'
 9414	write(*,'(A,$)')'            '
	if (ierr .eq. 1) write(*,*)' Version/shot reading error'
	if (ierr .eq. 2) write(*,*)' Data reading error'
	if (ierr .eq. 3) write(*,*)' Time dimension too small'
	if (ierr .eq. 4) write(*,*)' Radial dimension too small'
	goto	999
 942	write(*,*)'>>> Error in file "tmp/astra.log", line ',COMNAM
	goto	999
 999	call	a_stop
	end
C======================================================================|
	subroutine UF1RD(NCH,JDEV,DEVAR,FACTOR,ISHOT,IDEV,NAME,ERCODE)
C----------------------------------------------------------------------|
	implicit none
	integer	NCH,JDEV,ISHOT,IVAR,jdim,jtim,j3,jj,ERCODE,INDVAR
	logical	IFREE
	double precision	YY,VARDAT,FACTOR,DEVAR(*)
	include	'for/parameter.inc'
	include	'for/outcmn.inc'
	common	/EXPDAT/ VARDAT(3,NTVAR),INDVAR(NTVAR),IVAR
	character	STRI*32,STR1*32,IDEV*4,NAME*6
C----------------------------------------------------------------------|
	ERCODE = 0
	IFREE = .TRUE.
 1	continue
	read(NCH,32,ERR=63)STRI
C	if (STR1(32:32) .ne. ';')	goto
C Shot #, Device ID, U-file dimension,
	read(STRI(2:6),  fmt='(1I5)',ERR=63)ISHOT
	read(STRI(8:11), fmt='(1A4)',ERR=63)IDEV
	read(STRI(13:13),fmt='(1I1)',ERR=63)jdim
C	read(STRI(15:15),fmt='(1I1)',ERR=63)
	if(jdim.le.0 .or. jdim.ge.2)	goto	59
C Shot date
	read(NCH,32,ERR=63)STRI
C Number of associated scalar quantities
	read(NCH,32)STRI
	read(STRI,*,ERR=63)j3
	if (j3 .gt. 0)	then
	    do  jj=1,j3
		read(NCH,*,ERR=63)YY
		read(NCH,32,ERR=63)STRI
	    enddo
	endif
C Independent variable label
	read(NCH,32,ERR=63)STRI
C Dependent variable label
	read(NCH,32,ERR=63)STRI
	if (NAME .eq. ' ')	 goto	3
C	write(*,*)'"',STRI(2:7),'"',IFREE
C	write(*,*)'"',NAME,'"'
	if (STRI(2:7) .eq. NAME)	goto	3
 2	continue
	STR1 = STRI
	read(NCH,32,end=99)STRI
	if  ( STR1(1:10) .eq. STRI(1:10)
     +	.and. STR1(1:10) .eq. '**********')	goto	1
	goto	2

 3	continue
C Processing code
	read(NCH,32,ERR=63)STRI
	read(NCH,*,ERR=63)jtim
	if (jtim .lt. 1)		goto	67
	if (IVAR+jtim .gt. NTVAR)	goto	54
	if (jtim .eq. 1)	then
		IFDFVX(JDEV)=0
	else
		IFDFVX(JDEV)=1
	endif
C Read "time" array & function array
	if (.not.IFREE)	goto	83
C	write(*,*)"Floating format for ",FILENA(1:lname)
	read(NCH,*,ERR=82)(VARDAT(1,IVAR+jj),jj=1,jtim)
	read(NCH,*,ERR=82)(VARDAT(2,IVAR+jj),jj=1,jtim)
	goto	84
 82	IFREE = .FALSE.
	rewind(NCH)
	goto	1
 83	read(NCH,'(1X,6F13.1)',ERR=631)(VARDAT(1,IVAR+jj),jj=1,jtim)
	read(NCH,'(1X,6F13.1)',ERR=631)(VARDAT(2,IVAR+jj),jj=1,jtim)
C	write(*,*)"Fixed format for ",FILENA(1:lname)
 84	close(NCH)
C	write(*,'(1P,6E13.5)')(VARDAT(1,IVAR+jj),jj=1,jtim)
C	write(*,'(1P,6E13.5)')(VARDAT(2,IVAR+jj),jj=1,jtim)
	IFREE = .TRUE.
	DEVAR(JDEV) = factor*VARDAT(2,IVAR+1)
	do	jj = 1,jtim
	   IVAR = IVAR+1
	   INDVAR(IVAR) = JDEV
	   VARDAT(2,IVAR) = factor*VARDAT(2,IVAR)
	   VARDAT(3,IVAR)=0.
	enddo
	return
 54	ERCODE = 1
	call	a_stop
 59	ERCODE = 2
	call	a_stop
 63	ERCODE = 3
	call	a_stop
 67	ERCODE = 4
	call	a_stop
 99	ERCODE = 6			! Internal_name not found
	call	a_stop
 631	ERCODE = 5
	call	a_stop
 32	format(1A32)
	end
C=======================================================================
	subroutine	INTVAR
C-----------------------------------------------------------------------
C The time evolution of the input data is taken from
C	VARDAT(1,NTVAR)	- time			|
C	VARDAT(2,NTVAR)	- value			|  for variables
C	VARDAT(3,NTVAR)	- error (not used)	|
C
C For the current time, a value is stored in the array
C	DEVARX(NCONST)	- (description in the file for/const.inc)
C-----------------------------------------------------------------------
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
C----------------------------------------------------------------------|
C  NTVAR total number of the time slices for all variables in a data file
	common /EXPDAT/ VARDAT(3,NTVAR),INDVAR(NTVAR),IVAR
	integer	INDVAR,IVAR
	integer jj,N1,N2
	double precision	ydt,ydtr,ydtl,vardat
C----------------------------------------------------------------------|
C Input:
C	IVAR
C	INDVAR
C	IFDFVX
C	VARDAT
C	TIME
C Output:
C	DEVARX
C	DEVAR
C
C	jj  - ordinal number of the record in VARDAT(1-2-3,jj)
C N=INDVAR(jj) - ordinal number of the quantity in DEVAR(N) & DEVARX(N)
C IFDFVX(N) - type of variable
C IFDFVX = 0 - determined by the data file (independent on time),
C	= 1 - determined by the data file (dependent on time),
C	= 2 - determined by the MODEL,
C	= 3 - Keyboard
C	= 4 - any change of the variable is forbidden 
C		  eg. (AB, RTOR, ELONM, TRICH or set interactively)
	call	markloc("INTVAR"//char(0))
	jj = 0
	N1 = 0
	N2 = 0
 30	jj = jj+1
	if (jj .gt. IVAR)	return
	N2 = N1
	N1 = INDVAR(jj)
	if (IFDFVX(N1) .lt. 0)	goto	30
	if (IFDFVX(N1).eq.0 .or. N1.ne.N2)  DEVARX(N1) = VARDAT(2,jj)
	if (N1 .eq. N2)	    	then
	    if (TIME-VARDAT(1,jj-1)) 34, 31, 31
 31	    if (TIME-VARDAT(1,jj))   32, 33, 33
 32	    YDT = VARDAT(1,jj)-VARDAT(1,jj-1)
	    YDTR = (VARDAT(1,jj)-TIME)/YDT
	    YDTL = (TIME-VARDAT(1,jj-1))/YDT
	    DEVARX(N1) = VARDAT(2,jj)*YDTL+VARDAT(2,jj-1)*YDTR
	    goto	34
 33	    DEVARX(N1) = VARDAT(2,jj)
	    goto	34
	endif
 34	continue
	if (IFDFVX(N1) .le. 1)	DEVAR(N1) = DEVARX(N1)
	goto	30
	end
C=======================================================================
	integer	function	IFDEFX(XARNAM)
	implicit none
	character	XARNAM*(*)
	character	XARNAME*6
	integer		j
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/outcmn.inc'
C----------------------------------------------------------------------|
	XARNAME = XARNAM//"     "
	do	j=1,NARRX
	   if (XARNAME .eq. EXARNM(j))	goto	1
	enddo
	IFDEFX = 0			! False (not recognized)
	return
 1	continue
C	write(*,'(3A,$)')'"',XARNAME,'"'
C	if (IFDFAX(j) .eq. -1)	write(*,*)" is not defined"
C	if (IFDFAX(j) .ne. -1)	write(*,*)" is defined"
	if (IFDFAX(j) .eq. -1)	IFDEFX = 0 ! False (not defined)
	if (IFDFAX(j) .ne. -1)	IFDEFX = 1 ! True (X-array is defined)
	end
C=======================================================================
	blockdata ASTRA_DATA
	implicit none
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	integer	NCONST32,j
	parameter(NCONST32=NCONST-32)
C	integer	j1,A_Wtotal,A_W1total,A_CarTotal,A_Extotal
C	parameter(A_Wtotal=2*NRD*NRD,	A_W1total=(2*NRD+7)*NRD)
C	parameter(A_Cartotal=32*NRD,	A_Extotal=NRD*NARRX)
C----------------------------------------------------------------------|
C Default constants
	data	NSBR /0/	NXOUT /0/
C Time interval control names
	data	(DTNAME(j),j=1,24)/
     1		'dRout ','dTout ','dPout ','Time  ',
     2		'TAUmin','TAUmax','TAUinc','DELvar',
     3		'Iterex','N/used','Tinit ','Tscale',
     4		'NB1   ','NUF   ','Xaxis ','Xdeflt',
     5		'NB2EQL','NEQUIL','NBND  ','Xflag ',
     6		'DTeql ','MEQUIL','Tpause','Tend  '/
	data	(DTNAME(j),j=25,64)/
     4		'DTeq1 ','BEeq1 ','ENeq1 ',' Keq1 ',
     5		'DTeq2 ','BEeq2 ','ENeq2 ',' Keq2 ',
     6		'DTeq3 ','BEeq3 ','ENeq3 ',' Keq3 ',
     7		'DTeq4 ','BEeq4 ','ENeq4 ',' Keq4 ',
     8		'DTeq5 ','BEeq5 ','ENeq5 ',' Keq5 ',
     9		'DTeq6 ','BEeq6 ','ENeq6 ',' Keq6 ',
     1		'DTeq7 ','BEeq7 ','ENeq7 ',' Keq7 ',
     2		'DTeq8 ','BEeq8 ','ENeq8 ',' Keq8 ',
     3		'DTeq9 ','BEeq9 ','ENeq9 ',' Keq9 ',
     4		'DTeq10','BEeq10','ENeq10',' Keq10'/
	data	(DTNAME(j),j=65,104)/
     5		'DTeq11','BEeq11','ENeq11',' Keq11',
     6		'DTeq12','BEeq12','ENeq12',' Keq12',
     7		'DTeq13','BEeq13','ENeq13',' Keq13',
     8		'DTeq14','BEeq14','ENeq14',' Keq14',
     9		'DTeq15','BEeq15','ENeq15',' Keq15',
     1		'DTeq16','BEeq16','ENeq16',' Keq16',
     2		'DTeq17','BEeq17','ENeq17',' Keq17',
     3		'DTeq18','BEeq18','ENeq18',' Keq18',
     4		'DTeq19','BEeq19','ENeq19',' Keq19',
     5		'DTeq20','BEeq20','ENeq20',' Keq20'/
C Note: The order of arrays in this data statement must be the same
C	as in the COMMON block A_CARX file for/status.inc
	data    (EXARNM(j),j=1,32) /
     1		'CAR1X ','CAR2X ','CAR3X ','CAR4X ',
     2		'CAR5X ','CAR6X ','CAR7X ','CAR8X ',
     3		'CAR9X ','CAR10X','CAR11X','CAR12X',
     4		'CAR13X','CAR14X','CAR15X','CAR16X',
     5		'CAR17X','CAR18X','CAR19X','CAR20X',
     6		'CAR21X','CAR22X','CAR23X','CAR24X',
     7		'CAR25X','CAR26X','CAR27X','CAR28X',
     8		'CAR29X','CAR30X','CAR31X','CAR32X'/
	data    (EXARNM(j),j=33,NARRX) /
     9		'F0X   ','F1X   ','F2X   ','F3X   ','F4X   ',
     &		'F5X   ','F6X   ','F7X   ','F8X   ','F9X   ',
     1		'MUX   ','MVX   ','GNX   ','SNX   ','PEX   ',
     2		'PIX   ','PRADX ','TEX   ','TIX   ','NEX   ',
     3		'CUX   ','ZEFX  ','VRX   ','SHX   ','ELX   ',
     4		'TRX   ','G11X  ','G22X  ','G33X  ','DRODAX',
     5		'IPOLX ','NIX   ','VPOLX ','VTORX ','SLATX'/
	data	(ARXUSE(j),j=1,NARRX) /NARRX*0/
C Initial times
	data	TAU/.000001/  TEQ/NSBMX*-99999./  LTOUT/1/  IPOUT/1/
C Output and time step control values
	data	(DELOUT(j),j=1,24)/
C		DROUT,	DTOUT,	DPOUT,	TIME,
     1		.01,     .01,     .01,     0.,
C		TAUMIN,	TAUMAX,	TAUINC,	DELVAR,
     2		.000001, .05,     1.1,     .1,
C		ITEREX,	ITERIN,	TINIT,	TSCALE,
     3		1.,       1.,      0.,     1.,
C		NB1R,	NUFR,	XOUT,	XINPUT,
     4		41.,     41.,      1.,     1.,
C		NB2EQL,	NEQUIL,	NBNDR,	XFLAGR,
     5		1.,       0.,      0.,     0.,
C		DTEQL,  MEQUIL, TPAUSE,	TEND,
     6		0.,      1.,      1.e20,   1.e20/
C Default subroutine call intervals: DTEQ(4,20)
	data	(DELOUT(j),j=25,64)/
     1		.0,  -99999.,  99999.,    -1.,
     2		.0,  -99999.,  99999.,    -1.,
     3		.0,  -99999.,  99999.,    -1.,
     4		.0,  -99999.,  99999.,    -1.,
     5		.0,  -99999.,  99999.,    -1.,
     6		.0,  -99999.,  99999.,    -1.,
     7		.0,  -99999.,  99999.,    -1.,
     8		.0,  -99999.,  99999.,    -1.,
     9		.0,  -99999.,  99999.,    -1.,
     &		.0,  -99999.,  99999.,    -1./
	data	(DELOUT(j),j=65,104)/
     1		.0,  -99999.,  99999.,    -1.,
     2		.0,  -99999.,  99999.,    -1.,
     3		.0,  -99999.,  99999.,    -1.,
     4		.0,  -99999.,  99999.,    -1.,
     5		.0,  -99999.,  99999.,    -1.,
     6		.0,  -99999.,  99999.,    -1.,
     7		.0,  -99999.,  99999.,    -1.,
     8		.0,  -99999.,  99999.,    -1.,
     9		.0,  -99999.,  99999.,    -1.,
     &		.0,  -99999.,  99999.,    -1./
C Screen parameters
	data	XSCMAX/640/ YSCMAX/350/ DXLET/8/ DYLET/13/ IDT/5/
	data	XWX/470/    XWY/10/     XWW/660/ XWH/550/
C Astra colors: #0 - background, ##1-7 - plots 1-7
C		#7 - also color for REPORT (user's output)
C		#8 - not used
C		#9 - delete curve, #10 - standard text color (black)
C		#11 - Black
C		#12 - warning messages, axis upper marks in View
C		#13 -> #15 -reserved (black)
C	data	COLTAB/
C     1		0, 0,	50,0,	3, 0,	62,0,	37,0,
C     2		5,0,	26,0,	30,0,	50,0,	0, 0,
C     3		1,0,	1,0,	50,0,	1,0,	1,0,
C     4		1,0,	1,0,	1,0,	1,0,	1,0,
C     5		1,0,	1,0,	1,0,	1,0,	1,0,
C     6		1,0,	1,0,	1,0,	1,0,	1,0,
C     7		1,0,	1,0/
	data	NWINDX/NRW*0/      OSHIFT/NRW*0./    OSHIFR/NRW*0./
	data	MARKT/NRW*0/       MARKR/NRW*0/      GRAL/NRW*0./
	data	NAMEX/NRW*'      '/NAM7/'Tmin','Tmax','Tmark','Style'/
	data	IFDFVX /NCONST*-1/ TIM7/0,9999,9999,1/
	data	NSCR/9*0/	   MODEY/1/
C Boundary:
	data	NBND  /0/	   NBNT/0/
C General arrays:
	data	ELON  /NRD*1./     TRIA  /NRD*0./    SHIF  /NRD*0./
	data	SHIV  /NRD*0./	   VTOR /NRD*0./     VPOL  /NRD*0./
	data	TE    /NRD*.01/    TI    /NRD*.01/   NE    /NRD*.1/
	data	CU    /NRD*.1/     NN    /NRD*1.d-5/ TN    /NRD*.01/
	data	SNN   /NRD*0./	   PET   /NRD*0./    PIT   /NRD*0./
	data	SFF0  /NRD*0./
	data	SFF1  /NRD*0./	   SFF2  /NRD*0./    SFF3  /NRD*0./
	data	SFF4  /NRD*0./	   SFF5  /NRD*0./    SFF6  /NRD*0./
	data	SFF7  /NRD*0./	   SFF8  /NRD*0./    SFF9  /NRD*0./
	data	MU    /NRD*.3/     ZEF   /NRD*1./    AMAIN /NRD*1./
	data	MV    /NRD*.0/     CV    /NRD*0./    FV    /NRD*0./
	data	ZMAIN /NRD*1./     GP    /3.1415926/ GP2   /6.283185/
	data	G33   /NRD*1./	   IPOL  /NRD*1./    SQEPS /NRD*0.2/
	data	VOLUM /NRD*1./     NI    /NRD*.1/
	data	ER    /NRD*0./	   PBLON /NRD*0./    PBPER /NRD*0./
	data	NHYDR /NRD*0./     NDEUT /NRD*0./    NTRIT /NRD*0./
	data	NHE3  /NRD*0./     NALF  /NRD*0./
	data	NIZ1  /NRD*0./     NIZ2  /NRD*0./    NIZ3  /NRD*0./
	data	ZIM1  /NRD*1./     ZIM2  /NRD*1./    ZIM3  /NRD*1./
	data	B0DB2 /NRD*1./     BDB02 /NRD*1./
	data	PETOT /NRD*1.d-4/  PITOT /NRD*1.d-4/
	data	DVN   /NRD*0.d0/   DVE   /NRD*0.d0/  DVI   /NRD*0.d0/
	data	DSN   /NRD*0.d0/   DSE   /NRD*0.d0/  DSI   /NRD*0.d0/
	data	DVF0  /NRD*0.d0/   DVF1  /NRD*0.d0/  DVF2  /NRD*0.d0/
	data	DVF3  /NRD*0.d0/   DVF4  /NRD*0.d0/  DVF5  /NRD*0.d0/
	data	DVF6  /NRD*0.d0/   DVF7  /NRD*0.d0/  DVF8  /NRD*0.d0/
	data	DVF9  /NRD*0.d0/
	data	DSF0  /NRD*0.d0/   DSF1  /NRD*0.d0/  DSF2  /NRD*0.d0/
	data	DSF3  /NRD*0.d0/   DSF4  /NRD*0.d0/  DSF5  /NRD*0.d0/
	data	DSF6  /NRD*0.d0/   DSF7  /NRD*0.d0/  DSF8  /NRD*0.d0/
	data	DSF9  /NRD*0.d0/
	data	SDN  /NRD*0.d0/    PDE   /NRD*0.d0/  PDI   /NRD*0.d0/
	data	SD0  /NRD*0.d0/    SD1   /NRD*0.d0/  SD2   /NRD*0.d0/
	data	SD3  /NRD*0.d0/    SD4   /NRD*0.d0/  SD5   /NRD*0.d0/
	data	SD6  /NRD*0.d0/    SD7   /NRD*0.d0/  SD8   /NRD*0.d0/
	data	SD9  /NRD*0.d0/
C Set Xarrays to zero
	data	CAR1X  /NRD*0.d0/  CAR2X  /NRD*0.d0/  CAR3X  /NRD*0.d0/
	data	CAR4X  /NRD*0.d0/  CAR5X  /NRD*0.d0/  CAR6X  /NRD*0.d0/
	data	CAR7X  /NRD*0.d0/  CAR8X  /NRD*0.d0/  CAR9X  /NRD*0.d0/
	data	CAR10X /NRD*0.d0/  CAR11X /NRD*0.d0/  CAR12X /NRD*0.d0/
	data	CAR13X /NRD*0.d0/  CAR14X /NRD*0.d0/  CAR15X /NRD*0.d0/
	data	CAR16X /NRD*0.d0/  CAR17X /NRD*0.d0/  CAR18X /NRD*0.d0/
	data	CAR19X /NRD*0.d0/  CAR20X /NRD*0.d0/  CAR21X /NRD*0.d0/
	data	CAR22X /NRD*0.d0/  CAR23X /NRD*0.d0/  CAR24X /NRD*0.d0/
	data	CAR25X /NRD*0.d0/  CAR26X /NRD*0.d0/  CAR27X /NRD*0.d0/
	data	CAR28X /NRD*0.d0/  CAR29X /NRD*0.d0/  CAR30X /NRD*0.d0/
	data	CAR31X /NRD*0.d0/  CAR32X /NRD*0.d0/  DRODAX /NRD*0.d0/
	data	F0X /NRD*0.d0/  F1X /NRD*0.d0/  F2X /NRD*0.d0/
	data	F3X /NRD*0.d0/  F4X /NRD*0.d0/  F5X /NRD*0.d0/
	data	F6X /NRD*0.d0/  F7X /NRD*0.d0/  F8X /NRD*0.d0/
	data	F9X /NRD*0.d0/  MUX /NRD*0.d0/  MVX /NRD*0.d0/
	data	GNX /NRD*0.d0/  SNX /NRD*0.d0/  PRADX /NRD*0.d0/
	data	PEX /NRD*0.d0/  PIX /NRD*0.d0/  ZEFX  /NRD*0.d0/
	data	TEX /NRD*0.d0/  TIX /NRD*0.d0/  NEX   /NRD*0.d0/
	data	CUX /NRD*0.d0/  VRX /NRD*0.d0/  NIX   /NRD*0.d0/
	data	SHX /NRD*0.d0/  G11X/NRD*0.d0/  IPOLX /NRD*0.d0/
	data	ELX /NRD*0.d0/  G22X/NRD*0.d0/  VPOLX /NRD*0.d0/
	data	TRX /NRD*0.d0/  G33X/NRD*0.d0/  VTORX /NRD*0.d0/
	data	SLATX /NRD*0.d0/
C Global variables:
	data	AB    	/.3/	ABC	/.3/	AMJ	/2./
	data	AWALL	/.4/	BTOR	/3./	ELONG	/1./
	data	ELONM	/1./	ENCL	/.002/	ENWM	/.02/
	data	IPL	/.3/	RTOR	/1.5/   ROC	/.3/
	data	SHIFT	/0./	NNCL	/.001/	NNWM	/.0001/
	data	GN2E	/0./	GN2I	/0./    UEXT	/0./
	data	TRIAN	/.0/	TRICH	/.0/	ZMJ	/1./
	data	WNE	/.03/	WTE	/.03/	WTI	/.03/
	data	TSTART	/0./
	data	NB1	/41/	NA1	/41/	NNCX	/200/
	data	NAB	/41/	NUF	/41/	UPDWN	/0./
	data	NA	/40/	NSTEPS	/0/
C Constants:
	data	CPTOT/0./	CPTEQL/0./	CPTGRA/0./	CPT/0./
	data	(CPTSBR(j),j=1,NSBMX)	/NSBMX*0./
	data	(LEQ(j),j=1,19)		/19*-1/
	data	(CONSTF(j),j=1,16)	/16*1./
     1		(CONSTF(j),j=17,32)	/16*0./
     2		(CONSTF(j),j=33,NCONST)	/NCONST32*1./
	data	VERSION/'                                '/
	data	XLINE1 /'                                '/
	data	XLINE2 /'                                '/
	data	CNFILE/'***'/	NBFILE/'***'/	ECFILE/'***'/
	data	ICFILE/'***'/	MSFILE/'***'/
C	data	((WORK (j,j1),j=1,NRD),j1=1,2*NRD)	/A_Wtotal*0./
C	data	((WORK1(j,j1),j=1,NRD),j1=1,2*NRD+7)	/A_W1total*0./
C	data	((CAR(j,j1),j=1,NRD),j1=1,32)		/A_Cartotal*0./
C	data	((EXT(j,j1),j=1,NRD),j1=1,NARRX)	/A_Extotal*0./
	end
C=======================================================================
