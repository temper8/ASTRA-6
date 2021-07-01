C======================================================================|
C------------------------ MAIN, ENTEQU, NOCHANGE ----------------------|
C I/O unit usage:
C 1 read  "model.tmp", "status.inc", "const.inc", "fml.dir", "fnc.dir"
C 3 write "model.out" from the main
C 2 write  "equftn.tmp", "detvar.tmp", "radout.tmp", "timout.tmp"
C	   "inivar.tmp", "ininam.tmp", "inipsi.tmp", "xar_usage"
C 7 write  "model.txt"
C 8,9,... read include formulas from "astra/fml/" (sbrs APPFML and TMPFML)
C======================================================================|
C	implicit	UNDEFINED(A-Z)
	implicit	none
	include		'.srv/nambuf.inc'
	include		'.srv/tmpbuf.inc'
C NAMER	 -  4 maximal name length, NRW maximal number of channels
C NAMET	 -  4 maximal name length, NRW maximal number of channels
C NAMEX	 -  6 maximal variable name length, NRW maximal number of channels
	character*1	NAMARR(6,100),SCALER(6,NRW),SCALET(6,NRW)
	character*1	SBRBUF(2000), EQBUF(2000), PRFINI(78,NEQNS)
	character*1	LINTXT(133),  LINFOR(133), DISTR(131),  STR(80)
	character*1	FORMR(75,NRW),FORMT(75,NRW),CNSBUF(2000),NAME(6)
	character*1	STRI(132),    STRIN*132,   STRFOR*132,  STR80*80
	character*1	ch0, CHBS, KK, TEQ, SCNS*3, ALKEY*3
	character*4	NAMER(NRW),   NAMET(NRW)
C In PRFINI, EQNAM, COEFLN, LEQ, NCOEF "NEQNS" is a maximum number of eqs.
C In COEFLN(12,NEQNS)  "12" is a max number of recognizable lhs + eq name,
C		 eg.: TE, HE, DE, XE, CE, PE, PET, GN2E, QEB, QETB, ROE
C		 or,  F1, DF1, VF1, GF1, SF1, SFF1, QF1B, QFF1B, RO1
C			either TEB or (QEB and QETB) is used at a time
	character*78	COEFLN(12,NEQNS),  ARRSTR(100)
	character	MFNAME(3)*40,  INCNAM*80
	character*6	SP6,  NAME6,   NAMEX(NRW),  STSTEP(NSBMX,5)
	equivalence	(STRI,STRIN),(LINFOR(2),STRFOR),(STR,STR80)
	integer	iargc,J,J1,J2,JEQ,JSB,IFERR,jst,jlb,jrb,jout,jgt,L,jbs
	integer	LEN1,LEN2,LEN3,LEN4,LSTR,LOUT,LUR,LENC,LVAR,jsh,LCNS
	integer	LEQ(NEQNS),NCOEF(NEQNS),NARRS,LNAMAR(100),NWX(NRW),ICBUF
	integer	IFIPL,IFSYNT,IFTOUT(NRW),IFROUT(NRW),NSBR,LOCSBR(NSBMX)
	integer length,npos,rnpos,joneof,ioneof,in6eof,lengtf,ifnum    ! Fctns
	integer	NROUT,NXOUT,NTOUT,ICON,NPR,NPR1,ILNG,NB,NB1,NB2,NB3,NB4
	integer	IEBUF,ITYPE,lrb,jd,NCH,NCHI,NSC,ISYMB,NRW6,ERASYM,IFPEI
	double precision RR
	logical	JPEI
	parameter (NRW6=6*NRW)
	data	MFNAME /'tmp/model.txt','tmp/model.tmp','tmp/model.out'/
	data	STSTEP
     &/NSBMX*'Each_T',NSBMX*'-Infty',NSBMX*'+Infty',NSBMX*' ',NSBMX*' '/
	data	NAMER  /NRW*'    '/ NAMET/NRW*'    '/NAMEX/NRW*'      '/
	data	SCALER /NRW6*' '/   SCALET/NRW6*' '/ ARRSTR/100*' '/
	data	NARRS/0/ IFTOUT/NRW*0/ IFROUT/NRW*0/ TEQ/'='/ ICBUF/1/
	data	LEQ/NEQNS*-1/ NCOEF/NEQNS*0/ IEBUF/1/ IFPEI/0/ JSB/1/
	data	NCHI /0/  jsh/1/ LOCSBR/NSBMX*-4/
	data	NROUT/0/  NTOUT/0/ NXOUT/0/  IFIPL/0/ IFSYNT/0/ NSBR/0/
	data	SP6/'      '/
C	data	(STRI(j),j=1,132) /132*' '/
C----------------------------------------------------------------------|
	do j=1,2000
	   CNSBUF(j) = char(0)
	enddo
	do j=1,5000
	   DTVBUF(j) = char(0)
	enddo
	do j=1,10000
	   TMPBUF(j) = char(0)
	enddo
	do j=1,80
	   INCNAM(j:j) = char(0)
	enddo
	lvar = iargc()
	do j=0,lvar
	   call	getarg(j,STRIN)
	   j1=lengtf(STRIN)
C	   if (j .eq. 0) write(*,'(3a,1i6)')
C     &		'Task: "',STRIN(1:j1),'",   PID =',getpid()
C	   if (j .gt. 0) write(*,'(a,i2,3a)')
C     &		'      Argument #',j,'  "',STRIN(1:j1),'"'
	   if (STRIN(1:j1) .eq. 'T_flag')	Tflag = 1
	enddo
C	write(*,*)"-T option is enabled in the command line"
C----------------------------------------------------------------------|
C Read file ".exe/version"
C	call	OPENRD(1,'.exe/version',0,IERR)
C	if (IERR .gt. 1)	then
C	   write(*,*)'>>> Warning >>>  Unknown version, continue anyway'
C	endif
        J = 92
        CHBS = char(J)			! Define Backslash
        J = 0
        ch0 = char(J)			! Define null character
	write(EQBUF(IEBUF),132) char(0)
 	call	ENTEQU
C	write(*,*)(ARRNAM(j2),j2=1,NARR)
C	do	j2=1,NARR
C	   j1 = length(6,ARRNAM(j2))
C	   write(*,*)(ARRNAM(j2),j2=1,NARR)
C	   write(*,*)(length(6,ARRNAM(j2)),j2=1,NARR)
C	   if (ARRNAM(j2)(j1:j1) .eq. 'X')write(*,*)ARRNAM(j2)
C	enddo

C FML repeatition control (KEYOUT <= 0 - off, KEYOUT >= 1 - on)
	KEYOUT=0
Channel for "model.txt"
	call	OPENWT(7,MFNAME(1),0,IFERR)
	if (IFERR.gt.1)	then
	    write(*,*)'MODEL (MAIN): Cannot open file ',MFNAME(1)
	    call	exit(1)
	endif
Channel for "model.tmp"
	NCH = 1
	call	OPENRD(NCH,MFNAME(2),0,IFERR)
	if(IFERR.gt.1)	then
	    write(*,*)'MODEL analyser: Cannot open file ',MFNAME(2)
	    call	exit(1)
	endif
Channel for "model.out"
	call	OPENWT(3,MFNAME(3),0,IFERR)
	if (IFERR.gt.1)	then
	    write(*,*)'MODEL (MAIN): Cannot open file ',MFNAME(3)
	    call	exit(1)
	endif
 	write(7,*)'=====   Variables definition   ====='
 1	continue
C	write(*,*)"NCH =",NCH
	read(NCH,132,ERR=81,END=81)STRI
 132	format(132A1)
	j1 = length(132,STRI)
	write(3,132)(STRI(J),J=1,J1)
C '!' - 1st position - the entire string is ignored
	if(STRI(1).eq.'!')	then
	   write(7,132)(STRI(J),J=1,J1)
	   goto	1
	endif
	NB=1
	LEN1=132
	call BLEX(LEN1,STRI)
	if(LEN1.le.0)		goto	1
C	write(*,*)
C	write(*,*)'Complete string:   "',(STRI(J),J=1,LEN1),'"'

C ';' - control string delimiter
 2	LSTR=NPOS(LEN1-NB+1,STRI(NB),';')-1
C LEN1-NB+1 .eq. LSTR+1 - ";" is found 
C LEN1-NB+1 .eq. LSTR   - ";" is not found (end of the string)
C	write(*,*)
C	write(*,*)'String to analyze  "',(STRI(j),j=NB,NB+LSTR-1),'"'
C     >		  ,"   Length, Start, End ",LSTR,NB,NB+LSTR-1
	if (LSTR .le. 0)	goto	80
C '!' - control string comment
	ICON=NPOS(LSTR,STRI(NB),'!')-1
	if (ICON .le. 0)	goto	80
C '#' - pre-processor command
	if (STRIN(NB:NB+7) .eq. "#include")	then
	   if (NCHI .ne. 0)	goto	96
C	   write(*,*)NB+8,NB+LSTR-1,'  "',STRIN(NB+8:NB+LSTR-1),'"'
	   INCNAM = STRIN(NB+8:NB+LSTR-1)//char(0)
C	   write(*,*)'"',INCNAM(1:length(80,INCNAM(1:))),'"'
	   NCH = NCH+1
	   NCHI = NCH
	   call	OPENRD(NCH,"equ/"//INCNAM,0,IFERR)
	   if (IFERR .gt. 1)	goto	97
	   goto	80
	endif

C	NB=NB+LSTR+1
C	write(*,*)'Rest of the string "',(STRI(j),j=NB,LEN1),'"'
C     >		  ,"Length, Start, End ",LEN1-NB,NB,LEN1
C	if(LEN1-NB) 1,2,2

C ':' - equation or subroutine
	NPR=NPOS(ICON,STRI(NB),':')
	if (NPR .le. ICON)	then
	   if (STRI(NB+NPR) .eq. '=')	then		! ":="
	      if (STRIN(NB+NPR-3:NB+NPR-2) .eq. 'NI')	then
		 LEQ(6) = 0
	      else
		 write(*,'(A,$)')' >>> Illegal RHS of ":=" operator in'
		 write(*,'(A)')' line "'//STRIN(NB:NB+NPR)//'..."'
	      endif
C	      write(*,*)'"',STRI(NB+NPR-1:NB+LSTR-1),'"'
	      do j=NB+NPR-1,LEN1
		 STRI(j) = STRI(j+1)
	      enddo				! if ":=" found
	      LEN1 = LEN1-1			! remove ":"
	      ICON = ICON-1			! and ignore
	      LSTR = LSTR-1
C	      write(*,*)'"',STRI(NB+NPR-1:NB+LSTR-1),'"'
	      goto	3
	   endif
	   goto	40
	endif
 3	continue
C '&' - launch an independent process and continue Astra execution
C       (the symbol ':' should not appear in the command string)
	NPR=NPOS(ICON,STRI(NB),'&')
	if (NPR .le. ICON)	goto	60
C       
C '=' - equation coefficient or constant
	NPR=NPOS(ICON,STRI(NB),TEQ)
	if (NPR .le. ICON)	then
 	   if (STRI(NB+NPR) .eq. '>')	then	! "=>"
C	      write(*,*)'"',STRI(NB+NPR-1:NB+LSTR-1),'"'
	      do j=NB+NPR,LEN1
		 STRI(j) = STRI(j+1)
	      enddo				! if "=>" found
	      LEN1 = LEN1-1			! remove ">"
	      ICON = ICON-1			! and ignore
	      LSTR = LSTR-1
C	      write(*,*)'"',STRI(NB+NPR-1:NB+LSTR-1),'"'
	   endif
	   if (STRI(NB+NPR-1) .eq. '<')	goto	50	! if "=<"
	   goto	30
	endif
C '\' - radial output
	NPR=NPOS(ICON,STRI(NB),CHBS)
	if (NPR .le. ICON)	goto	20
C '_' - time output
	NPR=NPOS(ICON,STRI(NB),'_')
	if (NPR .le. ICON)	goto	10
				goto	80
C----------------------------------------------------------------------|
C Time output command line
 10	NTOUT=NTOUT+1
 104	format(' >>> Warning: More than',1I4,A)
	if (NTOUT .eq. NRW+1)	then
	   write(7,104)NRW,' time channels'
	   write(*,104)NRW,' time channels'
	endif
	if (NTOUT .ge. NRW+1)	goto	80

	jout = NPR+NPOS(ICON-NPR,STRI(NB+NPR),'>')-1	! Exclude the name
C	write(*,*)jout,ICON,NPR				! from analysis
	if (jout .lt. ICON)	then
	   IFTOUT(NTOUT) = 1
	endif
	ICON = min(ICON,jout)
	NPR1=NPR-1
	if(NPR.le.5)		goto	12
	NPR1=4
	write(7,103)NTOUT
 103	format(' >>> Warning: Too long name in channel',1I3,' <<<')
 12	if(NPR.le.1)		goto	16
	call COPY(NPR1,STRI(NB),NAMET(NTOUT))
 16	ILNG=ICON-NPR
	NB1=NB+NPR
	NSC=NPOS(ILNG,STRI(NB1),'_')-1
	if(NSC.ge.ILNG-1)	goto	13
	NB2=NB1+NSC+1
	LEN2=ILNG-NSC-1
	call COPY(LEN2,STRI(NB2),SCALET(1,NTOUT))
 13	call COPY(NSC,STRI(NB1),FORMT(2,NTOUT))
	FORMT(1,NTOUT)=char(NSC)
				goto	80
C----------------------------------------------------------------------|
C Radial output command line
C  NB   - position of the command line in the whole string,
C  ICON - the line length,
C  NPR  - position of the 1st '\',     (NPR-1) - the name length
C  NB1  - position after the 1st '\',	ILNG - remaining length

 20	NROUT=NROUT+1
	if (NROUT .eq. NRW+1)	then
	   write(7,104)NRW,' radial profiles'
	   write(*,104)NRW,' radial profiles'
	endif
	if (NROUT .ge. NRW+1)	goto	80
C	write(*,*)
C	write(*,*)'"',(STRI(j),j=1,LEN1),'"'		! Total string
C	write(*,*)'"',(STRI(j),j=NB,NB+ICON-1),'"'	! String to analyze

	jlb = NPR+NPOS(ICON-NPR,STRI(NB+NPR),'[')-1
	jrb = NPR+NPOS(ICON-NPR,STRI(NB+NPR),']')-1
	LEN2 = jrb-jlb
	if (LEN2 .lt. 0)	then
	   write(7,91)(STRIN(NB:NB+ICON-1))
	   write(*,91)(STRIN(NB:NB+ICON-1))
	   IFSYNT = 1
	   goto	99
	elseif (LEN2 .eq. 0)	then
	   goto	23
	endif
C	write(*,*)'"',(STRI(j),j=NB+jlb+1,NB+jlb+LEN2-1),'"'
	j1 = jlb+NPOS(LEN2-1,STRI(NB+jlb+1),',')
	j2 = j1-jlb-1
C	if (j2 .le. 0) write(*,*)
C     +		      "Left  boundary is set to the default value"
	if (j2 .le. 0)	goto	21
C	write(*,*)'"',(STRI(j),j=NB+jlb+1,NB+jlb+LEN2-1),'"',j2
C	write(*,*)'"',STRIN(NB+jlb+1:NB+j1-1),'"'
C     +		,IFNUM(STRI(NB+jlb+1),j2)
	call	NUM2ST(NROUT,SCNS,LCNS)
	if (IFNUM(STRI(NB+jlb+1),j2) .eq. 1)	then
	   LOUT = j2+LCNS+13
	   call	APPBUF(LOUT,SP6//'GRAL('//SCNS(1:LCNS)//')='//
     +		STRIN(NB+jlb+1:NB+j1-1),CNSBUF,ICBUF,2000)
	   if (LOUT .eq. 0)	goto	94
	   LOUT = j2
	   do	j=1,LOUT
	      LINTXT(j) = STRI(NB+jlb+j)
	   enddo	      
	else
	   KEYOUT=-1
	   STR80(1:) = char(j2)//STRIN(NB+jlb+1:NB+j1)
	   call	ANLFML(-1,STR,LINTXT,LINFOR) ! line '\'[...
C	   write(*,*)'for:  "',(LINFOR(1+j2),j2=1,ichar(LINFOR(1))),'"'
C	   write(*,*)'txt:  "',(LINTXT(1+j2),j2=1,ichar(LINTXT(1))),'"'
	   KEYOUT=0
	   LUR  = ichar(LINFOR(1))
	   LOUT = LUR+LCNS+13
	   call	APPBUF(LOUT,SP6//'GRAL('//SCNS(1:LCNS)//')='//
     +		STRFOR(1:LUR),DTVBUF,IVBUF,5000)
	   if (LOUT .eq. 0)	goto	93
	   LOUT=ichar(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
	endif
 21	continue
C	if (j1 .eq. jrb) write(*,*)"Comma not found"
C     +			"Right boundary is set to the default value"
	j2 = jrb-j1-1
	if (j2 .le. 0)	goto	22
C	write(*,*)'Right boundary: "',(STRI(NB+j1+j),j=1,j2),'"'
C	write(*,*)'"',STRIN(NB+j1+1:NB+jrb),'"',IFNUM(STRI(NB+j1+1),j2)
	call	NUM2ST(NROUT,SCNS,LCNS)
	if (IFNUM(STRI(NB+j1+1),j2) .eq. 1)	then
	   LOUT = j2+LCNS+13
	   call	APPBUF(LOUT,SP6//'GRAP('//SCNS(1:LCNS)//')='//
     +		STRIN(NB+j1+1:NB+jrb),CNSBUF,ICBUF,2000)
	   if (LOUT .eq. 0)	goto	94
	   LOUT = j2
	   do	j=1,LOUT
	      LINTXT(j) = STRI(NB+j1+j)
	   enddo	      
	else
	   KEYOUT = -1
	   STR80(1:) = char(j2)//STRIN(NB+j1+1:NB+jrb)
	   call	ANLFML(-1,STR,LINTXT,LINFOR) ! line '\' ...]
C	write(*,*)'LINFOR:   "',(LINFOR(1+j),j=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j),j=1,ichar(LINTXT(1))),'"'
	   KEYOUT = 0
	   LUR = ichar(LINFOR(1))
	   LOUT = LUR+LCNS+13
	   call	APPBUF(LOUT,SP6//'GRAP('//SCNS(1:LCNS)//')='//
     +		STRFOR(1:LUR),DTVBUF,IVBUF,5000)
	   if (LOUT .eq. 0)	goto	93
	   LOUT=ichar(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
	endif
 22	continue
	write(7,131)(NAME(J),J=1,NPR-1),TEQ,(LINTXT(J),J=1,LOUT)
	STRIN(1:) = STRIN(1:NB+jlb-1)//STRIN(NB+jrb+1:LEN1)
	LEN1 = LEN1-jrb+jlb-1
	LSTR = LSTR-jrb+jlb-1
	ICON = ICON-jrb+jlb-1

 23	continue
C	write(*,*)'"',(STRI(j),j=NB,NB+ICON-1),'"'	! Compressed string
C	write(*,*)'"',(STRI(j),j=1,LEN1),'"'
	jout = NPR+NPOS(ICON-NPR,STRI(NB+NPR),'>')-1
	if (jout .lt. ICON)	then
	   IFROUT(NROUT) = 1
	endif
	ICON = min(ICON,jout)
	isymb = -1
C	write(*,*)(STRI(j),j=NB,NB+ICON-1),'	Channel ',NROUT
C Set undefined channel to zero (will be overwritten if defined)
	call	COPY(2,"0.",FORMR(2,NROUT))
	FORMR(1,NROUT)=char(2)
	NPR1=NPR-1
	if(NPR.le.5)		goto	24
C Too long name in channel NROUT. Cut to 4.
	NPR1=4
	write(7,103)NROUT
 24	continue
C if name_length > 1 write channel name to NAMER
	if(NPR.gt.1)	call COPY(NPR1,STRI(NB),NAMER(NROUT))
	ILNG=ICON-NPR
	NB1=NB+NPR
 25	continue
	if (ILNG .le. 0)	goto	80
	NSC=NPOS(ILNG,STRI(NB1),CHBS)-1
C NSC = string-length between two '\'
	if (NSC .eq. -1)	goto	80	! No more '\' found
C if \\ encounered goto	28
	if (NSC .eq. 0)		goto	28	! Double backslash "\\" found
C if 2nd \ is not found do not write scaler (goto	26)
	if (isymb.lt.0)		goto	26
C	write(*,*)'   Writing "',(STRI(j),j=NB1,NB1+NSC-1),'" to SCALER'
	call COPY(NSC,STRI(NB1),SCALER(1,NROUT))
				goto	80
 26	continue
	isymb = isymb+1
	if (isymb.lt.1)	goto	27
	write(7,*)'>>> IFSYNT ERROR:  Radial output line "',
     >			(STRI(j),j=NB,NB+ICON-1),'"'
	write(*,*)'>>> IFSYNT ERROR:  Radial output line "',
     >			(STRI(j),j=NB,NB+ICON-1),'"'
	IFSYNT = 1
 27	continue
C	write(*,*)'   Writing "',(STRI(j),j=NB1,NB1+NSC-1),'" to FORMR'
	call COPY(NSC,STRI(NB1),FORMR(2,NROUT))
	FORMR(1,NROUT)=char(NSC)
	ILNG = ILNG-NSC-1
	NB1 = NB1+NSC+1
				goto	25
 28	continue
	ILNG = ILNG-1
	NB1 = NB1+1
	NSC=NPOS(ILNG,STRI(NB1),CHBS)-1
C skip empty channel
	if (NSC.eq.0)	goto	25
	isymb = isymb+1
	if (isymb.eq.0) isymb = 1
	NXOUT=NXOUT+1
	if (NXOUT .eq. NRW+1)	then
	   write(7,104)NRW,' data channels'
	   write(*,104)NRW,' data channels'
	endif
	if (NXOUT .ge. NRW+1)	goto	80

C	write(*,*)'   Writing "',(STRI(j),j=NB1,NB1+NSC-1),'" to NAMEX',
C     >	',    NXOUT =',NXOUT,',    NROUT =',NROUT,NSC,NAMEX(NXOUT)
	call COPY(NSC,STRI(NB1),NAMEX(NXOUT))
	call UPCASE(NSC,NAMEX(NXOUT))
	NWX(NXOUT)=NROUT
	ILNG = ILNG-NSC-1
	NB1 = NB1+NSC+1
				goto	25
C----------------------------------------------------------------------|
C Equation coefficient command line,  NPR - position of '='
 30	if(NPR.eq.1)		goto	39		!-> Syntax error
C	write(*,*)
C	write(*,*)'String to analyze  "',(STRI(j),j=NB,NB+LSTR),'"'
C     >		  ,"  Length =",LSTR-NB+1
C     >		  ,"  Start =",NB
C     >		  ,"  End =",NB+LSTR
	call	SETNAM(NPR-1,STRI(NB),NAME)
	call	UPCASE(6,NAME)
	L = IONEOF(NAME,NTMP,TMPNAM)
C lhs of '+' command:
	if (STRIN(NB:NB+2) .eq. 'NB1')	then
	   write(*,*)'>>> WARNING: The statement "',
     &		STRIN(NB:NB+LSTR),'" is ignored.'
	   write(*,*)'             Use data file to assign NB1.'
	   goto	80
	endif
C Enable NI:=... special treatment but keep back compatibility for NI=...
C if L==0 NI=... goes to NAMARR(,), otherwise to TMPBUF(NBEG(LNI))
C NBEG(l) > 0 results in defining TMPNAM(l)
	if (LEQ(6).lt.0 .and. L.eq.LNI)	L = 0
	if (L .gt. 0)	NBEG(L)=IBUF
	ILNG=ICON-NPR
	NB1=NB+NPR
	write(STRI(NB1-1),132)	char(ILNG)
C left: IPL, LEXT, UEXT (boundary conditions for flux Eq.)
	if( L.eq.LIPL .or. L.eq.LLEXT.or. L.eq.LUEXT)	goto	34
C left: NEB, TEB, TIB (boundary conditions for NE,... Eq.)
	if( L.eq.LNEB .or. L.eq.LTEB .or. L.eq.LTIB)	goto	33
C left: QNB, QNNB, QEB, QETB, QIB, QITB, FjB, QFjB, QFFjB
	if( L.eq.LQNB  .or.L.eq.LQEB  .or.L.eq.LQIB
     .	.or.L.eq.LQNNB .or.L.eq.LQETB .or.L.eq.LQITB)	goto	33
	do	j2=0,9
	   if (L.eq.LFB(j2))				goto	33
	   if (L.eq.LQFB(j2) .or. L.eq.LQFFB(j2))	goto	33
	enddo
C other TMPNAM (equation coeff. or TE,TI,NE,CU,Fj)
	if(L.ne.0)					goto	31
C left: Not from the TMPNAM list (38 - general array, 34 - constant)
	if(JONEOF(NAME,NARR,ARRNAM))			34,34,38
C left: Equation coefficient
 31	KEYOUT=1
C left: CU,MU,NE,TE,TI,Fj
C	if (L .eq. LNI)		then
	if (L.eq.LCU .or. L.eq.LMU .or. L.eq.LMV .or. 
     .	    L.eq.LNE .or. L.eq.LTE .or. L.eq.LTI)	KEYOUT=0
	do	j2=0,9
	   if (L.eq.LF(j2))	KEYOUT=0
	enddo

C appending fml to TMPBUF -> tmp/detvar.tmp
C	   write(*,*)"Call from model 1:",KEYOUT,0
	call ANLFML(0,STRI(NB1-1),LINTXT,LINFOR)	! line '='
C	write(*,*)'LINFOR:   "',(LINFOR(1+j2),j2=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j2),j2=1,ichar(LINTXT(1))),'"'
C	write(*,*)IBUF,NPR,ichar(LINFOR(1))
	if (L.eq.LPE.or.L.eq.LPET.or.L.eq.LPI.or.L.eq.LPIT)	then
	   if (KEYOUT .eq. 2) IFPEI = 1
	endif
	KEYOUT=0
	if(IBUF+NPR+11.gt.10000)		goto	92
C writing "	lhs=fml" to TMPBUF
C this record goes to tmp/detvar.tmp
	call COPY(6,SP6,TMPBUF(IBUF+1))
	call COPY(NPR-1,STRI(NB),TMPBUF(IBUF+7))
C	write(*,'(132A1)')'"',(TMPBUF(IBUF+6+j),j=1,NPR-1),'"'
 	call COPY(4,'(J)=',TMPBUF(IBUF+NPR+6))
	LUR=ichar(LINFOR(1))
	if(IBUF+NPR+LUR+11.gt.10000)	goto	92
	call COPY(LUR,LINFOR(2),TMPBUF(IBUF+NPR+10))
	LUR = NPR+LUR+9
	write(TMPBUF(IBUF),132)		char(LUR)
	write(TMPBUF(IBUF+LUR+1),132)	char(0)
C	call COPY(1,char(0),TMPBUF(IBUF-1))
	IBUF=IBUF+LUR+2

	LOUT=ichar(LINTXT(1))
	call COPY(NPR,NAME,DISTR)
	call COPY(3,'(r)',DISTR(NPR))
	DISTR(NPR+3)=TEQ
	call COPY(LOUT,LINTXT(2),DISTR(NPR+4))
 131	format(1X,130A1)
C initial distributions
	LENC=NPR+LOUT+3
C	write(*,*)l,':   "',TMPNAM(l),'"',LHE,':   "',TMPNAM(LHE),'"'
	if(L.eq.LCN.or.L.eq.LDN.or.L.eq.LHN.or.L.eq.LNE.or.
     +	   L.eq.LSN.or.L.eq.LSNN.or.L.eq.LXN)		JEQ=1
	if(L.eq.LCE.or.L.eq.LDE.or.L.eq.LHE.or.L.eq.LPE.or.
     +	   L.eq.LPET.or.L.eq.LTE.or.L.eq.LXE)		JEQ=2
	if(L.eq.LCI.or.L.eq.LDI.or.L.eq.LHI.or.L.eq.LPI.or.
     +	   L.eq.LPIT.or.L.eq.LTI.or.L.eq.LXI)		JEQ=3
	if(L.eq.LCC.or.L.eq.LCD.or.L.eq.LCU.or.L.eq.LCV.or.
     +	   L.eq.LMU.or.L.eq.LMV.or.L.eq.LDC.or.L.eq.LHC.or
     +	  .L.eq.LXC)					JEQ=4
	do	j2=0,9
	   if ( L.eq.LF(j2)  .or. L.eq.LDF(j2) .or.
     +		L.eq.LVF(j2) .or. L.eq.LGF(j2) .or.
     +		L.eq.LSF(j2) .or. L.eq.LSFF(j2) )	JEQ=10+j2
	enddo
	LENC=MIN0(LENC,78)
	if (JEQ .eq .0)	then
	   write(7,*)'>>> Coefficient of undefined equation is used'
	   write(7,*)(DISTR(jd),jd=1,LENC)
	   write(*,*)'>>> Coefficient of undefined equation is used'
	   write(*,*)(DISTR(jd),jd=1,LENC)
	   IFSYNT = 1
	   goto	99
	endif
	NCOEF(JEQ)=NCOEF(JEQ)+1
C write model.txt
C	write(*,*)(DISTR(jd),jd=1,LENC)
	call COPY(LENC,DISTR,COEFLN(NCOEF(JEQ),JEQ))
	if(L.eq.LNE .or. L.eq.LTE .or. L.eq.LTI .or. L.eq.LCU)	then
	   call COPY(LENC,DISTR,PRFINI(1,JEQ))
	   call COPY(1,char(0),PRFINI(LENC+1,JEQ))
	endif
	do	j2=0,9
	   if ( L.eq.LF(j2) )	then
	      call COPY(LENC,DISTR,PRFINI(1,JEQ))
	      call COPY(1,char(0),PRFINI(LENC+1,JEQ))
	   endif
	enddo
C	write(*,*)L,JEQ,LENGTF(PRFINI(1,JEQ))
C     .	,(PRFINI(J1,JEQ),J1=1,LENGTF(PRFINI(1,JEQ)))
				goto	80
 33	continue
C----------------------------------------------------------------------|
C Left hand side is one of: NEB,TEB,TIB,FjB,QNB,QEB,QIB,QNNB,
C			    QETB,QITB,QFjB,QFFjB
C RHS analysis:		LINTXT -> model.txt;	LINFOR  -> tmp/xxx.tmp
C----------------------------------------------------------------------|
	KEYOUT=-1
C appending boundary conditions to TMPBUF
C	   write(*,*)"Call from model3:",0,0
C RHS analysis:		LINTXT -> model.txt;	LINFOR  -> tmp/xxx.tmp
	call ANLFML(0,STRI(NB1-1),LINTXT,LINFOR)	! line "bnd="

C	write(*,*)'LINFOR:   "',(LINFOR(1+j2),j2=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j2),j2=1,ichar(LINTXT(1))),'"'
C	write(*,*)l,':   "',TMPNAM(l),'"'
	KEYOUT=0

	if(IBUF+13.gt.10000)		goto	92
C writing "	lhs=fml" to TMPBUF
	call COPY(6,SP6,TMPBUF(IBUF+1))
	if ( L.eq.LNEB .or. L.eq.LTEB .or. L.eq.LTIB )	then
	   J1=2
	   call COPY(J1,NAME,TMPBUF(IBUF+7))
	   call COPY(5,'(ND1)',TMPBUF(IBUF+J1+7))
	   J1 = J1+5
	endif
	do	j2=0,9
	   if (L.eq.LFB(j2))	then
	      J1=2
	      call COPY(J1,NAME,TMPBUF(IBUF+7))
	      call COPY(5,'(ND1)',TMPBUF(IBUF+J1+7))
	      J1 = J1+5
	   endif
	enddo
	if ( L.eq.LQNB .or. L.eq.LQEB .or. L.eq.LQIB)	then
	   J1=3
	   call COPY(J1,NAME,TMPBUF(IBUF+7))
	endif
	if ( L.eq.LQNNB .or. L.eq.LQETB .or. L.eq.LQITB)	then
	   J1=4
	   call COPY(J1,NAME,TMPBUF(IBUF+7))
	endif
	do	j2=0,9
	   if (L.eq.LQFB(j2))	then
	      J1=4
	      call COPY(J1,NAME,TMPBUF(IBUF+7))
	   endif
	enddo
	do	j2=0,9
	   if (L.eq.LQFFB(j2))	then
	      J1=5
	      call COPY(J1,NAME,TMPBUF(IBUF+7))
	   endif
	enddo
	if (IBUF+6+J1+8+LUR .gt. 10000)		goto	92

	TMPBUF(IBUF+J1+7)=TEQ
	LUR=ichar(LINFOR(1))

	if(IBUF+LUR+J1+9.gt.10000)	goto	92
	call COPY(LUR,LINFOR(2),TMPBUF(IBUF+J1+8))
	LUR = LUR+j1+7
	write(TMPBUF(IBUF),132)		char(LUR)
	write(TMPBUF(IBUF+LUR+1),'(1A1)')	ch0
C	write(*,*)'123456789012345678901234567890',LUR,LUR-j1-8,j1
C	write(*,132)'"',(TMPBUF(j2),j2=IBUF+1,IBUF+LUR),'"'
	IBUF=IBUF+LUR+2

	LOUT=ichar(LINTXT(1))
	if(L.eq.LNEB .or. L.eq.LQNB  .or. L.eq.LQNNB )	JEQ=1
	if(L.eq.LTEB .or. L.eq.LQEB  .or. L.eq.LQETB )	JEQ=2
	if(L.eq.LTIB .or. L.eq.LQIB  .or. L.eq.LQITB )	JEQ=3
	do	j2=0,9
	   if(L.eq.LFB(j2).or.L.eq.LQFB(j2).or.L.eq.LQFFB(j2)) JEQ=10+j2
	enddo
	NCOEF(JEQ)=NCOEF(JEQ)+1
	if (J1 .eq. 7)	then	! Boundary condition of the 1st kind
	   call	COPY(2,NAME,COEFLN(NCOEF(JEQ),JEQ))
C write (a)
	   call COPY(5,'(a_b)',COEFLN(NCOEF(JEQ),JEQ)(3:7))
	else
	   call	COPY(J1,NAME,COEFLN(NCOEF(JEQ),JEQ))
	endif
	COEFLN(NCOEF(JEQ),JEQ)(J1+1:J1+1)=TEQ
	call	COPY(LOUT,LINTXT(2),COEFLN(NCOEF(JEQ),JEQ)(J1+2:J1+2))
C	write(*,*)(COEFLN(NCOEF(JEQ),JEQ)(J3:J3),J3=1,J1+LOUT+1)
				goto	80
C----------------------------------------------------------------------|
C-- left hand side is treated as radially independent (NOT an array) --|
 34	continue
	LVAR = L
C Find NAME ordinal number on the list of variables in "const.inc":
	L=IONEOF(NAME,NMAIVA,MAIVAR)
C L=0 - NAME is not on the list:
	IF(L.eq.0)		goto	36
C left: variable
	if  (LVAR.eq.LIPL)   IFIPL = 1
	if( (LVAR.eq.LUEXT .or. LVAR.eq.LLEXT) .and. IFIPL.eq.0)   then
		J1 = -IN6EOF('IPL   ',NMAIVA,MAIVAR)
		IFIPL = 1
C forbid 'IPL' reading from a data file when assigned in a model	
		call	NOCHANGE(J1,DTVBUF,IVBUF)
	endif
C write strings preventing 'NAME' reading from a data file 
C	when it is assigned in a model	
C 		if(IFDFVX(L).le.2)	IFDFVX(L)=2
C 		if(IFDFVX(L).le.2)
C	     >	'NAME'='LINFOR'
	call	NOCHANGE(L,DTVBUF,IVBUF)

	KEYOUT=-1
C	   write(*,*)"Call from model4:",-1,-1
	call ANLFML(-1,STRI(NB1-1),LINTXT,LINFOR)	! line '=' 
	KEYOUT=0
	if(IVBUF+NPR+9.gt.5000)		goto	93
	call COPY(NPR-1,NAME,DTVBUF(IVBUF+9))
	DTVBUF(IVBUF+NPR+8)=TEQ
	LUR=ichar(LINFOR(1))
	if(IVBUF+NPR+LUR+9.gt.5000)	goto	93
	call COPY(LUR,LINFOR(2),DTVBUF(IVBUF+NPR+9))
	write(DTVBUF(IVBUF),132)	char(NPR+LUR+8)
	IVBUF=IVBUF+NPR+LUR+9
					goto	37
C left: unrecognised
C right: time dependent expression
 36	KEYOUT=-1
C	   write(*,*)"Call from model6:",-1,-1
	call ANLFML(-1,STRI(NB1-1),LINTXT,LINFOR)	! line '=' 
	KEYOUT=0
	LUR=ichar(LINFOR(1))
C	write(*,*)'LINFOR:   "',(LINFOR(1+j),j=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j),j=1,ichar(LINTXT(1))),'"'
C	if (STRFOR(1:7).eq.'RADIAL(' )	then
C	   write(*,*)STRFOR(1:10)
C	endif
C	write(STRFOR,'(7A1)')(LINFOR(1+j),j=1,min(7,LUR))
C2	if (STRFOR(1:7).eq.'RADIAL(' )	then
C	   write(*,*)'"',(LINFOR(1+j),j=1,LUR),'"',LUR,NPR
! Special treatment is needed in order to map RHO to AMETR
C Defining array variable at intermediate point
C This record goes to detvar.tmp before "jdetv" loop
C2	   if(IVBUF+NPR+LUR+36.gt.5000)	goto	93
C2	   J1 = index(STRFOR,',')+1
C2	   call COPY(9,'      YR=',DTVBUF(IVBUF+1))
C2	   call COPY(LUR-J1,LINFOR(J1+1),DTVBUF(IVBUF+10))
C2	   write(DTVBUF(IVBUF),132)	char(LUR-J1+9)
C2	   IVBUF=IVBUF+LUR-J1+10
C1	   call COPY(16,'      YR=RFA(YR)',DTVBUF(IVBUF+1))
C1	   write(DTVBUF(IVBUF),132)	char(16)
C1	   IVBUF=IVBUF+17
C2	   call COPY(6,SP6,DTVBUF(IVBUF+1))
C2	   call COPY(NPR-1,NAME,DTVBUF(IVBUF+7))
C2	   DTVBUF(IVBUF+NPR+6)=TEQ
C2	   call COPY(J1-1,LINFOR(2),DTVBUF(IVBUF+NPR+7))
C2	   call COPY(3,'YR)',DTVBUF(IVBUF+NPR+J1+6))
C2	   write(DTVBUF(IVBUF),132)	char(NPR+J1+8)
C2	   IVBUF=IVBUF+NPR+J1+9
C2	else
	   if(IVBUF+NPR+LUR+7.gt.5000)	goto	93
	   call COPY(6,SP6,DTVBUF(IVBUF+1))
	   call COPY(NPR-1,NAME,DTVBUF(IVBUF+7))
	   DTVBUF(IVBUF+NPR+6)=TEQ
	   call COPY(LUR,LINFOR(2),DTVBUF(IVBUF+NPR+7))
	   write(DTVBUF(IVBUF),132)	char(NPR+LUR+6)
	   IVBUF=IVBUF+NPR+LUR+7
C2	endif
 37	LOUT=ichar(LINTXT(1))
	do	j=1,LOUT
	   LINTXT(j) = LINTXT(j+1)
	enddo
	write(7,131)(NAME(J),J=1,NPR-1),TEQ,(LINTXT(J),J=1,LOUT)
	goto	80
C------------------- left hand side is an array -----------------------|
 38	continue
	call	SETNAM(6,NAME,NAME6)
	j1 = length(6,NAME6)
	if (NAME6(j1:j1) .ne. 'X')	goto	381
	JARX = JARX+1
	ARXNAM(JARX) = NARR+JONEOF(NAME,NARR,ARRNAM)
C	write(*,*)"lhs:    ",ARXNAM(JARX),' ',NAME6,JARX,NARR
C	if (JARX .eq. 1)write(*,*)JARX,(ARRNAM(j2),j2=1,JARX)
	if (JARX .eq. 1)	goto	381
	do	j1=1,JARX-1		! Exclude names repeated on lhs
	   if (ARXNAM(j1) .eq. ARXNAM(JARX))	then
	      ARXNAM(JARX) = -1
	      JARX = JARX-1
	      goto	381
	   endif
	enddo
 381	continue

	NARRS=NARRS+1
	call COPY(NPR-1,NAME,NAMARR(1,NARRS))
	LNAMAR(NARRS)=NPR-1
	call COPY(ILNG+1,STRI(NB1-1),ARRSTR(NARRS))
C	do	J1=1,NARRS
C	write(*,'(I3,60A)')J1,'	',(NAMARR(J,J1),J=1,LNAMAR(J1))
C	enddo
				goto	80
 39	write(7,134)(STRI(J),J=NB,NB+ICON-1)
 134	format(' >>> Coefficient string error <<< ',50A1)
				goto	80
C----------------------------------------------------------------------|
C Equation command string
C STRI(NB)	- first symbol
C STRI(NB+ICON)	- last symbol (';' or <Ret>)
C ICON 		- string length till ';'
C (NPR-1) - {name+formal_parameters} length before the 1st ':'
C	    (position of ':' in the command string)
C	NPR=NPOS(ICON,STRI(NB),':')

 40	continue
C	write(*,'(132A1)')'"',(STRI(j),j=NB,NB+NPR-1),'"',NPR
C	write(*,'(<LSTR+2>A,$)')
C     >		'String to analyze  "',(STRI(j),j=NB,NB+LSTR-1),'"'
C	write(*,'(A,4I5)')'   Length, Start, End ',LSTR,NB,NB+LSTR-1,NPR
	JST = NB
	if (NPR .le. 2)		goto	90	! Error
	if (NPR .eq. 3)		goto	410	! Main eqn
	if (NPR .gt. 4)		goto	46	! SBR line analysis

	NAME6 = STRIN(NB:NB+NPR-2)
	call	UPCASE(3,NAME6)
C	write(*,*)'"',NAME6,'"'
	LEQ(5) = -1			! Solver is set by NEQUIL (default)
	if (NAME6(1:3) .eq. "EQ0")	then
	   LEQ(5) = -2				! No equilibrium (forced)
	   goto	41
	endif
	if (NAME6(1:3) .eq. "EQX")	then
	   LEQ(5) = 0				! External equilibrium
	   goto	41
	endif
	if (NAME6(1:3) .eq. "E3M")	then
	   LEQ(5) = 1				! Electrodynamic moments
	   goto	41
	endif
	if (NAME6(1:3) .eq. "ESC")	then
	   LEQ(5) = 2				! ESC
	   goto	41
	endif
	if (NAME6(1:3).eq."ESP" .or. NAME6(1:6).eq."SPIDER")	then
	   LEQ(5) = 3				! SPIDER
	   goto	41
	endif
	if (NAME6(1:3).eq."EQS" .or. NAME6(1:6).eq."SCOPE")	then
	   LEQ(5) = 4				! SCoPE
	   goto	41
	endif
	if (NAME6(3:3) .eq. "*")	JST = NB+2
C	jst = NPOS(ICON,STRI(NB),'*')+NB-1
	goto	410

 41	continue
	JEQ = 5
C	jst = NPOS(LSTR,STRI(NB),'&')-1
C	if (jst .lt. ICON)	then		! Instruction '&' found
C	   LEQ(5) = LEQ(5)+11
C	endif
	goto	80

 410	continue
C set JEQ: 0 - sbr, 1 - NE, 2 - TE, 3 - TI, 4 - CU, 5 - EQ, 
C	   6<->9 - not used, 10 - F0, 11 - F1, 12 - F2, etc.
	JEQ = 0
	if (NPR.ne.3 .and. JST.eq.NB)	goto	46	! sbr line analysis
	call	UPCASE(2,STRI(NB))
	if (STRI(NB)//STRI(NB+1).eq.'NE') JEQ=1
	if (STRI(NB)//STRI(NB+1).eq.'TE') JEQ=2
	if (STRI(NB)//STRI(NB+1).eq.'TI') JEQ=3
	if (STRI(NB)//STRI(NB+1).eq.'CU') JEQ=4
	if (STRI(NB).eq.'F') then
	   do	j2=0,9
	      write(KK,'(I1)')j2
	      if (STRI(NB)//STRI(NB+1) .eq. 'F'//KK)  JEQ = NEQNS-9+j2
	   enddo
	endif
	if (JEQ .eq. 0)		goto	46	! sbr line analysis
C-----------------------ICON------------------------>;?
C<-(NB1-1)->:<----------ILNG------------------------>;?
C<-(NB1-1)->:<-(NSC+1)->:<--------(ILNG-NSC-1)------>;?
C<----------:------(NB2-1)->:
C ICON	  - full string length excluding ';' or <CR>
C ILNG	  - remaining string length after the 1st ':'
C (NB1-1) - absolute position of the 1st ':'	-> STRI(NB1-1)=':'
C (NB2-1) - absolute position of the 2nd ':'	-> STRI(NB2-1)=':'
C NPR	  - relative position of the 1st ':'	-> STRI(NB-1+NPR)=':'
C (NSC+1) - relative position of the 2nd ':', 0 if none
C (ILNG-NSC-1) - remaining string length after the 2nd ':'
C ITYPE: -2 - unrecognised or syntax error -> Stop)
C 	 -1   no equation & no assignment encountered
C	  0 - AS or "Name" appears on the l.h.s. of assignment command
C         1 - EQ, 2 - FU, 3 - Flux with a multiplicative factor
C	      Note: option Tj:AS:5.2; is not supported !!! Should be added
C	write(*,*)
C	write(*,*)'Eqn command line  "',(STRI(J),J=NB,NB+ICON-1),'"'
C	write(*,*)'Eqn command line  "',(STRI(NB+j),j=0,ICON-1),'"'
C	write(*,*)'Eqn command line  "',STRIN(NB:NB+ICON-1),'"'
	ILNG=ICON-NPR
	NB1=NB+NPR
	call	UPCASE(1,STRI(NB1))

	ITYPE=-2
	ITYPE=1
	NSC = NPOS(ILNG,STRI(NB1),':')
	if (NSC .lt. ILNG) ITYPE=2	! 2nd ":" is not the last symbol
	if (STRI(NB1) .eq. 'A') ITYPE=0			! :A...;
	if (JST .gt. NB)	ITYPE = ITYPE+8
C	if(STRI(NB1).eq.'E') ITYPE=1			! :E...
C	if(STRI(NB1).eq.'[') ITYPE=1			! :[
C	if(STRI(NB1).eq.':' .and. ILNG.eq.1) ITYPE=1	! ::;
C	if(STRI(NB1).eq.':' .and. ILNG.gt.1) ITYPE=2	! ::...;
C	if(    ILNG .eq. 0 ) ITYPE=1			! :;
C	if(STRI(NB1).eq.'F') ITYPE=2			! :F...;
C	write(*,'(3(A,I3))')"JEQ =",JEQ,'   ITYPE =',ITYPE
C     >		,'   ILNG =',ILNG
C	if (  ITYPE .eq. -2)		goto	90	! Error
C	write(*,'(/3A,I5)')'"',STRIN(NB:NB+ICON-1),'"    JEQ =',JEQ
	jlb = index(STRIN(NB:NB+ICON-1),'[')
	jrb = index(STRIN(NB:NB+ICON-1),']')
	jgt = index(STRIN(NB:NB+ICON-1),'>')		! Symbol is not used
	jbs = index(STRIN(NB:NB+ICON-1),CHBS)		! CHBS = '\'
	if (jlb.eq.0 .and. jrb.eq.0)	goto	43	! No "Tag"
	if (jlb.eq.0 .or.  jrb.eq.0)	goto	90	! Syntax error
	if (jlb .gt. jrb)		goto	90	! Syntax error
	if (jbs .eq. 0)			goto	42	! No '\'
	if (jbs .gt. jrb)		goto	90	! Syntax error
	if (jbs .lt. jlb)		goto	90	! Syntax error
	LINAPP(JEQ) = 1			! Linear extrapolation to the LCFS
	LRB = jrb-jbs-1
	if (jrb-jbs .gt. 1)	then			! value@LCFS is set
C	   write(*,*)LRB
C	   write(*,*)'"',STRIN(NB+jbs:NB+jbs+LRB-1),'"',NB+jbs,NB+jbs+LRB-1
	   STR80(1:)=char(LRB)//STRIN(NB+jbs:NB+jbs+LRB-1)
	   KEYOUT = 0
	   call	ANLFML(-2,STR,LINTXT,LINFOR)
C	   write(*,*)'LINFOR:   "',(LINFOR(1+j2),j2=1,ichar(LINFOR(1))),'"'
C     >			,ichar(LINFOR(1))
C	   write(*,*)'LINTXT:   "',(LINTXT(1+j2),j2=1,ichar(LINTXT(1))),'"'
	   j2=6+6+2+ichar(LINFOR(1))
	   STR80(1:)=char(j2)//SP6
	   if (JEQ .eq. 1)	STR80(8:) = 'NE'
	   if (JEQ .eq. 2)	STR80(8:) = 'TE'
	   if (JEQ .eq. 3)	STR80(8:) = 'TI'
	   if (JEQ .ge. 10)	then
	      write(KK,'(I1)')JEQ-10
	      STR80(8:) = 'F'//KK
	   endif
	   STR80(10:)='(NA1)='//STRFOR(1:ichar(LINFOR(1)))
	   LINAPP(JEQ) = 2			! Linear extrapolation between ND1 and NA1
	endif
C	jrb = jrb-1
C	ICON = ICON-1
C	ILNG = ILNG-1
C	LSTR = LSTR-1
C	write(*,*)'Old'
C	write(*,*)jlb,jbs,jrb,ICON,LRB
C	write(*,*)'"',STRIN(NB:NB+ICON-1),'"',jeq,LRB
	jrb = jrb-LRB-1
	ICON = ICON-LRB-1
	ILNG = ILNG-LRB-1
	LSTR = LSTR-LRB-1
C	write(*,*)NB+jbs-1,'   <-',NB+jbs+LRB,NB+ICON,jrb
	do j=NB+jbs-1,NB+ICON-1		! Remove "\" from the string
C	   write(*,*)j,'   <-',j+LRB+1
	   STRI(j) = STRI(j+LRB+1)
	enddo
C	write(*,*)'New','"',STRI(NB+ICON),'"'
C	write(*,*)jbs-1,ICON,jrb
C	write(*,*)'"',STRIN(NB:NB+ICON-1),'"',jeq,LINAPP(JEQ)

 42	j = index(STRIN(NB+jlb:NB+jrb-2),',')		! Position of ','
	jout = 0					! jout = "Tag"
	if (j .gt. 1)	read(STRIN(NB+jlb:NB+jlb+j-2),*,ERR=90)jout
C jout = 0,1,2 are processed as a/ABC, a, \rho/ROC, respectively.
C Other values are treated as \rho
C	write(*,*)'"',STRIN(NB+jlb-1:NB+jlb-1),'"',jlb
C	write(*,*)'"',STRIN(NB+jrb-1:NB+jrb-1),'"',jrb
C	write(*,*)'"',STRIN(NB+jlb:NB+jrb-2),
C     >		'"  Length =',jrb-jlb-1
C	write(*,*)'Type   =',jout,
C     >		'"    Rbndry ="',STRIN(NB+jlb+j:NB+jrb-2),'"'

C write the right boundary for each eqn to TMPBUF and then to equftn.tmp
	if (JEQ .eq. 1)  call SETNAM(3,'RON',NAME)
	if (JEQ .eq. 2)	 call SETNAM(3,'ROE',NAME)
	if (JEQ .eq. 3)	 call SETNAM(3,'ROI',NAME)
	if (JEQ .ge. 10)	 then
	   write(KK,'(I1)')JEQ-10
	   call SETNAM(2,'RO',NAME)
	   NAME(3)=KK
	   NAME(4)=' '
	   NAME(5)=' '
	   NAME(6)=' '
	endif
	J2=IONEOF(NAME,NTMP,TMPNAM)
	if(J2 .gt. 0)	NBEG(J2)=IBUF
	KEYOUT=0
	KK = STRI(NB+jlb+j-1)
	LRB = jrb-jlb-j-1
C	write(*,*)jlb,j,jbs,jrb,LRB
	write(STRI(NB+jlb+j-1),132) char(LRB)
C appending boundary conditions to TMPBUF
	call ANLFML(-2,STRI(NB+jlb+j-1),LINTXT,LINFOR)	! line ':' 
C	write(*,*)'LINFOR:   "',(LINFOR(1+j2),j2=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j2),j2=1,ichar(LINTXT(1))),'"'
	STRI(NB+jlb+j-1) = KK
	KEYOUT=0
	LUR=ichar(LINFOR(1))
C	write(*,*)'LINFOR:   "',(LINFOR(1+j2),j2=1,LUR),'"',LUR
C Length of RON,ROE,ROI,RO0,RO1,...,RO9
	J1=3
	if (IBUF+LUR+J1+29 .gt. 10000)		goto	92
	call COPY(6,SP6,TMPBUF(IBUF+1))
	call COPY(J1,NAME,TMPBUF(IBUF+7))
	TMPBUF(IBUF+J1+7)=TEQ
	call COPY(LUR,LINFOR(2),TMPBUF(IBUF+J1+8))
	LUR = LUR+J1+7
	write(TMPBUF(IBUF),132)		char(LUR)
	IBUF=IBUF+LUR+1
C	write(*,*)"jout =",jout
	call	SETNAM(3,NAME,NAME6)
	if (LINAPP(jeq) .eq. 2)	then
C	   write(*,*)jeq
C	   write(*,*)'"',STR80(2:2+ichar(STR80(1:1))),'"',ichar(STR80(1:1))
	   LUR = ichar(STR80(1:1))
	   call	COPY(LUR,STR80(2:),TMPBUF(IBUF+1))
	   write(TMPBUF(IBUF),132)	char(LUR)
	   IBUF = IBUF+LUR+1
	endif
	if (jout .eq. 0)	then
	   STR80 = SP6//NAME6(1:3)//'=RFA('//NAME6(1:3)//')'
	   LUR = 18
	endif
	if (jout .eq. 1)	then
	   STR80 = SP6//NAME6(1:3)//'=RFAN('//NAME6(1:3)//')'
	   LUR = 19
	endif
C	if (jout .eq. 1)	then
C	   STR80 = SP6//NAME6(1:3)//'=RFA('//NAME6(1:3)//')'
C	   LUR = 18
C	endif
	if (jout .eq. 2)	then
	   STR80 = SP6//NAME6(1:3)//'='//NAME6(1:3)//'*ROC'
	   LUR = 17
	endif
	call	COPY(LUR,STR80,TMPBUF(IBUF+1))
	write(TMPBUF(IBUF),132)	char(LUR)
	IBUF = IBUF+LUR+1
	LUR = 19
	STR80 = SP6//'ND1=NODE('//NAME6(1:3)//')'
	call	COPY(LUR,STR80,TMPBUF(IBUF+1))
	write(TMPBUF(IBUF),132)		char(LUR)
	write(TMPBUF(IBUF+LUR+1),132)	ch0
	IBUF = IBUF+LUR+2
	LOUT=ichar(LINTXT(1))
	do	j=1,LOUT
	   LINTXT(j) = LINTXT(j+1)
	enddo
C	write(*,*)"LINTXT:"
C	write(*,131)(NAME(J),J=1,NPR-1),TEQ,(LINTXT(J),J=1,LOUT)

 43	continue
	if (ILNG .eq. 0)	then
	   NSC = -1
	else   
	   NSC = NPOS(ILNG,STRI(NB1),':')-1
	endif
C	write(*,*)"ITYPE =",ITYPE,ILNG,NSC,'"',STRI(NB1+NSC),'"'
	if (ITYPE .eq. 0)		goto	44
	if (JEQ.ne.2 .and. JEQ.ne.3)	then
	   if (ITYPE.eq.2)	goto	90
	   goto	44
	endif
C	write(*,*)'"',STRI(NB1:NB1+ILNG-1),'"',NB
C	write(*,'(6(A,I5))')"    jlb =",jlb,"    NSC =",NSC
C     >	,"    NB1 =",NB1,"    NB1+NSC =",NB1+NSC
C	goto	90

 431	continue
	if (NPR+NSC .gt. ICON)		goto	90	! Error
	if (NPR+NSC .ge. ICON-1)	goto	44
	NB2=NB1+NSC+1
	J2 = NPOS(ILNG-NSC-1,STRI(NB2),':')-1
	if (J2 .lt. 0)	 goto	90
	if (J2 .eq. 0)	 goto	44

C Eq command line has a non-zero parameter after 2nd ':' (TE:EQ:...;)
	if (JEQ.ne.2 .and. JEQ.ne.3)	goto	90
	if (IVBUF+J2+12 .gt. 5000)	then
	   write(7,*)'"',(STRI(J),J=NB,NB+ICON-1),'"'
	   write(*,*)'"',(STRI(J),J=NB,NB+ICON-1),'"'
	   goto	93
	endif

!	write(*,*)'"',(STRI(j),j=NB2,NB2+j2-1),'"',j2
	if (IFNUM(STRI(NB2),j2) .eq. 1)	then
!	   LOUT = j2+12
!	   if (JEQ .eq. 2) call	APPBUF(LOUT,SP6//'GN2E='//
!     +		STRIN(NB2:NB2+j2-1),CNSBUF,ICBUF,2000)
!	   if (JEQ .eq. 3) call	APPBUF(LOUT,SP6//'GN2I='//
!     +		STRIN(NB2:NB2+j2-1),CNSBUF,ICBUF,2000)
!	   if (LOUT .eq. 0)	goto	94
	else
!	   KEYOUT=-1
!	   STR80(1:) = char(j2)//STRIN(NB2:NB2+j2-1)
!	   call	ANLFML(-1,STR,LINTXT,LINFOR) ! line '\'[...
!C	   write(*,*)'for:  "',(LINFOR(1+j2),j2=1,ichar(LINFOR(1))),'"'
!C	   write(*,*)'txt:  "',(LINTXT(1+j2),j2=1,ichar(LINTXT(1))),'"'
!	   KEYOUT=0
!	   LUR  = ichar(LINFOR(1))
!	   LOUT = LUR+j2+12
!	   if (JEQ .eq. 2) call	APPBUF(LOUT,SP6//'GN2E='//
!     +		STRFOR(1:LUR),DTVBUF,IVBUF,5000)
!	   if (JEQ .eq. 3) call	APPBUF(LOUT,SP6//'GN2I='//
!     +		STRFOR(1:LUR),DTVBUF,IVBUF,5000)
!	   if (LOUT .eq. 0)	goto	93
!	   LOUT=ichar(LINTXT(1))
!	   do	j=1,LOUT
!	      LINTXT(j) = LINTXT(j+1)
!	   enddo	   
	endif

	call COPY(6,SP6,DTVBUF(IVBUF+1))
	if (JEQ .eq. 2) call COPY(5,'GN2E=',DTVBUF(IVBUF+7))
	if (JEQ .eq. 3) call COPY(5,'GN2I=',DTVBUF(IVBUF+7))
	call COPY(J2,STRI(NB2),DTVBUF(IVBUF+12))
	write(DTVBUF(IVBUF),132)char(J2+11)
C write model.txt
	NCOEF(JEQ)=NCOEF(JEQ)+1
	call COPY(J2+5,DTVBUF(IVBUF+7),COEFLN(NCOEF(JEQ),JEQ))
C	write(*,*)" Equation",'  "',(STRI(J),J=NB,NB+1),
C     >		'":   Factor',' "', (STRI(J),J=NB2,NB2-1+J2),'"'
	IVBUF = IVBUF+J2+12
C	ITYPE=3

 44	continue
	if (ITYPE .lt.-1)	goto	90	! -> Error 
	if (ITYPE .ne. 2)	goto	45
C       if (JEQ .eq. 2) write(*,*)"Set GN2E=2.5"
C       if (JEQ .eq. 3) write(*,*)"Set GN2I=2.5"
	if (IVBUF+15 .gt. 5000)	then
	   write(7,*)'"',(STRI(J),J=NB,NB+ICON-1),'"'
	   write(*,*)'"',(STRI(J),J=NB,NB+ICON-1),'"'
	   goto	93
	endif
	call COPY(6,SP6,DTVBUF(IVBUF+1))
	if (JEQ .eq. 2) call COPY(8,'GN2E=2.5',DTVBUF(IVBUF+7))
	if (JEQ .eq. 3) call COPY(8,'GN2I=2.5',DTVBUF(IVBUF+7))
	write(DTVBUF(IVBUF),'(1A1)')char(14)
C write model.txt
	NCOEF(JEQ)=NCOEF(JEQ)+1
	call COPY(8,DTVBUF(IVBUF+7),COEFLN(NCOEF(JEQ),JEQ))
C       write(*,*)" Equation",'  "',(STRI(J),J=NB,NB+1),
C     >		'":   Factor',' "', (STRI(J),J=NB2,NB2-1+3),'"'
	IVBUF = IVBUF+15

 45	continue
C	write(*,*)"ITYPE =",ITYPE
	LEQ(JEQ)=ITYPE

C Write equation number to EQBUF(IEBUF) when the equation is not used
	if(IEBUF.gt.1999)	goto	95
	write(EQBUF(IEBUF),  132) char(129+4*JEQ)
	write(EQBUF(IEBUF+1),132) char(0)
	IEBUF=IEBUF+1
				goto	80
C--------------- Subroutine string analysis ---------------------------v
 46	continue
C External subroutine
C NSBR - ordinal subroutine number
C LEN1 - abs. pos. of the total string STRIN(NB:NB+LEN1-1)
C	 write(*,*)"---------------------------------------------------"
C	 write(*,*)"Input sbr string" 
C	 write(*,*)'"',STRIN(NB:NB+LEN1-1),'"'
C LSTR - meaningful string (till closing ";")
C ICON - meaningful string (till closing "!")
C	    STRIN(NB:NB+LSTR-1) or STRIN(NB:NB+ICON-1), respectively
C	    write(*,*)'"',(STRI(NB+J1-1),J1=1,LSTR),'"'
C	    write(*,*)'"',(STRI(NB+J1-1),J1=1,ICON),'"'
C NPR  - subroutine header (till 1st ":")	STRIN(NB:NB+NPR-2)
C	write(*,*)'"',(STRI(NB+J1-1),J1=1,NPR-1),'"',STRIN(NB:NB+LSTR-1),'"'
Command string analysis for a subroutine
	if (NSBR .ge. NSBMX)	then
	   NSC=NPOS(ICON-NPR,STRI(NB),':')-2
	   write(7,'(A,1I3)')'>>> Subroutine number >',NSBMX
	   write(*,'(A,1I3)')'>>> Subroutine number >',NSBMX
	   write(7,*)'>>> Subroutine "',
     >		(STRI(j),j=NB,NB+NSC),'" call ignored'
	   write(*,*)'>>> Subroutine "',
     >		(STRI(j),j=NB,NB+NSC),'" call ignored'
	   goto	80
	endif
	NAME6 = STRIN(NB:NB+NPR-2)
	call	UPCASE(6,NAME6)
	NSBR = NSBR+1

C LOCSBR() =-4 Default value (the sbr is never used in the model)
C LOCSBR() =-3 SBP, tag "<" used, to be called from detvar.tmp
C LOCSBR() =-2 SBP, no tag used,  to be called from init.inc and eqns.inc
C LOCSBR() =-1 SBR, tag "<" used, to be called from detvar.tmp
C LOCSBR() = 0 SBR, no tag used,  to be called  from init.inc and eqns.inc
C LOCSBR() > 0 SBR, tag ">" used. Ordinal number of a subroutine tagged
	LOCSBR(NSBR) = 0
	jst = NPOS(LSTR,STRI(NB),'&')-1
	if (jst .lt. ICON)	then		! Tag '&' found
	   do j=NB+jst,NB+ICON-1
	      STRI(j) = STRI(j+1)
	   enddo
	   jsh = jsh+1
	   NPR = NPR-1
	   ICON = ICON-1
	   LSTR = LSTR-1
	   LOCSBR(NSBR) = -2
	endif

	jlb = NPOS(LSTR,STRI(NB),'<')-1
C	if (jlb.lt.ICON .and. jst.lt.ICON)	then
C	   write(*,*)
C     >	  '>>> ERROR >>> Simultaneous use of "&" and "<" is not allowed'
C	   call exit(1)
C	endif
	if (jlb .lt. ICON)	then		! Tag '<' found
C	   write(*,*)'jlb =',jlb,NB,NPR,NPR-NB   
	   do j=NB+jlb,NB+ICON-1
	      STRI(j) = STRI(j+1)
	   enddo
	   jsh = jsh+1
	   if (jlb .lt. NPR-NB)	NPR = NPR-1
	   ICON = ICON-1
	   LSTR = LSTR-1
	   LOCSBR(NSBR) = LOCSBR(NSBR)-1
C	   if (jrb .gt. jlb)	jrb = jrb-1
	endif
	jrb = NPOS(LSTR,STRI(NB),'>')-1
	if (jrb.lt.ICON .and. jst.lt.ICON)	write(*,*)
     >	  '>>> ERROR >>> Simultaneous use of "&" and ">" is not allowed'
	if (jrb .lt. ICON)	then		! Tag '>' found
	   do j=NB+jrb,NB+ICON-1
	      STRI(j) = STRI(j+1)
	   enddo
	   jsh = jsh+1
	   if (jrb .lt. NPR-NB)	NPR = NPR-1
	   ICON = ICON-1
	   LSTR = LSTR-1
	   LOCSBR(NSBR) = LOCSBR(NSBR)+1
	endif
C	write(*,*)"string to process" 
C	write(*,*)'"',STRIN(NB:NB+LSTR-1),'"',LSTR,ICON

	if (jlb .le. ICON-1 .and. jrb .le. ICON-1)	then
	   write(7,*)'>>> Error in line:  "',STRIN(NB:NB+ICON-1),'"'
	   write(7,*)'    Contradictive pointers ignored'
	   write(*,*)'>>> Error in line:  "',STRIN(NB:NB+ICON-1),'"'
	   write(*,*)'    Contradictive pointers ignored'
	   LOCSBR(NSBR) = 0
	endif
	if (IEBUF+NPR.gt.1999)	goto	95

C	write(*,*)'Input:  "',STRIN(NB:NB+NPR-1),'"'
	J = NPOS(NPR-1,STRI(NB),'(')-1		! Store SBR name
	J1= NPOS(NPR-1,STRI(NB),':')-1		! SBR without parameters
	j = min(j,j1)				! j - SBR/SBP name length
	J2 = NPOS(NPR-1,STRI(NB),'/')
	if (J2 .gt. j)	then			! Ignore "/"
	   J2 = 0				! in the SBR
	   goto	461				! parameter list
	endif
	J2 = RNPOS(NPR-1,STRI(NB),'/')
	if (J2 .lt. NPR-1)	then
	   if (LOCSBR(NSBR) .gt. -2)	goto	89	! "/" in sbr
C	   write(*,*)'SBRBUF: "',STRIN(NB+J2:NB+NPR-2),'"',j2,NPR-1
C	   write(*,*)'SBPATH: "',STRIN(NB:NB+J2-1),'"'
	   ! Cut path and write it to a separate place
C	   if (IPBUF+J2 .gt. 1000)	goto	88
C	   SBPATH(IPBUF)=char(J2)
C	   call COPY(J2,STRI(NB),SBPATH(IPBUF+1))
C	   IPBUF = IPBUF+J2 
	endif
	if (J2 .eq. NPR-1)	goto	89		! "/" at end
	if (J2 .eq. NPR)	J2 = 0

 461	continue
	SBRBUF(JSB)=char(J)
	if (j .ge. 25)	goto	88
	call COPY(J,STRI(NB),SBRBUF(JSB+1))
	JSB=JSB+J+1
	call STOBUF(NPR-1-J2,STRI(NB+J2),IEBUF,EQBUF(IEBUF))
	EQBUF(IEBUF) = char(0)
	if (NAME6.eq."MIXINT")	LOCSBR(NSBR) = 1
	if (NAME6.eq."MIXEXT")  LOCSBR(NSBR) = 1
	if (NAME6.eq."TSCTRL")  LOCSBR(NSBR) = 1
C	write(*,*)NAME6,NSBR,'  ',STRIN(NB:NB+j-1),LOCSBR(NSBR)

	ILNG=ICON-NPR
 	if (ILNG .le. 0)	goto	80
C	write(SCNS,'(1I2)')NSBR
	NB1=NB+NPR
	J  = NPOS(ILNG,STRIN(NB1:),',')-1
	J1 = NPOS(ILNG,STRIN(NB1:),':')-1
	NSC = min(J,J1)
C	NSC = NPOS(ILNG,STRIN(NB1:),':')-1

C Subroutine time step DTEQ(1, )
	if (NSC .eq. 0)			goto 47		! Default value
C	write(*,*)(STRI(j2),j2=NB1,NB1+NSC-1),NB1,NSC
C	write(*,*)STRIN(NB1:NB1+NSC-1),NB1,NSC,IFNUM(STRI(NB1),NSC)
	call	NUM2ST(NSBR,SCNS,LCNS)
	if (IFNUM(STRI(NB1),NSC) .eq. 1)	then
	   LOUT = NSC+LCNS+15
	   call	APPBUF(LOUT,SP6//'DTEQ(1,'//SCNS(1:LCNS)//')='//
     +		STRIN(NB1:NB1+NSC-1),CNSBUF,ICBUF,2000)
	   if (LOUT .eq. 0)	goto	94
	   STSTEP(NSBR,1)=' '
	   call	COPY(NSC,STRI(NB1),STSTEP(NSBR,1))
	else
	   KEYOUT = -1
	   STR80(1:) = char(NSC)//STRIN(NB1:NB1+NSC)
	   call	ANLFML(-1,STR,LINTXT,LINFOR) ! line '\' ...]
C	write(*,*)'LINFOR:   "',(LINFOR(1+j),j=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j),j=1,ichar(LINTXT(1))),'"'
	   KEYOUT = 0
	   LUR = ichar(LINFOR(1))
	   LOUT = LUR+LCNS+15
	   call	APPBUF(LOUT,SP6//'DTEQ(1,'//SCNS(1:LCNS)//')='//
     +		STRFOR(1:LUR),DTVBUF,IVBUF,5000)
	   if (LOUT .eq. 0)	goto	93
	   LOUT=ichar(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
	endif

 47	if (NSC .lt. J1)	then
	   write(*,*)(STRI(NB1+NSC+1+J2),J2=0,J1-NSC)
	   write(*,*)NB1
     >    ,'   "',STRIN(NB1+NSC+1:NB1+J1-1),'"',STRI(NB1+NSC+1)
	   write(*,*)STSTEP(NSBR,5)
	   NSC = J1
	endif
	LEN2=ILNG-NSC-1
 	if(LEN2.le.0)		goto	80
C Subroutine start time
	NB2=NB1+NSC+1
	NSC=NPOS(LEN2,STRI(NB2),':')-1
	if(NSC.eq.0)		goto	48
C	if(ICBUF+NSC+18.gt.2000)	goto	94
C	call COPY(13,'      DTEQ(2,',CNSBUF(ICBUF+1))
C	call COPY(2,SCNS,CNSBUF(ICBUF+14))
C	call COPY(2,')=',CNSBUF(ICBUF+16))
C	call COPY(NSC,STRI(NB2),CNSBUF(ICBUF+18))
C	STSTEP(NSBR,2)=' '
C	call COPY(NSC,STRI(NB2),STSTEP(NSBR,2))
C	write(CNSBUF(ICBUF),132)	char(NSC+17)
C	ICBUF=ICBUF+NSC+18
	call	NUM2ST(NSBR,SCNS,LCNS)
	if (IFNUM(STRI(NB2),NSC) .eq. 1)	then
	   LOUT = NSC+LCNS+15
	   call	APPBUF(LOUT,SP6//'DTEQ(2,'//SCNS(1:LCNS)//')='//
     +		STRIN(NB2:NB2+NSC-1),CNSBUF,ICBUF,2000)
	   if (LOUT .eq. 0)	goto	94
	   STSTEP(NSBR,2)=' '
	   call	COPY(NSC,STRI(NB1),STSTEP(NSBR,2))
	else
	   KEYOUT = -1
	   STR80(1:) = char(NSC)//STRIN(NB2:NB2+NSC)
	   call	ANLFML(-1,STR,LINTXT,LINFOR) ! line '\' ...]
C	write(*,*)'LINFOR:   "',(LINFOR(1+j),j=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j),j=1,ichar(LINTXT(1))),'"'
	   KEYOUT = 0
	   LUR = ichar(LINFOR(1))
	   LOUT = LUR+LCNS+15
	   call	APPBUF(LOUT,SP6//'DTEQ(2,'//SCNS(1:LCNS)//')='//
     +		STRFOR(1:LUR),DTVBUF,IVBUF,5000)
	   if (LOUT .eq. 0)	goto	93
	   LOUT=ichar(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
	endif
 48	LEN3=LEN2-NSC-1
 	if(LEN3.le.0)		goto	80
C Subroutine end time
	NB3=NB2+NSC+1
	NSC=NPOS(LEN3,STRI(NB3),':')-1
	if(NSC.eq.0)		goto	49
C	if(ICBUF+NSC+18.gt.2000)	goto	94
C	call COPY(13,'      DTEQ(3,',CNSBUF(ICBUF+1))
C	call COPY(2,SCNS,CNSBUF(ICBUF+14))
C	call COPY(2,')=',CNSBUF(ICBUF+16))
C	call COPY(NSC,STRI(NB3),CNSBUF(ICBUF+18))
C	STSTEP(NSBR,3)=' '
C	call COPY(NSC,STRI(NB3),STSTEP(NSBR,3))
C	write(CNSBUF(ICBUF),132)	char(NSC+17)
C	ICBUF=ICBUF+NSC+18
	call	NUM2ST(NSBR,SCNS,LCNS)
	if (IFNUM(STRI(NB3),NSC) .eq. 1)	then
	   LOUT = NSC+LCNS+15
	   call	APPBUF(LOUT,SP6//'DTEQ(3,'//SCNS(1:LCNS)//')='//
     +		STRIN(NB3:NB3+NSC-1),CNSBUF,ICBUF,2000)
	   if (LOUT .eq. 0)	goto	94
	   STSTEP(NSBR,3)=' '
	   call	COPY(NSC,STRI(NB1),STSTEP(NSBR,3))
	else
	   KEYOUT = -1
	   STR80(1:) = char(NSC)//STRIN(NB3:NB3+NSC)
	   call	ANLFML(-1,STR,LINTXT,LINFOR) ! line '\' ...]
C	write(*,*)'LINFOR:   "',(LINFOR(1+j),j=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j),j=1,ichar(LINTXT(1))),'"'
	   KEYOUT = 0
	   LUR = ichar(LINFOR(1))
	   LOUT = LUR+LCNS+15
	   call	APPBUF(LOUT,SP6//'DTEQ(3,'//SCNS(1:LCNS)//')='//
     +		STRFOR(1:LUR),DTVBUF,IVBUF,5000)
	   if (LOUT .eq. 0)	goto	93
	   LOUT=ichar(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
	endif
 49	LEN4=LEN3-NSC-1
 	if(LEN4.le.0)		goto	80
	NB4=NB3+NSC+1
	NSC=NPOS(LEN4,STRI(NB4),':')-1
	if(NSC.ne.1)		goto	80
C Subroutine call ASCII code
	call	UPCASE(1,STRI(NB4))
	KK=STRI(NB4)
	if (KK.eq.'C')	then
	write(7,*)'>>> Warning: Illegal subroutine call code ',KK
	write(*,*)'>>> Warning: Illegal subroutine call code ',KK
				goto	80
			endif
	if(ICBUF+21.gt.2000)	goto	94
	call COPY(13,'      DTEQ(4,',CNSBUF(ICBUF+1))
	call COPY(2,SCNS,CNSBUF(ICBUF+14))
	call COPY(2,')=',CNSBUF(ICBUF+16))
	write(ALKEY,'(1I3)')ichar(STRI(NB4))-64
	call COPY(3,ALKEY,CNSBUF(ICBUF+18))
	STSTEP(NSBR,4)=STRI(NB4)
	write(CNSBUF(ICBUF),132)	char(20)
	ICBUF=ICBUF+21
	goto	80
 50	continue
	write(*,*)'"',STRIN(NB:NB+LSTR),'"'  ! Symbol "=<" or ":=" found
	goto	80
 60	continue
	write(*,*)'"',STRIN(NB:NB+LSTR),'"'  ! Symbol "&" found
	goto	80
 70	continue
	write(*,*)'"',STRIN(NB:NB+LSTR),'"'  ! Symbol "=<" or ":=" found
	goto	80
C----------------------------------------------------------------------|
 80	NB=NB+LSTR+jsh
	jsh = 1
C	if (LEN1-NB .ge. 0) write(*,'(3A,3I5)')'Rest of the string "',
C     >	   STRIN(NB:LEN1),'"   Length, Start, End ',LEN1-NB,NB,LEN1
	if(LEN1-NB) 1,2,2
C End of file "model.tmp"
 81	if (NCHI .ne. 0) then
	   close(NCHI)
	   NCH = NCHI-1
	   NCHI = 0
	   goto	1
	endif
	close(NCH)
	close(3)

	if (NBEG(LCV).gt.0 .and. NBEG(LMV).gt.0)	then
		write(*,*)">>> WARNING: CV assignment is overridden,",
     >			  "  MV assignment is used"
		write(7,*)">>> WARNING: CV assignment is overridden,",
     >			  "  MV assignment is used"
C Suppress writing "CU(J)=..." into "tmp/detvar.tmp"
		NBEG(LCV)=0
	endif
	goto	99

 90	write(7,91)(STRIN(NB:NB+LSTR-1))
	write(*,91)(STRIN(NB:NB+LSTR-1))
 91	format(' >>> Error in the command line  "',A,'"')
	IFSYNT = 1
	goto	80
 88	write(7,*)'>>> Error in line:  "',STRIN(NB:NB+ICON-1),'"'
	write(7,*)'    The subroutine name is too long'
	write(*,*)'>>> Error in line:  "',STRIN(NB:NB+ICON-1),'"'
	write(*,*)'    The subroutine name is too long'
	IFSYNT = 1
	goto	99
 89	write(7,*)'>>> Error in line:  "',STRIN(NB:NB+ICON-1),'"'
	write(7,*)'    Symbol "/" found in the subroutine header'
	write(*,*)'>>> Error in line:  "',STRIN(NB:NB+ICON-1),'"'
	write(*,*)'    Symbol "/" found in the subroutine header'
	IFSYNT = 1
	goto	99
 92	write(7,*)'>>> Temporary buffer length',IBUF,' > 10000'
	write(*,*)'>>> Temporary buffer length > 10000'
	IFSYNT = 1
	goto	99
 93	write(7,*)'>>> Constant buffer length',IVBUF,' > 5000'
	write(*,*)'>>> Constant buffer length > 5000'
	IFSYNT = 1
	goto	99
 94	write(7,*)'>>> Constant buffer length',ICBUF,' > 2000'
	write(*,*)'>>> Constant buffer length > 2000'
	IFSYNT = 1
	goto	99
 95	write(7,*)'>>> Equation buffer length',IEBUF,' > 2000'
	write(*,*)'>>> Equation/subroutine buffer length > 2000'
	IFSYNT = 1
	goto	99
 96	write(7,*)'>>> Model analyser: model ',MFNAME(2)
	write(7,*)'    Recurrent "#include" is not allowed'
	write(7,*)'    Include file ',INCNAM
	write(*,*)'>>> Model analyser: model ',MFNAME(2)
	write(*,*)'    Recurrent "#include" is not allowed'
	write(*,*)'    Include file ',INCNAM
	write(*,*)
	IFSYNT = 1
	goto	99
 97	write(7,*)'>>> Model analyser: Cannot open file ',
     >			"equ/"//INCNAM(1:length(80,INCNAM(1:)))
	write(*,*)'>>> Model analyser: Cannot open file ',
     >			"equ/"//INCNAM(1:length(80,INCNAM(1:)))
	write(*,*)
	IFSYNT = 1
	goto	99
 98	write(7,*)'>>> Model analyser: Incompatible description of ',
     >			'e-i equipartition' 
	write(*,*)'>>> Model analyser: Incompatible description of ',
     >			'e-i equipartition' 
	write(*,*)
	IFSYNT = 1
	goto	99
 99	if (IFSYNT .ne. 0)	then
		write(*,*)'>>> Model analyzer STOP: Syntax error'
		call	exit(1)
	endif
C----------------------------------------------------------------------|
C Print buffer example:
C	write(*,*)"DTVBUF:"
C	do	j=1,NTMP
C	if (NBEG(j).ne.0)	then
C		write(*,*)j,TMPNAM(j),NBEG(j)
C		call	showbu(TMPBUF(NBEG(j)))
C		call	showbu(DTVBUF(NBEG(j)))
C	endif
C	enddo
C Print entire buffers
C	call	SHOWBU(TMPBUF)
C	call	SHOWBU(DTVBUF)
C	write(*,'(//20("-"),"SBRBUF",20("-")/)')
C	call	SHOWBU(SBRBUF)
C	write(*,*)(LOCSBR(j),j=1,NSBR)
C	write(*,'(//20("-"),"EQBUF",21("-")/)')
C	call	SHOWBU(EQBUF)		! eqn ID & sbr line
C	do	j=1,7
C	   write(*,*)"Eqn ID:",j
C	   do	j1=1,NCOEF(j)		! No. of coefs + in. + bnd.
C	      write(*,*)COEFLN(j1,j)	! eqn coefficients for model.txt
C	   enddo
C	enddo
C	write(*,*)(ARRNAM(j2),j2=1,NARR)
C	
C----------------------------------------------------------------------|
	if (NROUT .gt. NRW)	NROUT = NRW
	if (NXOUT .gt. NRW)	NXOUT = NRW
	if (NTOUT .gt. NRW)	NTOUT = NRW
C	write(*,'(3A,2I5)')('"',TMPNAM(j),'"',NBEG(j),j,j=31,35)
C	write(*,'(3A,2I5)')('"',TMPNAM(j),'"',NBEG(j),j,j=41,44)
	do  51	j=1,NEQNS
	   if (LEQ(j) .ge. 0)	goto	51	! Eqn is defined
Check if any of NE,TE,TI,CU,MU,F(*) is defined as initial value 
	   if (j .eq. 1 .and. NBEG(LNE) .eq. 0)		goto	51
	   if (j .eq. 2 .and. NBEG(LTE) .eq. 0)		goto	51
	   if (j .eq. 3 .and. NBEG(LTI) .eq. 0)		goto	51
	   if (j .eq. 6 .and. NBEG(LNI) .eq. 0)		goto	51
	   if (j .ge. 5 .and. j .le. 9)			goto	51
	   if (j .gt. 9 .and. NBEG(LF(j-10)) .eq. 0)	goto	51
	   if (j .eq. 4 .and. NBEG(LCU)+NBEG(LMU).eq.0)	goto	51
	   LEQ(j) = 0
	   if(IEBUF.gt.1999)	goto	95
	   write(EQBUF(IEBUF),  132) char(129+4*j)
	   write(EQBUF(IEBUF+1),132) char(0)
	   IEBUF=IEBUF+1
 51	continue
	if (LEQ(2).gt.8 .and. LEQ(3).le.8)	goto	98
	if (LEQ(3).gt.8 .and. LEQ(2).le.8)	goto	98

C	J = IONEOF(NAME,NTMP,TMPNAM)
C	if( J.eq.LPE.or.J.eq.LPET.or.J.eq.LPI.or.J.eq.LPIT)	then
C	   write(*,*)J,LPE,LPET,LPI,LPIT
C	endif
	if (IFPEI.eq.1 .and. LEQ(2).gt.0 .and. LEQ(3).gt.0)	then
	   inquire(file='tmp/impei.msg',exist=JPEI)
	   if (JPEI) then
	      call system('tput bold')
	      write(*,'(A,$)') "               "
	      call system('tput smso')
	      write(*,*) ">>>>>>>>>>>>>> Warning <<<<<<<<<<<<<< "
	      call system('tput sgr0')
	      call system('cat tmp/impei.msg')
	   endif
	endif

C	do	J1=1,NARRS
C	   write(*,'(I3,60A)')j1,'  ',(NAMARR(J,J1),J=1,LNAMAR(J1))
C	enddo
C	call	SHOWBU(TMPBUF)	       ! Detailed print of the entire buffer
C	call	APPTMP(TMPBUF,6)	! Print of the entire buffer
C	call	APPTMP(TMPBUF(NBEG(LROE)),6)
C	call	APPTMP(TMPBUF(NBEG(LHE)),6)	! Print of the specified line
C	write(*,*)"CNSBUF"
C	call	APPTMP(CNSBUF,6)
C	write(*,*)"DTVBUF"
C	call	APPTMP(DTVBUF,6)
	NCH = 2
C	write(*,*)"SETVAR:"
	call	SETVAR(CNSBUF,NCH)
C	write(*,*)"DETVAR:"
	call	DETVAR(NCH,NARRS,LNAMAR,NAMARR,ARRSTR,EQBUF,LOCSBR)
C	write(*,*)"INIVAR:"
	call	INIVAR(NCH,LEQ)
C	write(*,*)"APPTXT:"
	call	APPTXT(PRFINI,LEQ)
C	write(*,*)"INIPSI"
C	call	INIPSI(NCH,NBEG(LCU),NBEG(LMU))
C	write(*,*)"EQUFTN:"
C	call	EQUFTN(NCH,LEQ,NCOEF,EQBUF,COEFLN,STSTEP,LOCSBR)
	call	ALLDEF(NCH,LEQ,EQBUF,LOCSBR)
C	write(*,*)"init.inc", "eqns.inc"
C
C	write(*,*)
C	write(*,*)"Before writing XARNAM",JARX
C	do	j=1,JARX
C	   if (ARXNAM(j) .gt. 0 .and. ARXNAM(j) .le. NARR)
C     &			write(*,*)ARRNAM(ARXNAM(j)),ARXNAM(j)
C	   if (ARXNAM(j) .ge. NARR)	then
C	      write(*,*)ARRNAM(ARXNAM(j)-NARR),ARXNAM(j),ARXNAM(j)-NARR
C	     endif
C	enddo
C	write(*,*)"ININAM:"
Contents of the array ARXNAM which is used in XARNAM is changed
C		by calling RADOUT, TIMOUT. Present location of the call
C		means that the warning is suppressed if an X-array 
C       	is used for output only
	call	ININAM(NCH,NSBR,SBRBUF,NBEG(LMV),NTOUT,NROUT,NXOUT,NWX,
     +		NAMET,NAMER,NAMEX,SCALET,SCALER,IFTOUT,IFROUT,LOCSBR,
     +		LEQ,NEQNS)
C	write(*,*)"RADOUT:"
	call	RADOUT(NCH,NROUT,FORMR,SCALER,NAMER)
C	write(*,*)"TIMOUT:"
	call	TIMOUT(NCH,NTOUT,FORMT,SCALET,NAMET)
C	write(*,*)"Done"
	close(7)
C	write(*,*)NBEG(LNI),LEQ(6),NBEG(LNE),LEQ(1)
C	write(*,*)NBEG(LTE),NBEG(LTI)
C	call	APPTMP(TMPBUF(NBEG(LNE)),6)
C	call	APPTMP(TMPBUF(NBEG(LNI)),6)
C	call	APPTMP(TMPBUF(NBEG(LTE)),6)
C	call	APPTMP(TMPBUF(NBEG(LTI)),6)
	end
C======================================================================|
C End of the model analyzer run!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C The following buffers are defined:
C For all parts:
C	NBEG(L_var_name)=Address_in_TMPBUF
C	NBEG(LTE) = 0 if "TE:..." is absent 
C		  > 0 if "TE:..." is found
C	MFNAME()
C	include nambuf.inc
C	KEYOUT
C For equftn.tmp: (equations)
C	EQBUF(IEBUF <= 1000) One position is used to transfer Eq. ID
C	LEQ(JEQ)=ITYPE(-1<=ITYPE<=3)
C	LEQ(JEQ) =-1 no equation
C		 = 0 type: AS
C		 = 1 type: EQ (default)
C		 = 2 type: FU
C		 = 3 type: heat conductivity flux + heat convection
C	LEQ(6) Definition of NI
C	LEQ(6) =-1 default. No definition for NI
C	LEQ(6) = 0 Old type definition for NI
C	LEQ(6) = 1 NI is assigned and has to be placed in eqns.inc
C	LEQ(5) separate treatment 
C	LEQ(7-9) not used
C
C	EQBUF(IEBUF <= 1000) Subroutine_call_string
C For detvar.tmp:
C	DTVBUF(IVBUF < 5000)
C	NARRS
C	ARRSTR(NARRS<
C	NAMARR(J,J1<LUR)
C	LNAMAR(J1<LUR)
C	NTMP
C	include nambuf.inc
C	TMPBUF(IBUF <= 10000)	- Contents of the do loop in detvar.tmp
C		Includes rhs for the command lines where lhs are eqn coef-ts
C		all formula proccessing results according to the ASTRA rules
C				  
C For inivar.tmp:
C	TMPBUF(IBUF <= 10000)	- Initial profiles for TE, TI, NE, Fj
C For inipsi.tmp:	
C	NCH, NBEG(LCU) and NBEG(LMU) are used only
C For ininam.tmp:
C	NSBR	- number of subroutines called form a model
C	CNSBUF(ICBUF<=2000)
C	NROUT
C	SCALER
C	NAMER
C	NTOUT
C	SCALET
C	NAMET
C	NXOUT
C	NAMEX
C	NWX
C	SBRBUF	- buffer for subroutine names
C For model.txt:
C	COEFLN(NCOEF(J),J)) rhs, tr.coefs, initial conditions
C	PRFINI initial conditions or assignments for TE, TI, NE, Fj, CU, MU
C	NCOEF
C ENTEQU ===================== Names preparing ========================|
	subroutine ENTEQU
C----------------------------------------------------------------------|
C Fill arrays
C   ARRNAM	!  Arrays   from for/status.inc between AMAIN and ZMAIN
C   MAIVAR	! Variables from for/const.inc between AB and ZRD9
C   CONNAM	! Variables from for/const.inc between CF1 and CSOL4
C   SRVNAM	! Variables from for/const.inc between DROUT and ...
C   FMLNAM	! 
C   FNCNAM 	! 
C   TMPNAM is defined in blockdata
C----------------------------------------------------------------------|
	implicit none
	include '.srv/nambuf.inc'
	integer	i,j,IERR,NCH,IARR,IRET
	character LFNAME(4)*40,NAME*6,STRING*132
	data	LFNAME/	'for/status.inc','for/const.inc',
     &			'tmp/fml.dir','tmp/fnc.dir'/
C----------------------------------------------------------------------|
	NCH = 1
C Array list
	IARR=0
	NARR=0
	call	OPENRD(NCH,LFNAME(1),0,IERR)
	if(IERR.gt.0)	then
	   write(*,*)'MODEL (ENTEQU): Can not open file ',LFNAME(1)
	   pause
	endif
 21	NARR=NARR+IARR
	if(NARR.gt.NCVA)		then
		write(7,*)'>>> Error:   Total array number >',NCVA
		write(*,*)'>>> Error:   Total array number >',NCVA
		pause
		return
				endif
	read(NCH,102,ERR=41,END=41)NAME
 102	format(2X,1A6)
	if(NAME.eq.'End ar')	goto 20
	if(IARR.gt.0)	ARRNAM(NARR)=NAME
C	write(*,*)NARR,'"',NAME,'"',ARRNAM(NARR),'"'
	if(NAME.eq.'Array ')	IARR=1
			goto 21
 20	close(NCH)
	NARR=NARR-1
C Main variable list
	call	OPENRD(NCH,LFNAME(2),0,IERR)
	if(IERR.gt.0)	then
		write(*,*)'MODEL (ENTEQU): Can not open file ',LFNAME(2)
		pause
			endif
	IARR=0
	NMAIVA=0
 22	NMAIVA=NMAIVA+IARR
	if(NMAIVA.gt.NCVA)	then
		write(7,*)'>>> ERROR:   Total variable number >',NCVA
		write(*,*)'>>> ERROR:   Total variable number >',NCVA
		pause
		return
				endif
	read(NCH,102,ERR=42,END=42)NAME
	if(NAME.eq.'End va')	goto 23
	if(IARR.gt.0)	MAIVAR(NMAIVA)=NAME
	if(NAME.eq.'Variab ')	IARR=1
			goto 22
C Constant list
 23	NMAIVA=NMAIVA-1
	IARR=0
	NCONS=0
 24	NCONS=NCONS+IARR
	if(NCONS.gt.NCVA)	then
		write(7,*)'>>> ERROR:   Total constant number >',NCVA
		write(*,*)'>>> ERROR:   Total constant number >',NCVA
		pause
		return
				endif
	read(NCH,102,ERR=43,END=43)NAME
	if(NAME.eq.'End co')	goto 25
	if(IARR.gt.0)	CONNAM(NCONS)=NAME
	if(NAME.eq.'Consta')	IARR=1
			goto 24
C Internal variable list
 25	NCONS=NCONS-1
	IARR=0
	NSRV=0
 26	NSRV=NSRV+IARR
	if(NSRV.gt.NCVA)	then
	    write(7,*)'>>> ERROR:   Internal variable number >',NCVA
	    write(*,*)'>>> ERROR:   Internal variable number >',NCVA
	    call	exit(1)
	    return
				endif
	read(NCH,102,ERR=44,END=44)NAME
	if(NAME.eq.'End in')	goto 27
	if(IARR.gt.0)	SRVNAM(NSRV)=NAME
	if(NAME.eq.'Intern')	IARR=1
			goto 26
 27	NSRV=NSRV-1
	close(NCH)

C Standard function list
C	goto	30	! Disable check 
	call	OPENRD(NCH,"for/stdfun.f",0,IERR)
	if(IERR.gt.0)	then
	   write(*,*)' >>> Error: Can not open file "for/stdfun.f"'
	   stop
	endif
	NSTDF=0
 28	continue
C	if (NSTDF .gt. NCVA)	then
C	   write(7,*)'>>> ERROR:   Internal variable number >',NCVA
C	   write(*,*)'>>> ERROR:   Internal variable number >',NCVA
C	   call	exit(1)
C	   return
C	endif
	read(NCH,'(A)',ERR=45,END=29)STRING
	j = index(STRING,"double")
	if (j .eq. 0)	goto	28
	i = j
	j = index(STRING(i:),"precision")
	if (j .eq. 0)	goto	28
	i = i+j-1
	j = index(STRING(i:),"function")
	if (j .eq. 0)	goto	28
	i = i+j-1
	j = index(STRING(i:),' ')
	if (j .eq. 0)	goto	28
	i = i+j
	j = index(STRING(i:),'(')
	if (j .eq. 0)	goto	28
C	write(*,*)'"',STRING(i:i+j-2),'"'
	NSTDF=NSTDF+1
	STDF(NSTDF)=STRING(i:i+j-2)
	call	UPCASE(6,STDF(NSTDF))
	if (j .eq. 0)	goto	45
	goto	28
 29	close(NCH)
C	write(*,'(I5/1(15A)))')NSTDF,('"',STDF(j),'"',j=1,NSTDF)
 30	continue

C 500 is a maximal number of FMLNAM elements
C 6   is a maximal length of FMLNAM element 
C nfml equal to the total formula number is returned
	IERR = 0
	NFML = 500
	I = 6
	call	getfml(IERR,NFML,FMLNAM,I)
C	write(*,'(3a)')('"',FMLNAM(j),'"',j=1,NFML)
C	write(*,*)(FMLNAM(j),j=1,NFML)
	if(IERR .lt. 0)	  call	exit(1)
	do	j=1,NFML
	    call	UPCASE(I,FMLNAM(j))
	enddo
	IERR = 0
	NFNC = 200
	I = 6
	call	getfnc(IERR,NFNC,FNCNAM,I)
C	write(*,'(3a)')('"',FNCNAM(j),'"',j=1,NFNC)
C	write(*,*)(FNCNAM(j),j=1,NFNC)
	if(IERR .lt. 0)	  call	exit(1)
	do	j=1,NFNC
	    call	UPCASE(I,FNCNAM(j))
	enddo
C----------------------------------------------------------------------|
C Formula list: read file tmp/fml.dir
C	NFML=0
C	call	OPENRD(NCH,LFNAME(3),0,IERR)
C	if(IERR.gt.0)	then
C		write(*,*)'MODEL (ENTEQU): Can not open file ',LFNAME(3)
C		pause
C			endif
C 1	NFML=NFML+1
C	if(NFML.gt.500)		then
C		write(7,*)'>>> ERROR:   Total FML number > 500'
C		write(*,*)'>>> ERROR:   Total FML number > 500'
C		pause
C		return
C				endif
C	read(NCH,101,ERR=3,END=3)FMLNAM(NFML)
C	call	UPCASE(6,FMLNAM(NFML))
C				goto 1
C 3	continue
C	write(*,*)'"',FMLNAM(NFML-1),'"'
C	close(NCH)
C	NFML=NFML-1
C----------------------------------------------------------------------|
C Function list: read file tmp/fnc.dir
C For IBM:
C	STRING='ls -1 /afs/ipp/home/g/grp/Alpha/astra/fnc'//
C     >		' | grep ".f"$ | cut -d"." -f1 > ~/temp/fnc.dir'
C For SUN:
C	STRING='ls -1 /afs/ipp/home/g/grp/Alpha/astra/fnc'//
C     >		' | cut -d"." -f1 | fgrep -v "%" | fgrep -v "~"'//
C     > 		' > ~/temp/fnc.dir'
C The next line gives
C	j = system(STRING)
C Note: IEEE floating-point exception flags raised: 
C    Invalid Operation; 
C See the Numerical Computation Guide, ieee_flags(3M) 
C	NFNC=0
C	call	OPENRD(NCH,LFNAME(4),0,IERR)
C	if(IERR.gt.0)	then
C		write(*,*)'MODEL (ENTEQU): Can not open file ',LFNAME(4)
C		pause
C			endif
C 4	NFNC=NFNC+1
C	if(NFNC.gt.200)		then
C		write(7,*)'>>> ERROR:   Total FNC number > 200'
C		write(*,*)'>>> ERROR:   Total FNC number > 200'
C		pause
C		return
C				endif
C 101	format(13X,1A6)
C	read(NCH,101,ERR=5,END=5)FNCNAM(NFNC)
C	call	UPCASE(6,FNCNAM(NFNC))
C				goto 4
C 5	close(NCH)
C	NFNC=NFNC-1
C----------------------------------------------------------------------|
C Function & formula list
 	I=1
	J=1
	NFIN=0
 6	if(I.gt.NFML.or.J.gt.NFNC)	goto 9
	if(FMLNAM(I).ne.FNCNAM(J))	goto 7
	NFIN=NFIN+1
	if (NFIN .gt. NFINMX)		then
	 write(7,*)'>>> ERROR:   FNC and FML same names number >',NFINMX
	 write(*,*)'>>> ERROR:   FNC and FML same names number >',NFINMX
	   pause
	   return
	endif
	FINNAM(NFIN)=FMLNAM(I)
	I=I+1
	J=J+1
	goto 6
 7	if(FMLNAM(I).gt.FNCNAM(J))	goto 8
	I=I+1
	goto 6
 8	J=J+1
	goto 6
 9	call AORDER(NARR,ARRNAM,IRET)
	if(IRET.eq.0)	goto 32
 	write(*,112)ARRNAM(IRET),ARRNAM(IRET+1)
 	write(7,112)ARRNAM(IRET),ARRNAM(IRET+1)
 112	format('  >>> Warning: Array list out of alphabet order: '/
     >	'  >>> ',1A6,' before ',1A6/'  >>> Check file status.inc')
 32	call AORDER(NFML,FMLNAM,IRET)
	if(IRET.eq.0)	goto 33
 	write(*,113)FMLNAM(IRET),FMLNAM(IRET+1)
 	write(7,113)FMLNAM(IRET),FMLNAM(IRET+1)
 113	format('  >>> Warning: Formula list out of alphabet order: '/
     >	'  >>> ',1A6,' before ',1A6/'  >>> Check file tmp/declar.fml')
 33	call AORDER(NFNC,FNCNAM,IRET)
	if(IRET.eq.0)	goto 34
 	write(*,114)FNCNAM(IRET),FNCNAM(IRET+1)
 	write(7,114)FNCNAM(IRET),FNCNAM(IRET+1)
 114	format('  >>> Warning: Function list out of alphabet order: '/
     >	'  >>> ',1A6,' before ',1A6/'  >>> Check file tmp/declar.fnc')
 34	goto	35
C	call AORDER(NTMP,TMPNAM,IRET)
C	if(IRET.eq.0)	goto 35
C	write(*,115)TMPNAM(IRET),TMPNAM(IRET+1)
C	write(7,115)TMPNAM(IRET),TMPNAM(IRET+1)
 115  format('  >>> Warning: Coefficients list out of alphabet order:'
     >	/'  >>> ',1A6,' before ',1A6/
     >	'  >>> Check BLOCKDATA in the file .srv/model2.for')
 35	return
 41	write(*,*)' >>> Error in the file for/status.inc'
	write(*,*)'     Check array list'
	goto	45
 42	write(*,*)' >>> Error in the file for/const.inc'
	write(*,*)'     Check global variable list'
	goto	45
 43	write(*,*)' >>> Error in the file for/const.inc'
	write(*,*)'     Check constant list'
	goto	45
 44	write(*,*)' >>> Error in file for/status.inc'
	write(*,*)'     Check internal variable list'
 45	stop
	end
C======================================================================|
	subroutine	NOCHANGE(J,DTVBUF,IVBUF)
	implicit none
	integer		j,jj,IVBUF
	character	SCNS*2,DTVBUF(*)*1
C  Strings are appended to "detvar.tmp" which forbid assignment
C  of the main Astra variable NAME from the data file 
C  if the variable is defined by a model
C	J=IONEOF(NAME,NMAIVA,MAIVAR)
C	NAME=MAIVAR(J)
C 	if(IFDFVX(J).le.2)	IFDFVX(J)=2
C 	if(IFDFVX(J).le.2)
C     >  IPL=.6
C
C left: variable
	JJ = 1
	if(J.lt.0)	JJ = -1
	J = abs(J)
	write(SCNS,'(1I2)')J
C string:	if(IFDFVX(J).le.2)	IFDFVX(J)=2
	call COPY(16,'      if(IFDFVX(',DTVBUF(IVBUF+1))
	call COPY(2,SCNS,DTVBUF(IVBUF+17))
	call COPY(8,').lt.2)	',DTVBUF(IVBUF+19))
	call COPY(7,'IFDFVX(',DTVBUF(IVBUF+27))
	call COPY(2,SCNS,DTVBUF(IVBUF+34))
	call COPY(3,')=2',DTVBUF(IVBUF+36))
	write(DTVBUF(IVBUF),'(1A1)')	char(38)
	IVBUF=IVBUF+39
	if(IVBUF.gt.4900)	then
		write(7,*)'>>> Constant buffer length',IVBUF,'>4900'
		pause '>>> Constant buffer length >4900'
				endif
	IVBUF=MIN0(IVBUF,4900)
	if (JJ.lt.0)	return
C strings:	if(IFDFVX(J).le.2)
C	     >	'NAME'='LINFOR'
	call COPY(16,'      if(IFDFVX(',DTVBUF(IVBUF+1))
	call COPY(2,SCNS,DTVBUF(IVBUF+17))
	call COPY(7,').le.2)',DTVBUF(IVBUF+19))
	write(DTVBUF(IVBUF),'(1A1)')	char(25)
	IVBUF=IVBUF+26
	if(IVBUF.gt.4900)	then
		write(7,*)'>>> Constant buffer length',IVBUF,'>4900'
		pause ' >>> Constant buffer length >4900'
				endif
	IVBUF=MIN0(IVBUF,4900)
	call COPY(8,'     >  ',DTVBUF(IVBUF+1))
	end
C======================================================================|
	subroutine	APPBUF (LSTRI,STRI,XBUF,IBUF,IBMAX)
C Example:	call	APPBUF(LSTRI,STRI,CNSBUF,ICBUF,5000)
c		call	APPBUF(NSC,STRI,DTVBUF,IVBUF,5000)
	implicit none
	integer		LSTRI,IBUF,IBMAX,LNUM,j
	character*1	SNUM*3,STRI(*),XBUF(*)
	j = LSTRI+1
	if (IBUF+j .gt. IBMAX)	then
	   LSTRI = 0
	   return
	endif
	call COPY(LSTRI,STRI(1),XBUF(IBUF+1))
	write(XBUF(IBUF),132)	char(LSTRI)
	IBUF=IBUF+j
	return
 132	format(132A1)
	end
C======================================================================|
	subroutine	NUM2ST (NUM,SNUM,LNUM)
C Input:   NUM  - integer
C Output:  SNUM - string
C          LNUM - length of SNUM
	implicit none
	integer		NUM,LNUM
	character*3	SNUM
	if     (NUM .gt. 999)	then
	   write(*,*)" >>> NUM2ST >>> ERROR: too large number"
	   LNUM = 0
	elseif (NUM .gt. 99)	then
	   LNUM = 3
	   write(SNUM,'(1I3)')NUM
	elseif (NUM .gt. 9)	then
	   LNUM = 2
	   write(SNUM,'(1I2)')NUM
	elseif (NUM .lt. 1)	then
	   write(*,*)" >>> NUM2ST >>> ERROR: too small number"
	   LNUM = 0
	else
	   LNUM = 1
	   write(SNUM,'(1I1)')NUM
	endif
	end
C======================================================================|
