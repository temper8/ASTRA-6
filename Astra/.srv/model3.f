C ================ Service subroutines ====================
C ASSTR: subroutine ASSTR(LEN,STRIN,STROUT)
C	Concatenates string: 
C		STROUT(1:LEN)=STRIN(1)//STRIN(2)//...
C				 NAME[LEN+1:6]=' '
C SETNAM: subroutine SETNAM(LEN,STRI,NAME)
C	Appends NAME with spaces NAME[1:LEN]=STRI[1:LEN]
C				 NAME[LEN+1:6]=' '
C BLEX:   subroutine BLEX(LEN,STRING)
C	Exclude BLANCS and TAB from STRING[1:LEN], LEN=LENGTH(STRING)
C COPY:   subroutine COPY(LEN,STFROM,STTO)
C	Sets STTO[1:LEN]=STFROM[1:LEN]
C IONEOF: integer function IONEOF(NAME,N,BUF)
C	BUF in arbitrary order
C JONEOF: integer function JONEOF(NAME,N,BUF)
C	BUF in alphabetical order
C 	JONEOF(or IONEOF)=Position of NAME in BUF
C			 =0 if NAME is not in BUF
C NPOS:   integer function NPOS(LEN,STRING,SYM)
C	NPOS=First occurence of SYM position in STRING(1:LEN), else LEN+1
C RNPOS:  integer function RNPOS(LEN,STRING,SYM)
C	RNPOS=Last occurence of SYM position in STRING(1:LEN), else LEN+1
C ERASYM:  integer function ERASYM(LEN,STRING,SYM)
C       Find and remove the 1st ocurrence of the symbol SYM in STRING(1:LEN).
C	ERASYM=Number of SYM position in STRING(1:LEN), else ERASYM=0
C LINSUM: subroutine LINSUM(LINE,STRI,LEN)
C	Assigns LEN charcters of STRI[1:LEN] to the end of LINE
C		i.e.	LINE[LINE(1)+1:LINE(1)+LEN]=STRI[1:LEN]
C	and changes the length stored in LINE(1)=LINE(1)+LEN
C NDEL:   integer function NDEL(LEN,STRING)
C	NDEL=Number of delimiter position  + - * : ) ( ,  in STRING,
C 		else NDEL=LEN+1
C LENGTH: integer function LENGTH(LEN,STRING)
C	LENGTH=Length of STRING without ' ' & 0 at the end
C LENGTF: integer function LENGTF(STRI)
C	LENGTF=Position No. of the first appearance ' ' or <Tab> or \0 in STRI
C APPFML: subroutine APPFML(NCH,LENS,STRI)
C	Appends formula 'STRI(1:LENS)' to channel NCH
C APPTMP: subroutine APPTMP(BUF,NCH)
C	Appends BUF to channel NCH
C Example: 	writing "klm" to standard output (the same as SHOWBU(BUF))
C	call APPTMP(char(3)//'klm',6)
C TMPFML: subroutine TMPFML(LENS,STRI,IB,BUF)
C	Appends formula 'STRI(1:LENS)' to BUF(1:IB)
C STOBUF: subroutine STOBUF(LEN,STR,IB,BUF)
C	Appends STR to BUF (BUF(1)=LEN, BUF[2,LEN+1]=STR[1:LEN])
C	IB is increased by LEN+1
C AORDER: subroutine AORDER(NEN,ARR,IRET)
C	Checks if ARR[6,1:NEN] in alphabetical order
C IFNUM:   integer function IFNUM1(STR,N)
C	Check if STR[1:N] is a number
C UPCASE: subroutine UPCASE(LEN,STRI)
C DNCASE: subroutine DNCASE(LEN,STRI)
C OPENRD:
C OPENWT: subroutines	OPENRD,OPENWT(NCHAN,FILNAM,ITYPE,IERR)
C	input:	NCHAN - channel number
C		FILNAM - file name (character variable or string CALL 
C		opepator). Extention (name.ext) is necessary for string.
C		ITYPE - format status (0 - default, > 0 - 'UNFORMATTED')
C	output:	IERR - output status (0 - O'kay, 1 - warning, 2 - error)
C ASSTR ============================================
	subroutine ASSTR(LEN,STRIN,STROUT)
	implicit	none
	character	STROUT*(*)
	character*1	STRIN(*)
	integer		LEN,J
	do 1	J=1,LEN
	STROUT(J:J)=STRIN(J)
 1	continue
	STROUT(LEN+1:LEN+1)=char(0)
	end
C SETNM6 ============================================
	subroutine SETNM6(LEN,STRI,NAME)
	implicit	none
	character*1 NAME*6,STRI(*)
	integer LEN,J
	if (LEN .le. 0)	return
	do 1	J=1,6
	if(J.gt.LEN)		goto 2
	NAME(J:J)=STRI(J)
				goto 1
 2	NAME(J:J)=' '
 1	continue
	end
C SETNAM ============================================
	subroutine SETNAM(LEN,STRI,NAME)
	implicit	none
	character*1 NAME(6),STRI(*)
	integer LEN,J
	do 1	J=1,6
	if(J.gt.LEN)		goto 2
	NAME(J)=STRI(J)
				goto 1
 2	NAME(J)=' '
 1	continue
	end
C BLEX	=============================================
	subroutine BLEX(LEN,STRI)
	implicit	none
	character*1 STRI(*),S,TAB,TZER
	integer KEY,J,L,LEN
	J = 0
	TZER= char(J)
	J = 9
	TAB  = char(J)
	L=0
	do 1	J=1,LEN
	S=STRI(J)
	KEY=ICHAR(STRI(J))
	if(ICHAR(S).eq.0)		goto 2
	if(S.eq.' '.or.S.eq.TAB)	goto 1
	L=L+1
C	if(KEY.ge.97.and.KEY.le.122)	KEY=KEY-32
C	STRI(L)=CHAR(KEY)
	STRI(L)=STRI(J)
 1	continue
 2	LEN=L
	write(STRI(L+1),'(1A1)')TZER
	end
C COPY	=============================================
	subroutine COPY(LEN,STFROM,STTO)
	implicit	none
	character*1 STFROM(*),STTO(*)
	integer J,LEN
	do 1	J=1,LEN
	STTO(J)=STFROM(J)
 1	continue
	end
C IN6EOF ============================================
C The same as IONEOF with 1st argument char*6
	integer function IN6EOF(NAME,N,BUF)
	implicit	none
	character*6 NAME*6,BUF(*)
	integer J,N
	do	1	J=1,N
	if(NAME.eq.BUF(J))	then
		IN6EOF=J
		return
				endif
1	continue
	IN6EOF=0
	end
C IONEOF ============================================
	integer function IONEOF(NAME,N,BUF)
	implicit	none
	character*6 NAME6*6,BUF(*),NAME(6)*1
	integer J,N
	write(NAME6,'(6A1)')(NAME(j),j=1,6)
	do	1	J=1,N
	if(NAME6.eq.BUF(J))	then
		IONEOF=J
		return
				endif
1	continue
	IONEOF=0
	end
C JN6EOF ============================================
C The same as JONEOF with 1st argument char*6
	integer function JN6EOF(NAME,N,BUF)
	implicit	none
	character*6 NAME*6,BUF(*)
	integer JE,JJ,JB,N
	JB=1
	JE=N
	JN6EOF=0
 1	if(JE.lt.JB)		return
	JJ=(JB+JE)/2
	if(NAME.gt.BUF(JJ))	goto 4
	if(NAME.lt.BUF(JJ))	goto 3
 	JN6EOF=JJ
	return
 3	JE=JJ-1
				goto 1
 4	JB=JJ+1
				goto 1
	end
C JONEOF ============================================
	integer function JONEOF(NAME,N,BUF)
	implicit	none
	character*6 NAME6*6,BUF(*),NAME(6)*1
	integer JE,JJ,JB,N
	write(NAME6,'(6A1)')(NAME(jj),jj=1,6)
	JB=1
	JE=N
	JONEOF=0
 1	if(JE.lt.JB)		return
	JJ=(JB+JE)/2
	if(NAME6.gt.BUF(JJ))	goto 4
	if(NAME6.lt.BUF(JJ))	goto 3
 	JONEOF=JJ
	return
 3	JE=JJ-1
				goto 1
 4	JB=JJ+1
				goto 1
	end
C NPOS	=============================================
	integer function NPOS(LEN,STRI,SYM)
	implicit	none
	character*1 STRI(*),SYM
	integer J,LEN
	do 1	J=1,LEN
	if(STRI(J).ne.SYM)	goto 1
	NPOS=J
	return
 1	continue
	NPOS=LEN+1
	end
C RNPOS	=============================================
	integer function RNPOS(LEN,STRI,SYM)
	implicit	none
	character*1 STRI(*),SYM
	integer J,LEN
	do 1	J=LEN,1,-1
	if(STRI(J).ne.SYM)	goto 1
	RNPOS=J
	return
 1	continue
	RNPOS=LEN+1
	end
C ERASYM ============================================
	integer function ERASYM(LEN,STRI,SYM)
	implicit	none
	character*1 STRI(*),SYM
	integer J,LEN
	ERASYM = 0
	if (LEN .le. 0)	return
	do  1	J=1,LEN
	   if(STRI(J) .ne. SYM)	goto 1
	   ERASYM = J
	   goto	2
 1	continue
	return
 2	continue
	if (ERASYM .eq. LEN)	return
	do  3	J=ERASYM,LEN-1
	   STRI(J) = STRI(J+1)
 3	continue
	end
C LINSUM ============================================
C	Assigns LEN charcters of STRI[1:LEN] to the end of LINE
C		i.e.	LINE[LINE(1)+1:LINE(1)+LEN]=STRI[1:LEN]
C	and changes the length stored in LINE(1)=LINE(1)+LEN
	subroutine LINSUM(LINE,STR,LEN)
	implicit	none
	character*1 LINE(133),STR(*)
	integer J,LEN,N
	if (LEN .le. 0)	return
	N=ICHAR(LINE(1))
	do 1	J=1,LEN
	LINE(N+J+1)=STR(J)
 1	continue
	write(LINE(1),'(1A1)')	char(N+LEN)
	end
C NDEL	=============================================
	integer function NDEL(LEN,STRI)
	implicit	none
	character*1 STRI(*),SYM(7)
	integer J,JJ,LEN
	data SYM/'+','-','(',')','*','/',','/
	do 1	J=1,LEN
	do 1	JJ=1,7
	if(STRI(J).ne.SYM(JJ))	goto 1
	NDEL=J
	return
 1	continue
	NDEL=LEN+1
	end
C LENGTH ============================================
	integer function LENGTH(LEN,STRI)
	implicit	none
	character*1 STRI(*)
	integer J,LEN,NLIN
	do 4	J=1,LEN
	NLIN=LEN-J+1
	if(ICHAR(STRI(NLIN)).ne.0.and.STRI(NLIN).ne.' ') goto 5
 4	continue
 5	LENGTH=NLIN
	end
C LENGTF ============================================
C----------------------------------------------------------------------|
	integer	function LENGTF(STRI)
	implicit	none
	character*1	STRI(*),TZER,TAB
	integer	j
	j = 0
	TZER= char(j)
	j = 9
	TAB = char(j)
	LENGTF	= 0
	do	j=1,132
	if(STRI(j).eq.' '.or.STRI(j).eq.TZER.or.STRI(j).eq.TAB) goto 5
	LENGTF	= j
	enddo
 5	continue
	end
C SHOWBU ============================================
	subroutine SHOWBU(BUF)
C----------------------------------------------------------------------|
Called as "call SHOWBU(BUF)" prints out the entire buffer
C----------------------------------------------------------------------|
	implicit	none
	character*1	BUF(*)
	integer 	NB,M,J
	NB=1
 1	if (ichar(BUF(NB)).gt.129) then		! When printing EQBUF
	   write(*,*)"Equation No.",(ichar(BUF(NB))-129)/4,"is used"
	   NB=NB+1				! skip transport eqns
	   goto	1				! ID records: 
	endif					! character*1 each
	M=ICHAR(BUF(NB))+NB
C print the buffer until a zero length word encountered
	if(M.eq.NB)write(*,*)
     >	"                   Buffer ends at the position = ",NB
	if(M.eq.NB)		return
	write(*,'(20X,7(A,I4))')"Field = ",NB," ->",M,
     >		",      Length = ",M-NB,"  String: (",NB+1," ->",M,")"
	write(*,101)'"',(BUF(J),J=NB+1,M),'"'
 101	format(79A1)
	NB=M+1
	goto 1
	end
C APPTMP ==============================================================|
	subroutine  APPTMP(BUF,NCH)
	implicit    none
	integer	    NCH,M,NB,j,M1,LENS,NPOS
	character*1 BUF(*),STR*72,STR1*72,STRI*7
	NB = 1
 1	M  = ICHAR(BUF(NB))
	if (M.eq.0)		return

	M1 = NPOS(M,BUF(NB+1),'/')
	if (M1.lt.12 .or. M1.ge.M .or. NCH.eq.6)	goto	2
	if (BUF(NB+M1-1).ne.'l' .and. BUF(NB+M1-1).ne.'L') goto 2
	do	j=1,M
		STR(j:j)  = BUF(NB+j)
		STR1(j:j) = BUF(NB+j)
	enddo
	j = M1
	call	BLEX  (j,STR1)
	if (j.lt.12)	goto	2
	call	DNCASE(j,STR1)
	if (STR1(1:12) .eq. "include'fml/")	then
	   LENS = M-M1-1
	   if (LENS .gt. 6)	goto	99
	   STRI = STR(M1+1:M-1)//char(0)
	   call	APPFML(NCH,LENS,STRI)
	   NB = NB+M+1
	   goto 1
	endif

 2	if (M.gt.61)	then
	   write(NCH,'(61A1)')(BUF(J),J=NB+1,NB+61)
	   NB = NB+61
	   M = M-61
	else
	   write(NCH,'(61A1)')(BUF(J),J=NB+1,NB+M)
	   NB = NB+M+1
	   goto	1
	endif
 3	continue
	if (M.gt.60)	then
	   write(NCH,'(5X,62A1)')'+',' ',(BUF(J),J=NB+1,NB+60)
	   NB = NB+60	   
	   M = M-60
	   goto	3
	else
	   write(NCH,'(5X,62A1)')'+',' ',(BUF(J),J=NB+1,NB+M)
	   NB = NB+M+1
	   goto	1
	endif
 99	write(*,*)
     >	   ' >>> Warning: formula name "',STR(M1+1:M-1),'" too long'
	end
C======================================================================|
	subroutine APPFML(NCH,LENS,STRI)
C----------------------------------------------------------------------|
C NCH > 0 - write formula to the logical unit NCH (must be open)
C LENS - length of STRI
	implicit	none
	character*6 STRI,STR1,ST(79)*1,NREAD*40,ST79*79,MARKER*79,
     >		    OUTNAM(500)
	integer LENGTH,LENGTF,noblan,IONEOF
	integer LENS,NCH,J,J1,NL,NCHFML,MINCH,IERR
	include	'.srv/nambuf.inc'	! KEYOUT,NOUT,Tflag are used
	data 	OUTNAM/500*'      '/	NREAD/'fml/'/	MINCH/8/
	save	OUTNAM,NREAD,MINCH
C	write(*,*)"APPFML Input:",LENS,'   "',STRI,'"',NCH
	if (Tflag .ne. 0)	write(MARKER,*)
     >	'      call add2loc("Formula ',STRI(1:LENS),'"//char(0))'
C	j = index(MARKER,'char(0)')+7
C	write(*,*)'"',MARKER(1:40+LENS),'"',j
C     >,40+LENS,'"',MARKER(j-2:j),'"'
	call	DNCASE(LENS,STRI)
	NREAD(5:) = STRI(1:LENS)
	STR1 = STRI(1:LENS)
	NCHFML = MINCH
C	write(*,*)"NOUT =",NOUT
 1	NL = LENGTH(40,NREAD)
	if (KEYOUT .le. 0)	goto	2
C By-pass repeated formula 
	J1 = IONEOF(STR1,NOUT,OUTNAM)
	if (J1 .gt .0)	then
C	   write(*,*)"Requested ",NREAD(1:NL),",	include level"
C     +		,NCHFML-MINCH,"   ID:",J1,"   Skipped"
	   if (NCHFML.eq.MINCH)	return
	   if (NCHFML.gt.MINCH)	goto	7
	else
	   NOUT=NOUT+1
	   OUTNAM(NOUT)=STR1
C	   write(*,*)"Requested ",NREAD(1:NL),",	include level"
C     +		,NCHFML-MINCH,"   ID:",NOUT,"	Read unit ",NCHFML
	endif
 2	continue
	call	OPENRD(NCHFML,NREAD,0,IERR)
	if(IERR .gt. 0)	then
		write(*,*)'MODEL (APPFML):  Cannot open file: ',NREAD
		call	exit(1)
			endif
C	if (NCHFML.eq.MINCH .and. Tflag.ne.0)
C     +   write(*,*)'"',MARKER(1:index(MARKER,'char(0)')+7),'"'
	if (NCHFML.eq.MINCH .and. Tflag.ne.0)
     +  	write(NCH,*)MARKER(1:index(MARKER,'char(0)')+7)
 3	read(NCHFML,'(79A1)',ERR=6,END=6)ST
	if(ST(1).eq.'C'.or.ST(1).eq.'c'.or.ST(1).eq.'!')  goto	3
	NL = LENGTH(79,ST)
	if (NL.le.2)	goto	3
C write included formula
C	write(*,'(1I5,79A1)')NL,(ST(J),J=1,NL)
	write(ST79,'(79A1)')(ST(J),J=1,NL)
	call	DNCASE(NL,ST79)
	J = noblan(ST79)
	if (j .eq. 0)			goto	5
	if (ST79(J:J+6) .ne. 'include')	goto	5
	J = J+7+noblan(ST79(J+7:))
	if (ST79(J:J+3) .ne. 'fml/')	goto	5
	J1 = J+2+LENGTF(ST79(j+4:))
	NREAD(5:) = ST79(j+4:J1)
	STR1 = ST79(j+4:J1)
	NCHFML = NCHFML+1
					goto	1
 5	write(NCH,'(79A1)')(ST(J),J=1,NL)
	goto	3
 6	close(NCHFML)
C	write(*,*)"				Closed:",NCHFML
	goto	(8,7)	NCHFML-MINCH+1
 7 	NCHFML = NCHFML-1
 	goto	3
 8	if (NCHFML .lt. MINCH)
     .		write(*,*)">>> APPFML include formula error"
	return
	end
C TMPFML ============================================
	subroutine TMPFML(LENS,STRI,IB,BUF)
	implicit	none
	character*1 STRI*6,BUF(*),NREAD*40,ST(79)
	integer LENS,IB,J,NL,NCHFML,MINCH,NB,LENGTH,IERR
	data	NREAD/'fml/'/	MINCH/8/
C	NREAD(5:4+LENS)=STRI(1:LENS)
	call	DNCASE(LENS,STRI)
	NREAD(5:) = STRI(1:LENS)
	NCHFML = MINCH
	call	OPENRD(NCHFML,NREAD,0,IERR)
	if(IERR.gt.0)	then
		write(*,*)'>> MODEL (TMPFML): Cannot open file ',NREAD
		pause 
			endif
	NB=1
 3	read(NCHFML,'(79A1)',ERR=6,END=6)ST
	if(ST(1).eq.'C'.or.ST(1).eq.'c'.or.ST(1).eq.'!')  goto	3
	NL=LENGTH(79,ST)
C	if (NL.le.1)	goto	3
C append formula to buffer
C	write(*,'(1I5,79A1)')NL,(ST(J),J=1,NL)
 5	write(BUF(NB),'(79A1)')	char(NL)
	do	J=1,NL
	    BUF(J+NB)=ST(J)
	enddo
	NB=NB+NL+1
	IB=IB+NL+1
	goto	3
 6	close(NCHFML)
	end
C INCFML ============================================
	subroutine INCFML(LENS,STRI,IB,BUF)
	implicit	none
	character BUF(*)
	character STRI*6,NREAD*40,ST(79)*1,ST79*79
	integer LENS,IB,J,NL,NCHFML,IERR
	equivalence	(ST,ST79)
	data	NREAD/'fml/'/	NCHFML/8/
C	NREAD(5:4+LENS)=STRI(1:LENS)
	call	DNCASE(LENS,STRI)
	NREAD(5:) = STRI(1:LENS)
	call	OPENRD(NCHFML,NREAD,0,IERR)
	if(IERR.gt.0)	then
		write(*,*)'>> MODEL (INCFML): Cannot open file ',NREAD
		pause 
			endif
	close(NCHFML)
	ST79(1:) ="      include '"//NREAD(1:4+LENS)//"'"
	NL = 20+LENS
 5	write(BUF(1),'(79A1)')	char(NL)
	do	J=1,NL
	    BUF(J+1)=ST(J)
	enddo
	IB=IB+NL+1
	end
C=======================================================================
C	The function returns the position of the 1st non-blanck (space 
C		or tabulation) in the string STRI
C		from the beginning till '\0' or <CR>
C	or  0	if no blancks found
C-----------------------------------------------------------------------
	integer	function noblan(STRI)
	implicit	none
	character*1 STRI(*),TZER,TCR,TAB
	integer	j
	j = 0
	TZER=char(j)
	j = 9
	TAB=char(j)
	j = 13
	TCR=char(j)
	noblan	= 0
	do	j=1,132
	    if (STRI(j).eq.TZER .or. STRI(j).eq.TCR)	return
	    if (STRI(j).ne.' ' .and. STRI(j).ne.TAB)	then
		noblan = j
		return
	    endif
	enddo
	end
C STOBUF ============================================
	subroutine STOBUF(LEN,STR,IB,BUF)
	implicit	none
	character*1 STR(*),BUF(*)
	integer J,LEN,IB
	write(BUF(1),'(1A1)')	char(LEN)
	do 1	J=1,LEN
 1	BUF(J+1)=STR(J)
	IB=IB+LEN+1
	end
C AORDER ============================================
	subroutine AORDER(NEN,ARR,IRET)
	implicit	none
	character*6 ARR(*)
	integer JJ,NEN,IRET
	IRET=0
	do 1	JJ=1,NEN-1
C	   write(*,*)ARR(jj)(1:6),jj,'  "',ARR(JJ)(1:1),'"'
C	   if(ARR(JJ)(1:1).eq.'T'.or.ARR(JJ)(1:1).eq.'S')	then
C	      if(ARR(JJ).lt.ARR(JJ+1))	write(*,*)ARR(jj)(1:6),jj
C	   endif
	if(ARR(JJ).lt.ARR(JJ+1)) goto 1
	IRET=JJ
 1	continue
	end
C IFNUM ==============================================
	integer function IFNUM(STR,N)
Check if str(1:N) is a decimal number
C returns 1 true
C         0 false
	implicit	none
	character*1 STR(*),NU(17)
	integer I,J,N,NDOT,NFLO
	data NU/'0','1','2','3','4','5','6','7','8','9','+','-'
     +	,'E','e','D','d','.'/
	NDOT = 0
	NFLO = 0
	i = 0
	do	j=11,12
	   if (STR(1) .eq. NU(j))	i = 1
	enddo

 10	continue
	i = i+1
	if (i .gt. N)		goto	20
	do	j=1,10
	   if (STR(i) .eq. NU(j))	goto 10
	enddo
C is not a digit:
	if (STR(i) .eq. NU(17))	then
	   if (NDOT .ne. 0)	goto	77
	   NDOT = i
	   goto	10
	endif
C is not a digit, not a dot:
	if (NFLO .ne. 0)	goto	77
	do	j=13,16
	   if (STR(i) .eq. NU(j))	NFLO = i
	enddo
	if (NFLO .eq. 1)	goto	77
	if (NFLO .eq. N)	goto	77
	if (STR(i+1).eq.NU(11) .or. STR(i+1).eq.NU(12))	i = i+1
	if    (i .eq. N)	goto	77
	if (NFLO .ne. 0)	goto	10
	goto	77
 20	continue
	if (NFLO.gt.0 .and. NFLO.lt.NDOT)	goto	77
C the next line requires usage of dot if the exponent is present
C omitting this line allows numbers as 1d3
	if (NFLO.gt.0 .and.    0.eq.NDOT)	goto	77
	if (NDOT .eq. 1)	then
	   if (NFLO .eq. 2)	goto	77
	   if (N    .eq. 1)	goto	77
	endif
	IFNUM = 1
	return
 77	IFNUM = 0
	end
C UPCASE ===========================================================
	subroutine UPCASE(LEN,STRI)
	implicit	none
	integer KEY,J,LEN
	character*1 STRI(*)
	do	1	J=1,LEN
	KEY=ICHAR(STRI(J))
	if(KEY.ge.97.and.KEY.le.122)	KEY=KEY-32
	STRI(J)=CHAR(KEY)
1	continue
	end
C DNCASE ===========================================================
	subroutine DNCASE(LEN,STRI)
	implicit	none
	integer KEY,J,LEN
	character*1 STRI(*)
	do	1	J=1,LEN
	KEY=ICHAR(STRI(J))
	if (KEY.lt.65 .or. KEY.gt.90)	goto	1
C	if(KEY.ge.65.and.KEY.le.90)	KEY=KEY+32
	KEY=KEY+32
	STRI(J) = CHAR(KEY)
1	continue
	end
C LINDEX ===========================================================
C	integer function INDX(STRI,STR,L)
	subroutine INIDX(STRI,STR,L)
	implicit	none
	integer J,I,L,LEN,indx
	character*1 STRI(*),STR(L)
	LEN = ichar(STRI(1))
C	write(*,*)LEN,L
C	write(*,*)(STRI(j),j=2,LEN-1)
	indx = LEN
	do	J=0,LEN-L
	   do 1	i=1,L
	      if (STR(i) .ne. STRI(1+i+j))	goto	2
C	      write(*,'(5A1)')'"',STR(i),'"',STRI(1+i+j),'"'
 1	   continue
C	   goto	2
	   write(*,'(i3,132A1)')2+j,'  ',(STRI(1+i+j),i=1,L)
C	      indx = 2+j
C	      return
 2	   continue
	enddo
	end
c OPENRD ===============================================================
	subroutine	OPENRD(NCHAN,FILNAM,ITYPE,IERR)
	implicit	none
	integer 	NCHAN,ITYPE,IERR
	character*(*)	FILNAM
	if(ITYPE.eq.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1)
	rewind(NCHAN)
	if(ITYPE.gt.0)	open(UNIT=NCHAN,FILE=FILNAM,STATUS='OLD',ERR=1,
     *			FORM='UNFORMATTED')
	rewind(NCHAN)
	IERR=0
	return
1	IERR=2
	end
C OPENWT ===============================================================
	subroutine	OPENWT(NCHAN,FILNAM,ITYPE,IERR)
	implicit	none
	integer 	NCHAN,ITYPE,IERR
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
C======================================================================|
	block data
	implicit	none
	include 	'.srv/nambuf.inc'
	include 	'.srv/tmpbuf.inc'
	integer		j
	data	NFML/0/ NFNC/0/ IBUF/1/ IVBUF/1/ LINAPP/NEQNS*0/
C Equation coefficients
C The order of data in this DATA statement should be consistent with
C the COMMON block description in the file "tmpbuf.inc"
	data	(TMPNAM(j),j=1,100) /
     &	'NE', 'DN', 'HN', 'XN', 'CN', 'SN', 'SNN', 'NEB', 'QNB', 'QNNB',
     1	'TE', 'DE', 'HE', 'XE', 'CE', 'PE', 'PET', 'TEB', 'QEB', 'QETB',
     2	'TI', 'DI', 'HI', 'XI', 'CI', 'PI', 'PIT', 'TIB', 'QIB', 'QITB',
     3	'CU', 'DC', 'HC', 'XC', 'CV', 'CD', 'CC',  'IPL', 'LEXT','UEXT',
     4	'MU', 'MV','CUBS','RON','ROE','ROI','DVN', 'DSN', 'DVE', 'DSE',
     5	'DVI','DSI','NI', 'DU2','DU3','DU4','DU5', 'DU6', 'DU7', 'DU8',
     6	'RO0','RO1','RO2','RO3','RO4','RO5','RO6', 'RO7', 'RO8','RO9',
     7	'F0', 'F1', 'F2', 'F3', 'F4', 'F5', 'F6',  'F7',  'F8', 'F9',
     8	'DF0','DF1','DF2','DF3','DF4','DF5','DF6', 'DF7', 'DF8','DF9',
     9	'VF0','VF1','VF2','VF3','VF4','VF5','VF6', 'VF7', 'VF8','VF9'/
      	data	(TMPNAM(j),j=101,NTMP) / 
     &	'SF0','SF1','SF2','SF3','SF4','SF5','SF6', 'SF7', 'SF8','SF9',
     1	'SFF0','SFF1','SFF2','SFF3','SFF4','SFF5','SFF6','SFF7','SFF8',
     1	'SFF9',
     2	'GF0','GF1','GF2','GF3','GF4','GF5','GF6','GF7','GF8','GF9',
     3	'F0B','F1B','F2B','F3B','F4B','F5B','F6B','F7B','F8B','F9B',
     4	'QF0B','QF1B','QF2B','QF3B','QF4B','QF5B','QF6B','QF7B','QF8B',
     4	'QF9B',
     5	'QFF0B','QFF1B','QFF2B','QFF3B','QFF4B','QFF5B','QFF6B','QFF7B',
     5	'QFF8B','QFF9B',
     6	'DVF0','DVF1','DVF2','DVF3','DVF4','DVF5','DVF6','DVF7','DVF8',
     6	'DVF9',
     7	'DSF0','DSF1','DSF2','DSF3','DSF4','DSF5','DSF6','DSF7','DSF8',
     7	'DSF9'/

C LTMPBUF overlays the COMMON block "TMPFIL" and is used to set
C	  its elements
	data	(LTMPBUF(j),j=1,100) /
     &		 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
     1		11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
     2		21, 22, 23, 24, 25, 26, 27, 28, 29, 30,
     3		31, 32, 33, 34, 35, 36, 37, 38, 39, 40,
     4		41, 42, 43, 44, 45, 46, 47, 48, 49, 50,
     5		51, 52, 53, 54, 55, 56, 57, 58, 59, 60,
     6		61, 62, 63, 64, 65, 66, 67, 68, 69, 70,
     7		71, 72, 73, 74, 75, 76, 77, 78, 79, 80,
     8		81, 82, 83, 84, 85, 86, 87, 88, 89, 90,
     9		91, 92, 93, 94, 95, 96, 97, 98, 99,100/
	data	(LTMPBUF(j),j=101,NTMP) /
     &	       101,102,103,104,105,106,107,108,109,110,
     1	       111,112,113,114,115,116,117,118,119,120,
     2	       121,122,123,124,125,126,127,128,129,130,
     3	       131,132,133,134,135,136,137,138,139,140,
     4	       141,142,143,144,145,146,147,148,149,150,
     5	       151,152,153,154,155,156,157,158,159,160,
     6	       161,162,163,164,165,166,167,168,169,170,
     7	       171,172,173,174,175,176,177,178,179,180/
	data	(ARXNAM(j),j=1,NARRX) /NARRX*-1/ JARX/0/
	data	Tflag/0/
	end
C======================================================================|
