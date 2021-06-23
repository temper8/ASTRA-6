C======================================================================|
	implicit none
	character*40	revnam,FILENA,RDNAME,EQNAME
	character*1	TZER,STR4*4,WDATE*132,DEVICE*132
	integer		IDAY,IMON,IYR,IHOUR,IMIN,LM,j,i,k,leng,lonlen
	data		j/0/k/0/
	TZER=char(j)
1	j = j+1
	read(*,'(1A40)',ERR=21,END=24)revnam
	LM	=leng(revnam(1:1))
C	write(*,*)'"',revnam(1:LM),'"'
	if ( revnam(1:10) .eq. "_Retrieve_" )	then
	   FILENA=revnam(11:LM)//TZER
	   i = 1
	else
	   FILENA=revnam(1:LM)//TZER
	   i = 0
	endif
	if (j.eq.1 .and. i.eq.0)	then
	   write(*,'(/A//2A)')
     +	   '  >>>>>>>>>>   Available "review" files   <<<<<<<<<<<',
     +	   '  Name                  ',
     +	   'Date      Time         Data_file           Model'
	endif
	if (LM .le. 0)	goto	1
C	if (LM .gt. 20)	goto	21
C	write(*,*)revnam
C	write(*,*)LM,i,"  ",FILENA
C Initial infofmation, profile number, names and scales
	open(3,FILE=FILENA,STATUS='OLD',ERR=23,FORM='UNFORMATTED')
	call	RESTORE(3,i,k)	! Skip model, model.log; define type
C	write(*,*)k
	if (k.ne.0)	then
	   read(3,ERR=22)RDNAME,EQNAME,DEVICE(1:32),DEVICE
     +		,IYR,IMON,IDAY,IHOUR,IMIN
	else
	   read(3,ERR=22)RDNAME(1:8),EQNAME(1:8),DEVICE(2:17)
     +		,IYR,IMON,IDAY,IHOUR,IMIN
	endif
C	write(*,*)IYR,IMON,IDAY,IHOUR,IMIN
	write(STR4,'(1I2.2)')IMON
	write(STR4(3:4),'(1I2.2)')IMIN
	if (k.ne.0)	then
	   write(WDATE,100)FILENA,IDAY,IMON,1900+IYR,IHOUR,IMIN,
     +		RDNAME(1:leng(RDNAME)),EQNAME(1:leng(EQNAME))
	else
	   write(WDATE,100)FILENA,IDAY,IMON,1900+IYR,IHOUR,IMIN,
     +		RDNAME(1:8),EQNAME(1:8)
	endif
 100	format(A20,I3,'-',I2.2,'-',I4,I5,':',I2.2,4X,A14,A16)
	if (i.eq.0)	write(*,*)WDATE(1:lonlen(WDATE))
	close (3)
	go to 1
 21	stop	'Review: >>>> File name error  <<<<'
 22	write(*,*)'Review:  Illegal review file "',FILENA,'"'
C	i=system('rm -i '//FILENA)
	stop
 23	if ( revnam(1:10) .eq. "_Retrieve_" )	LM = LM-10
	write(*,*)
	write(*,*)
     >	    'Review: >>>> File "',FILENA(1:LM),'" does not exist  <<<<'
	write(*,*)
	stop
 24	write(*,*)
	end
C======================================================================|
	subroutine RESTORE(NCHR,IFW,IRET)
C NCHR channel No.
C IFW =/= 0	Restore model & model.log file
C IFW  = 0	Skip model & model.log records
C   Output:
C IRET = 0	old type record
C IRET = 1	new type record
C----------------------------------------------------------------------|
	implicit none
	integer NCHR,IFW,IRET,j,n
	character STRI*132,STR*32,CH1*1
	read(NCHR) STR
C	write(*,*)"IFW = ",IFW
	if (STR.ne."^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^")	then
	   rewind(NCHR)
	   IRET = 0
	   return
	endif
	IRET = 1
	j = 0
	n = 1
	STR = "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"
	if (IFW .ne. 0) open (IFW,file='model.tmp',STATUS='UNKNOWN')
 10	read(NCHR) CH1,STRI(1:ichar(CH1))
	if (ichar(CH1).eq.32 .and. STR.eq.STRI(1:32)) then
	   j = j+1
	   n = n+1
	   if (IFW .ne. 0) close(IFW)
	   if (j.eq.2)	return
	   if (n.eq.2 .and. IFW.ne.0)
     >			open(IFW,file='model.log',STATUS='UNKNOWN')
	else
	   j = 0
	endif
	if (j.eq.0 .and. IFW.ne.0) write(IFW,'(A)')STRI(1:ichar(CH1))
C	if (j.eq.0) write(*,'(A)')STRI(1:ichar(CH1))
	goto	10
	end
C======================================================================|
	integer	function	leng(string)
C----------------------------------------------------------------------|
	character*1	string*(*),ch1,TZER,TAB,TCR
	integer		j
	j = 0
	TZER=char(j)
	j = 9
	TAB=char(j)
	j = 13
	TCR=char(j)
	leng	=0
	do	j=1,120
	ch1	=string(j:j)
	if(ch1.eq.TZER.or.ch1.eq.' '.or.ch1.eq.TAB.or.ch1.eq.TCR)
     &			then
		leng=j-1
		return
			endif
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
C	   if (STRI(j:j).eq.TZER .or. STRI(j:j).eq.TCR)	return
C	   if (STRI(j:j).eq.TZER)	return
C	   if (STRI(j:j).eq.TCR)	return
	   if (STRI(j:j).ne.' ' .and. STRI(j:j).ne.TAB)	lonlen = j
	enddo
	end
C======================================================================|
