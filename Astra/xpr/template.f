C======================================================================|
C ./template /home/grp/TAstra/.tsk/Intel/neut.exe 14638 772309392 1
C----------------------------------------------------------------------|
C Make as
C $> cc -c xpr/template.c -o xpr/temp.o
C $> cc -c xpr/ipc.c -o xpr/ipc.o
C $> f95f xpr/template.f xpr/temp.o xpr/ipc.o -o xpr/template
C----------------------------------------------------------------------|
	program   main
	implicit  none
	integer	  j,iargc,getpid,mampid,mamkey,eignr
	character*132 STRING,eigpath,mampath
	if (iargc() .ne. 4)	then
	   write(*,*)"Error"
	   call a_stop
	endif
	call	getarg(0,STRING)
	eigpath = STRING(1:len_trim(STRING))//char(0)
	call	getarg(1,STRING)
	mampath = STRING(1:len_trim(STRING))//char(0)
	call	getarg(2,STRING)
	read(STRING,*)mampid
	call	getarg(3,STRING)
	read(STRING,*)mamkey
	call	getarg(4,STRING)
	read(STRING,*)eignr
	write(*,*)
	call sbp2shm(eigpath, mampath, mampid, mamkey, eignr)
	end
C======================================================================|
