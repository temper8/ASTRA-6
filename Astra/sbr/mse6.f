C======================================================================|
	subroutine MSE6(YNOC)
C----------------------------------------------------------------------|
C Subroutine description:
C Note!	This subroutine is specialized to set exclusively 
C	the array "CAR6X"
C----------------------------------------------------------------------|
C Input:
C       YNOC Number of channels. ! Note! YNOC should be real for
C       			 ! compatibility with C-parameters)
C
C----------------------------------------------------------------------|
C Output:
C	CAR6X
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	integer		j,jf,length,JN,JMSE,IFIELD,ICHAN,ERCODE
        double precision YNOC,r,z,br,bz,bt,gamma
	logical		EXI
	save		JMSE
	data		JMSE/0/
C----------------------------------------------------------------------|
	jf = length(MSFILE)
	inquire(file=MSFILE(1:jf),exist=EXI)
C	write(*,*)'>>> MSE: Input file "',MSFILE(1:jf),'" will be used'
C----------------------------------------------------------------------|
        JN = abs(int(YNOC+0.5))
        if (YNOC .lt. .0)	JN = JN+1
        if (   JN .eq. 0)	goto	99
        if (JMSE .eq. 0)	JMSE = JN	! 1st call
	if (JMSE .ne. JN)	then		! No. of sources changed
           JMSE = JN
           goto	1
        endif
	if (EXI .and. YNOC.gt.0 )	goto	2
C----------------------------------------------------------------------|
 1      continue
C This piece is not ready for the time being.
C It is supposed to be used for interactive change of the MSE input file
C It should be activated if
C   (1) the subroutine is called with YNOC < 0, then abs(YNOC) gives 
C       the number of channels to be set,
C   (2) YNOC is changed during the Astra run.
C	write(*,*)'>>> MSE: Entering interactive parameter setting'
C Interactive viewer/editor of the beam source parameter list
C	open(2,file=FILNAM(1:lnbinp),status='UNKNOWN',err=95)
C	call	MQUERY(2,FILNAM(1:lnbinp),JMSE,ERCODE)
C       close(2)
	write(*,*)
     >	     '>>> MSE: Sorry, interactive data setting is not activated'
	return
C----------------------------------------------------------------------|
 2	continue
        JN = 1
	open(2,file=MSFILE(1:jf),status='OLD')
C	write(*,*)'>>> MSE: Reading parameter list'
C Read a record for one MSE channel
 3	continue
C Input to the subroutine STREAD:
C   2  is alogical unit connected to the input file file
C   8  is a number of fields in one group
C  WORK1 work array is be used for the data exchange
C  ERCODE error code
        call	STREAD(2,8,WORK1(1,JN),ERCODE)
C	write(*,*)"After STREAD: ercode =",ERCODE
        goto	(77,96,97,98,95),ERCODE+1
C----------------------------------------------------------------------|
 77	continue
        JN = JN+1
 	if (JN .le. JMSE)	goto	3
	close(2)
C----------------------------------------------------------------------|
	do	JN =1,JMSE
C	   write(*,'(2I4)')JN,JMSE
C	   write(*,'(1P,5G12.4)')(WORK1(j,JN),j=1,8)
	   r=work1(1,JN)
	   z=work1(2,JN)
C	   write(*,*)"r=",r
C	   write(*,*)"z=",z
	   call brz(r,z,br,bz)
C	   write(*,*)"br=",br
C	   write(*,*)"bz=",bz
C	   bz=-1*bz
C	   br=-1*br
	   bt=-1*BTOR*RTOR/r
C	   write(*,*)"bt=",bt
	   gamma=atan((work1(3,JN)*br+work1(4,JN)*bt+work1(5,JN)*bz)
     >     /(work1(6,JN)*br+work1(7,JN)*bt+work1(8,JN)*bz))*180/3.14159
C	   gamma=atan(work1(4,JN)/work1(7,JN))*180/3.14159
C	   write(*,*)gamma
	   CAR6X(JN)=gamma
	enddo
	call	MSET("CAR6X")
C The 2D array WORK1(1:NRD,2*NRD+7) is described in "for/status.inc" 
C     and can be used as a working array inside user subroutines.
C The subroutine A2MSE is an interface between Astra and MSE code
C	call	A2MSE(WORK1)
	return
 94     write(*,*)'>>> MSE >>> Input file name is too long'
        stop
 95     write(*,*)'>>> MSE calling STREAD: array out of limits'
        stop
 96     write(*,*)'>>> MSE >>> Error in file "',MSFILE(1:jf),'"'
     >		 , ': unrecognized variable name'
        stop
 97     write(*,*)'>>> MSE >>> File "',MSFILE(1:jf),'" read error'
        stop
 98     write(*,*)'>>> MSE >>> Wrong configuration file format. '
        write(*,*)'            More records expected than available.'
        stop
 99     write(*,*)'>>> MSE >>> Zero number of beams'
        write(*,*)"            Don't know what to do"
        stop
	end
C======================================================================|
	subroutine	MSET(CHAR6)
C----------------------------------------------------------------------|
C Note! DATARR is described as real*4
C----------------------------------------------------------------------|
	implicit 	none
	include 'for/parameter.inc'
	include 'for/const.inc'			! NAB is used
	include 'for/status.inc'
	include 'for/outcmn.inc'
	character	CHAR6*6,STRI6*6
        integer		jj,jarr,jtyp,jpnt,js,jn,length
	double precision RZ2A
C----------------------------------------------------------------------|
	STRI6 = CHAR6
	js = length(STRI6)
	if (js .lt. 6) write(STRI6(js+1:),'(5A1)')(' ',jj=js+1,6)
	jn = 0
	do	jj=1,NARRX
C	      write(*,*)jj,'"',STRI6,'"',EXARNM(jj)(1:js),'"',js
	   if (EXARNM(jj) .eq. STRI6)	jn = jj
	enddo
	if (jn .eq. 0)		return
	jarr = IFDFAX(jn)
	if (jarr .le. 0)	return
	jtyp = NTYPEX(jarr)
	if (jtyp .lt. 18)	return
	jpnt = NGRIDX(jarr)
	if (jpnt .le. 0)	return
	js = GDEX(jarr)

	if (jtyp .eq. 18 .or. jtyp .eq. 19)	then
C              write(*,*)EXARNM(jn),jn,jpnt,js,DATARR(js)
C              write(*,*)(DATARR(jj),jj=js+1,js+jpnt)
C              write(*,*)(DATARR(jj),jj=js+1,js+jpnt)
	elseif (jtyp .eq. 20)	then
	   do	jj=jpnt,1,-1
	      DATARR(js+2*jpnt+jj-1) = CAR6X(jj)
	   enddo
C              write(*,*)'Points:',jpnt,'  Address:',js
C              write(*,*)(DATARR(jj),jj=js,js+jpnt-1)
C              write(*,*)(DATARR(jj),jj=js+jpnt,js+2*jpnt-1)
C              write(*,*)(DATARR(jj),jj=js+2*jpnt,js+3*jpnt-1)
C              write(*,'(A,10F6.2)')"CAR6X  ",(CAR6X(jj),jj=10,1,-1)
	endif
	end
C======================================================================|
