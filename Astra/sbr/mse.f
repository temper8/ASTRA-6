C======================================================================|
	subroutine MSE(YNOC,YG)
C----------------------------------------------------------------------|
C Subroutine description:
C
C----------------------------------------------------------------------|
C Input:
C       YNOC Number of channels. ! Note! YNOC should be double precision
C       			 ! for compatibility with C-parameters)
C
C----------------------------------------------------------------------|
C Output:
C
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	include 'for/outcmn.inc'
	integer		j,jf,length,JN,JMSE,IFIELD,ICHAN,ERCODE
        double precision YNOC
	logical		EXI
	save		JMSE
	data		JMSE/0/
	double precision r,z,br,bz,bt,gamma,YG(NA1)
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
	   YG(JN)=gamma
	enddo
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
