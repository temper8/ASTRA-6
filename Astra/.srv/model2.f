C ====== SETVAR, DETVAR, INIVAR, ALLDEF, APPSBN, RADOUT, TIMOUT =======|
C =================== ININAM, ANLFML, BLOCKDATA =======================|
C======================================================================|
C======================================================================|
	subroutine  SETVAR(CNSBUF,NCH)
C----------------------------------------------------------------------|
C Write "setvar.tmp" (formerly a part of "detvar.tmp")
C    Initialize arrays
C	ZEF=ZEFX, ZMAIN=ZMJ, AMAIN=AMJ, NI=NE/ZMJ
C    Check boundary consistence
C
C----------------------------------------------------------------------|
	implicit	none
	include		'for/parameter.inc'
	include		'.srv/tmpbuf.inc'
	integer		NCH,IERR,j
	character*14	FNAME
	character*1	CNSBUF(*)
	data		FNAME/'tmp/setvar.tmp'/
	call	OPENWT(NCH,FNAME,0,IERR)
	if(IERR.gt.1)	then
	   write(*,*)'MODEL (MAIN): Cannot open file ',FNAME
	   call	exit(1)
	endif
	write(NCH,50)
 50	format(
     1	'      call	markloc("setvar.tmp"//char(0))'/
     2	'      do	j = 1,NA1'/
     3	'      ZEF(J)=max(1.d0,ZEFX(J))'/
     4	'      ZMAIN(J)=ZMJ'/
     5	'      AMAIN(J)=AMJ'/
     6	'      NI(J)=NE(J)/ZMJ'/
     7	'      enddo')
	write(NCH,'(10(A/),A)')
     >   '      if (ABC+abs(SHIFT).gt.AB) then'
     >	,'         write(*,*)char(7),'
     >	,'     >   ">>> Warning >>> Inconsistent boundary setting."'
     >	,'         if (NEQUIL.ne.42) then'
     >	,'            write(*,*)'
     >	,'     >"    Plasma beyond the vacuum vessel has been cut off"'
     >	,'         else'
     >	,'            write(*,*)'
     >	,'     >"    Plasma boundary intersects the vacuum vessel"'
     >	,'         endif'
     >	,'      endif'
	if (NBEG(LRON).eq.0)write(NCH,'(A)')'      RON=ROC'
	if (NBEG(LROE).eq.0)write(NCH,'(A)')'      ROE=ROC'
	if (NBEG(LROI).eq.0)write(NCH,'(A)')'      ROI=ROC'
	if (NBEG(LRO0).eq.0)write(NCH,'(A)')'      RO0=ROC'
	if (NBEG(LRO1).eq.0)write(NCH,'(A)')'      RO1=ROC'
	if (NBEG(LRO2).eq.0)write(NCH,'(A)')'      RO2=ROC'
	if (NBEG(LRO3).eq.0)write(NCH,'(A)')'      RO3=ROC'
	if (NBEG(LRO4).eq.0)write(NCH,'(A)')'      RO4=ROC'
	if (NBEG(LRO5).eq.0)write(NCH,'(A)')'      RO5=ROC'
	if (NBEG(LRO6).eq.0)write(NCH,'(A)')'      RO6=ROC'
	if (NBEG(LRO7).eq.0)write(NCH,'(A)')'      RO7=ROC'
	if (NBEG(LRO8).eq.0)write(NCH,'(A)')'      RO8=ROC'
	if (NBEG(LRO9).eq.0)write(NCH,'(A)')'      RO9=ROC'
	call	APPTMP(CNSBUF,NCH)
	close(NCH)
	end
C======================================================================|
C Write "detvar.tmp"
C	IVBUF
C	DTVBUF(1:IVBUF<5000) 	descibed in "tmpbuf.inc"
C		includes quantities to be written out of "jdetv" loop
C	NTMP			parameter descibed in "tmpbuf.inc"
C	NARRS
C	ARRSTR(1:NARRS)		astra arrays on rhs of "=" operator
C	LNAMAR(1:NARRS)
C	NAMARR(1:LNAMAR(1:NARRS),1:NARRS) array on lhs of "=" operator
C	TMPBUF(IBUF <= 10000)	- Contents of the do loop in detvar.tmp
C				  Includes all formula proccessing results
C				  according to the ASTRA rules
C======================================================================|
	subroutine  DETVAR(NCH,NARRS,LNAMAR,NAMARR,ARRSTR,EQBUF,LOCSBR)
	implicit	none
	include		'.srv/nambuf.inc'
	include		'.srv/tmpbuf.inc'
	character*1	EQBUF(*),NAMARR(6,100),ARRSTR(100)*78,FNAME*14,
     +			STR6*6,LINTXT(133),LINFOR(133)
	integer	LEQ(NEQNS),NCH,NARRS,LNAMAR(100),J,J1,J2,J3,J4
	integer	JJ,JL,LUR,NBEGB(NTMP),LOUT,NTMMIN,LOCSBR(*),JSBP,NSBP
	data	FNAME/'tmp/detvar.tmp'/NBEGB/NTMP*0/

	call	OPENWT(NCH,FNAME,0,JJ)
	if(JJ.gt.1)	then
	   write(*,*)'MODEL (MAIN): Cannot open file ',FNAME
	   call	exit(1)
	endif
	write(NCH,'(2A)')'      call	markloc("detvar.tmp',
     +	' (time signals)"//char(0))'

C	write(NCH,'(A)')'      JEX = max(1,int(ITEREX+0.5))-JIT'
C Print the buffer DTVBUF to stdout 
C	write(*,*)'>>> DETVAR >>> DTVBUF:'
C	write(*,*)'DTVBUF start'
C	call	APPTMP(DTVBUF,6)
C	write(*,*)'DTVBUF end'
C	write(*,*)'TMPBUF start'
C	call	APPTMP(TMPBUF,6)
C	write(*,*)'TMPBUF end'

C Add analysis for the case (IVBUF .eq. 1)	i.e. no arrays specified

C Constants: write variable definitions before the main loop
C DTVBUF is created in main and ANLFML and is used here only
	call	APPTMP(DTVBUF,NCH)

	do	j=1,NEQNS
	   LEQ(j) = -1
	enddo
C Find non-defined main variables
	JJ=1
 48	if (ichar(EQBUF(JJ)).gt.129) then
	   J1 = (ichar(EQBUF(JJ))-129)/4
	   LEQ(J1) = 0
	   JJ=JJ+1
	   goto	48
	endif
	J=ICHAR(EQBUF(JJ))+JJ
	if(J.eq.JJ)	goto	49
	JJ=J+1
	go to 48
 49	continue

C	do	J1=1,NARRS
C	   write(*,'(I3,60A)')J1,'	',(NAMARR(J,J1),J=1,LNAMAR(J1))
C	   call	COPY(6,NAMARR(1,J1),STR6)
C	   if (STR6(1:LNAMAR(J1)) .eq. "NI")write(*,*)"Yes"
C	enddo

	write(NCH,'(2A)')'C **** Radial profile computation'
	write(NCH,'(2A)')'      call	markloc("detvar.tmp',
     +	' (profiles)"//char(0))'
 100	format(132A)
	write(NCH,100)'      do	jdetv = 1,NA1'
	write(NCH,100)'      J = jdetv'
	J1 = 0			! Write default definitions 
	J2 = 0			! for these four arrays 
	J3 = 0			! only if they are not defined
	J4 = 0			! in a model
	if (NBEG(LNI) .gt. 0)	J4 = 1
	if (NARRS .lt. 1)	goto	50	! No arrays defined
	do	J=1,NARRS
	   write(STR6,100)(NAMARR(JJ,J),JJ=1,LNAMAR(J))
	   if (STR6 .eq. "ZEF   ")	J1 = 1
	   if (STR6 .eq. "ZMAIN ")	J2 = 1
	   if (STR6 .eq. "AMAIN ")	J3 = 1
	   if (STR6 .eq. "NI    ")	J4 = 1
	enddo
 50	if (J1 .eq. 0)	write(NCH,100)'      ZEF(J)=MAX(1.d0,ZEFX(J))'
	if (J2 .eq. 0)	write(NCH,100)'      ZMAIN(J)=ZMJ'
	if (J3 .eq. 0)	write(NCH,100)'      AMAIN(J)=AMJ'
	if (J4 .eq. 0)	write(NCH,100)'      NI(J)=NE(J)/ZMJ'
	if (NARRS .lt. 1)	goto	55	! No arrays defined

C Exclude repeated formulae
	KEYOUT=1
	NOUT=0
C----------------------------------------------------------------------|
C NARRS, NAMARR, LNAMAR; RHS - ARRSTR, NBEG
	do	54	J1=1,NARRS
	   if (Tflag .ne. 0)	write(NCH,'(60A)')
     +		'      call add2loc("Radial array ',
     +     	(NAMARR(J,J1),J=1,LNAMAR(J1)),'(J)"//char(0))'
C	   write(*,100)'"',(NAMARR(JJ,J1),JJ=1,LNAMAR(J1)),'"'
	   write(STR6,100)(NAMARR(JJ,J1),JJ=1,LNAMAR(J1))
	   if(STR6 .ne. "CUBS  ")	goto	51
	      if(NBEG(LHC)+NBEG(LXC)+NBEG(LDC).ne.0)	then
		 NBEG(LCUBS) = 0
		 write(*,'(A,1X/22X,132A)')
     +			" >>> Warning >>> CUBS definition in the model "
     +			,(NAMARR(J,J1),J=1,LNAMAR(J1)),'(J)=',
     +			(ARRSTR(J1)(j+1:j+1),j=1,ichar(ARRSTR(J1)(1:1)))
		 write(*,'(A)')"                 will be overwritten"
	      endif
 51	   continue
C	LUR = ICHAR(ARRSTR(J1)(1:1))
C	write(*,'(A,I3,A,I3,132A1)')
C     +		"No.",J1,"   RHS_length",LUR,'		',
C     +		(ARRSTR(J1)(j+1:j+1),j=1,LUR)
C Form rhs
C	   write(*,*)"Call from DETVAR:",KEYOUT,NCH
	   call ANLFML(NCH,ARRSTR(J1),LINTXT,LINFOR)

	   LUR = ICHAR(LINFOR(1))
	   do	j=1,LUR
	      LINFOR(j) = LINFOR(j+1)
	   enddo
C Write output line to stdout
C	write(*,'(A$)')'Output LINE: "'
C	write(*,'(60(A1$))')(LINFOR(J),J=1,LUR)
C	write(*,'(A,1I4)')'",  length = ',LUR
C Write resulting line to stdout
C	write(*,'(60A1)')'	',(NAMARR(J,J1),J=1,LNAMAR(J1)),
C     +		'(','J',')','=',(LINFOR(J),J=1,LUR)
	   JL = LUR
	   if (LUR .gt. 56-LNAMAR(J1))	   LUR = 56-LNAMAR(J1)
	   write(NCH,'(6X,61A1)')(NAMARR(J,J1),J=1,LNAMAR(J1)),
     +				'(','J',')','=',(LINFOR(J),J=1,LUR)
	   if (LUR .eq. JL)	goto	53
	   LUR = JL-LUR
	   JL = 56-LNAMAR(J1)
 52	   continue
	   if (LUR.gt.60)	then
	      write(NCH,'(5X,1A2,62A1)')'> ',(LINFOR(J),J=JL+1,JL+60)
	      JL = JL+60
	      LUR = LUR-60
	      goto	52
	   else
	      write(NCH,'(5X,1A2,62A1)')'> ',(LINFOR(J),J=JL+1,JL+LUR)
	   endif
 53	   LOUT=ICHAR(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
C	write(*,131)(NAMARR(J,J1),J=1,LNAMAR(J1)),'(','r',')','=',
C     1			(LINTXT(J),J=1,LOUT)
 131	   format(1X,130A1)
	   write(7,131)(NAMARR(J,J1),J=1,LNAMAR(J1)),'(','r',')','=',
     1			(LINTXT(J),J=1,LOUT)
 54	continue
 55	continue

C------------------------------
C Definition of all the quantitiies from LTMPBUF(1:NTMP) is moved to eqns.inc
C except for the case when the correspondent main variable is not mentioned.
C CC must be always defined
	if (NBEG(LCC) .eq. 0)	then
	   call APPFML(NCH,4,"ccsp")
	   write(NCH,'(A)')'      CC(J)=CCSP'
	else
C	   call	APPTMP(TMPBUF(NBEG(LCC)),NCH)
	endif
	do	58	J1=1,NTMP
	   NTMMIN = 10000
	   JJ = 0
	   do	57	J=1,NTMP
	      if(NBEG(J).le.0)	goto	57
	      if(NBEGB(J).gt.0)	goto	57
C Suppress writing scalars to "detvar.tmp"
C left: IPL,LEXT,UEXT,NEB,TEB,TIB,CU,MU,NE,TE,TI,QNB,QEB,QIB,
C	QNNB,QETB,QITB,QFjB,QFFjB
	      if ( J.eq.LIPL .or. J.eq.LLEXT.or. J.eq.LUEXT
     .	      .or. J.eq.LNEB .or. J.eq.LTEB .or. J.eq.LTIB
     .	      .or. J.eq.LRON .or. J.eq.LROE .or. J.eq.LROI
     .	      .or. J.eq.LCU  .or. J.eq.LMU  .or. J.eq.LMV  .or. J.eq.LCV
     .	      .or. J.eq.LNE  .or. J.eq.LQNB .or. J.eq.LQNNB
C     .	      .or. J.eq.LNI  ! add this line to move NI=... to "eqns.inc"
     .	      .or. J.eq.LTE  .or. J.eq.LQEB .or. J.eq.LQETB
     .	      .or. J.eq.LTI  .or. J.eq.LQIB .or. J.eq.LQITB)	goto  57
C Vector definitions will be placed in eqns.inc
	      if ((J.eq.LDN .or. J.eq.LDE .or. J.eq.LDI .or.
     .		   J.eq.LDC .or. J.eq.LSN .or. j.eq.LSNN.or.
     .		   J.eq.LDVN.or. J.eq.LDSN )
     .					.and. LEQ(1).ge.0)	goto  57
	      if ((J.eq.LHN .or. J.eq.LHE .or. J.eq.LHI .or.
     .		   J.eq.LHC .or. J.eq.LPE .or. j.eq.LPET.or.
     .		   j.eq.LDVE.or. j.eq.LDSE )
     .					.and. LEQ(2).ge.0)	goto  57
	      if ((J.eq.LXN .or. J.eq.LXE .or. J.eq.LXI .or. 
     .		   J.eq.LXC .or. J.eq.LPI .or. j.eq.LPIT.or.
     .		   J.eq.LDVI.or. J.eq.LDSI )
     .					.and. LEQ(3).ge.0)	goto  57
	      if ((J.eq.LCN .or. J.eq.LCE .or. J.eq.LCI .or. 
     .		   J.eq.LCD)		.and. LEQ(4).ge.0 )	goto  57
	      do	j2=0,9
		 if ( J.eq.LRO(j2) .or. J.eq.LFB(j2) .or. J.eq.LF(j2)
     .		 .or. J.eq.LQFB(j2).or. J.eq.LQFFB(j2) )	goto  57
		 if ((J.eq.LDF(j2) .or. J.eq.LVF(j2) .or.
     .		      J.eq.LSF(j2) .or. J.eq.LSFF(j2).or.
     .		      j.eq.LDVF(j2).or. J.eq.LDSF(j2) )
     .					.and. LEQ(10+j2).ge.0)	goto  57
	      enddo
	      if (NTMMIN .le. NBEG(J))				goto  57
	      NTMMIN = NBEG(J)
	      JJ = J
 57	   continue
	   if (JJ .le. 0)	goto	58
	   J2 = ICHAR(TMPBUF(NTMMIN))
	   NBEGB(JJ) = 1
	   call	APPTMP(TMPBUF(NTMMIN),NCH)
C print TMPBUF:	call APPTMP(TMPBUF(NTMMIN),6)
C	   call	APPTMP(TMPBUF(NTMMIN),6)
 58	continue
C------------------------------
 60	continue
	write(NCH,'(A)')'      VRO(J)=VR(J)'
	write(NCH,'(A)')'      enddo'
	write(NCH,'(A)')'      if (NEQUIL.lt.42) then'
	write(NCH,'(A)')'      SHIFT=min(SHIFT,0.9*AB)'
	write(NCH,'(A)')'      ABC = min(ABC,AB-abs(SHIFT))'
	write(NCH,'(A)')'      ABC = max(ABC,0.1*AB)'
	write(NCH,'(A)')'      endif'

	NSBP = 0
	do	j=1,NSBMX
	   if (LOCSBR(j).eq.-2 .or. LOCSBR(j).eq.-3)	then
	      NSBP = NSBP+1			! Count number of SBPs
	   endif
	   if (LOCSBR(j) .eq. -1)	then	! put SBR tagged with "<"
	      J2 = 0
	      call	APPSBN(NCH,EQBUF,j,j2)
C	      write(*,*)"SBR No.",j," will be called from detvar.tmp"
	      write(7,*)"SBR No.",j," will be called from detvar.tmp"
	   endif
	enddo

	if (NSBP .gt. 0)	then
	   write(NCH,'(A)')'C **** Fill shared memory segments'
	   write(NCH,'(A)')'      call	markloc("setvars"//char(0))'
	   write(NCH,'(A)')'      call	setvars(AB,TIME,HRO,NRD)'
	   write(NCH,'(A)')'      call	markloc("setarrs"//char(0))'
	   write(NCH,'(A)')'      call	setarrs(TE,AMETR,NN,NRD)'
	endif
	J1 = 0
	JSBP = 0
	JJ = 1
	do	j=1,NSBMX
	   J1 = J1+1
	   J3 = ichar(EQBUF(JJ))
	   if (LOCSBR(j) .eq. -2)	then
	      JSBP = JSBP+1
	   endif
	   if (LOCSBR(j) .eq. -3)	then	! put SBP tagged with "<"
	      JSBP = JSBP+1
	      J2 = 0
	      call SUBPROC(J1,EQBUF(JJ),NCH,JSBP,j2)
	   endif
	   JJ=JJ+J3+1
	enddo
	close(NCH)
	end		! subroutine DETVAR
C======================================================================|
	subroutine	INIVAR(NCH,LEQ)
	implicit	none
	character*1	FNAME*14,CJ
	integer		NCH,IERR,LEQ(*),J,j2,LENC,NB
	include 	'.srv/nambuf.inc'
	include 	'.srv/tmpbuf.inc'
	data		FNAME/'tmp/inivar.tmp'/
	KEYOUT=0			! used in APPTMP -> APPFML
	call	OPENWT(NCH,FNAME,0,IERR)
	if (IERR.gt.1)	then
	   write(*,*)'MODEL (MAIN): Cannot open file ',FNAME
	   call	exit(1)
	endif

	write(NCH,'(A)')
     +	'      call	markloc("inivar.tmp"//char(0))'
	write(NCH,'(A)')'      DTEQL=0.d0'
	write(NCH,'(A)')'      TAU=TAUMIN'
	if(NBEG(LTE).eq.0)	then
	   write(NCH,'(A)')'      j1 = 0'
	   write(NCH,'(A)')'      do	j=1,NARRX'
	   write(NCH,'(A)')
     >	'         if(EXARNM(j).eq."TEX   ".and.IFDFAX(j).lt.0) j1=j'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')
     >	'      if (NA1E.eq.NA1 .and. j1.ne.0 .and. JIT.eq.0) then'
	   write(NCH,'(A/2A)')
     >	'        write(*,*)" >>> Warning: TE and TEX are not defined"',
     >	'        write(*,*)"              Default value TE=10eV ',
     >	'will be used"'
	   write(NCH,'(A)')'      endif'
	   write(NCH,'(A)')'      if (IFDFAX(j1) .gt. 0)	then'
	   write(NCH,'(A)')'         do	J=1,NA1'
	   write(NCH,'(A)')'            TE(J)=TEX(J)'
	   write(NCH,'(A)')'         enddo'
	   write(NCH,'(A)')'      endif'
	else
	   write(NCH,'(A)')'      do	J=1,NA1'
	   call APPTMP(TMPBUF(NBEG(LTE)),NCH)
	   write(NCH,'(A)')'      enddo'
	endif
	if (LEQ(2).ne.0 .and. LINAPP(2).ne.0)	then
	   if (NBEG(LTEB).gt.0)	then
	      write(NCH,'(A)')'      ND1 = NA1'
	      call APPTMP(TMPBUF(NBEG(LTEB)),NCH)
	   endif
	   if (NBEG(LROE) .eq. 0)	then
	      write(NCH,'(A)')'      ND1 = NA1'
	   else
	      call APPTMP(TMPBUF(NBEG(LROE)),NCH)
	   endif
	   write(NCH,'(A)')'      NA1E = ND1'
	   write(NCH,'(A)')'      if (ND1 .lt. NA)	then'
	   write(NCH,'(A)')'      YA = (TE(NA1)-TE(ND1))/(ROC-ROE)'
	   write(NCH,'(A)')'      YB = TE(ND1)*ROC-TE(NA1)*ROE'
	   write(NCH,'(A)')'      YB = YB/(ROC-ROE)'
	   write(NCH,'(A)')'      do j=ND1+1,NA'
	   write(NCH,'(A)')'      TE(j) = YB+RHO(j)*YA'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	endif

	if(NBEG(LTI).eq.0)	then
	   write(NCH,'(A)')'      j1 = 0'
	   write(NCH,'(A)')'      do	j=1,NARRX'
	   write(NCH,'(A)')
     >	'         if(EXARNM(j).eq."TIX   ".and.IFDFAX(j).lt.0) j1=j'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')
     >	'      if (NA1E.eq.NA1 .and. j1.ne.0 .and. JIT.eq.0) then'
	   write(NCH,'(A/2A)')
     >	'        write(*,*)" >>> Warning: TI and TIX are not defined"',
     >	'        write(*,*)"              Default value TI=10eV ',
     >	'will be used"'
	   write(NCH,'(A)')'      endif'
	   write(NCH,'(A)')'      if (IFDFAX(j1) .gt. 0)	then'
	   write(NCH,'(A)')'         do	J=1,NA1'
	   write(NCH,'(A)')'            TI(J)=TIX(J)'
	   write(NCH,'(A)')'         enddo'
	   write(NCH,'(A)')'      endif'
	else
	   write(NCH,'(A)')'      do	J=1,NA1'
	   call APPTMP(TMPBUF(NBEG(LTI)),NCH)
	   write(NCH,'(A)')'      enddo'
	endif
	if (LEQ(3).ne.0 .and. LINAPP(3).ne.0)	then
	   if (NBEG(LTIB).gt.0)	then
	      write(NCH,'(A)')'      ND1 = NA1'
	      call APPTMP(TMPBUF(NBEG(LTIB)),NCH)
	   endif
	   if (NBEG(LROI) .eq. 0)	then
	      write(NCH,'(A)')'      ND1 = NA1'
	   else
	      call APPTMP(TMPBUF(NBEG(LROI)),NCH)
	   endif
	   write(NCH,'(A)')'      NA1I = ND1'
	   write(NCH,'(A)')'      if (ND1 .lt. NA)	then'
	   write(NCH,'(A)')'      YA = (TI(NA1)-TI(ND1))/(ROC-ROI)'
	   write(NCH,'(A)')'      YB = TI(ND1)*ROC-TI(NA1)*ROI'
	   write(NCH,'(A)')'      YB = YB/(ROC-ROI)'
	   write(NCH,'(A)')'      do j=ND1+1,NA'
	   write(NCH,'(A)')'      TI(j) = YB+RHO(j)*YA'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	endif

	if(NBEG(LNE).eq.0)	then
	   write(NCH,'(A)')'      if (NA1N .eq. NB1)	then'
	   write(NCH,'(A)')'         j1 = 0'
	   write(NCH,'(A)')'         do	j=1,NARRX'
	   write(NCH,'(A)')
     >	'            if(EXARNM(j).eq."NEX   ".and.IFDFAX(j).lt.0) j1=j'
	   write(NCH,'(A)')'         enddo'
	   write(NCH,'(A)')'         if (j1 .ne. 0) write(*,*)'
	   write(NCH,'(A)')'     >" >>> Warning: NEX is not defined"'
	   write(NCH,'(A)')'      endif'
	   write(NCH,'(A)')'      do	J=1,NA1'
	   write(NCH,'(A)')'      NE(J)=NEX(J)'
	   write(NCH,'(A)')'      enddo'
	else
	   write(NCH,'(A)')'      do	J=1,NA1'
	   call APPTMP(TMPBUF(NBEG(LNE)),NCH)
	   write(NCH,'(A)')'      enddo'
	endif
	if (LEQ(1).ne.0 .and. LINAPP(1).ne.0)	then
	   if (NBEG(LNEB).gt.0)	then
	      write(NCH,'(A)')'      ND1 = NA1'
	      call APPTMP(TMPBUF(NBEG(LNEB)),NCH)
	   endif
	   if (NBEG(LRON) .eq. 0)	then
	      write(NCH,'(A)')'      ND1 = NA1'
	   else
	      call APPTMP(TMPBUF(NBEG(LRON)),NCH)
	   endif
	   write(NCH,'(A)')'      NA1N = ND1'
	   write(NCH,'(A)')'      if (ND1 .lt. NA)	then'
	   write(NCH,'(A)')'      YA = (NE(NA1)-NE(ND1))/(ROC-RON)'
	   write(NCH,'(A)')'      YB = NE(ND1)*ROC-NE(NA1)*RON'
	   write(NCH,'(A)')'      YB = YB/(ROC-RON)'
	   write(NCH,'(A)')'      do j=ND1+1,NA'
	   write(NCH,'(A)')'      NE(j) = YB+RHO(j)*YA'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	endif

	write(NCH,'(A)')'      do	J=1,NA1'
	write(NCH,'(A)')
     >	'      NI(j)=min(NI(j),NE(j)*ZEF(j)/ZMAIN(j)**2)'
	write(NCH,'(A)')'      enddo'

 	write(NCH,'(A)')'      do	J=1,NA1'
	if(NBEG(LCU).eq.0 .and. NBEG(LMU).ne.0)	
     >		call APPTMP(TMPBUF(NBEG(LMU)),NCH)
	if(NBEG(LCU).eq.0)	then
	   if ( LEQ(4) .eq. 0)	then
	      write(NCH,'(A)')'      CU(J)=CC(J)'
	   endif
	else					! NBEG(LCU) .ne. 0
C          call APPTMP(TMPBUF(NBEG(LCU)),NCH)
	endif
	write(NCH,'(A)')'      enddo'

	do	j2=1,9
	   write(CJ,'(I1)')j2
	   write(NCH,'(A)')'      do	J=1,NA1'
	   if (NBEG(LF(j2)) .eq. 0)	then
	      write(NCH,'(3A)')'      F',CJ,'(J)=1.'
	   else
	      call APPTMP(TMPBUF(NBEG(LF(j2))),NCH)
	   endif
	   write(NCH,'(A)')'      enddo'

	   if (LEQ(j2).ne.0 .and. LINAPP(j2+10).ne.0)	then
	      if (NBEG(LFB(j2)).gt.0)	then
		 write(NCH,'(A)')'      ND1 = NA1'
		 call APPTMP(TMPBUF(NBEG(LFB(j2))),NCH)
	      endif
	      if (NBEG(LRO(j2)) .eq. 0)	then
		 write(NCH,'(A)')'      ND1 = NA1'
	      else
		 call APPTMP(TMPBUF(NBEG(LRO(j2))),NCH)
	      endif
	      write(NCH,'(3A)') '      NA1',CJ,' = ND1'
	      write(NCH,'(A)')  '      if (ND1 .lt. NA)	then'
	      write(NCH,'(7A)')
     >	      '      YA = (F',CJ,'(NA1)-F',CJ,'(ND1))/(ROC-RO',CJ,')'
	      write(NCH,'(6A)')
     >	      '      YB = F',CJ,'(ND1)*ROC-F',CJ,'(NA1)*RO',CJ
	      write(NCH,'(3A)') '      YB = YB/(ROC-RO',CJ,')'
	      write(NCH,'(A)')  '      do j=ND1+1,NA'
	      write(NCH,'(3A)') '      F',CJ,'(j) = YB+RHO(j)*YA'
	      write(NCH,'(A)')  '      enddo'
	      write(NCH,'(A)')  '      endif'
	   endif
	enddo

	close(NCH)
 165	format(1X,78A1)
	end
C======================================================================|
	subroutine	APPTXT(DISEQ,LEQ)
C DISEQ	(now obsolete) was used for writing model.txt exclusively.
C	It contains information about the initial conditions
C	or expressions assigning main variables 
	implicit	none
	character*1	DISEQ(78,7),KK
	integer		LEQ(*),LENGTF,j,j1,j2
	include 	'.srv/nambuf.inc'
	include 	'.srv/tmpbuf.inc'
	write(7,*)'=====   Initial distributions   ====='
	do	j1=1,4
	   j2 = 1+10*(j1-1)
C	   write(*,*)j2,NBEG(j2)
	   if(NBEG(j2).eq.0)	then
	      if (j1.eq.4)	then
		 if (LEQ(4) .eq. 0) write(7,*)'CU(r)=CC(r)*const'
	      else
		 write(7,'(4A)')TMPNAM(j2),'(r)=',TMPNAM(j2),'X(r)'
	      endif
	   else
C	      write(*,165) (DISEQ(j,j1),j=1,LENGTF(DISEQ(1,j1)))
	      write(7,165) (DISEQ(j,j1),j=1,LENGTF(DISEQ(1,j1)))
	   endif
	enddo
	do	j2=0,9
	   if (NBEG(LF(j2)) .eq. 0)	then
	      write(KK,'(I1)')j2
	      write(7,'(3A)')' F',KK,'(r)=1'
	   else
C	      write(*,165) (DISEQ(J,j2+10),J=1,LENGTF(DISEQ(1,j2+10)))
	      write(7,165) (DISEQ(J,j2+10),J=1,LENGTF(DISEQ(1,j2+10)))
	   endif
	enddo
 165	format(1X,78A1)
	end
C======================================================================|
C Files: "init.inc", "eqns.inc"
C Input:
C	NCH	- channel No.
C	NBEG(LCU),NBEG(LMU),NBEG(LCC),NBEG(LMU),NBEG(LPET)
C	EQBUF(IEBUF <= 1000) Subroutine_call_string
C			     One position is used to transfer Eq. ID
C	LEQ(JEQ)=ITYPE(-1<=ITYPE<=3)
C	LEQ(1)  LEQ(2)  LEQ(3)  LEQ(4)  LEQ(5)  LEQ(6) LEQ(7-9)
C	  NE	  TE	  TI	  CU	EqSolv	  NI   Not used
C       LEQ(10)  LEQ(11)  LEQ(12) etc. correspond to
C	  F0	  F1	  F2,  respectively
C       NSBMX (include for/parameter.inc)
C Below:    j=1 for NE; =2 for TE; =3 for TI; =4 for CU;
C	    j=5 for F0; =6 for F1; =7 for F2;
C	LEQ(j)    =-1 no equation (default)
C		  = 0 type: AS
C		  = 1 type: EQ
C		  = 2 type: FU
C		  = 3 type: heat conductivity flux + heat convection
C Note: option Tj:AS:5.2; is not supported !!! Should be added
C----------------------------------------------------------------------|
	subroutine	ALLDEF(NCH,LEQ,EQBUF,LOCSBR)
	implicit	none
	include		'for/parameter.inc'
	include		'.srv/tmpbuf.inc'
	character	EQBUF(*)*1,EQNAM(NEQNS)*3
	character	FNAME(2)*12,SBRNAM*132
	integer		LEQ(*),NPOS,LOCSBR(*),LOCSBP(NSBMX),EQNPOS
	integer		IMPEI,JEQ,JSB,J,J1,J2,J3,jj,M1,M2,M3,RNPOS
	integer		NCH,ICH,IERR,IEQ,ISB,ITYPE,NB,TCTR,MIX,MEX
	integer		NSBP,NSBR,MSBP,MSBR,JSBP
	data	IMPEI/0/ MIX/0/ MEX/0/ TCTR/0/ EQNPOS/0/ LOCSBP/NSBMX*0/
	data	( FNAME(j), j=1,2 )  /'tmp/eqns.inc','tmp/init.inc'/
	data	( EQNAM(j), j=1,9 )  /'NE ','TE ','TI ','CU ',
     >				      '   ','   ','   ','   ','   '/
	data	( EQNAM(j), j=10,14) /'F0 ','F1 ','F2 ','F3 ','F4 '/
	data	( EQNAM(j), j=15,19) /'F5 ','F6 ','F7 ','F8 ','F9 '/
C	call	SHOWBU(EQBUF)		! eqn ID & sbr line

	NB = 1
	j1 = 0
	ISB = 0
	NSBP = 0	! Total SBP
	NSBR = 0	! Total SBR
	MSBP = 0	! Total SBP number in eqns[init].inc
	MSBR = 0	! Total SBR number in eqns[init].inc
	do  55	j=1,NSBMX+14
	   JEQ = ichar(EQBUF(NB))
C	   write(*,*)NB,JEQ
	   if (JEQ .lt. 0)	JEQ = JEQ+256
	   if (JEQ .eq. 0)	goto	56	! end
	   if (JEQ .ge. 132)	then
	      if (EQNPOS .eq. 0)	EQNPOS = NB	! 1st TrEq position
	      NB = NB+1
	      goto	55
	   endif
C  sbr/sbp name analyser
	   M1 = NPOS(JEQ,EQBUF(NB+1),'(')-1
	   M1 = min(6,M1,JEQ)
	   write(SBRNAM,'(6A1)')(EQBUF(NB+j1),j1=1,M1)
	   call	UPCASE(M1,SBRNAM)
	   ISB = ISB+1
	   if (ISB .gt. NSBMX)	goto	56
C	   write(*,*)SBRNAM,ISB,NB,NB+JEQ
	   if	  (SBRNAM.eq."TSCTRL")  then
	      TCTR = ISB
	   elseif (SBRNAM.eq."MIXINT")  then
	      MIX = ISB
	   elseif (SBRNAM.eq."MIXEXT")  then
	      MEX = ISB
	   endif
	   if (LOCSBR(ISB) .eq. -3 .or. LOCSBR(ISB) .eq. -2)
     >					LOCSBP(ISB) = NB
	   if (LOCSBR(ISB) .eq. -4)	write(*,*)"???",ISB,NB
	   if (LOCSBR(ISB) .eq. -3)	NSBP = NSBP+1
	   if (LOCSBR(ISB) .eq. -2)	NSBP = NSBP+1
	   if (LOCSBR(ISB) .eq. -2)	MSBP = MSBP+1
	   if (LOCSBR(ISB) .ge.  0)	NSBR = NSBR+1
	   if (LOCSBR(ISB) .eq.  0)	MSBR = MSBR+1

	   NB = NB+JEQ+1
 55	continue
 56	continue

	do  57	j=1,NSBMX
	   if (LOCSBP(j).gt.EQNPOS .and. EQNPOS.ne.0)	then
	      M1 = ichar(EQBUF(locsbp(j)))
	      write(*,'(66A)')' >>> The subprocess "',
     >		(EQBUF(LOCSBP(J)+j1),j1=1,M1),'" is ordered'
	      write(*,'(32A)')
     >		'     after the main transport equation'
	      write(*,'(32A)')
     >		'     This calling sequence is not supported'
	      call	exit(1)
	   endif
 57	continue
 59	continue
C	call	exit(1)

	IEQ = 10
	do  70	ICH=NCH,NCH+1
	   
	call	OPENWT(ICH,FNAME(ICH-NCH+1),0,IERR)
	if (IERR .gt. 1)	then
	   write(*,*)'MODEL (MAIN): Cannot open file ',FNAME(ICH-NCH+1)
	   call	exit(1)
	endif

	if (NSBP+NSBR .eq. 0)	then
	   write(ICH,'(A)')'C **** No external subroutines'
	endif

	if (NSBP .ne. 0)	then
	   if (ICH .eq. NCH+1)	then		! writing init.inc only
	    write(ICH,'(A)')'      call	markloc("initipc"//char(0))'
	    write(ICH,'(A)')'      call	initipc(NB1)'
	   endif
	   write(ICH,'(A)')'C **** Fill shared memory segments'
	   write(ICH,'(A)')'      call	markloc("setvars"//char(0))'
	   write(ICH,'(A)')'      call	setvars(AB,TIME,HRO,NRD)'
	   write(ICH,'(A)')'      call	markloc("setarrs"//char(0))'
	   write(ICH,'(A)')'      call	setarrs(TE,AMETR,NN,NRD)'
	   if (ICH .eq. NCH+1)	then		! writing init.inc only
	    write(ICH,'(A)')'      call	markloc("inikids"//char(0))'
	    write(ICH,'(A)')'      call	inikids(NSBP,64,LISTSB)'
	   endif
	endif
	write(ICH,'(3A)')'      call	markloc("',FNAME(ICH-NCH+1),
     +				'"//char(0))'
	if (ICH.eq.NCH) write(NCH,'(A)')' 2100 continue'

	do	J=1,3
	   if (LEQ(J).lt.0 .or. ICH.ne.NCH)	then	! writing init.inc
	      IEQ = IEQ+1
	      if (J .eq. 1)	call	NEEQN(IEQ,-1,ICH)
	      if (J .eq. 2)	call	TEEQN(IEQ,-1,ICH)
	      if (J .eq. 3)	call	TIEQN(IEQ,-1,ICH)
	   endif
	enddo
	if (LEQ(4).lt.0 .or. ICH.ne.NCH)	then	! writing init.inc
	   IEQ = IEQ+1
	   if (NBEG(LCU).gt.0)	then
				call	CUASN(IEQ,-1,ICH)
	   elseif (NBEG(LMU).gt.0)	then
				call	CUASN(IEQ,0,ICH)
	   else
				call	CUASN(IEQ,1,ICH)
	   endif
	endif
	do	J=10,19
	   if (LEQ(J).lt.0 .or. ICH.ne.NCH)	then	! writing init.inc
	      IEQ = IEQ+1
	      call	FJEQN(IEQ,-1,J-10,ICH)
	   endif
	enddo

	NB = 1		! Start position in EQBUF
	JSBP = 0
	JSB = 0
	ISB = 0		! Ordinal SBP & SBR number

	if (LEQ(2).ge.9 .and. LEQ(3).ge.9)	IMPEI = 1
 60	continue
	JEQ = ichar(EQBUF(NB))
	if (JEQ.lt.0)	JEQ = JEQ+256
	if (JEQ.eq.0 .and. ICH.eq.NCH+1)	goto	71	! -> Exit
	if (JEQ.eq.0)	goto	70		! -> NCH+1
	if (JEQ.ge.132)	goto	68		! transport equation

C Subroutine and Xprocess treatment
	if (ISB .ge. NSBMX)	goto	60
	ISB = ISB+1		! Global ordinal number of : or &:
	JEQ = ichar(EQBUF(NB))
	M1 = NPOS(JEQ,EQBUF(NB+1),'(')-1
	M1 = min(6,M1,JEQ)

	write(SBRNAM,'(6A1)')(EQBUF(NB+j),j=1,M1)
	call	UPCASE(M1,SBRNAM)
C	write(*,*)ISB,LOCSBR(ISB),NB,ichar(EQBUF(NB)),'  "',SBRNAM,'"'
C Enable a special treatment for some sbrs
C	write(*,*)"  JEQ =",JEQ,'   "',(EQBUF(NB+j),j=1,JEQ)
C	write(*,*)"   M1 =",M1, '   "',(EQBUF(NB+j),j=1,M1)

C	if (ICH.eq.NCH)write(*,'(3I5,80(A))')ISB,LOCSBR(ISB),LOCSBP(ISB)
C     >	,'   ',(EQBUF(LOCSBP(ISB)+j),j=1,ichar(EQBUF(LOCSBP(ISB))))

	if (LOCSBR(ISB) .eq. 0)	then
	   JSB = JSB+1
	   j2 = 1
	   call SUBROUT(ISB,EQBUF(NB),ICH,j2)
C	   write(*,*)"Adding SBR No.",ISB," '",SBRNAM,"'  to eqns.inc"
	endif

	if (LOCSBR(ISB) .eq. -3) then
	   JSBP = JSBP+1
	endif
C JSBP is the global ordinal number of SBP
	if (LOCSBR(ISB) .eq. -2) then
	   JSBP = JSBP+1
	   JSB = JSB+1
	   j2 = 1
	   call SUBPROC(ISB,EQBUF(NB),ICH,JSBP,j2)
	endif

C ISB common ordinal number
C JSB sbr ordinal number
C JSBP xbr ordinal number
	NB=NB+JEQ+1
	if (NSBP .eq. 0)	goto	60	! No SBPs
C	if (JSB .lt. MSBR+MSBP)	goto	60	! List not finished
	if (NB  .ne. EQNPOS)	goto	60	! Write only once
	write(ICH ,'(A)')'C **** Synchronisation point'
	write(ICH ,'(A)')'      call	wait4all'
	write(ICH ,'(A)')'C **** Collect data from ShMem'
	j1 = 1				! NB
	j3 = 0				! ISB
	JSBP = 0
	do 65	j=1,NSBMX+14
	   J2 = ichar(EQBUF(j1))	! JEQ
	   if (J2 .lt. 0)	J2 = J2+256
	   if (J2 .eq. 0)	goto	66	! end
	   if (J2 .ge. 132)	then
	      j1 = j1+1
	      goto	65
	   endif
	   j3 = j3+1
	   M1 = NPOS(J2,EQBUF(j1+1),'(')-1
	   M2 = RNPOS(J2,EQBUF(j1+1),')')-1
	   M2 = min(M2,J2)
	   write(SBRNAM,'(132A)')(EQBUF(j1+jj),jj=1,M2),char(0)
C	   write(*,*)(EQBUF(j1+jj),jj=1,M1),j3
C	   write(*,*)(EQBUF(j1+jj),jj=1,M1),ISB,j1,j1+ichar(EQBUF(j1))
C	   write(*,'(3A)')'"',SBRNAM(1:M1),'"'
	   if (LOCSBR(j3).eq.-3 .or. LOCSBR(j3).eq.-2) 	then
	      JSBP = JSBP+1
	      write(ICH,101)JSBP
 101	      format('      if (IFSBP(',1I2,').ne.0) then')
	      if (M1 .eq. M2)	then
		 write(ICH,102)SBRNAM(1:M2),JSBP,JSBP
 102		 format('         call ot',A,'(',1I2,',IFSBP(',1I2,'))')
	      elseif (M2-M1 .le. 50)	then
		 write(ICH,103)SBRNAM(1:M1),SBRNAM(M1+1:M2),JSBP,JSBP
 103		 format('         call ot',A/
     &		 '     >',A,',',1I2,',IFSBP(',1I2,'))')
	      else
		 M3 = RNPOS(50,SBRNAM(M1+2:M2),',')
		 write(ICH,104)SBRNAM(1:M1+1),SBRNAM(M1+2:M1+M3+1),
     &				SBRNAM(M1+M3+2:M2),JSBP,JSBP
 104		 format('         call ot',A/'     >   ',A/
     &                  '     >   ',A,',',1I2,',IFSBP(',1I2,'))')
	      endif
	      write(ICH,105)JSBP
 105	      format('         IFSBP(',1I2,') = 0'/
     &	             '      endif')
	   endif
	   j1 = j1+j2+1
 65	continue
 66	continue
	goto	60

C Transport equations:
 68	continue
	JEQ = (JEQ-128)/4
	ITYPE = LEQ(JEQ)
	if(ITYPE.lt.0)	then
Control print: should never occur if the equation command line is present
	   write(*,*)'---->  ',EQNAM(JEQ),' string error?'
	   call	exit(1)
	endif

	if (ICH .eq. NCH+1)	goto	69
	IEQ = IEQ+1
	if (JEQ.eq.2 .and. IMPEI.eq.2)	IEQ = IEQ-1
	if (JEQ.eq.3 .and. IMPEI.eq.2)	IEQ = IEQ-1
C	write(*,*)JEQ,IEQ

C	if (LEQ(6).lt.0 .or. ICH.ne.NCH)call APPTMP(TMPBUF(NBEG(LNI)),6)

	if (JEQ.eq.1)			call NEEQN(IEQ,ITYPE,NCH)
	if (JEQ.eq.2 .and. IMPEI.eq.0)	call TEEQN(IEQ,ITYPE,NCH)
	if (JEQ.eq.3 .and. IMPEI.eq.0)	call TIEQN(IEQ,ITYPE,NCH)
	if (JEQ.eq.2 .and. IMPEI.eq.1)	then
	   write(NCH,'(A)')'C **** Electron temperature equation '
	   write(NCH,'(A)')'      call	markloc("TE equation"//char(0))'
	   if (NBEG(LDVI)+NBEG(LDVE)+NBEG(LDSI)+NBEG(LDSE) .gt. 0)
     >		write(NCH,123)IEQ
	   call TEEQN(IEQ,ITYPE,NCH)
	   IMPEI = IMPEI+1
	   ITYPE = LEQ(3)
	   IEQ = IEQ+1
	   write(NCH,'(A)')'C **** Ion temperature equation '
	   write(NCH,'(A)')'      call	markloc("TI equation"//char(0))'
	   call TIEQN(IEQ,ITYPE,NCH)
	   write(NCH,124)
	   if (NBEG(LDVI)+NBEG(LDVE)+NBEG(LDSI)+NBEG(LDSE) .gt. 0) then
	      write(NCH,126)IEQ-1
	   endif
	endif
 123	format(
     >'      NDV = 0'/
     >'      KDV =-1'/
     >' 2',1I2,'0 continue'/
     >'      KDV = KDV+1')
 124	format(
     >'      call NURTT(TE,TI,QE,QI,PETOT,PITOT,ND,WORK1)'/
     >'      if (ND1.lt.NA1) then'/
     >'         do j=ND1+1,NA1'/
     >'            QE(j)=QE(ND1)'/
     >'            QI(j)=QI(ND1)'/
     >'         enddo'/
     >'      endif')
	if (JEQ.eq.3 .and. IMPEI.eq.1)	then
	   write(NCH,'(A)')'C **** Ion temperature equation '
	   write(NCH,'(A)')'      call	markloc("TI equation"//char(0))'
	   if (NBEG(LDVI)+NBEG(LDVE)+NBEG(LDSI)+NBEG(LDSE) .gt. 0)
     >		write(NCH,123)IEQ
	   call TIEQN(IEQ,ITYPE,NCH)
	   IMPEI = IMPEI+1
	   ITYPE = LEQ(2)
	   IEQ = IEQ+1
	   write(NCH,'(A)')'C **** Electron temperature equation '
	   write(NCH,'(A)')'      call	markloc("TE equation"//char(0))'
	   call TEEQN(IEQ,ITYPE,NCH)
	   write(NCH,125)
	   if (NBEG(LDVI)+NBEG(LDVE)+NBEG(LDSI)+NBEG(LDSE) .gt. 0) then
	      write(NCH,126)IEQ-1
	   endif
	endif
 125	format(
     >'      call NURTT(TI,TE,QI,QE,PITOT,PETOT,NA,WORK1)'/
     >'      if (ND1.lt.NA1) then'/
     >'         do j=ND1+1,NA1'/
     >'            QE(j)=QE(ND1)'/
     >'            QI(j)=QI(ND1)'/
     >'         enddo'/
     >'      endif')
 126	format(
     >'      YH = HRO'/
     >'      do	j=1,NA'/
     >'         if (j.eq.NA) YH = HROA'/
     >'         YA = -1.6d-3*G11(j)'/
     >'         PDE(j) = YA*PDE(j)*(TEO(j)*TE(j+1)-TEO(j+1)*TE(j))'/
     >'     +	        +YA*DSE(j)*0.5*(NE(j)+NE(j+1))/YH'/
     >'     +	        *(TE(j+1)-TE(j)-TEO(j+1)+TEO(j))'/
     >'         PDI(j) = YA*PDI(j)*(TIO(j)*TI(j+1)-TIO(j+1)*TI(j))'/
     >'     +	        +YA*DSI(j)*0.5*(NI(j)+NI(j+1))/YH'/
     >'     +	        *(TI(j+1)-TI(j)-TIO(j+1)+TIO(j))'/
     >'      enddo'/
     >'      do	j=NA,2,-1'/
     >'         YA = 1./(HRO*VR(j))'/
     >'         PDE(j) = YA*(PDE(j-1)-PDE(j))'/
     >'         PDI(j) = YA*(PDI(j-1)-PDI(j))'/
     >'      enddo'/
     >'      YA = 1./(HRO*VR(1))'/
     >'      PDE(1) =-YA*PDE(1)'/
     >'      PDI(1) =-YA*PDI(1)'/
     >'      if (KDV.lt.NDV) goto 2',1I2,'0')
	if (JEQ.eq.4)	then
	   if (ITYPE.eq.1)	then
C	      write(*,*)"Current equation is solved"
	      call CUEQN(IEQ,ITYPE,NCH)
	   else
	      if     (NBEG(LCU).gt.0)	then
		 call	CUASN(IEQ,-1,NCH)	! prescribed CU = ...
	      elseif (NBEG(LMU).gt.0)	then
		 call	CUASN(IEQ,0,NCH)	! prescribed MU = ...
	      else
		 call	CUASN(IEQ,1,NCH)	! prescribed U = Const
	      endif
	    endif
	endif
C	write(*,*)EQNAM(JEQ),LEQ(JEQ),JEQ,JEQ-10
	if(JEQ.ge.10)	call FJEQN(IEQ,ITYPE,JEQ-10,NCH)
 69	continue
	NB = NB+1
	goto	60
 70	continue		! End ICH loop
 71	continue

 	write(NCH+1,'(A)')'      NA1N = NA1'
 	write(NCH+1,'(A)')'      NA1E = NA1'
 	write(NCH+1,'(A)')'      NA1I = NA1'
 	write(NCH+1,'(A)')'      NA10 = NA1'
 	write(NCH+1,'(A)')'      NA11 = NA1'
 	write(NCH+1,'(A)')'      NA12 = NA1'
 	write(NCH+1,'(A)')'      NA13 = NA1'
 	write(NCH+1,'(A)')'      NA15 = NA1'
 	write(NCH+1,'(A)')'      NA16 = NA1'
 	write(NCH+1,'(A)')'      NA17 = NA1'


 	write(NCH+1,'(A)')'      NA18 = NA1'
 	write(NCH+1,'(A)')'      NA19 = NA1'
 	write(NCH+1,'(A)')'      do	j=1,NA1'
	write(NCH+1,'(A)')'         NEO(j) = NE(j)'
	write(NCH+1,'(A)')'         NIO(j) = NI(j)'
	write(NCH+1,'(A)')'         TEO(j) = TE(j)'
	write(NCH+1,'(A)')'         TIO(j) = TI(j)'
	write(NCH+1,'(A)')'         do	jj=0,9'
	write(NCH+1,'(A)')'            FJO(j,jj) = FJ(j,jj)'
	write(NCH+1,'(A)')'         enddo'
	write(NCH+1,'(A)')'         VRO(j)=VR(j)'
 	write(NCH+1,'(A)')'      enddo'
	close(NCH+1)

	if (TCTR .ne. 0)	then
	   if (LOCSBR(TCTR) .eq. 0) call APPSBN(NCH,EQBUF,TCTR,1)
	endif
	write(NCH,'(A)')'      Jcall = 1'
	write(NCH,'(A)')'      if (JIT .eq. JEX) Jcall = IFSTEP(Jcall)'
	write(NCH,'(A)')'      if (Jcall .eq. 0) goto 2100'
	if (MIX .ne. 0)	then
	   if (LOCSBR(MIX) .eq. 0) call APPSBN(NCH,EQBUF,MIX,1)
	endif
	if (MEX .ne. 0)	then
	   if (LOCSBR(MEX) .eq. 0) call APPSBN(NCH,EQBUF,MEX,1)
	endif
	do	j=1,20
	   if (LOCSBR(j) .eq. 1)	then
	      call APPSBN(NCH,EQBUF,j,1)
	   endif
	enddo
	close(NCH)
C	call	exit(1)
	end
C======================================================================|
C Append a subroutine No. NSB to the channel NCH
C Input:
C	NCH	- channel No.
C	SBNAM	- the subroutine name
C	EQBUF(IEBUF <= 1000) Subroutine_call_string
C
C	EQBUF(IEBUF <= 1000) One position is used to transfer Eq. ID
C----------------------------------------------------------------------|
	subroutine	APPSBN(NCH,EQBUF,NSB,JPART)
	implicit	none
	character	EQBUF(*)*1
	integer		NPOS,NCH,ISB,NB,JEQ,NSB,JPART
C	call	APPTMP(EQBUF,6)		! Print buffer
	ISB=0
	NB=1
 1	JEQ=ICHAR(EQBUF(NB))
	if(JEQ.lt.0)	JEQ=JEQ+256
	if(JEQ.eq.0)		goto	69
	if(JEQ.lt.132)		goto	68
	NB=NB+1
	goto	1

 68	if (ISB .ge. 20)	goto	1
	ISB = ISB+1
	JEQ = ichar(EQBUF(NB))
	if (ISB.eq.NSB)	  call	SUBROUT(ISB,EQBUF(NB),NCH,JPART)
	NB=NB+JEQ+1
				goto	1
 69	return
	end
C======================================================================|
C Write "radout.tmp"
C Input:
C	NCH	- 
C	FNAME	- 
C	NROUT
C	FORMR
C	SCALER(,)
C	NAMER*4(96)
C----------------------------------------------------------------------|
	subroutine	RADOUT(NCH,NROUT,FORMR,SCALER,NAMER)
	implicit	none
C	include		'for/parameter.inc'
	include		'.srv/nambuf.inc'
	integer		NCH,NROUT,IERR,J,JJ,J1,JL,LUR,LOUT
	character*1	FORMR(75,NRW),SCALER(6,NRW),NAMER(NRW)*4,
     +			FNAME*14,LINTXT(133),LINFOR(133),STRFOR*132,
     +			CH3*3
	equivalence	(LINFOR(1),STRFOR)
	data		FNAME/'tmp/radout.tmp'/
	call	OPENWT(NCH,FNAME,0,IERR)
	if(IERR.gt.1)	then
	   write(*,*)'MODEL (MAIN): Cannot open file ',FNAME
	   call	exit(1)
	endif
	write(NCH,'(A)')
     +	'      call	markloc("radout.tmp"//char(0))'
	write(NCH,171)
 171	format (
C     1	'      if (int(XOUT+.49) .eq. 0)	then'/
C     2	'         NRADO = NAB'/
C     3	'      else'/
C     4	'         NRADO = NA1'/
C     5	'      endif'/
C     6	'      if (MOD10 .eq. 3) then'/
C     7	'         NRADO = NA1'/
C     8	'      endif'/
C     9		'      do 1 irado=1,NRADO'/
     +	        '      do 1 irado=1,NAB'/
     +	        '      J = irado')
	write(7,*)'=====   Radial profiles output   ====='
	write(7,*)' #  Scale  Name  Output expression'
C Exclude repeating formulae
	NOUT=0
	KEYOUT=1
 	do  56	J1=1,NROUT
C	LUR = ICHAR(FORMR(1,J1)(1:1))
C	write(*,'(A,I3,A,I3,132A1)')
C     +		"No.",J1,"   RHS_length",LUR,'		',
C     +		(FORMR(1,J1)(j+1:j+1),j=1,LUR)
C	   write(*,*)"Call from RADOUT:",KEYOUT,NCH
	   write(CH3,'(1I3)')J1
	   if (Tflag .ne. 0)	write(NCH,'(3A)')
     +	'      call add2loc("Radial channel ',CH3,'"//char(0))'
	   call ANLFML(NCH,FORMR(1,J1),LINTXT,LINFOR)
	   LUR=ICHAR(LINFOR(1))
C	   write(*,*)(LINFOR(1+j),j=1,LUR)
	   do	j=1,LUR
	      LINFOR(j) = LINFOR(j+1)
	   enddo
C Defining array variable at intermediate point
C	   if (STRFOR(1:7).eq.'RADIAL(' )  then !Special treatment needed
C	      JJ = index(STRFOR,',')+1
C	      write(NCH,'(2A)')'      YR=',STRFOR(JJ:LUR-1)
C	      write(STRFOR,'(132A1)')(LINFOR(J),J=1,JJ-1),'Y','R',')'
C	      LUR = JJ+2
C	   endif
	   JL = LUR
	   if (LUR .gt. 60)	   LUR = 60
	   write(NCH,107)J1,(LINFOR(J),J=1,LUR)
 107	   format('      ROUT(J,',1I3,')='/'     =	',60A1)
	   if (LUR .eq. JL)	goto	55
	   LUR = JL-60
	   JL = 60
 54	   continue
	   if (LUR.gt.60)	then
	      write(NCH,'(5X,1A2,62A1)')'> ',(LINFOR(J),J=1,60)
	      JL = JL+60
	      LUR = LUR-60
	      goto	54
	   else
	      write(NCH,'(5X,1A2,62A1)')'> ',(LINFOR(J),J=JL+1,JL+LUR)
	   endif
 55	   LOUT=ICHAR(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
	   if(SCALER(1,J1).eq.' ')	then
C		write(7,106)J1,(NAMER(J,J1),J=1,4),
	      write(7,106)J1,NAMER(J1),
     1		(LINTXT(J),J=1,LOUT)
	   else
	      write(7,104)J1,(SCALER(J,J1),J=1,6),
C     +		(NAMER(J,J1),J=1,4),(LINTXT(J),J=1,LOUT)
     +		NAMER(J1),(LINTXT(J),J=1,LOUT)
	   endif
 106	   format(1I3,9X,1A4,2X,114A1)
 104	   format(1I3,2X,6A1,1X,1A4,2X,114A1)
 56	continue
	write(NCH,108)
 108	format(' 1    continue')
	close(NCH)
	end
C======================================================================|
C Write timout.tmp
C Input:
C	NCH	- 
C	FNAME	- 
C	NTOUT	- Number of time channels
C	FORMT	- rhs as given in the model
C	SCALET	- Scale
C	NAMET	- Name
C----------------------------------------------------------------------|
	subroutine	TIMOUT(NCH,NTOUT,FORMT,SCALET,NAMET)
	implicit	none
C	include		'for/parameter.inc'
	include		'.srv/nambuf.inc'
	integer		NCH,NTOUT,IERR,J,JJ,J1,JL,LUR,LOUT
	character*1	FORMT(75,NRW),SCALET(6,NRW),NAMET(NRW)*4,CH3*3,
     +			FNAME*14,LINTXT(133),LINFOR(133),STRFOR*132
	equivalence	(LINFOR(1),STRFOR)
	data		FNAME/'tmp/timout.tmp'/
C	do   j1 =1,NTOUT
C 165	format(1X,78A1)
C	write(*,165)(FORMT(j,j1),j=1,75)
C	write(*,165)(scalet(j,j1),j=1,6)
C	write(*,165)(namet(j,j1),j=1,4)
C	enddo
	write(7,*)'=====   Time dependent values output   ====='
	write(7,*)' #  Scale  Name  Output expression'
C Check radial dependences
	KEYOUT=-1
	NOUT=0
	call	OPENWT(NCH,FNAME,0,IERR)
	if(IERR.gt.1)	then
	   write(*,*)'MODEL (MAIN): Cannot open file ',FNAME
	   call	exit(1)
	endif
	write(NCH,'(A)')
     +	'      call	markloc("timout.tmp"//char(0))'
	do  56	J1=1,NTOUT
	   write(CH3,'(1I3)')J1
	   if (Tflag .ne. 0)	write(NCH,'(3A)')
     +	'      call add2loc("Time channel ',CH3,'"//char(0))'
C	   JL=ICHAR(FORMT(1,J1))
C	   write(*,*)FORMT(2,J1)(1:JL),NCH,JL
C	   write(*,*)"Call from TIMOUT:",KEYOUT,NCH
	   call ANLFML(NCH,FORMT(1,J1),LINTXT,LINFOR)
	   LUR=ICHAR(LINFOR(1))
C	   write(*,*)(LINFOR(1+j),j=1,LUR)
	   do	j=1,LUR
	      LINFOR(j) = LINFOR(j+1)
	   enddo
C Write LINFOR(1:LUR) to "timout.tmp"
C	   write(NCH,107)J1,(LINFOR(J),J=1,LUR)
C	   write(*,  107)J1,(LINFOR(J),J=1,LUR)
C Defining array variable at intermediate point
C	   if (STRFOR(1:7).eq.'RADIAL(' )  then
! Special treatment is needed when
C	write(*,'(132(A1))')'"',(LINFOR(J),J=2,ichar(LINFOR(1))+1),'"'
C	      JJ = index(STRFOR,',')+1
C	      write(NCH,'(2A)')'      YR=',STRFOR(JJ:LUR-1)
C	      write(STRFOR,'(132A1)')(LINFOR(J),J=1,JJ-1),'Y','R',')'
C	      LUR = JJ+2
C	   endif
	   JL = LUR
	   if (LUR .gt. 60)	   LUR = 60
 107	   format('      TOUT(LTOUT,',1I3,')='/'     =	',60A1)
	   write(NCH,107)J1,(LINFOR(J),J=1,LUR)
	   if (LUR .eq. JL)	goto	55
	   LUR = JL-60
	   JL = 60
 54	   continue
	   if (LUR.gt.60)	then
	      write(NCH,'(5X,1A2,62A1)')'> ',(LINFOR(J),J=1,60)
	      JL = JL+60
	      LUR = LUR-60
	      goto	54
	   else
	      write(NCH,'(5X,1A2,62A1)')'> ',(LINFOR(J),J=JL+1,JL+LUR)
	   endif
 55	   LOUT=ICHAR(LINTXT(1))
	   do	j=1,LOUT
	      LINTXT(j) = LINTXT(j+1)
	   enddo
	   if(SCALET(1,J1).eq.' ')	then
	      write(7,106)J1,NAMET(J1),
     1		(LINTXT(J),J=1,LOUT)
	   else
	      write(7,104)J1,(SCALET(J,J1),J=1,6),
     +		NAMET(J1),(LINTXT(J),J=1,LOUT)
	   endif
 104	   format(1I3,2X,6A1,1X,1A4,2X,114A1)
 106	   format(1I3,9X,1A4,2X,114A1)
 56	continue
	close(NCH)
	end
C======================================================================|
C Write "ininam.tmp"
C	NCH
C	SBRBUF	- subroutine names
C	NSBR	- number of subroutines called from a model
C	NBMV	= NBEG(LMV)
C	NROUT	- number of radial channels
C	NTOUT	- number of time channels
C	NXOUT	- number of experimental channels
C	NWX(NRW)- window number for exp channel
C	NAMET
C	NAMER
C	SCALET
C	SCALER
C	NAMEX
C	ARXNAM(NARRX) the array is filled successively 
C		  = 0	 if X-array not used
C		Ord. No. if X-array is used on the lhs
C		The content of ARXNAM is changed when RADOUT, TIMOUT
C		are called
C----------------------------------------------------------------------|
	subroutine	ININAM(NCH,NSBR,SBRBUF,NBMV,NTOUT,NROUT,NXOUT,
     +	      NWX,NAMET,NAMER,NAMEX,SCALET,SCALER,IFTOUT,IFROUT,LOCSBR,
     +	      LEQ,NEQ)
C	subroutine	XARNAM(NCH,
C	implicit	none
C	integer		NCH,LEQ(*),NEQ,NARRX,ARXNAM(*),NARR,IERR,J

	implicit	none
	include		'.srv/nambuf.inc'
	integer		LEQ(*),NEQ
	character*1	SBRBUF(*),SCALER(6,NRW),SCALET(6,NRW),SBNAME*32,
     +			NAMER(NRW)*4,NAMET(NRW)*4,NAMEX(NRW)*6,FNAME*14
	integer		NCH,NSBR,NBMV,NTOUT,NROUT,NXOUT,NWX(*),IERR,JSB,
     +			RNPOS,J,J1,J2,JJ,IFTOUT(*),IFROUT(*),LOCSBR(*)
	data		FNAME/'tmp/ininam.tmp'/
C----------------------------------------------------------------------|
	call	OPENWT(NCH,FNAME,0,IERR)
	if(IERR.gt.1)	then
	   write(*,*)'MODEL (MAIN): Cannot open file ',FNAME
	   call	exit(1)
	endif

	write(NCH,'(A)')
     +	'      call	markloc("xar_usage"//char(0))'
	do	j=1,NARRX
	   if (ARXNAM(j).gt.0 .and. ARXNAM(j).le.NARR)
     t			write(NCH,111)j,ARXNAM(j)
C	   if (ARXNAM(j).gt.0 .and. ARXNAM(j).le.NARR)
C     t	   		write(*,*)ARRNAM(ARXNAM(j)),j,JARX
C	   if (ARXNAM(j) .gt. NARR)	write(*,*)ARRNAM(ARXNAM(j)-NARR)
C	   if (ARXNAM(j) .gt. 0)	write(*,*)ARRNAM(ARXNAM(j))
	enddo
 111	format('      ARXUSE(',I2,')=',I3)
 	write(NCH,112)(j,LEQ(j),j=1,NEQ)
 112	format(1('      LEQ(',I2,')=',I2))
 	write(NCH,113)
 113	format('      do j=1,NSBMX'/
     +         '         IFSBP(j) = 0'/
     +         '         LISTSB(j)=char(0)'/
     +         '         IFSBX(j) = 0'/
     +         '      enddo')
	j1 = 0
	do	j=1,NSBMX
	   if (LOCSBR(j).eq.-2 .or. LOCSBR(j).eq.-3)	then
	      j1 = j1+1				! Count number of SBPs
	      write(NCH,'(A,I2,A,I2)')'      IFSBX(',j1,') = ',j 
	   endif
	enddo

	write(NCH,'(A)')
     +	'      call	markloc("ininam.tmp"//char(0))'

 	write(NCH,'("      NSBR =",I3)')NSBR
	if (NBMV.gt.0)	write(NCH,101)
 101	format(
     t	'      if (NEQUIL.gt.0) then'/
     t	'      write(*,*)">>> Inconsistency warning:"'/
     t	'      write(*,*)"    Equilibrium solver is not compatible"'/
     t	'      write(*,*)"    with a vacuum rotational transform MV"'/
     t	'      write(*,*)'/
     t	'      endif')
	write(NCH,102)NTOUT,NROUT,NXOUT
 102	format('      NTOUT =',1I3/
     +	       '      NROUT =',1I3/
     +         '      NXOUT =',1I3)
	do	J1=1,NROUT
	   if(IFROUT(J1) .ne. 0)	write(NCH,103)J1
 103	   format('      MARKR(',1I3,')=1')
	   if(SCALER(1,J1).ne.' ')
     +		write(NCH,104)J1,(SCALER(J,J1),J=1,6)
 104	   format('      SCALER(',1I3,')=',6A1)
	   write(NCH,105)J1,NAMER(J1)
 105	   format('      NAMER (',1I3,')="',1A4,'"')
C	   write(NCH,105)J1,(NAMER(J,J1),J=1,4)
	enddo
	do	J1=1,NXOUT
	   write(NCH,106)J1,NAMEX(J1)
 106	   format('      NAMEX (',1I3,')="',1A6,'"')
	   write(NCH,107)J1,NWX(J1)
 107	   format('      NWINDX(',1I3,')=',1I3)
	enddo

	do	J1=1,NTOUT
	   if(IFTOUT(J1) .ne. 0)	write(NCH,108)J1
 108	   format('      MARKT(',1I3,')=1')
	   if(SCALET(1,J1).ne.' ')
     +		write(NCH,109)J1,(SCALET(J,J1),J=1,6)
 109	   format('      SCALET(',1I3,')=',6A1)
	   write(NCH,110)J1,NAMET(J1)
 110	   format('      NAMET (',1I3,')="',1A4,'"')
	enddo

C	call APPTMP(SBRBUF,NCH)
	JJ = 1
	do JSB=1,NSBR
	   J2 = ICHAR(SBRBUF(JJ))
C	   write(*,*)JJ,J2
C	   write(*,*)'SBRBUF: "',SBRBUF(JJ+1:JJ+J2),'"',J2,JJ
	   J1 = 0
	   if (LOCSBR(JSB).eq.-2 .or. LOCSBR(JSB).eq.-3)	then
	      J1 = RNPOS(J2,SBRBUF(JJ+1),'/')
	      if (J1 .eq. J2+1)	J1 = 0
C	      write(*,*)'SBpath: "',SBRBUF(JJ+1:JJ+J1),'"',J1
C	      write(*,*)'SBname: "',SBRBUF(JJ+J1+1:JJ+J2),'"'
C	      write(SBNAME,'(32A1)')SBRBUF(JJ+J1+1:JJ+J2),char(0)
C	      write(*,*)'SBNAME: "',SBNAME(1:J2-J1),'"',J2-J1
	   endif
 121	   format('      DTNAME(',1I2,'*4+24)="',8A)
	   write(NCH,121)JSB,
     +		(SBRBUF(JJ+J1+J),J=1,min(6,J2-J1)),'"//char(0)'
	   JJ = JJ+1+J2
	enddo

	JJ = 1
	J2 = 0
	do JSB=1,NSBR
	   J = ICHAR(SBRBUF(JJ))
C	   write(*,*)"Length =",j
	   write(NCH,122)JSB,LOCSBR(JSB)
 122	   format('      SIGNSB(',1I2,')=',1I2)
	   if (LOCSBR(JSB).eq.-2 .or. LOCSBR(JSB).eq.-3)	then
	      J2 = J2+1
 123	      format('      LISTSB(',1I2,')="',64A)
	   write(NCH,123)J2,(SBRBUF(JJ+J1),J1=1,min(32,J)),'"//char(0)'
	   endif
	   JJ = JJ+1+J
	enddo
	write(NCH,124)J2
 124	format('      NSBP =',1I3)
	if (J2 .ne. 0)	   write(NCH,125)
 125	format('      call	checkexec(NSBP,64,LISTSB)')
	close(NCH)
	return
C 88	write(7,*)'>>> Temporary buffer SBPATH length',IPBUF,' > 1000'
C	write(*,*)'>>> Temporary buffer SBPATH length',IPBUF,' > 1000'
C	IFSYNT = 1
C	goto	99
 99	write(7,*)'>>> Error in line:  "',(SBRBUF(JJ+J1),J1=1,J2),'"'
	write(7,*)'    Symbol "/" found in the subroutine header'
	write(*,*)'>>> Error in line:  "',(SBRBUF(JJ+J1),J1=1,J2),'"'
	write(*,*)'    Symbol "/" found in the subroutine header'
	write(*,*)'>>> Model analyzer STOP: Syntax error'
	call	exit(1)
	end
C======================================================================|
C ANLFML =============== Formula analisis =============================|
	subroutine	ANLFML(NCH,LINP,LINTXT,LINFOR)
C----------------------------------------------------------------------|
C Converts LINP into LINTXT (model.txt) & LINFOR (fortran text) 
C    according to the ASTRA interpretation rules.
C Output is written to the buffer 
C    DTVBUF (NCH<0, call from MAIN),
C    TMPBUF (NCH=0, call from MAIN)
C    or to a real file (NCH>0, call from DETVAR, RADOUT, TIMOUT)
C
C KEYOUT - type of external request;	    NCH -> 
C  	=-1,0,1    called from MAIN	    =-1 write fml to DTVBUF
C  					    = 0 write fml to TMPBUF
C  					    > 0 write fml to file NCH
C 	=-1   called from TIMOUT	    = 2   write to file
C	      suppresses some diagnostics
C  	= 1   called from DETVAR, RADOUT    = 2   write to file
C	      never used
C NFML,FNLNAM,NARR,ARRNAM,NFIN,FINNAM,NFNC,FNCNAM,NTMP,TMPNAM
C NSRV,SRVNAM,NCONS,CONNAM,NMAIVA,MAIVAR,IBUF,TMPBUF,IVBUF,DTVBUF
C
C OUTNAM(500) - list of formulas used in detvar.tmp, radout.tmp
C		and timout.tmp (subsequent usage destroys previous)
C		The buffer is used for excluding repeated formulas
C----------------------------------------------------------------------|
	implicit	none
	include 	'.srv/nambuf.inc'
	include 	'.srv/tmpbuf.inc'
	character*6	NAME,NAMEO,NAMEX
	character*1	LINE(133),LINP(80),LINTXT(133),LINFOR(133),CHAR1
	character*80	LINE80
	integer		NCH,NDEL,NPOS,JONEOF,IGAUSS,INDEX
	integer		LENS,LEN,IEND,STFUN,IRP,NB,NB1,NE,JO,j1,j,j2
	integer		RFLAG,IFNUM,JN6EOF,IN6EOF,LENGTF,LENGTH
	double precision YY
	save		IGAUSS
	data 		IGAUSS/0/
C----------------------------------------------------------------------|
C IEND:	=1 - at radial position, =2 - radially dependent
C	=3 - at the center,	 =4 - at the boundary
C STFUN:
C  -1 - .not.'standard function', 0 - GRAD[S], 1 - Vint,   2 - Iint
C   3 - STEP, XSTEP, ASTEP, RSTEP,   4 - GAUSS,   5 - FRMIN,  6 - FRMAX,
C   7 - RFMIN, 8 - RFMAX, 9 - RFVAL, AFVAL, RFVEX, AFVEX, RFVIN, AFVIN,
C  10 - SMOOTH (not used)
C  20 - ATR, ATX
C----------------------------------------------------------------------|
	IRP = 0
        CHAR1 = char(IRP)
	RFLAG = 0		! Default=False: don't add right ")"
	STFUN = -1
	write(LINTXT,'(1A1)')(CHAR1,j=1,133)
	write(LINFOR,'(1A1)')(CHAR1,j=1,133)
	NB=1
	LEN=ICHAR(LINP(1))
	if(LEN.eq.0)		go to 74
	do	j=1,LEN+1
	   LINE(j) = LINP(j)
	enddo
	call	UPCASE(LEN,LINE(2))
C	write(*,*)
C	write(*,*)'             123456789 123456789 123456789 123456789'
C	write(*,*)'Input LINE: "',(Line(j),j=2,len+1),'"'
C     >		,',	length =',LEN
C     >		,',	KEYOUT =',KEYOUT,',	NCH =',NCH

C----------------------------------------------------------------------|
 1	continue				! Starting the main loop
	NB1 = NB+1
	if (ichar(LINFOR(1)) .ne. 0)	then
C	write(*,*)"			Intermediate output"
C	write(*,'(132(A1))')'"',(LINFOR(J),J=2,ichar(LINFOR(1))+1),'"-'
C	write(*,*)'LINFOR:   "',(LINFOR(1+j2),j2=1,ichar(LINFOR(1))),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j2),j2=1,ichar(LINTXT(1))),'"'
C	write(*,*)
	endif

C Literal transfer:
	if (LINE(NB1) .eq. '"')	then
	   j2=NPOS(LEN-NB,LINE(NB1+1),'"')
C	   write(*,*)'Quoted string: "',(Line(j),j=2,len+1)
C     >	   ,'",  length =',LEN,j2,(Line(j),j=NB1+1,NB1+j2-1)
	   call LINSUM(LINTXT,LINE(NB1+1),j2-1)
	   call LINSUM(LINFOR,LINE(NB1+1),j2-1)
	   NE = NB1+j2
	   goto		70
	endif
C----------------------------------------------------------------------|
	if (IRP.eq.NB .and. RFLAG.eq.1)	then
	   call	LINSUM(LINFOR,')',1)
	   IRP = 0
	   RFLAG = 0
	endif

C	write(*,*)'Rest of the line: "',(LINE(NB+j),j=1,LEN+1-NB),'"'
C     &		,',  Rest length',LEN+1-NB
	LENS = NDEL(LEN+1-NB,LINE(NB1))-1
C	write(*,*)'Delimiter: "',LINE(NB1+LENS),'"',LEN+1-NB,LENS
	if (LENS.eq.0 .or. LENS.eq.LEN+2-NB)	goto	52
C (1) Find the next delimiter i.e. one of { + - * / ) ( , } and 
C (2) check if the term is a float number of the type ddd.E+dd
C     If yes, do not consider '+' or '-' as a delimiter
C	write(*,*)'Search next delimiter in "'
C     &		,(LINE(NB1+LENS+j),j=1,LEN-NB-LENS),'"'
	j2 = NDEL(LEN-NB-LENS,LINE(NB1+LENS+1))-1	! Next delimeter
	if (j2 .eq. LEN+2-NB-LENS)		goto	52
	CHAR1 = LINE(NB1+LENS)				! Current delimiter
	if (CHAR1.ne.'+' .and. CHAR1.ne.'-')	goto	52
	CHAR1 = LINE(NB1+LENS-1)
	if (CHAR1.ne.'E' .and. CHAR1.ne.'D')	goto	52
	CHAR1 = LINE(NB1+LENS-2)
	if (ichar(CHAR1) .eq. 47 )		goto	52
	if (ichar(CHAR1) .gt. 57 .or.		! Not a digit, not a "."
     >	    ichar(CHAR1) .lt. 46 )		goto	52
C	write(*,*)LENS,j2
	LENS = LENS+j2+1	! Extend the term till the next delimiter
C	LENS = LENS+j2		! Extend the term till the next delimiter

 52	continue
C	write(*,*)LENS
	if (LENS .le. 0)	goto	56
C	write(*,*)'      Term: "',(Line(j),j=nb1,nb+lens),'"'
C     >		,',	NB =',NB,',	lens =',LENS,',	stfun =',STFUN
	call	ASSTR(LENS,LINE(NB1),LINE80)

C	write(*,*)'      Line: "',(Line(j),j=nb1,nb+lens),'"'
	if (IFNUM(LINE(NB1),LENS) .eq. 0)	goto	56
	read(LINE80(1:LENS),*,ERR=56)YY
	write(LINE80,'(1P,E15.5,A1)')YY,char(0)
	j1 = 1
	j2 = 15
 53	continue
	if (LINE80(1:1) .eq. ' ')	then	! Remove leading blanks
	   LINE80(1:) = LINE80(2:)
	   j1 = j1+1
	   j2 = j2-1
	   goto	53
	endif
	j = index(LINE80,'E')
	if (j .ne. 0)	LINE80(j:j) = 'd'
 54	continue
	if (LINE80(j-1:j) .eq. "0d")	then	! Remove trailing zeros
	   LINE80(j-1:) = LINE80(j:)
	   j = j-1
	   j2 = j2-1
	   goto	54
	endif
C 55	continue
	if (LINE80(j2-1:j2-1).eq."0" .and. LINE80(j:j).eq."d")	then
	   LINE80(j2-1:) = LINE80(j2:)		! Remove the 1st (of two) zero
	   j2 = j2-1
	endif
	if (LINE80(j:j+1) .eq. "d+")	then
	   LINE80(j+1:) = LINE80(j+2:)		! Remove '+' in exponent
	   j2 = j2-1
	endif

C	write(*,'(133A)')'Old term: "',(LINE(NB+j),j=1,LENS),'"'
C	write(*,'(133A)')'New term: "',LINE80(1:j2),'"'
C	write(*,'(133A)')'Old line: "',(LINE(j+1),j=1,LEN),'"'
C F2D

C	if (LEN+1 .ge. NB1+LENS)	then
	if (LENS .gt. j2)	then
	   do	j=NB1+LENS,LEN+1
	      LINE(j+j2-LENS) = LINE(j)
	   enddo
	else
	   do	j=LEN+1,NB1+LENS,-1
	      LINE(j+j2-LENS) = LINE(j)
	   enddo
	endif
	do	j=1,j2
	   LINE(NB+j) = LINE80(j:j)
	enddo
C	write(*,'(133A)')'New line: "',(LINE(j+1),j=1,LEN+j2-LENS),'"'
C	NE = NE+j2-LENS
	if (IRP .ne. 0)	IRP = IRP+j2-LENS
	LEN = LEN+j2-LENS
	LINE(1) = char(LEN)
	LENS = j2
 56	continue
	NE  = NB1+LENS

C NB - old delimiter position, NE - new delimiter position
C LENS - string length between the delimiters
	IEND=2
	if(LINE(NE-1).eq.'C' .and. STFUN.lt.3)		IEND=3
	if(LINE(NE-1).eq.'B' .and. STFUN.lt.3)		IEND=4
	if(  LINE(NE).eq.'(' .and. NE.le.LEN )		IEND=1
	if(IEND.eq.1 .and. RFLAG.eq.0)		then
C	   write(*,*)'Rest string:   "',(Line(j),j=NE,len+1),'"'
	   j  = NE
	   j1 = NE-1
 57	   j1 = j1+NPOS(LEN-j1,LINE(j1+1),')')
	   j  = j+NPOS(LEN-j+1,LINE(j+1),'(')
	   if (j .lt. j1)	goto	57
	   IRP = j1
C j > LEN+1  - no "(" found
C j < IRP
	   if (IRP .eq. LEN+2)	write(*,*)"Syntax error"
C	   write(*,*)'Parsed clause: "',(Line(j),j=NE,IRP),'"'
C	   write(*,*)'Position of "(" ',NE
C     >		 ,'   Position of ")" ',IRP
C     >		 ,NE,IRP,'    Check: "',LINE(NE),LINE(IRP),'"'
Check if ',' is present in the clause:	j = NPOS(IRP-NE,LINE(NE+1),',')
	endif
C IEND:	=1 - at radial position, =2 - radially dependent
C	=3 - at the center,	 =4 - at the boundary

	if (LENS.eq.0)		goto	60
C	write(*,*)'Rest string:   "',(Line(j),j=NE,len+1),'"',IEND

C	if (IEND.ne.2)		goto	2
C Preliminary:
C        to be used for prescribing type of radial argument
C Define type of the radial argument (Default rho_tor)
C	if (LINE(NE-1)//LINE(NE).eq.'M)')	then
C	   call SETNAM(LENS,LINE(NB1),NAME)
C	   if(JONEOF(NAME,NARR,ARRNAM).gt.0) goto 2
C	   write(*,*)'    Before: "',(Line(j),j=2,len+1),'"'
C     >	    ,',	length =',LEN,',  Argument "',(LINE(j),j=NB1,NE-1),'"'
C	   do j=NE-1,LEN
C	      LINE(j) = LINE(j+1)
C	   enddo
C	   NE = NE-1
C	   LEN = LEN-1
C	   LENS = LENS-1
C	endif

	call SETNM6(LENS,LINE(NB1),NAME)
	if (IEND.gt.2)	then
Cut off final "C" or "B"
		NAMEO=NAME
		LENS=LENS-1
		call SETNM6(LENS,LINE(NB1),NAME)
C Override cutting if both "name" and "nameC"/"nameB" exist
		if ( JONEOF(NAME, NFML,FMLNAM).gt.0 .and.
     +		     JONEOF(NAMEO,NFML,FMLNAM).gt.0 )	goto	60
	endif

 3	continue
C 4	continue
C 5	continue
C Few syntax checks:
	if(STFUN.ge.3 .and. STFUN.le.4 .and. LINE(NE).eq.')')	then
C the argument cannot be an array,   otherwise -> 83
	   if(JONEOF(NAME,NARR,ARRNAM).gt.0) goto 83
C 2nd argument will be added -> 15
	   goto	15
	endif
	if(STFUN.ge.5 .and. STFUN.le.9)	then
C 		     1st argument must be ASTRA array,   otherwise -> 84
	   if(LINE(NB).eq.'('.and.JONEOF(NAME,NARR,ARRNAM).eq.0) goto 84
	   if(LINE(NE).eq.')')	goto  15
	endif
	if(STFUN.ge.3)	goto	82

C	write(*,*)JONEOF(NAME,NFIN,FINNAM),JONEOF(NAME,NARR,ARRNAM)
C     >	,JONEOF(NAME,NFML,FMLNAM),JONEOF(NAME,NFNC,FNCNAM)

	if(STFUN.ge.0..and.JONEOF(NAME,NARR,ARRNAM).le.0)go to 83
	if(JN6EOF(NAME,NFIN,FINNAM).gt.0)	go to 10
	if(JN6EOF(NAME,NARR,ARRNAM).gt.0)	go to 15
	if(JN6EOF(NAME,NFML,FMLNAM).gt.0)	go to 30
	if(JN6EOF(NAME,NFNC,FNCNAM).gt.0)	go to 40
	if(NAME.eq.'GRAD')			STFUN=0
	if(NAME.eq.'GRADS')			STFUN=0
	if(NAME.eq.'VINT')			STFUN=1
	if(NAME.eq.'IINT')			STFUN=2
C All *STEP are treated similarly and can be linked to STFUN=3
C	if(NAME.eq.'STEP')			STFUN=3
	if(NAME.eq.'ASTEP')			STFUN=3
	if(NAME.eq.'XSTEP')			STFUN=3
	if(NAME.eq.'RSTEP')			STFUN=3
	if(NAME.eq.'FABOX')			STFUN=3
	if(NAME.eq.'FRBOX')			STFUN=3
	if(NAME.eq.'FXBOX')			STFUN=3
	if(NAME.eq.'GAUSS'.and.IGAUSS.eq.0)	STFUN=4
	if(NAME.eq.'GAUSS'.and.IGAUSS.ne.0)	goto 85
	if(NAME.eq.'GAUSS'.and.IGAUSS.eq.0)    IGAUSS=1
	if(NAME.eq.'FRMIN')			STFUN=5
	if(NAME.eq.'FRMAX')			STFUN=6
	if(NAME.eq.'RFMIN')			STFUN=7
	if(NAME.eq.'RFMAX')			STFUN=8
	if(NAME.eq.'RFVAL')			STFUN=9
	if(NAME.eq.'AFVAL')			STFUN=9
	if(NAME.eq.'RFVEX')			STFUN=9
	if(NAME.eq.'AFVEX')			STFUN=9
	if(NAME.eq.'RFVIN')			STFUN=9
	if(NAME.eq.'AFVIN')			STFUN=9
C	if(NAME.eq.'SMOOTH')			STFUN=9
	if(NAME.eq.'ATR' .or. NAME.eq.'ATX')	STFUN=20
	if(NAME.eq.'QETB' .or. NAME.eq.'QITB' .or. NAME.eq.'QNNB')
     .							goto	60
	if (NAME(1:3).eq."QFF" .and. NAME(5:5).eq."B")	goto	60
	if(STFUN.ge.0)					goto	80
	JO=IN6EOF(NAME,NTMP,TMPNAM)
C	write(*,*)'The name: "',name,'"'
C     >	,',     Place in buffer: ',JO,',   Used in function? ',STFUN
C	write(*,*)'Defined? ',NBEG(JO)	! ,' [==0 - no, !=0 - yes]'
C     >	,'             Radial type: ',IEND
	if(JO) 60,60,50			! Standard name? if (yes) -> 50
C If NAME - one of FML and one of FNC
 10	go to (11,12,13,14),IEND
 11	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,'R(RFA',5)
	RFLAG = 1
				go to 25
 12	continue
C Formula appending when called from DETVAR or RADOUT 
	if(KEYOUT.lt.0
C Exclude some special names from warning:
     +	  .and. NAME .ne. 'TITER '
     +	  .and. NAME .ne. 'THQ99 '
     +	  .and. NAME .ne. 'TLQ97 '
     +	  .and. NAME .ne. 'TAUNA '
     +	  .and. NAME .ne. 'TAU89 '
     +				   ) then
		write(*,105)NAME
		write(7,105)NAME
 105		format(' >>> Formula ',1A6,
     .			' is used in a radially independent expression')
	endif
	if (KEYOUT.gt.0 .and. NCH.lt.0)		then	! Never happens
	    write(*,*)NCH,LENS,NFML,KEYOUT,
     >		JN6EOF(NAME,NFML,FMLNAM),'	"',NAME,'"',NOUT
	endif
	if (KEYOUT.eq.1 .and. NCH.eq.0)	then
C	   write(*,*)"KEYOUT =",KEYOUT,'	"',NAME,'"',
C     >		JN6EOF(NAME,NFML,FMLNAM)
	   if (NAME.eq."PEI   " .or. NAME.eq."PEICL ")	KEYOUT = 2
	endif
 	if (NCH .gt. 0)	call APPFML(NCH,LENS,NAME)
	if (NCH .eq. 0)	call INCFML(LENS,NAME,IBUF,TMPBUF(IBUF))
	if (NCH .eq.-1)	call TMPFML(LENS,NAME,IVBUF,DTVBUF(IVBUF))
 122	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,'(r)',3)
	call LINSUM(LINFOR,LINE(NB1),LENS)
				go to 25
 13	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,'(0)',3)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,'R(0.d0)',7)
				go to 25
 14	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,'(a)',3)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,'R(ROC)',6)
				go to 25
C If NAME - one of ARR
 15	if(STFUN.lt.0)		go to 20
C after standard functions GRAD[S],*STEP,GAUSS,FRMIN,FRMAX,RFMIN,RFMAX
	if(STFUN.eq.0 .or. STFUN.ge.3 .and. STFUN.le.9)
     >					go to 16
C after standard functions Vint or Iint
	if(IEND.eq.1.or.IEND.eq.3)	go to 81
	if(IEND.eq.4)		go to 19
	if(LINE(NE).eq.')')	go to 17
	if(LINE(NE).eq.',')	go to 18
				go to 81
 16	if(IEND.eq.1)		go to 81
C Add ",j)" to functions GAUSS, ASTEP, etc.
	if(IEND.eq.2)	then
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,')(r)',4)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	if(STFUN.le.4)	call LINSUM(LINFOR,',J)',3)
C	if(STFUN.ge.5 .and. STFUN.le.9)	call LINSUM(LINFOR,'(1))',4)
	if(STFUN.ge.5 .and. STFUN.le.9) call LINSUM(LINFOR,')',1)
				go to 70
			endif
	if(IEND.eq.3)	then
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,')(0)',4)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,',0)',3)
				go to 70
			endif
	if(IEND.eq.4)	then
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,')(a)',4)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,',NA1)',5)
				go to 70
			endif
C radial profile
 17	continue
C VINT is attributed to a shifted grid
C	if (NCH .gt. 0)	write(NCH,'(      YR=J*HRO)')
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,')(r)',4)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,',J*HRO)',7)
				go to 70
C at radial position
 18	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,')(',2)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,',RFA(',5)
	RFLAG = 1
C	write(*,*)"??"
C	write(*,'(132(A1))')(LINFOR(J),J=2,ichar(LINFOR(1))+1)
				go to 70
C total integral
 19	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,')(a)',4)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,',ROC)',5)
				go to 70
 20	continue
	if (NAME(LENS:LENS).ne.'X')	goto	26	! Store X-array names
	JARX = JARX+1			! requested in the model
	j2 = JONEOF(NAME,NARR,ARRNAM)
	if (JARX .eq. 1)	then
	   ARXNAM(JARX) = j2
C	   write(*,*)"ANLFML: ",j2,' ',NAME,JARX
	   goto	26
	endif

	do	j1=1,JARX-1			! Exclude repeated LHS names
	   if (ARXNAM(j1) .eq. j2+NARR)	then
	      JARX = JARX-1			! The X-array has already
	      goto	26			! appeared on a LHS
	   endif
	enddo
 201	ARXNAM(JARX) = j2
C	write(*,*)"ANLFML: ",j2,' ',NAME,JARX
	do	j1=1,JARX-1			! Exclude repeated RHS names
	   if (ARXNAM(j1) .eq. j2)	then
	      ARXNAM(JARX) = -1
	      JARX = JARX-1
	      goto	26
	   endif
	enddo

 26	continue
	go to (21,22,23,24),IEND
C array at radial position
 21	continue
C	if(NAME(LENS:LENS).eq.'X')	write(*,*)LENS,' atr ',NAME
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,LINE(NE),1)
	call LINSUM(LINFOR,'RADIAL(',7)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,',RFA(',5)
	RFLAG = 1
C This piece is used when calling from MAIN (call No.6) 
C				and from TIMOUT
C	write(*,*)'"',(LINFOR(1+j),j=1,ichar(LINFOR(1))),'"'
C     >			,ICHAR(LINFOR(1))
C	write(*,*)"NE,LEN",NE,LEN,NB1,LENS
C	write(*,*)LINE(NE),LINE(NB1)
				go to 70
C radially dependent
 22	if(KEYOUT.lt.0)		then
	write(*,106)NAME
	write(7,106)NAME
 106	format(
     1	' >>> Variable ',1A6,
     2	' is used in a radially independent expression')
				endif
C	if(NAME(LENS:LENS).eq.'X')	write(*,*)LENS,' rad ',NAME
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,'(r)',3)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,'(J)',3)
				go to 25
C evaluate at the center
 23	continue
C	if(NAME(LENS:LENS).eq.'X')	write(*,*)LENS,' ctr ',NAME
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,'(0)',3)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,'(1)',3)
				go to 25
C evaluate at the boundary
 24	continue
C	if(NAME(LENS:LENS).eq.'X')	write(*,*)LENS,' bnd ',NAME
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,'(a)',3)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,'(NA1)',5)
				go to 25
C If NAME - one of FML
 30	go to (31,12,31,31),IEND
 31	write(7,102)NAME(1:INDEX(NAME,' ')-1)
	write(*,*)'>>> Error in line  "',(LINE(J),J=2,LEN+1),'"'
	write(*,102)NAME(1:INDEX(NAME,' ')-1)
 102	format(' >>> This use of "',A,'" is not permitted now')
	go to 77
C If NAME - one of FNC
 40	continue
	goto (11,42,13,14),IEND
 42	if(KEYOUT.lt.0)		then
	write(*,107)NAME
	write(7,107)NAME
 107	format(
     1	' >>> Function ',1A6,
     2	' is used in a radially independent expression')
				endif
C	if (NCH .gt. 0)	write(NCH,103)
C 103	format('      YR=RHO(j)')
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINTXT,'(r)',3)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINFOR,'R(RHO(j))',9)
				go to 25
C If NAME - one of TMP
 50	continue
C	write(*,*)'The name: "',name,'"'
C     >	,',     Place in buffer: ',JO,',   Used in function? ',STFUN
C	write(*,*)'Defined? ',NBEG(JO)	! ,' [==0 - no, !=0 - yes]'
C     >	,'             Radial type: ',IEND
	if(NBEG(JO).gt.0)	go to 59	! If defined -> 59
C Processing boundary conditions:
	if(     JO.eq.LIPL  .or.JO.eq.LLEXT .or.JO.eq.LUEXT
     +	   .or. JO.eq.LNEB  .or.JO.eq.LTEB  .or.JO.eq.LTIB
     +	   .or. JO.eq.LQEB  .or.JO.eq.LQIB  .or.JO.eq.LQNB
     +	   .or. JO.eq.LROE  .or.JO.eq.LROI  .or.JO.eq.LRON
     +	   .or. JO.eq.LQETB .or.JO.eq.LQITB .or.JO.eq.LQNNB) goto 60
	do	j2=0,9
	   if ( JO.eq.LRO(j2)  .or. JO.eq.LFB(j2) .or.
     +		JO.eq.LQFB(j2) .or. JO.eq.LQFFB(j2) )	goto	60
	enddo
 104	format(' >>> Warning: Variable "',A,'" is used but not defined')
	write(*,104)NAME(1:lengtf(NAME))
	write(7,104)NAME(1:lengtf(NAME))
 59	if (   JO.ne.LIPL.and.JO.ne.LLEXT.and.JO.ne.LUEXT
     +	  .and.JO.ne.LRON.and.JO.ne.LROE .and.JO.ne.LROI
     +	  .and.JO.ne.LRO1.and.JO.ne.LRO2 .and.JO.ne.LRO3 )	goto 15
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINFOR,LINE(NB1),LENS)
				go to 25
C If NAME - Vint or Iint
C IEND:	=1 - at radial position, =2 - radial profile
C	=3 - at the center,	 =4 - at the boundary
 80	continue
	go to (82,81,81,81),IEND
 81	write(*,*)'>>> ERROR in ',(LINE(J),J=2,LEN+1)
	write(*,*)'>>> ERROR: Incorrect use of the function ',NAME
	write(7,*)'>>> ERROR: Incorrect use of the function ',NAME
	go to 77
 83	write(*,*)'>>> ERROR in ',(LINE(J),J=2,LEN+1)
	if(STFUN.eq.1)	then
	write(*,*)'>>> ERROR: Incorrect use of the function Vint'
	write(7,*)'>>> ERROR: Incorrect use of the function Vint'
			endif
	if(STFUN.eq.2)	then
	write(*,*)'>>> ERROR: Incorrect use of the function Iint'
	write(7,*)'>>> ERROR: Incorrect use of the function Iint'
			endif
	if(STFUN.eq.3)	then
	write(*,*)'>>> Incorrect use of the radial function STEP'
	write(7,*)'>>> Incorrect use of the radial function STEP'
			endif
	if(STFUN.eq.4)	then
	write(*,*)'>>> Incorrect use of the radial function GAUSS'
	write(7,*)'>>> Incorrect use of the radial function GAUSS'
			endif
	go to 77
 84	write(*,*)'>>> ERROR in ',(LINE(J),J=2,LEN+1)
C argument must be ASTRA array
	if(STFUN.eq.5)	then
	write(*,*)'>>> Incorrect argument of the function FRMIN'
	write(7,*)'>>> Incorrect argument of the function FRMIN'
			endif
	if(STFUN.eq.6)	then
	write(*,*)'>>> Incorrect argument of the function FRMAX'
	write(7,*)'>>> Incorrect argument of the function FRMAX'
			endif
	if(STFUN.eq.7)	then
	write(*,*)'>>> Incorrect argument of the function RFMIN'
	write(7,*)'>>> Incorrect argument of the function RFMIN'
			endif
	if(STFUN.eq.8)	then
	write(*,*)'>>> Incorrect argument of the function RFMAX'
	write(7,*)'>>> Incorrect argument of the function RFMAX'
			endif
	if(STFUN.eq.9)	then
	write(*,*)'>>> Incorrect argument of the function R[A]FVAL'
	write(7,*)'>>> Incorrect argument of the function R[A]FVAL'
			endif
	go to 77
 85	continue
	write(*,*)'>>> Multiple use of the radial function GAUSS'
     >			,' is not allowed'
	write(7,*)'>>> Multiple use of the radial function GAUSS'
     >			,' is not allowed'
	go to 77
 82	continue
	call LINSUM(LINTXT,LINE(NB1),LENS)
	call LINSUM(LINFOR,LINE(NB1),LENS)
	call LINSUM(LINTXT,LINE(NE),1)
	call LINSUM(LINFOR,LINE(NE),1)
C	write(*,'(132(A1))')(LINFOR(J),J=2,ichar(LINFOR(1))+1),'"'
	if(IEND.eq.2 .and. STFUN.eq.20)	STFUN = -1
				go to 75
C If NAME - unknown

C right: number

 60	if(IEND.le.2)		go to 61
C A formula name "nameC" or "nameB" was used but "name" does not exist
	LENS=LENS+1
	IEND=2
	NAME=NAMEO
				go to 3
 61	continue
	if(IN6EOF(NAME,NSRV,SRVNAM).gt.0)	go to 63
	if(IN6EOF(NAME,NCONS,CONNAM).gt.0)	go to 63
	if(IN6EOF(NAME,NMAIVA,MAIVAR).gt.0)	go to 63
	do	j=1,5
	    if (NAME(j:j+1).eq.'X ')	then
		NAMEX(j:j) = ' '
	    else
		NAMEX(j:j) = NAME(j:j)
	    endif
	enddo
	if (NAME(6:6).eq.'X')	then
	    NAMEX(6:6) = ' '
	else
	    NAMEX(6:6) = NAME(6:6)
	endif
	if(IN6EOF(NAMEX,NMAIVA,MAIVAR).gt.0)	go to 63
	if(IFNUM(LINE(NB1),LENS).eq.0)	then
C Warning!        it is often LENS==1
	write(*,*)'>>> Warning: unrecognized variable ',
     >		(LINE(j),j=NB1,NE-1)
	write(7,*)'>>> Warning: unrecognized variable ',
     >		(LINE(j),j=NB1,NE-1)
					endif
 63	call LINSUM(LINTXT,LINE(NB1),NE-NB)
	call LINSUM(LINFOR,LINE(NB1),NE-NB)
C	write(*,*)'LINFOR:   "',(LINFOR(1+j2),j2=1,NE-NB-1),'"'
C	write(*,*)'LINTXT:   "',(LINTXT(1+j2),j2=1,NE-NB-1),'"'
				go to 70
 25	continue
C write delimeter at the end
C	write(*,*)LINE(NE)
	call LINSUM(LINTXT,LINE(NE),1)
	call LINSUM(LINFOR,LINE(NE),1)
 70	STFUN=-1
 75	continue
	if(NE-LEN-1)		71,72,73
 71	NB=NE
				go to 1
 73	continue
	write(LINTXT(1),'(1A1)')	char(ICHAR(LINTXT(1))-1)
	write(LINFOR(1),'(1A1)')	char(ICHAR(LINFOR(1))-1)
 72	continue
	if (IRP.eq.NE .and. RFLAG.eq.1)	call	LINSUM(LINFOR,')',1)
C	write(*,*)'   IRP =',IRP,'     RFLAG =',RFLAG
C	write(*,'(A$)')'Output LINE: '
C	write(*,'(132(A1))')'"',(LINFOR(J),J=2,ichar(LINFOR(1))+1),'"'
C	write(*,'(132(A1))')'"',(LINTXT(J),J=2,ichar(LINFOR(1))+1),'"'
	return
C Empty string
 74	write(LINFOR(1),'(1A1)')	char(2)
	call COPY(2,'0.',LINFOR(2))
	return
 77	call	exit(1)
	end
C======================================================================|
