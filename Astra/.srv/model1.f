C SUBROUT ================ External subroutine call ====================|
	subroutine	SUBROUT(NEQ,LINE,NCH,JPART)
	implicit none
	integer		NEQ,NPOS,J,NCH,M,M1,JPART
	character	LINE(132)*1,SBRNAM*16
	M  = ICHAR(LINE(1))+1
C Enable a special treatment for some subroutines
	M1 = NPOS(M-1,LINE(2),'(')-1
	M1 = min(6,M1,M-1)
	write(SBRNAM,'(6A1)')(LINE(J),J=2,M1+1)
	call	UPCASE(M1,SBRNAM)
C	write(*,*)"M =",M,'   "',(LINE(j),j=2,M)	! Full input line
C	write(*,*)"M1 =",M1,'   "',SBRNAM(1:M1)		! Sbr name only

	if (SBRNAM(1:M1).eq."MIXINT" .or. SBRNAM(1:M1).eq."MIXEXT") then
	   write(NCH,100)NEQ,NEQ,NEQ,NEQ,NEQ,').gt.0.0.and.JIT.eq.JEX'
	elseif (SBRNAM(1:M1).eq."TSCTRL")	then
! Special treatment for TSCTRL
	   write(NCH,100)NEQ,NEQ,NEQ,NEQ,NEQ,').gt.0.0.and.JIT.eq.JEX'	
	else
	   write(NCH,100)NEQ,NEQ,NEQ,NEQ,NEQ,').gt.0.0'
	endif

 100	format('C **** External subroutine'/
     1	'      Jcall=0'/
     2	'      if(KEY.ne.0.and.ABS(KEY-DTEQ(4,',1I2,')).lt.0.1)',
     3	' Jcall=1'/
     4	'      if(TIME.ge.DTEQ(2,',1I2,').and.TIME.le.DTEQ(3,',1I2,')',
     5	'.and.'/
     6	'     > TIME-TEQ(',1I2,')+.1E-7-DTEQ(1,',1I2,A,
     7	') Jcall=1')
	if (JPART.eq.0)	write(NCH,'(A)')
     1	'      if (IPART.eq.0) Jcall = 0'

	if (SBRNAM(1:M1) .eq. "STRAHL")	then
		write(NCH,102)NEQ,NEQ,(LINE(J),J=9,M)
 102		format(
     1	'      if (Jcall.eq.1) then'/
     2	'      TEQ(',1I2,')=TIME'/
     3	'      call	ADDTIME(CPT)'/
     4	'      call	markloc('/
     5	'     >"subroutine STRAHL"//char(0))'/
     6	'      call'/
     7	'     >STRAHL(DTEQ(1,',1I2,'),',50A1)
		write(NCH,'(A,I2,A)')
     &		'      call	ADDTIME(CPTSBR(',NEQ,'))'
		write(NCH,'(A)')'      endif'
		return
	endif

	write(NCH,101)NEQ,SBRNAM(1:M1)
 101	format(
     1	'      if(Jcall.eq.1) then'/
     2	'      TEQ(',1I2,')=TIME'/
     3	'      call	ADDTIME(CPT)'/
     4	'      call	markloc("subroutine ',A,'"//char(0))'/
     5	'      call')
C Splitting long lines
	M1 = 2
 10	if(M-M1 .le. 66)	then
	   write(NCH,110)(LINE(J),J=M1,M)
 110	   format('     >',67A1)
	else
	   write(NCH,110)(LINE(J),J=M1,67)
	   M1 = M1+66
	   goto	10
	endif
	write(NCH,'(A,I2,A)')
     &		'      call	ADDTIME(CPTSBR(',NEQ,'))'
	write(NCH,'(A)')'      endif'
	end
C SUBPROC ================ External subroutine call ====================|
	subroutine	SUBPROC(NEQ,LINE,NCH,JSBP,JPART)
	implicit none
C JPART == 0 the record goes to detvar.tmp
C JPART == 1 the record goes to init.inc or eqns.inc
	integer		NEQ,NCH,M,JSBP,JPART,J,J1,NPOS,RNPOS
	character	LINE(132)*1,SBRNAM*132
	M = ICHAR(LINE(1))
	J1 = NPOS(M,LINE(2),'(')
	if (J1 .eq. M+1) J1 = M
	write(SBRNAM,'(132A1)')(LINE(J),J=2,M+1)
C	call	UPCASE(JSBP,SBRNAM)
C	write(*,*)"M =",M,'   "',(LINE(j),j=2,M+1)  ! Full input line

	write(NCH,100)NEQ,NEQ,NEQ,NEQ,NEQ,').gt.0.0'
 100	format('C **** External subprocess'/
     1	'      Jcall=0'/
     2	'      if(KEY.ne.0.and.ABS(KEY-DTEQ(4,',1I2,')).lt.0.1)',
     3	' Jcall=1'/
     4	'      if(TIME.ge.DTEQ(2,',1I2,').and.TIME.le.DTEQ(3,',1I2,')',
     5	'.and.'/
     6	'     > TIME-TEQ(',1I2,')+.1E-7-DTEQ(1,',1I2,A,') Jcall=1'/
     7	'      if (Jcall.eq.1) Jcall = ifipc()')
	if (JPART.eq.0)	write(NCH,'(A)')
     1	'      if (IPART.eq.0) Jcall = 0'
	write(NCH,101)NEQ,JSBP,NEQ,SBRNAM(1:J1)
 101	format(
     1	'      if (Jcall.eq.1) then'/
     2	'      TEQ(',1I2,')=TIME'/
     3	'      IFSBP(',1I2,') =',1I2/
     4	'      call	ADDTIME(CPT)'/
     5	'      call	markloc("subroutine ',A,'"//char(0))'/
     6	'      call')
	if (J1 .eq. M)	then		! No parameters
	   write(NCH,102)SBRNAM(1:M),JSBP
 102	   format(
     1	   '     >to',A,'(',1I2,')')
	elseif (M .le. 60)	then
	   write(NCH,103)SBRNAM(1:M-1),JSBP
 103	   format('     >to',A,',',1I2,')')
	else				! Long parameter list
	   J = RNPOS(60,LINE(2),',')
	   write(NCH,104)SBRNAM(1:J1),SBRNAM(J1+1:J),
     &		SBRNAM(J+1:M-1),JSBP
 104	   format('     >to',A,/,'     >',A,/,'     >',A,',',1I2,')')
	endif
	write(NCH,105)JSBP
 105	format(
     2	'      call	letsbp(',1I2,')'/
     3	'      endif')
	end
C NEEQN ============= Density equation ================================|
	subroutine NEEQN(NEQ,ITYPE,NCH)
C ITYPE < 0   no equation, no time dependence
C ITYPE = 0   NE:assigned
C ITYPE = 1   NE:equation
	include 'for/parameter.inc'
	include '.srv/tmpbuf.inc'
	integer	j,jj
	jj = 0
	do	j=2,7		! 'DN', 'HN', 'XN', 'CN', 'SN', 'SNN'
	   jj = jj+NBEG(j)
	enddo
	if (jj .eq. 0)	then
	   do	j=8,10		! 'NEB', 'QNB', 'QNNB'
	      if (NBEG(j).gt.0)	then
		 write(*,'(2A)')
     >		   '  >>> Warning >>> Equation for "NE" is requested',
     >			' but not defined.'
		 write(*,'(18X,A)')'Boundary condition ignored.'
		 call APPTMP(TMPBUF(NBEG(j)),6)
		 NBEG(j) = 0			! Kill definition
	      endif
	   enddo
	   if (ITYPE .lt. 0)	return
	endif

	if (ITYPE .eq. 1) write(NCH,'(A)')'C **** Density equation '
	if (ITYPE .eq. 1) write(NCH,'(A)')
     >	'      call	markloc("NE equation"//char(0))'
	if (ITYPE .le. 0) write(NCH,'(A)')'C **** Density assignment'
	if (ITYPE .le. 0) write(NCH,'(A)')
     >	'      call	markloc("NE assignment"//char(0))'
	write(NCH,'(A)')'      YHRO = HRO'
	write(NCH,101)NEQ
 101 	format('      do 2',1I2,'2 J=1,NA1')
	do	j=2,7		! 'DN', 'HN', 'XN', 'CN', 'SN', 'SNN'
	   jj = NBEG(j)
	   if (jj .gt. 0)	then
	      call APPTMP(TMPBUF(jj),NCH)
C	      call APPTMP(TMPBUF(jj),6)
	   endif
	enddo
	if (NBEG(LDVN) .gt. 0)	then		! 'DVN'
	   call APPTMP(TMPBUF(NBEG(LDVN)),NCH)
C	   call APPTMP(TMPBUF(NBEG(LDVN)),6)
	endif
	if (NBEG(LDSN) .gt. 0)	then		! 'DSN'
	   call APPTMP(TMPBUF(NBEG(LDSN)),NCH)
C	   call APPTMP(TMPBUF(NBEG(LDSN)),6)
	endif

	if(NBEG(LSN).eq.0)  write(NCH,'(A)')'      SN(J)=0.'
	if(NBEG(LSNN).eq.0) write(NCH,'(A)')'      SNN(J)=0.'
			    write(NCH,'(A)')'      SNTOT(J)=SN(J)'
	if(NBEG(LDN)+NBEG(LDVN) .eq. 0)	then
	   write(NCH,'(A,$)')'      YWA(J)=0.'
	else
	   write(NCH,'(A,$)')'      YWA(J)='
	endif
	if(NBEG(LDN) .gt.0)	write(NCH,'(A,$)')'+DN(J)'
	if(NBEG(LDVN).gt.0)	write(NCH,'(A,$)')'+DVN(J)'
C	if(NBEG(LDSN).gt.0)	write(NCH,'(A,$)')'+DSN(J)'
				write(NCH,*)
C	if(NBEG(LDN).eq.0)  write(NCH,'(A)')'      YWA(J)=0.'
C	if(NBEG(LDN).gt.0)  write(NCH,'(A)')'      YWA(J)=DN(J)'
C	if(NBEG(LDVN).gt.0) write(NCH,'(A)')'      YWA(J)=YWA(J)+DVN(J)'
C	if(NBEG(LDSN).gt.0) write(NCH,'(A)')'      YWA(J)=YWA(J)+DSN(J)'
C
C	if(NBEG(LDVN).eq.0)	then
C	   if(NBEG(LDN).gt.0)	write(NCH,'(A)')'      YWA(J)=DN(J)'
C	   if(NBEG(LDN).eq.0)	write(NCH,'(A)')'      YWA(J)=0.'
C	else
C	   if(NBEG(LDN).gt.0)write(NCH,'(A)')'      YWA(J)=DN(J)+DVN(J)'
C	   if(NBEG(LDN).eq.0)	write(NCH,'(A)')'      YWA(J)=DVN(j)'
C	endif
	write(NCH,'(A,1I2,A)')'      if (j.gt.NA)	goto 2',NEQ,'2'
	write(NCH,'(A)')      '      if (j.eq.NA)	YHRO = HROA'
	if(NBEG(LDVN).eq.0)	then
	   write(NCH,'(A)')'      YWB(J)=0.'
	else
	   write(NCH,'(A)')'      YWB(J)=DVN(j)*log(NE(j)/NE(j+1))/YHRO'
	endif
	if(NBEG(LHN).gt.0)	write(NCH,102)
	if(NBEG(LXN).gt.0)	write(NCH,103)
	if(NBEG(LCN).gt.0)	write(NCH,104)
 102	format('     > +2.*HN(J)/(TE(J+1)+TE(J))*(TE(J+1)-TE(J))/YHRO')
 103	format('     > +2.*XN(J)/(TI(J+1)+TI(J))*(TI(J+1)-TI(J))/YHRO')
 104	format('     > -CN(J)')
	write(NCH,108)NEQ
 108	format(' 2',1I2,'2 continue')
	if(ITYPE.lt.0)					goto	8
 121	format('  >>> Warning: Initial condition for NE is not defined'
     >	      /'               NE=NEX(TSTART) will be used')
	if(NBEG(LNE).eq.0)	then
C	   if (ITYPE.eq.0)	goto	8
	   write(*,121)
	   write(7,121)
	endif
	if (ITYPE.eq.0 .and. NBEG(LNEB).gt.0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	   call APPTMP(TMPBUF(NBEG(LNEB)),NCH)
C	if(ITYPE.eq.0)		then
C	   write(NCH,'(A)')'      do j=1,NA1'
C	   call APPTMP(TMPBUF(NBEG(LNE)),NCH)
C	   write(NCH,'(A)')'      enddo'
C	   goto	8
	endif
	if (NBEG(LRON) .eq. 0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	else
	   call APPTMP(TMPBUF(NBEG(LRON)),NCH)
	endif
	write(NCH,'(A)')'      NA1N = ND1'
	write(NCH,'(A)')'      ND = ND1-1'
C	write(NCH,'(A)')'      YWC(1) = RHO(ND1)-RHO(ND)'
	if(ITYPE.gt.0)write(NCH,'(A)')'      YWC(1) = RHO(ND1)-RHO(ND)'
	if(NBEG(LNE).gt.0 .and. LINAPP(1).eq.0)		then
	   write(NCH,'(A)')'      if (ND1.lt.NA1) then'
C	   write(NCH,'(A)')'      do j=ND1,NA1'
	   if (ITYPE.eq.0 .and. NBEG(LNEB).gt.0)	then
	      write(NCH,'(A)')'      do j=ND1,NA'
	   else
	      write(NCH,'(A)')'      do j=ND1,NA1'
	   endif
	   call APPTMP(TMPBUF(NBEG(LNE)),NCH)
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	endif
	if(ITYPE.eq.0)		then
	   write(NCH,'(A)')'      do j=1,ND1'
	   call APPTMP(TMPBUF(NBEG(LNE)),NCH)
	   write(NCH,'(A)')'      enddo'
	   goto	6
	endif
	if (NBEG(LQNNB)+NBEG(LQNB)+NBEG(LNEB) .eq. 0)	then
	   write(*,*)
     >	      ' >>> Warning: Boundary condition for NE is not given.'
	   if (NBEG(LRON).eq.0)write(*,*)
     >		'              It is set to NEX(t0)'
	   if (NBEG(LRON).gt.0 .and. NBEG(LNE).eq.0) write(*,*)
     >	    '              NEX(t0) will be used at the shifted boundary'
	   if (NBEG(LRON).gt.0 .and. NBEG(LNE).gt.0) write(*,*)
     >	    '              NEX(t) will be used at the shifted boundary'
	   write(7,*)
     >	      ' >>> Warning: Boundary condition for NE is not set.'
	   if (NBEG(LRON).eq.0) write(7,'(A)')
     >		'              It is set to NEX(t0)'
	   if (NBEG(LRON).gt.0 .and. NBEG(LNE).eq.0) write(7,'(A)')
     >	    '              NEX(t0) will be used at the shifted boundary'
	   if (NBEG(LRON).gt.0 .and. NBEG(LNE).gt.0) write(7,'(A)')
     >	    '              NEX(t) will be used at the shifted boundary'
	   goto	6
	endif
C Two options:	-> 5  NEB = ...  is used for the boundary condition
C		-> 6  NE  = ...  is used for the boundary condition
C	if (NBEG(LRON) .gt. 0)	goto	5
C	if (NBEG(LRON) .gt. 0)	goto	6
	if (NBEG(LNEB) .gt. 0)			goto	5
	if (NBEG(LQNB)+NBEG(LQNNB) .eq. 0)	goto	6
	if (NBEG(LQNB) .gt. 0)	then
				call APPTMP(TMPBUF(NBEG(LQNB)),NCH)
				write(NCH,'(A)')'      YWC(2)=QNB'
	else
				write(NCH,'(A)')'      YWC(2)=0.'
	endif
	if (NBEG(LQNNB) .gt. 0)	then
				call APPTMP(TMPBUF(NBEG(LQNNB)),NCH)
				write(NCH,'(A)')'      YWC(3)=QNNB'
	else
				write(NCH,'(A)')'      YWC(3)=0.'
	endif
				write(NCH,'(A)')'      YWC(4)=-1.'
	goto	7
 5	if (NBEG(LNEB) .gt. 0)	call APPTMP(TMPBUF(NBEG(LNEB)),NCH)
 6	continue
	if(LINAPP(1) .ne. 0)		then
	   write(NCH,'(A)')'      if (ND1 .lt. NA)	then'
	   write(NCH,'(A)')'      YA = (NE(NA1)-NE(ND1))/(ROC-RON)'
	   write(NCH,'(A)')'      YB = NE(ND1)*ROC-NE(NA1)*RON'
	   write(NCH,'(A)')'      YB = YB/(ROC-RON)'
	   write(NCH,'(A)')'      do j=ND1+1,NA'
	   write(NCH,'(A)')'      NE(j) = YB+RHO(j)*YA'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	   write(7,*)'NE is linearly interpolated between RON and ROC'
	endif
	if (ITYPE .eq. 0)	goto	8	
				write(NCH,'(A)')'      NEO(ND1)=NE(ND1)'
				write(NCH,'(A)')'      YWC(4)=1.'
 7				write(NCH,124)
 124	format('      call RUNN(YWA,YWB,DSN,SNN,SN,',
     >	'NEO,ND,TAU,HRO,YWC,NE,VRO,VR,G11)')
 8	continue
	if(ITYPE.gt.0)	then
	   write(NCH,125)NEQ
	else
	   write(NCH,126)NEQ
	endif
 125	format(
     >	'      do 2',1I2,'4 J=1,NA'/
     >	'      QN(J)=-G11(J)*(YWA(J)*NE(J+1)-YWB(J)*NE(J))')
 126	format(
     >	'      do 2',1I2,'4 J=1,NA-1'/
     >	'      QN(J)=G11(J)*(-YWA(J)*(NE(J+1)-NE(J))/HRO-'/
     >	'     > 0.5*YWB(J)*(NE(J+1)+NE(J)))')
	write(NCH,'(A)')'      QN(J)=QN(J)+SLAT(J)*GNX(J)'
	write(NCH,'(A)')'      GN(J)=QN(J)/SLAT(J)'
	if(NBEG(LSNN).ne.0)
     >	write(NCH,'(A)')'      SNTOT(J)=SNTOT(J)+SNN(J)*NE(J)'
	write(NCH,127)NEQ
 127	format(' 2',1I2,'4 continue')
	if(ITYPE.le.0)	then
	   write(NCH,'(A)')
     >	'      QN(NA)=G11(NA)*(-YWA(NA)*(NE(NA1)-NE(NA))/HROA-'
	   write(NCH,'(A)')
     >	'     > 0.5*YWB(NA)*(NE(NA1)+NE(NA)))+SLAT(J)*GNX(J)'
	   write(NCH,'(A)')'      GN(NA)=QN(NA)/SLAT(NA)'
	   if(NBEG(LSNN).ne.0)
     >	   write(NCH,'(A)')'      SNTOT(NA)=SNTOT(NA)+SNN(NA)*NE(NA)'
	endif
	if( NBEG(LQNB).gt.0 .and. NBEG(LQNNB).eq.0
     >			    .and. NBEG(LRON) .eq.0)	then
	   write(NCH,'(A)')'      QN(ND1)=QNB'
	else
	   write(NCH,'(A)')'      QN(NA1)=QN(NA)'
	endif
	write(NCH,'(A)')'      GN(NA1)=QN(NA1)/SLAT(NA1)'
	write(NCH,'(A)')'      SNTOT(NA1)=SNTOT(NA)'
	end
C NIEQN ============= Ion density assignment ==========================|
	subroutine NIEQN(NEQ,ITYPE,NCH)
C NEQ	      not used
C ITYPE < 0   NI is assigned in "detvar.inc"
C ITYPE = 0   NI is assigned in "eqns.inc"
	include 'for/parameter.inc'
	include '.srv/tmpbuf.inc'
	integer	j,jj
	if (ITYPE .le. 0)write(NCH,'(A)')'C **** Ion density assignment'
	if (ITYPE .le. 0)write(NCH,'(A)')
     >	'      call	markloc("NI assignment"//char(0))'

	if(ITYPE.eq.0)		then
	   write(NCH,'(A)')'      do j=1,ND1'
	   call APPTMP(TMPBUF(NBEG(LNI)),NCH)
	   write(NCH,'(A)')'      enddo'
	endif
	end
C FJEQN ============ j-th dummy equation ==============================|
	subroutine FJEQN(NEQ,ITYPE,INUM,NCH)
C----------------------------------------------------------------------|
C Input:
C   NEQ 	ordinal number of the equation according to its position
C		in a model (is used to construct loop labels)
C   ITYPE < 0   no equation, no time dependence
C   ITYPE = 0   Fj:assigned
C   ITYPE = 1   Fj:equation
C   INUM	number of the equation: INUM=0 - F0; 1 - F1; ...; 9 - F9.
C   NCH 	logical unit to write include file
C----------------------------------------------------------------------|
	implicit none
	integer	NEQ,ITYPE,INUM,NCH
	character   J*1
C	character*3 FJX,DFJ,VFJ,GFJ,SFJ,FJB
C	character   J*1,FJ*2,SFFJ*4,QFJB*4,QFFJB*5,SFJTOT*5
	include 'for/parameter.inc'
	include '.srv/tmpbuf.inc'
	integer	j1,jj
	save	jj
	data	jj/0/
C----------------------------------------------------------------------|
	if (INUM.lt.0 .or. INUM.gt.9) then
	   write(*,'(A,I2,A)')'>>>  "F',INUM,'"  equation error'
	   stop
	endif
	write(J,'(1I1)')INUM
C	FJ = 'F'//J
C	DFJ = 'DF'//J
C	VFJ = 'VF'//J
C	GFJ = 'GF'//J
C	SFJ = 'SF'//J
C	SFFJ = 'SFF'//J
C	FJB = 'F'//J//'B'
C	QFJB = 'QF'//J//'B'
C	QFFJB = 'QFF'//J//'B'
C	SFJTOT = 'SF'//J//'TOT'
	jj = 0
C Search through DF?, VF?, SF?, SFF?:
C 81  ->  90:  DF0, DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8, DF9
C 91  -> 100:  VF0, VF1, VF2, VF3, VF4, VF5, VF6, VF7, VF8, VF9
C 101 -> 110:  SF0, SF1, SF2, SF3, SF4, SF5, SF6, SF7, SF8, SF9
C 111 -> 120:  SFF0,SFF1,SFF2,SFF3,SFF4,SFF5,SFF6,SFF7,SFF8,SFF9
	do	j1=81+INUM,111+INUM,10
	   jj = jj+NBEG(j1)
	enddo
	if (jj .eq. 0)	then
C Scan over F?B, QF?B, QFF?B:
C 131 -> 140:  F0B, F1B, F2B, F3B, F4B, F5B, F6B, F7B, F8B, F9B
C 141 -> 150:  QF0B,QF1B,QF2B,QF3B,QF4B,QF5B,QF6B,QF7B,QF8B,QF9B
C 151 -> 160:  QFF0B,QFF1B,QFF2B,QFF3B,QFF4B,QFF5B,QFF6B,QFF7B,QFF8B,QFF9B
	   do	j1=131+INUM,151+INUM,10
	      if (NBEG(j1).gt.0)	then
		 write(*,'(4A)')
     >		 '  >>> Warning >>> Equation for "F',J,'" is requested',
     >			' but not defined.'
		 write(*,'(18X,A)')'Boundary condition ignored.'
		 call APPTMP(TMPBUF(NBEG(j1)),6)
		 NBEG(j1) = 0			! Kill definition
	      endif
	   enddo
	   if (ITYPE .lt. 0)	return
	endif

	if (ITYPE .eq. 1) write(NCH,'(3A)')'C **** F',J,' equation '
	if (ITYPE .eq. 1) write(NCH,'(3A)')
     >	'      call	markloc("F',J,' equation"//char(0))'
	if (ITYPE .le. 0) write(NCH,'(3A)')'C **** F',J,' assignment'
	if (ITYPE .le. 0) write(NCH,'(3A)')
     >	'      call	markloc("F',J,' assignment"//char(0))'
	write(NCH,'(A)')'      YHRO = HRO'
	write(NCH,101)NEQ
 101 	format('      do 2',1I2,'2 J=1,NA1')
C
C 81  ->  90:  DF0, DF1, DF2, DF3, DF4, DF5, DF6, DF7, DF8, DF9
C 91  -> 100:  VF0, VF1, VF2, VF3, VF4, VF5, VF6, VF7, VF8, VF9
C 101 -> 110:  SF0, SF1, SF2, SF3, SF4, SF5, SF6, SF7, SF8, SF9
C 111 -> 120:  SFF0,SFF1,SFF2,SFF3,SFF4,SFF5,SFF6,SFF7,SFF8,SFF9
	do	j1=81+INUM,111+INUM,10
	   jj = NBEG(j1)
	   if (jj .gt. 0)	then
	      call APPTMP(TMPBUF(jj),NCH)
	   endif
	enddo
C 161 -> 170:  DVF0, DVF1, DVF2, DVF3, DVF4, DVF5, DVF6, DVF7, DVF8, DVF9
	jj = NBEG(161+INUM)
	if (jj .gt. 0)	then
	   call APPTMP(TMPBUF(jj),NCH)
C	   call APPTMP(TMPBUF(jj),6)
	endif
C 171 -> 180:  DSF0, DSF1, DSF2, DSF3, DSF4, DSF5, DSF6, DSF7, DSF8, DSF9
	jj = NBEG(171+INUM)
	if (jj .gt. 0)	then
	   call APPTMP(TMPBUF(jj),NCH)
C	   call APPTMP(TMPBUF(jj),6)
	endif
	if(NBEG(LSF(INUM)).eq.0) write(NCH,'(3A)')'      SF',J,'(J)=0.'
	if(NBEG(LSFF(INUM)).eq.0)write(NCH,'(3A)')'      SFF',J,'(J)=0.'
	   write(NCH,'(5A)')'      SF',J,'TOT(J)=SF',J,'(J)'
	if(NBEG(LDF(INUM))+NBEG(LDVF(INUM)) .eq. 0) then
	   write(NCH,'(A,$)')'      YWA(J)=0.'
	else
	   write(NCH,'(A,$)')'      YWA(J)='
	endif
	if(NBEG(LDF(INUM)) .gt.0) write(NCH,'(3A,$)')'+DF',J,'(J)'
	if(NBEG(LDVF(INUM)).gt.0) write(NCH,'(3A,$)')'+DVF',J,'(J)'
				  write(NCH,*)
C	if(NBEG(LDF(INUM)).eq.0)	then
C	   write(NCH,'(A)')'      YWA(J)=0.'
C	else
C	   write(NCH,'(3A)')'      YWA(J)=DF',J,'(J)'
C	endif
C	if(NBEG(LDVF(INUM)).gt.0)
C     >	   write(NCH,'(5A)')'      YWA(J)=DF',J,'(J)+DVF',J,'(J)'
C	if(NBEG(LDSF(INUM)).gt.0)
C     >	   write(NCH,'(5A)')'      YWA(J)=DF',J,'(J)+DSF',J,'(J)'
C
C	if(NBEG(LDVF(INUM)).eq.0)	then
C	   if(NBEG(LDF(INUM)).gt.0)
C     >			write(NCH,'(3A)')'      YWA(J)=DF',J,'(J)'
C	   if(NBEG(LDF(INUM)).eq.0)
C     >			write(NCH,'(A)')'      YWA(J)=0.'
C	else
C	   if(NBEG(LDF(INUM)).gt.0)
C     >		write(NCH,'(5A)')'      YWA(J)=DF',J,'(J)+DVF',J,'(J)'
C	   if(NBEG(LDF(INUM)).eq.0)
C     >		write(NCH,'(3A)')'      YWA(J)=DVF',J,'(j)'
C	endif
	write(NCH,'(A,1I2,A)')'      if (j.gt.NA)	goto 2',NEQ,'2'
	write(NCH,'(A)')      '      if (j.eq.NA)	YHRO = HROA'
	if(NBEG(LDVF(INUM)).eq.0)	then
	   write(NCH,'(A)')'      YWB(J)=0.'
	else
	   write(NCH,'(7A)')
     >	      '      YWB(J)=DVF',J,'(j)*log(F',J,'(j)/F',J,'(j+1))/YHRO'
	endif
	if(NBEG(LVF(INUM)).gt.0) write(NCH,'(3A)')'     > -VF',J,'(J)'
	write(NCH,108)NEQ
 108	format(' 2',1I2,'2 continue')
	if(ITYPE.lt.0)					goto	8
 121	format('  >>> Warning: Initial condition for F',A,' is not set'
     >	      /'               F',A,'=F',A,'X(TSTART) will be used')
	if(NBEG(LF(INUM)).eq.0)	then
C	   if (ITYPE.eq.0)	goto	8
	   write(*,121)J,J,J
	   write(7,121)J,J,J
	endif
	if (ITYPE.eq.0 .and. NBEG(LFB(INUM)).gt.0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	   call APPTMP(TMPBUF(NBEG(LFB(INUM))),NCH)
C	if(ITYPE.eq.0)		then
C	   write(NCH,'(A)')'      do j=1,NA1'
C	   call APPTMP(TMPBUF(NBEG(LF(INUM))),NCH)
C	   write(NCH,'(A)')'      enddo'
C	   goto	8
	endif
	if (NBEG(LRO(INUM)) .eq. 0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	else
	   call APPTMP(TMPBUF(NBEG(LRO(INUM))),NCH)
	endif
	write(NCH,'(3A)')'      NA1',J,' = ND1'
	write(NCH,'(A)')'      ND = ND1-1'
C	write(NCH,'(A)')'      YWC(1) = RHO(ND1)-RHO(ND)'
	if(ITYPE.gt.0)write(NCH,'(A)')'      YWC(1) = RHO(ND1)-RHO(ND)'
	if(NBEG(LF(INUM)).gt.0 .and. LINAPP(INUM+10).eq.0)	then
	   write(NCH,'(A)')'      if (ND1.lt.NA1) then'
C	   write(NCH,'(A)')'      do j=ND1,NA1'
	   if (ITYPE.eq.0 .and. NBEG(LFB(INUM)).gt.0)	then
	      write(NCH,'(A)')'      do j=ND1,NA'
	   else
	      write(NCH,'(A)')'      do j=ND1,NA1'
	   endif
	   call APPTMP(TMPBUF(NBEG(LF(INUM))),NCH)
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	endif
	if(ITYPE.eq.0)		then
	   write(NCH,'(A)')'      do j=1,ND1'
	   call APPTMP(TMPBUF(NBEG(LF(INUM))),NCH)
	   write(NCH,'(A)')'      enddo'
	   goto	6
	endif
	if (NBEG(LQFFB(INUM))+NBEG(LQFB(INUM))+NBEG(LFB(INUM)) .eq. 0)
     >	then
	   write(*,*)
     >	    ' >>> Warning: Boundary condition for F',J,' is not given.'
	   if (NBEG(LRO(INUM)).eq.0)write(*,'(3A)')
     >		'               It is set to F',J,'X(t0)'
	   if (NBEG(LRO(INUM)).gt.0 .and. NBEG(LF(INUM)).eq.0)
     >			write(*,'(3A)')
     > '               F',J,'X(t0) will be used at the shifted boundary'
	   if (NBEG(LRO(INUM)).gt.0 .and. NBEG(LF(INUM)).gt.0)
     >			write(*,'(3A)')
     >	'               F',J,'X(t) will be used at the shifted boundary'
	   write(7,*)
     >	      ' >>> Warning: Boundary condition for F',J,' is not set.'
	   if (NBEG(LRO(INUM)).eq.0) write(7,'(3A)')
     >		'               It is set to F',J,'X(t0)'
	   if (NBEG(LRO(INUM)).gt.0 .and. NBEG(LF(INUM)).eq.0)
     >			write(7,'(3A)')
     > '               F',J,'X(t0) will be used at the shifted boundary'
	   if (NBEG(LRO(INUM)).gt.0 .and. NBEG(LF(INUM)).gt.0)
     >			write(7,'(3A)')
     >	'               F',J,'X(t) will be used at the shifted boundary'
	   goto	6
	endif
	if (NBEG(LFB(INUM)) .gt. 0)			goto	5
	if (NBEG(LQFB(INUM))+NBEG(LQFFB(INUM)) .eq. 0)	goto	6
	if (NBEG(LQFB(INUM)).gt.0)	then
			call APPTMP(TMPBUF(NBEG(LQFB(INUM))),NCH)
			write(NCH,'(3A)')'      YWC(2)=QF',J,'B'
	else
			write(NCH,'(A)')'      YWC(2)=0.'
	endif
	if(NBEG(LQFFB(INUM)).gt.0)	then
			call APPTMP(TMPBUF(NBEG(LQFFB(INUM))),NCH)
			write(NCH,'(3A)')'      YWC(3)=QFF',J,'B'
	else
				write(NCH,'(A)')'      YWC(3)=0.'
	endif
				write(NCH,'(A)')'      YWC(4)=-1.'
	goto	7
 5	if(NBEG(LFB(INUM)) .gt. 0)
     >			call APPTMP(TMPBUF(NBEG(LFB(INUM))),NCH)
 6	continue
	if(LINAPP(INUM+10) .ne. 0)	then
	   write(NCH,'(A)')'      if (ND1 .lt. NA)	then'
	   write(NCH,'(7A)')
     >		'      YA = (F',J,'(NA1)-F',J,'(ND1))/(ROC-RO',J,')'
	   write(NCH,'(6A)')
     >		'      YB = F',J,'(ND1)*ROC-F',J,'(NA1)*RO',J
	   write(NCH,'(3A)')'      YB = YB/(ROC-RO',J,')'
	   write(NCH,'(A)')'      do j=ND1+1,NA'
	   write(NCH,'(3A)')'      F',J,'(j) = YB+RHO(j)*YA'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	   write(7,'(5A)')
     >	      'F',J,' is linearly interpolated between RO',J,' and ROC'
	endif
	if (ITYPE .eq. 0)	goto	8
		write(NCH,'(5A)')'      F',J,'O(ND1)=F',J,'(ND1)'
		write(NCH,'(A)')'      YWC(4)=1.'
 7		write(NCH,123)J,J,J,J,J
 123		format(
     >'      call RUNN(YWA,YWB,DSF',A,',SFF',A,',SF',A,',F',A,'O,',
     >'ND,TAU,HRO',',YWC,F',A,',VRO,VR,G11)')
 8	continue
	if(ITYPE.gt.0)	then
	   write(NCH,124)NEQ,J,J,J
	else
	   write(NCH,125)NEQ,J,J,J,J,J
	   if(NBEG(LGF(INUM)).gt.0)    write(NCH,'(7A)')
     >	'      QF',J,'(J)=QF',J,'(J)+SLAT(J)*GF',J,'(J)'
	endif
 124	format(
     >	'      do 2',1I2,'4 J=1,NA'/
     >	'      QF',A,'(J)=-G11(J)*(YWA(J)*F',A,'(J+1)-YWB(J)*F',A,'(J))')
 125	format(
     >	'      do 2',1I2,'4 J=1,NA-1'/
     >	'      QF',A,'(J)=-G11(J)*(YWA(J)*(F',A,'(J+1)-F',A,'(J))/HRO+'/
     >	'     > 0.5*YWB(J)*(F',A,'(J+1)+F',A,'(J)))')
	if(NBEG(LGF(INUM)).eq.0)
     >	write(NCH,'(5A)')'      GF',J,'(J)=QF',J,'(J)/SLAT(J)'
	if(NBEG(LSFF(INUM)).ne.0)	write(NCH,'(9A)')
     >	'      SF',J,'TOT(J)=SF',J,'TOT(J)+SFF',J,'(J)*F',J,'(J)'
	write(NCH,126)NEQ
 126	format(' 2',1I2,'4 continue')
	if(ITYPE.le.0)	then
	   write(NCH,'(7A)')'      QF',J,
     >	    '(NA)=-G11(NA)*(YWA(NA)*(F',J,'(NA1)-F',J,'(NA))/HROA+'
	   write(NCH,'(5A)')
     >	'     > 0.5*YWB(NA)*(F',J,'(NA1)+F',J,'(NA)))'
	if(NBEG(LGF(INUM)).eq.0)
     >	   write(NCH,'(5A)')'      GF',J,'(NA)=QF',J,'(NA)/SLAT(NA)'
	if(NBEG(LSFF(INUM)).ne.0)	write(NCH,'(9A)')
     >	'      SF',J,'TOT(NA)=SF',J,'TOT(NA)+SFF',J,'(NA)*F',J,'(NA)'
	endif
	if( NBEG(LQFB(INUM)).gt.0 .and. NBEG(LQFFB(INUM)).eq.0
     >			     .and. NBEG(LRO(INUM))  .eq.0 )	then
		   write(NCH,'(5A)')'      QF',J,'(ND1)=QF',J,'B'
		else
		   write(NCH,'(5A)')'      QF',J,'(NA1)=QF',J,'(NA)'
		endif
	if(NBEG(LGF(INUM)).eq.0)
     >	write(NCH,'(5A)')'      GF',J,'(NA1)=QF',J,'(NA1)/SLAT(NA1)'
	write(NCH,'(5A)')'      SF',J,'TOT(NA1)=SF',J,'TOT(NA)'
	end
C TEEQN ============= Electron temperature equation ===================|
C ITYPE < 0   fixed
C ITYPE = 0   assigned
C ITYPE = 1   equation is solved without particle flux
C ITYPE = 2   equation is solved with particle flux, factor 5/2
C ITYPE = 3   equation is solved with particle flux, factor defined in a model
	subroutine TEEQN(NEQ,ITYPE,NCH)
	include 'for/parameter.inc'
	include	'.srv/tmpbuf.inc'
	integer	j,jj
	jj = 0
	do	j=12,17		! 'DE', 'HE', 'XE', 'CE', 'PE', 'PET'
	   jj = jj+NBEG(j)
	enddo
	if (jj .eq. 0)	then
	   do	j=18,20		! 'TEB', 'QEB', 'QETB'
	      if (NBEG(j).gt.0)	then
		 write(*,'(2A)')
     >		   '  >>> Warning >>> Equation for "TE" is requested',
     >			' but not defined.'
		 write(*,'(18X,A)')'Boundary condition ignored.'
		 call APPTMP(TMPBUF(NBEG(j)),6)
		 NBEG(j) = 0			! Kill definition
	      endif
	   enddo
	   if (ITYPE .lt. 0)	return
	endif

	IMPEI = 0
	if (ITYPE .gt. 8)	then
	   IMPEI = 1
	   ITYPE = ITYPE-8
	endif
	if (ITYPE .ge. 1)	then
	   if (IMPEI .eq. 0)	then
	      write(NCH,'(A)')'C **** Electron temperature equation '
	      write(NCH,'(A)')
     >	      '      call	markloc("TE equation"//char(0))'
	   endif
	else
	   write(NCH,'(A)')'C **** Electron temperature assignment'
	   write(NCH,'(A)')
     >	'      call	markloc("TE assignment"//char(0))'
	endif
	write(NCH,102)NEQ
 102		format('      do 2',1I2,'1 J=1,NA1')
	do	j=12,17		! 'DE', 'HE', 'XE', 'CE', 'PE', 'PET'
	   jj = NBEG(j)
	   if (jj .gt. 0)	then
	      call APPTMP(TMPBUF(jj),NCH)
C	      write(*,*)"TE",jj,j
C	      call APPTMP(TMPBUF(jj),6)
	   endif
	enddo
	jj = NBEG(LDVE)
	if (jj .gt. 0)	then
	   call APPTMP(TMPBUF(jj),NCH)
C	   call APPTMP(TMPBUF(jj),6)
	endif
	jj = NBEG(LDSE)
	if (jj .gt. 0)	then
	   call APPTMP(TMPBUF(jj),NCH)
C	   call APPTMP(TMPBUF(jj),6)
	endif
	write(NCH,'(A,1I2,A)')'      if (j.gt.NA)	goto 2',NEQ,'1'
	write(NCH,'(A)')      '      YWA(J)=0.'
	if (NBEG(LHE).gt.0)	write(NCH,'(A)')'     > +HE(J)'
	if (NBEG(LHN).gt.0 .and. ITYPE.ge.2)
     >				write(NCH,'(A)')'     > +GN2E*HN(J)'
	write(NCH,'(A)')'      YWA(J)=YWA(J)*(NE(J+1)+NE(J))*0.5'

	LB=1
	if(ITYPE.ge.2.and.NBEG(LDN).gt.0) go to 3
	if(NBEG(LDE).gt.0)	write(NCH,'(A)')
     >	'      YWB(J)=DE(J)'
	if(NBEG(LDE).eq.0)	LB=0
					go to 4
 3	if(NBEG(LDE).gt.0)	write(NCH,'(A)')
     >	'      YWB(J)=(DE(J)+GN2E*DN(J))'
	if(NBEG(LDE).eq.0)	write(NCH,'(A)')
     >	'      YWB(J)=GN2E*DN(J)'
 4	LC=1
	if(ITYPE.ge.2.and.NBEG(LXN).gt.0) go to 5
	if(NBEG(LXE).gt.0)	write(NCH,'(A)')
     >	'      YWC(J)=XE(J)*(NE(J+1)+NE(J))/(TI(J+1)+TI(J))'
	if(NBEG(LXE).eq.0)	LC=0
					go to 6
 5	if(NBEG(LXE).gt.0)	write(NCH,'(A)')
     >	'      YWC(J)=(XE(J)+GN2E*XN(J))*(NE(J+1)+NE(J))/(TI(J+1)+TI(J))'
	if(NBEG(LXE).eq.0)	write(NCH,'(A)')
     >	'      YWC(J)=GN2E*XN(J)*(NE(J+1)+NE(J))/(TI(J+1)+TI(J))'
 6	continue
				write(NCH,'(A)')'      YWD(J)=0.'
	if(NBEG(LCE).gt.0)	write(NCH,'(A)')'     >	+CE(J)'
	if(NBEG(LCN).gt.0 .and. ITYPE.ge.2)
     >				write(NCH,'(A)')'     >	+GN2E*CN(J)'
	write(NCH,'(2A)')'      YWD(J)=YWD(J)*(NE(J+1)+NE(J))*0.5',
     >					'+GN2E*GNX(J)*SLAT(J)/G11(J)'
	write(NCH,116)NEQ
 116			format(' 2',1I2,'1 continue')
	write(NCH,'(A,1I2,A)')' 2',NEQ,'2 YHRO = HRO'
	write(NCH,'(A,1I2,A)')'      do 2',NEQ,'3 J=1,NA1'
	if(NBEG(LPET).eq.0)	write(NCH,'(A)')'      PET(J)=0.'
	if(NBEG(LPE).eq.0)	write(NCH,'(A)')'      PE(J)=0.'
				write(NCH,'(A)')'      PETOT(J)=PE(J)'
	write(NCH,'(A,1I2,A)')'      if (j.gt.NA)	goto 2',NEQ,'3'
	if(LB+LC.gt.0)	write(NCH,'(A)')'      if (j.eq.NA) YHRO = HROA'
	if(ITYPE.lt.0)			go to 10
				write(NCH,'(A)')'      YWB(J)=-YWD(J)'
	if(LB.gt.0)	write(NCH,103)
	if(LC.gt.0)	write(NCH,104)
 103	format('     > +YWB(J)*(NE(J+1)-NE(J))/YHRO')
 104	format('     > +YWC(J)*(TI(J+1)-TI(J))/YHRO')
	write(NCH,'(A)')'      YWC(J)=PE(j)'
	if(NBEG(LDVE)+NBEG(LDSE).gt.0) write(NCH,'(A)')
     >	'      if (KDV.gt.0) YWC(j)=YWC(j)-PDE(j)'

	if(NBEG(LDVE).gt.0)	write(NCH,'(A)')'      YWD(J)=DVE(j)'

C	if(NBEG(LDN).eq.0)  write(NCH,'(A)')'      YWA(J)=0.'
C	if(NBEG(LDN).gt.0)  write(NCH,'(A)')'      YWA(J)=DN(J)'
C	if(NBEG(LDVN).gt.0) write(NCH,'(A)')'      YWA(J)=YWA(J)+DVN(J)'
C	if(NBEG(LDSE).gt.0) write(NCH,'(A)')'      YWA(J)=YWA(J)+DSN(J)'
C
C	if(NBEG(LDN)+NBEG(LDVN)+NBEG(LDSN).eq.0)	then
C	   write(NCH,'(A,$)')'      YWA(J)=0.'
C	else
C	   write(NCH,'(A,$)')'      YWA(J)='
C	endif
C	if(NBEG(LDN) .gt.0)	write(NCH,'(A,$)')'+DN(J)'
C	if(NBEG(LDVN).gt.0)	write(NCH,'(A,$)')'+DVN(J)'
C	if(NBEG(LDSN).gt.0)	write(NCH,'(A,$)')'+DSN(J)'
C				write(NCH,*)


 10	continue
	write(NCH,133)NEQ
 133	format(' 2',1I2,'3 continue')
	if(ITYPE.lt.0)					goto	14
 121	format('  >>> Warning: Initial condition for TE is not defined'
     >	      /'               TE=TEX(TSTART) will be used')
	if(NBEG(LTE).eq.0)	then
	   write(*,121)
	   write(7,121)
	endif
	if (ITYPE.eq.0 .and. NBEG(LTEB).gt.0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	   call APPTMP(TMPBUF(NBEG(LTEB)),NCH)
	endif
	if (NBEG(LROE) .eq. 0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	else
	   call APPTMP(TMPBUF(NBEG(LROE)),NCH)
	endif
	write(NCH,'(A)')'      NA1E = ND1'
	write(NCH,'(A)')'      ND = ND1-1'
	if(ITYPE.gt.0)write(NCH,'(A)')'      QE(1) = RHO(ND1)-RHO(ND)'
	if(NBEG(LTE).gt.0 .and. LINAPP(2).eq.0)		then
	   write(NCH,'(A)')'      if (ND1.lt.NA1) then'
	   if (ITYPE.eq.0 .and. NBEG(LTEB).gt.0)	then
	      write(NCH,'(A)')'      do j=ND1,NA'
	   else
	      write(NCH,'(A)')'      do j=ND1,NA1'
	   endif
	   call APPTMP(TMPBUF(NBEG(LTE)),NCH)
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	endif
	if(ITYPE.eq.0)		then
	   write(NCH,'(A)')'      do j=1,ND1'
	   call APPTMP(TMPBUF(NBEG(LTE)),NCH)
	   write(NCH,'(A)')'      enddo'
	   goto	12
	endif
	if (NBEG(LQETB)+NBEG(LQEB)+NBEG(LTEB) .eq. 0)	then
	   write(*,*)
     >		' >>> Warning: Boundary condition for TE is not given.'
	   if (NBEG(LROE).eq.0)write(*,*)
     >		'              It is set to TEX(t0)'
	   if (NBEG(LROE).gt.0 .and. NBEG(LTE).eq.0) write(*,*)
     >	    '              TEX(t0) will be used at the shifted boundary'
	   if (NBEG(LROE).gt.0 .and. NBEG(LTE).gt.0) write(*,*)
     >	    '              TEX(t) will be used at the shifted boundary'
	   write(7,*)
     >		' >>> Warning: Boundary condition for TE is not set.'
	   if (NBEG(LROE).eq.0) write(7,'(A)')
     >		'              It is set to TEX(t0)'
	   if (NBEG(LROE).gt.0 .and. NBEG(LTE).eq.0) write(7,'(A)')
     >	    '              TEX(t0) will be used at the shifted boundary'
	   if (NBEG(LROE).gt.0 .and. NBEG(LTE).gt.0) write(7,'(A)')
     >	    '              TEX(t) will be used at the shifted boundary'
	   goto	12
	endif
	if (NBEG(LTEB) .gt. 0)			goto	11
	if (NBEG(LQEB)+NBEG(LQETB) .eq. 0)	goto	12
	if (NBEG(LQEB) .gt. 0)	then
				call APPTMP(TMPBUF(NBEG(LQEB)),NCH)
				write(NCH,'(A)')'      QE(2)=QEB'
	else
				write(NCH,'(A)')'      QE(2)=0.'
	endif
	if(NBEG(LQETB).gt.0)	then
				call APPTMP(TMPBUF(NBEG(LQETB)),NCH)
				write(NCH,'(A)')'      QE(3)=QETB'
	else
				write(NCH,'(A)')'      QE(3)=0.'
	endif
				write(NCH,'(A)')'      QE(4)=-1.'
	goto	13
 11	if(NBEG(LTEB) .gt. 0)	call APPTMP(TMPBUF(NBEG(LTEB)),NCH)
 12	continue
	if(LINAPP(2) .ne. 0)		then
	   write(NCH,'(A)')'      if (ND1 .lt. NA)	then'
	   write(NCH,'(A)')'      YA = (TE(NA1)-TE(ND1))/(ROC-ROE)'
	   write(NCH,'(A)')'      YB = TE(ND1)*ROC-TE(NA1)*ROE'
	   write(NCH,'(A)')'      YB = YB/(ROC-ROE)'
	   write(NCH,'(A)')'      do j=ND1+1,NA'
	   write(NCH,'(A)')'      TE(j) = YB+RHO(j)*YA'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	   write(7,*)'TE is linearly interpolated between ROE and ROC'
	endif
	if (ITYPE .eq. 0)	goto	14
				write(NCH,'(A)')'      TEO(ND1)=TE(ND1)'
				write(NCH,'(A)')'      QE(4)=1.'
 13	continue
	if(NBEG(LDVE).eq.0)	then
				write(NCH,'(A)')'      YWD(ND1)=0.'
	else
				write(NCH,'(A)')'      YWD(ND1)=1.'
	endif
	if(NBEG(LDSE).eq.0)	then
				write(NCH,'(A)')'      DSE(ND1)=0.'
	else
				write(NCH,'(A)')'      DSE(ND1)=1.'
	endif
	if (IMPEI .eq. 0)	then
	   write(NCH,135)
 135	   format(
     >	'      call RUNT(YWA,YWB,PET,PETOT,YWD,NEO,NE,TEO,'/
     >	'     > ND,TAU,HRO,QE,TE,VRO,VR,G11)')
	else
	   write(NCH,136)
 136	   format('      call RUNTT(YWA,YWB,PET,YWC,NEO,NE,TEO,'/
     >	   '     > ND,TAU,HRO,QE(1),YWD,DSE,VRO,VR,G11,WORK1)'/
     >	   '      do	J=ND1,NB1'/
     >	   '         PDE(j) = 0.'/
     >	   '      enddo')
	   write(NCH,'(A)')'      do	J=1,ND'
	   if (NBEG(LDVE) .eq. 0)	then
	      write(NCH,'(A)')'         PDE(j) = 0.'
	   else
	      write(NCH,'(A)')'         PDE(j) = YWD(j)'
	   endif
C	   if (NBEG(LDVE)+NBEG(LDSE) .eq. 0)	then
C	      write(NCH,'(A)')'         PDE(j) = 0.'
C	   elseif (NBEG(LDVE) .eq. 0)	then
C	      write(NCH,'(A)')'         PDE(j) = YWA(j)'
C	   elseif (NBEG(LDSE) .eq. 0)	then
C	      write(NCH,'(A)')'         PDE(j) = YWD(j)'
C	   else
C	      write(NCH,'(A)')'         PDE(j) = YWA(j)+YWD(j)'
C	   endif
	   write(NCH,'(A)')'      enddo'
	endif
	if (ITYPE .gt. 0)	return
 14	continue
	write(NCH,'(A)')'      do j=1,NA-1'
	write(NCH,'(A)')
     >		'      QE(J)=-G11(J)*(YWA(J)*(TE(J+1)-TE(J))/HRO+'
	write(NCH,'(A)')'     > 0.5*YWB(J)*(TE(J+1)+TE(J)))*0.0016'
	write(NCH,'(A)')'      enddo'
	write(NCH,'(A)')
     >		'      QE(NA)=-G11(NA)*(YWA(NA)*(TE(NA1)-TE(NA))/HROA+'
	write(NCH,'(A)')'     > 0.5*YWB(NA)*(TE(NA1)+TE(NA)))*.0016'
	if( NBEG(LQEB).gt.0 .and. NBEG(LQETB).eq.0
     >			    .and. NBEG(LROE) .eq.0 )	then
			write(NCH,'(A)')'      QE(ND1)=QEB'
		else
			write(NCH,'(A)')'      QE(NA1)=QE(NA)'
		endif
	end
C TIEQN ============== Ion temperature equation =======================|
C ITYPE < 0   fixed
C ITYPE = 0   assigned
C ITYPE = 1   equation is solved without particle flux
C ITYPE = 2   equation is solved with particle flux, factor 5/2
C ITYPE = 3   equation is solved with particle flux, factor defined in a model
	subroutine TIEQN(NEQ,ITYPE,NCH)
	include 'for/parameter.inc'
	include '.srv/tmpbuf.inc'
	integer	j,jj
	jj = 0
	do	j=22,27		! 'DI', 'HI', 'XI', 'CI', 'PI', 'PIT'
	   jj = jj+NBEG(j)
	enddo
	if (jj .eq. 0)	then
	   do	j=28,30		! 'TIB', 'QIB', 'QITB'
	      if (NBEG(j).gt.0)	then
		 write(*,'(2A)')
     >		   '  >>> Warning >>> Equation for "TI" is requested',
     >			' but not defined.'
		 write(*,'(18X,A)')'Boundary condition ignored.'
		 call APPTMP(TMPBUF(NBEG(j)),6)
		 NBEG(j) = 0			! Kill definition
	      endif
	   enddo
	   if (ITYPE .lt. 0)	return
	endif

	IMPEI = 0
	if (ITYPE .gt. 8)	then
	   IMPEI = 1
	   ITYPE = ITYPE-8
	endif
	if (ITYPE .ge. 1)	then
	   if (IMPEI .eq. 0)	then
	      write(NCH,'(A)')'C **** Ion temperature equation '
	      write(NCH,'(A)')
     >	      '      call	markloc("TI equation"//char(0))'
	   endif
	else
	   write(NCH,'(A)')'C **** Ion temperature assignment'
	   write(NCH,'(A)')
     >	'      call	markloc("TI assignment"//char(0))'
	endif
	write(NCH,102)NEQ
 102		format('      do 2',1I2,'1 J=1,NA1')
	do	j=22,27		! 'DI', 'HI', 'XI', 'CI', 'PI', 'PIT'
	   jj = NBEG(j)
	   if (jj .gt. 0)	then
	      call APPTMP(TMPBUF(jj),NCH)
C	      write(*,*)"TI",jj,j
C	      call APPTMP(TMPBUF(jj),6)
	   endif
	enddo
	jj = NBEG(LDVI)
	if (jj .gt. 0)	then
	   call APPTMP(TMPBUF(jj),NCH)
C	   call APPTMP(TMPBUF(jj),6)
	endif
	jj = NBEG(LDSI)
	if (jj .gt. 0)	then
	   call APPTMP(TMPBUF(jj),NCH)
C	   call APPTMP(TMPBUF(jj),6)
	endif
	write(NCH,'(A,1I2,A)')'      if (j.gt.NA)	goto 2',NEQ,'1'
	write(NCH,'(A)')      '      YWA(J)=0.'
	if (NBEG(LXI).gt.0)	write(NCH,'(A)')'     > +XI(J)'
	if (NBEG(LXN).gt.0 .and. ITYPE.ge.2)
     >				write(NCH,'(A)')'     > +GN2I*XN(J)'
	write(NCH,'(A)')'      YWA(J)=YWA(J)*(NI(J+1)+NI(J))*0.5'

	LB=1
	if(ITYPE.ge.2.and.NBEG(LDN).gt.0) go to 3
	if(NBEG(LDI).gt.0)	write(NCH,'(A)')
     >	'      YWB(J)=DI(J)*(NI(J+1)+NI(J))/(NE(J+1)+NE(J))'
	if(NBEG(LDI).eq.0)	LB=0
					go to 4
 3	if(NBEG(LDI).gt.0)	write(NCH,'(A)')
     >'      YWB(J)=(DI(J)+GN2I*DN(J))*(NI(J+1)+NI(J))/(NE(J+1)+NE(J))'
	if(NBEG(LDI).eq.0)	write(NCH,'(A)')
     >	'      YWB(J)=GN2I*DN(J)*(NI(J+1)+NI(J))/(NE(J+1)+NE(J))'
 4	LC=1
	if(ITYPE.ge.2.and.NBEG(LHN).gt.0) go to 5
	if(NBEG(LHI).gt.0)	write(NCH,'(A)')
     >	'      YWC(J)=HI(J)*(NI(J+1)+NI(J))/(TE(J+1)+TE(J))'
	if(NBEG(LHI).eq.0)	LC=0
					go to 6
 5	if(NBEG(LHI).gt.0)	write(NCH,'(A)')
     >'      YWC(J)=(HI(J)+GN2I*HN(J))*(NI(J+1)+NI(J))/(TE(J+1)+TE(J))'
	if(NBEG(LHI).eq.0)	write(NCH,'(A)')
     >	'      YWC(J)=GN2I*HN(J)*(NI(J+1)+NI(J))/(TE(J+1)+TE(J))'
 6	continue

	write(NCH,'(A)')
     >	'      YWD(J)=2.*GN2I*GNX(J)*SLAT(J)/G11(J)/(NE(J+1)+NE(J))'
	if(NBEG(LCI).gt.0)	write(NCH,'(A)')'     >	+CI(J)'
	if(NBEG(LCN).gt.0 .and. ITYPE.ge.2)
     >				write(NCH,'(A)')'     >	+GN2I*CN(J)'
	write(NCH,'(A)')'      YWD(J)=0.5*YWD(J)*(NI(J+1)+NI(J))'
	write(NCH,116)NEQ
 116			format(' 2',1I2,'1 continue')
	write(NCH,'(A,1I2,A)')' 2',NEQ,'2 YHRO = HRO'
	write(NCH,'(A,1I2,A)')'      do 2',NEQ,'3 J=1,NA1'
	if(NBEG(LPIT).eq.0)	write(NCH,'(A)')'      PIT(J)=0.'
	if(NBEG(LPI).eq.0)	write(NCH,'(A)')'      PI(J)=0.'
				write(NCH,'(A)')'      PITOT(J)=PI(J)'
	write(NCH,'(A,1I2,A)')'      if (j.gt.NA)	goto 2',NEQ,'3'
	if(LB+LC.gt.0)	write(NCH,'(A)')'      if (j.eq.NA) YHRO = HROA'
	if(ITYPE.lt.0)			go to 10
			write(NCH,'(A)')'      YWB(J)=-YWD(J)'
	if(LB.gt.0)	write(NCH,103)
	if(LC.gt.0)	write(NCH,104)
 103	format('     > +YWB(J)*(NE(J+1)-NE(J))/YHRO')
 104	format('     > +YWC(J)*(TE(J+1)-TE(J))/YHRO')
	write(NCH,'(A)')'      YWC(J)=PI(j)'
	if(NBEG(LDVI)+NBEG(LDSI).gt.0) write(NCH,'(A)')
     >	'      if (KDV.gt.0) YWC(j)=YWC(j)-PDI(j)'
	if(NBEG(LDVI).gt.0)	write(NCH,'(A)')'      YWD(J)=DVI(j)'
 10	continue
	write(NCH,133)NEQ
 133	format(' 2',1I2,'3 continue')
	if(ITYPE.lt.0)					goto	14
 121	format('  >>> Warning: Initial condition for TI is not defined'
     >	      /'               TI=TIX(TSTART) will be used')
	if(NBEG(LTI).eq.0)	then
	   write(*,121)
	   write(7,121)
	endif
	if (ITYPE.eq.0 .and. NBEG(LTIB).gt.0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	   call APPTMP(TMPBUF(NBEG(LTIB)),NCH)
	endif
	if (NBEG(LROI) .eq. 0)	then
	   write(NCH,'(A)')'      ND1 = NA1'
	else
	   call APPTMP(TMPBUF(NBEG(LROI)),NCH)
	endif
	write(NCH,'(A)')'      NA1I = ND1'
	write(NCH,'(A)')'      ND = ND1-1'
	if(ITYPE.gt.0)write(NCH,'(A)')'      QI(1) = RHO(ND1)-RHO(ND)'
	if(NBEG(LTI).gt.0 .and. LINAPP(3).eq.0)		then
	   write(NCH,'(A)')'      if (ND1.lt.NA1) then'
	   if (ITYPE.eq.0 .and. NBEG(LTIB).gt.0)	then
	      write(NCH,'(A)')'      do j=ND1,NA'
	   else
	      write(NCH,'(A)')'      do j=ND1,NA1'
	   endif
	   call APPTMP(TMPBUF(NBEG(LTI)),NCH)
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	endif
	if(ITYPE.eq.0)		then
	   write(NCH,'(A)')'      do j=1,ND1'
	   call APPTMP(TMPBUF(NBEG(LTI)),NCH)
	   write(NCH,'(A)')'      enddo'
	   goto	12
	endif
	if (NBEG(LQITB)+NBEG(LQIB)+NBEG(LTIB) .eq. 0)	then
	   write(*,*)
     >		' >>> Warning: Boundary condition for TI is not given.'
	   if (NBEG(LROI).eq.0)write(*,*)
     >		'              It is set to TIX(t0)'
	   if (NBEG(LROI).gt.0 .and. NBEG(LTI).eq.0) write(*,*)
     >	    '              TIX(t0) will be used at the shifted boundary'
	   if (NBEG(LROI).gt.0 .and. NBEG(LTI).gt.0) write(*,*)
     >	    '              TIX(t) will be used at the shifted boundary'
	   write(7,*)
     >		' >>> Warning: Boundary condition for TI is not set.'
	   if (NBEG(LROI).eq.0) write(7,'(A)')
     >		'              It is set to TIX(t0)'
	   if (NBEG(LROI).gt.0 .and. NBEG(LTI).eq.0) write(7,'(A)')
     >	    '              TIX(t0) will be used at the shifted boundary'
	   if (NBEG(LROI).gt.0 .and. NBEG(LTI).gt.0) write(7,'(A)')
     >	    '              TIX(t) will be used at the shifted boundary'
	   goto	12
	endif
C Two options:	-> 11  TIB = ...  is used for the boundary condition
C		-> 12  TI  = ...  is used for the boundary condition
C	if (NBEG(LROI) .gt. 0)	goto	11
C	if (NBEG(LROI) .gt. 0)	goto	12
	if (NBEG(LTIB) .gt. 0)			goto	11
	if (NBEG(LQIB)+NBEG(LQITB).eq.0)	goto	12
	if (NBEG(LQIB) .gt. 0)	then
				call APPTMP(TMPBUF(NBEG(LQIB)),NCH)
				write(NCH,'(A)')'      QI(2)=QIB'
	else
				write(NCH,'(A)')'      QI(2)=0.'
	endif
	if(NBEG(LQITB).gt.0)	then
				call APPTMP(TMPBUF(NBEG(LQITB)),NCH)
				write(NCH,'(A)')'      QI(3)=QITB'
	else
				write(NCH,'(A)')'      QI(3)=0.'
	endif
				write(NCH,'(A)')'      QI(4)=-1.'
	goto	13
 11	if(NBEG(LTIB) .gt. 0)	call APPTMP(TMPBUF(NBEG(LTIB)),NCH)
 12	continue
	if(LINAPP(3) .ne. 0)		then
	   write(NCH,'(A)')'      if (ND1 .lt. NA)	then'
	   write(NCH,'(A)')'      YA = (TI(NA1)-TI(ND1))/(ROC-ROI)'
	   write(NCH,'(A)')'      YB = TI(ND1)*ROC-TI(NA1)*ROI'
	   write(NCH,'(A)')'      YB = YB/(ROC-ROI)'
	   write(NCH,'(A)')'      do j=ND1+1,NA'
	   write(NCH,'(A)')'      TI(j) = YB+RHO(j)*YA'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	   write(7,*)'TI is linearly interpolated between ROI and ROC'
	endif
	if (ITYPE .eq. 0)	goto	14
				write(NCH,'(A)')'      TIO(ND1)=TI(ND1)'
				write(NCH,'(A)')'      QI(4)=1.'
 13	continue
	if(NBEG(LDVI).eq.0)	then
				write(NCH,'(A)')'      YWD(ND1)=0.'
	else
				write(NCH,'(A)')'      YWD(ND1)=1.'
	endif
	if(NBEG(LDSI).eq.0)	then
				write(NCH,'(A)')'      DSI(ND1)=0.'
	else
				write(NCH,'(A)')'      DSI(ND1)=1.'
	endif
	if (IMPEI .eq. 0)	then
	   write(NCH,135)
 135	   format(
     >	'      call RUNT'/
     >	'     > (YWA,YWB,PIT,PITOT,YWD,NIO,NI,TIO,ND,TAU,HRO,QI,TI,',
     >						'VRO,VR,G11)')
	else
	   write(NCH,136)
 136	   format('      call RUNTT(YWA,YWB,PIT,YWC,NIO,NI,TIO,'/
     >	   '     > ND,TAU,HRO,QI(1),YWD,DSI,VRO,VR,G11,WORK1)'/
     >	   '      do	J=ND1,NB1'/
     >	   '         PDI(j) = 0.'/
     >	   '      enddo')
	   write(NCH,'(A)')'      do	J=1,ND'
	   if (NBEG(LDVI) .eq. 0)	then
	      write(NCH,'(A)')'         PDI(j) = 0.'
	   else
	      write(NCH,'(A)')'         PDI(j) = YWD(j)'
	   endif
C	   if (NBEG(LDVI)+NBEG(LDSI) .eq. 0)	then
C	      write(NCH,'(A)')'         PDI(j) = 0.'
C	   elseif (NBEG(LDVI) .eq. 0)	then
C	      write(NCH,'(A)')'         PDI(j) = YWA(j)'
C	   elseif (NBEG(LDSI) .eq. 0)	then
C	      write(NCH,'(A)')'         PDI(j) = YWD(j)'
C	   else
C	      write(NCH,'(A)')'         PDI(j) = YWA(j)+YWD(j)'
C	   endif
	   write(NCH,'(A)')'      enddo'
	endif
	if (ITYPE .gt. 0)	return
 14	continue
	write(NCH,'(A)')'      do j=1,NA-1'
	write(NCH,'(A)')
     >		'      QI(J)=-G11(J)*(YWA(J)*(TI(J+1)-TI(J))/HRO+'
	write(NCH,'(A)')'     > 0.5*YWB(J)*(TI(J+1)+TI(J)))*0.0016'
	write(NCH,'(A)')'      enddo'
	write(NCH,'(A)')
     >		'      QI(NA)=-G11(NA)*(YWA(NA)*(TI(NA1)-TI(NA))/HROA+'
	write(NCH,'(A)')'     > 0.5*YWB(NA)*(TI(NA1)+TI(NA)))*.0016'
	if( NBEG(LQIB).gt.0 .and. NBEG(LQITB).eq.0
     >			    .and. NBEG(LROI) .eq.0 )	then
			write(NCH,'(A)')'      QI(ND1)=QIB'
		else
			write(NCH,'(A)')'      QI(NA1)=QI(NA)'
		endif
	end
C CUEQN ============= Current equation ================================|
	subroutine CUEQN(NEQ,ITYPE,NCH)
C ITYPE = 1  Current equation is solved
C ITYPE # 1  Not used
C-----------------------------------------------------------------------
	integer	j,jj,NEQ,ITYPE,NCH
	character PLATFORM*79
	include 'for/parameter.inc'
	include '.srv/tmpbuf.inc'
	if(ITYPE.ne.1)	then
		write(*,*)' >>> Warning <<< Unknown CU:?? type',ITYPE
		return
	endif
	write(NCH,'(A)')'C **** Current equation'
	write(NCH,'(A)')
     >	'      call	markloc("Current equation"//char(0))'
	write(NCH,'(A)')'      YC=.4*GP*RTOR/BTOR'
	write(NCH,'(A)')'      YD=-0.8*GP*GP*RTOR'
	write(NCH,'(A)')'      YA=2./(HRO**2*YD)'
C stellarator:---------------------------------------------------------|
	if (NBEG(LMV).gt.0)	then
	write(NCH,'(A)')'      YFV=GP2*HRO*HRO*BTOR'
	write(NCH,'(A)')'      FV(1)=0.'
	write(NCH,'(A)')'      YM=0.'
	write(NCH,"('      do 2',1I2.2,'1 J=1,NA1')")NEQ
C MU vacuum distribution is prescribed, CV is calculated
		call APPTMP(TMPBUF(NBEG(LMV)),NCH)
	write(NCH,'(A)')'      YM1=MV(j)*j'
	write(NCH,'(A)')'      YM2=YM1*G22(j)'
	write(NCH,*)'      CV(j)=(YM2-YM)*G33(j)*IPOL(j)**3/YC/RHO(j)'
C	write(NCH,*)'      CV(j+1)=(YM2-YM)*G33(j)*IPOL(j)**3/YC/RHO(j)'
	write(NCH,'(A)')'      if (j.lt.NA)	then'
	write(NCH,'(A)')'         YDF=YFV*YM1'
	write(NCH,'(A)')'      elseif(j.eq.NA)	then'
	write(NCH,'(A)')'         YDF=GP2*HROA*BTOR*YM1*HRO'
	write(NCH,'(A)')'      endif'
	write(NCH,'(A)')'      if (j.le.NA)	then'
	write(NCH,'(A)')'         FV(j+1)=FV(j)+YDF'
	write(NCH,'(A)')'      else'
	write(NCH,"('         CV(NA1)=CV(NA)')")
	write(NCH,'(A)')'      endif'
	write(NCH,'(A)')'      YM=YM2'
	write(NCH,"(' 2',1I2.2,'1 continue')")NEQ
C	elseif (NBEG(LCV).gt.0)	then
C CV vacuum distribution is prescribed, MV is calculated
C		write(NCH,"('      do 2',1I2.2,'1 J=1,NA')")NEQ
C		call APPTMP(TMPBUF(NBEG(LCV)),NCH)
C		write(NCH,"(' 2',1I2.2,'1 continue')")NEQ
	endif
C stellarator:^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^|
	write(NCH,"('      do 2',1I2,'4 J=1,NA1')")NEQ
	do	j=32,34		! 'DC', 'HC', 'XC'
	   jj = NBEG(j)
	   if (jj .gt. 0)	then
	      call APPTMP(TMPBUF(jj),NCH)
C	      call APPTMP(TMPBUF(jj),6)
	      jj = 0
	   endif
	enddo
	do	j=36,37		! 'CD', 'CC'
	   jj = NBEG(j)
	   if (jj .gt. 0)	then
	      call APPTMP(TMPBUF(jj),NCH)
C	      call APPTMP(TMPBUF(jj),6)
	      jj = 0
	   endif
	enddo
	if(NBEG(LHC)+NBEG(LXC)+NBEG(LDC).eq.0)	then
C		if (ICUBS .eq. 0) write(NCH,'(A)')'      CUBS(J)=0.'
	   if (NBEG(LCUBS) .eq. 0) write(NCH,'(A)')'      CUBS(J)=0.'
	   if (NBEG(LCUBS) .gt. 0) call APPTMP(TMPBUF(NBEG(LCUBS)),NCH)
	else
C For assigned NE, TE, TI the following line should read:
C		write(NCH,*)'      if (j .eq. NA) YA = 2./(YD*HROA*HRO)'
C     or
C		write(NCH,*)'      if (j .eq. NA) YA = 2./(YD*HRO*HRO)'
C The next is only good when NE, TE and TI are calculated:
	   write(NCH,'(A,1I2,A)')'      if (j .eq. NA1) goto 2',NEQ,'2'
	   write(NCH,'(A)')'      if (j .eq. NA) YA = 2./(YD*HROA**2)'
	   write(NCH,'(A)')'      CUBS(J)=YA*(FP(J+1)-FP(J))*(0.'
 105	   format('     > +HC(J)*(TE(J+1)-TE(J))/(TE(J+1)+TE(J))')
 106	   format('     > +XC(J)*(TI(J+1)-TI(J))/(TI(J+1)+TI(J))')
 107	   format('     > +DC(J)*(NE(J+1)-NE(J))/(NE(J+1)+NE(J))')
 108	   format('     > )')
	   if(NBEG(LHC).gt.0)	write(NCH,105)
	   if(NBEG(LXC).gt.0)	write(NCH,106)
	   if(NBEG(LDC).gt.0)	write(NCH,107)
	   write(NCH,108)
	   write(NCH,'(A,1I2,A)')'      goto 2',NEQ,'3'
	   write(NCH,'(A,1I2,A)')' 2',NEQ,'2 continue'
	   write(NCH,'(A)')
     &	'      CUBS(NA1)=CUBS(NA)+(CUBS(NA)-CUBS(NA-1))*HROA/HRO'
	endif
	write(NCH,'(A,1I2,A)')' 2',NEQ,'3 continue'
	if (NBEG(LCD).eq.0)	then
	   write(NCH,'(A)')'      YYD=CUBS(J)'
	else
	   call APPTMP(TMPBUF(NBEG(LCD)),NCH)
	   write(NCH,'(A)')'      YYD=CUBS(J)+CD(J)'
	endif
	write(NCH,'(A)')'      YWD(J)=YYD*YD/(IPOL(J)**3*G33(J))'
	write(NCH,"(' 2',1I2,'4 continue')")NEQ
C keeping this loop separately is needed for appropriate use of CCMHD
	write(NCH,'(A)')'      YC=.4*GP'
	write(NCH,"('      do 2',1I2,'5 J=1,NA1')")NEQ
	write(NCH,'(A)')'      YWB(J)=CC(J)*YC/IPOL(J)**2'
	write(NCH,"(' 2',1I2,'5 continue')")NEQ
	if(NBEG(LUEXT)+NBEG(LLEXT).eq.0.or.NBEG(LIPL).gt.0)	then
C IPL is set by exp data or in a model
 102 	format('C Prescribed plasma current:')
	write(NCH,102)
	write(NCH,'(A)')
     &		'      FP(NA)=.4*GP*HROA*RTOR/(G22(NA)*IPOL(NA1))'
	write(NCH,'(A)')
     &	'      FP(NA1)=(HROA-.5*HRO)*ROC*CC(NA1)/RTOR/IPOL(NA1)/TAU'
	write(NCH,'(A)')'      FP(NA1)=0.'	! Suppress extrapolation
	write(NCH,'(A)')'      FP(NA-1)=1.+FP(NA)*FP(NA1)'
C FP(NA1) - HROA*(d\Psi)/(d\rho)|_{at_the_edge}
	if (NBEG(LMV).eq.0)	then
	  write(NCH,'(A)')'      FP(NA1)=FP(NA)*(IPL+FPO(NA1)*FP(NA1))'
	else
	  write(NCH,'(A)')
     &		'      FP(NA1)=FP(NA)*(IPL+FPO(NA1)*FP(NA1))+YDF'
	  write(NCH,'(A/A)')'      CV(NA1) = CV(NA)+YWD(NA1)+',
     &		'     +((CV(NA)+YWD(NA1))-(CV(NA-1)+YWD(NA)))/HRO*HROA'
	endif
	write(NCH,'(A)')'      FP(NA)=-1.'
C----------------------------------------------------------------------|
C----------------------------------------------------------------------|
C otherwise: 			(NBEG(LUEXT).gt.0 .or. NBEG(LLEXT).gt.0)
		elseif(NBEG(LLEXT).eq.0)	then
C UEXT is set, LEXT isn't
 103 	format(		'C Prescribed loop voltage:')
	write(NCH,103)
C	write(NCH,'(A)')'      YWD(NA1)=1.'
C	write(NCH,'(A)')'      YWB(NA1)=0.'
	write(NCH,'(A)')'      FP(NA-1)=1.'
	write(NCH,'(A)')'      FP(NA)=0.'
	write(NCH,'(A)')'      FP(NA1)=TAU*UEXT+FPO(NA1)'
		else
C LEXT is set 
 104 	format(		'C Circuit equation:')
	write(NCH,104)
C	write(NCH,'(A)')'      YL=5.*LEXT*G22(NA)*IPOL(NA1)/(GP*RTOR)'
C	write(NCH,'(A)')'      YWB(NA1)=HROA-YL'
C	write(NCH,'(A)')'      YWD(NA1)=HROA+YL'
C	write(NCH,'(A)')'      FP(NA1)=YWB(NA1)*FP(NA)+YWD(NA1)*FP(NA1)'
C----------------------------------------------------------------------|
	write(NCH,'(A)')'      YL=5.*LEXT*G22(NA)*IPOL(NA1)/(GP*RTOR)'
	write(NCH,'(A)')
     &		'      FP(NA1)=(HROA-YL)*FP(NA)+(HROA+YL)*FP(NA1)'
	if(NBEG(LUEXT).gt.0)
     &	write(NCH,'(A)')'     > +2.*HROA*TAU*UEXT'
C	if (NBEG(LMV).gt.0)	write(NCH,100)
C 100	format('     > +2.*HROA*(FV(NA1)-FVO(NA1))')
	write(NCH,'(A)')'      FP(NA)=HROA-YL'
	write(NCH,'(A)')'      FP(NA-1)=HROA+YL'
		endif
	write(NCH,'(A)')
     >	'      call RUNF(G22,YWB,YWC,YWD,FPO,NA,TAU,HRO,HROA,FP,FV)'
Ca	write(NCH,'(A)')'      YCU=1.25/(GP*GP*RTOR)'
	write(NCH,"('      do 2',1I2,'6 J=1,NA1')")NEQ
Ca	write(NCH,'(A)')'      CU(J)=YCU*G33(J)*IPOL(J)**3*YWC(j)'
	if (NBEG(LMV).gt.0) write(NCH,'(A)')'      CU(J)=CU(J)-CV(j)'
Ca	write(NCH,'(A)')'      MU(J)=YWB(j)/(GP2*BTOR)'
	write(NCH,'(A)')'      UPL(J)=YWD(J)'
	write(NCH,'(A)')'      ULON(J)=IPOL(J)*G33(J)*UPL(J)'
	write(NCH,"(' 2',1I2,'6 continue')")NEQ
	write(NCH,'(A)')'      UPL(NA1)=ARRNA1(UPL(NA-3:NA),HROA/HRO)'
	write(NCH,'(A)')'      ULON(NA1)=IPOL(NA1)*G33(NA1)*UPL(NA1)'
C	write(NCH,'(A)')'      ULON(NA1)=ULON(NA)'
	write(NCH,'(A)')'      call CUOFP'
	if (NBEG(LMV).gt.0)	then
	write(NCH,"('      do 2',1I2,'7 J=1,NA1')")NEQ
	write(NCH,'(A)')'      CU(J)=CU(J)-CV(j)'
	write(NCH,"(' 2',1I2,'7 continue')")NEQ
	endif
Ca	write(NCH,'(A)')
Ca     &	'      CU(NA)=CUBS(NA)+CD(NA)+CC(NA)*ULON(NA)/(RTOR*GP2)'
Ca	write(NCH,'(A)')
Ca     &	'      CU(NA1)=CUBS(NA1)+CD(NA1)+CC(NA1)*ULON(NA1)/(RTOR*GP2)'
	if(NBEG(LUEXT).eq.0 .and. NBEG(LLEXT).eq.0)	return
	write(NCH,'(A)')
     &		'      YIPL=(FP(NA1)-FP(NA)-(FV(NA1)-FV(NA)))/HROA'
	write(NCH,'(A)')'      IPL=5.*IPOL(NA1)*G22(NA)*YIPL/GP2/RTOR'
C	write(NCH,'(A)')'      write(*,*)"IPL =",IPL,"   Upl =",UPL(NA1)'
	end
C CUASN ============= Current iterations ==============================|
	subroutine CUASN(NEQ,ITYPE,NCH)
C-----------------------------------------------------------------------
C NEQ = 0 - write file "inipsi.tmp"
C 	ITYPE = -1 CU distribution is prescribed
C 	ITYPE = 0  MU distribution is prescribed
C 	ITYPE = 1  current is adjusted to have U = Const
C NEQ > 0 - write file "equftn.tmp"
C-----------------------------------------------------------------------
	integer		NEQ,ITYPE,NCH
	character	string*80
	include 'for/parameter.inc'
	include '.srv/tmpbuf.inc'
	write(NCH,'(A)')'C **** Current profile adjustment'
	write(NCH,'(A)')
     >	'      call	markloc("CU adjustment"//char(0))'
	write(NCH,'(A)')'      YB=.4*GP'
	write(NCH,'(A)')'      YC=YB*RTOR/BTOR'
	write(NCH,'(A)')'      YD=-0.8*GP*GP*RTOR'
	write(NCH,'(A)')'      YA=2./(HRO**2*YD)'
C stellarator:v--------------------------------------------------------v
	if (NBEG(LMV).gt.0)	then
	write(NCH,'(A)')'      YFV=GP2*HRO*HRO*BTOR'
	write(NCH,'(A)')'      FV(1)=0.'
	write(NCH,'(A)')'      YM=0.'
	write(NCH,"('      do 2',1I2.2,'4 J=1,NA1')")NEQ
C MU vacuum distribution is prescribed, CV is calculated
		call APPTMP(TMPBUF(NBEG(LMV)),NCH)
	write(NCH,'(A)')'      YM1=MV(j)*j'
	write(NCH,'(A)')'      YM2=YM1*G22(j)'
	write(NCH,'(A)')
     &		'      CV(j)=(YM2-YM)*G33(j)*IPOL(j)**3/YC/RHO(j)'
C	write(NCH,'(A)')
C     &		'      CV(j+1)=(YM2-YM)*G33(j)*IPOL(j)**3/YC/RHO(j)'
	write(NCH,'(A)')'      if (j.lt.NA)	then'
	write(NCH,'(A)')'         YDF=YFV*YM1'
	write(NCH,'(A)')'      elseif(j.eq.NA)	then'
	write(NCH,'(A)')'         YDF=GP2*HROA*BTOR*YM1*HRO'
	write(NCH,'(A)')'      endif'
	write(NCH,'(A)')'      if (j.le.NA)	then'
	write(NCH,'(A)')'         FV(j+1)=FV(j)+YDF'
	write(NCH,'(A)')'      else'
	write(NCH,"('         CV(NA1)=CV(NA)')")
	write(NCH,'(A)')'      endif'
	write(NCH,'(A)')'      YM=YM2'
	write(NCH,"(' 2',1I2.2,'4 continue')")NEQ
	else
	write(NCH,'(A)')'      do  J=1,NA1'
	write(NCH,'(A)')'         FV(j)=0.'
	write(NCH,'(A)')'      enddo'
	endif
C stellarator:^--------------------------------------------------------^
	if (ITYPE.eq.1)		then
	   write(NCH,'(A)')'      if (JIT .eq. 1)	then'
	   write(NCH,'(A)')'      do  j=1,NA1'
	   write(NCH,'(A)')
     &		'      FPO(j)=FV(j)+0.2*GP*RTOR*IPL*(RHO(j)/ROC)**2'
	   write(NCH,'(A)')'      enddo'
	   write(NCH,'(A)')'      endif'
	   write(NCH,'(A)')'      JCALL = 10'
	   write(NCH,"(' 2',1I2.2,'0 continue')")NEQ
	   write(NCH,'(A)')'      JCALL = JCALL-1'
	   write(NCH,'(A)')'      YF=1.E3'
	endif
	write(NCH,"('      do 2',1I2.2,'1 J=1,NA1')")NEQ
	do	j=32,34		! 'DC', 'HC', 'XC'
	   jj = NBEG(j)
	   if (jj .gt. 0)	then
	      call APPTMP(TMPBUF(jj),NCH)
C	      call APPTMP(TMPBUF(jj),6)
	      jj = 0
	   endif
	enddo
	if (NBEG(LCD) .gt. 0)	then
	   call APPTMP(TMPBUF(NBEG(LCD)),NCH)
C	   call APPTMP(TMPBUF(NBEG(LCD)),6)
	endif
	if (NBEG(LCC) .gt. 0)	then
	   call	APPTMP(TMPBUF(NBEG(LCC)),NCH)
C	else
C Default CC definition is written to detvar.tmp
C	   call APPFML(NCH,4,"ccsp")
C	   write(NCH,'(A)')'      CC(J)=CCSP'
	endif
	if(NBEG(LHC)+NBEG(LXC)+NBEG(LDC).eq.0)	then
	   if (NBEG(LCUBS) .eq. 0) write(NCH,'(A)')'      CUBS(J)=0.'
	   if (NBEG(LCUBS) .gt. 0) call APPTMP(TMPBUF(NBEG(LCUBS)),NCH)
	else
	   write(NCH,'(A,I2.2,A)')'      if (j .eq. NA1) goto 2',NEQ,'5'
	   write(NCH,'(A)')'      if (j .eq. NA) YA = 2./(YD*HROA*HROA)'
	   write(NCH,'(A)')'      CUBS(J)=YA*(FP(J+1)-FP(J))*(0.'
 102	   format('     > +HC(J)*(TE(J+1)-TE(J))/(TE(J+1)+TE(J))')
 103	   format('     > +XC(J)*(TI(J+1)-TI(J))/(TI(J+1)+TI(J))')
 104	   format('     > +DC(J)*(NE(J+1)-NE(J))/(NE(J+1)+NE(J))')
 105	   format('     > )')
	   if(NBEG(LHC).gt.0)	write(NCH,102)
	   if(NBEG(LXC).gt.0)	write(NCH,103)
	   if(NBEG(LDC).gt.0)	write(NCH,104)
	   write(NCH,105)
	   write(NCH,'(A,1I2.2,A)')'      goto 2',NEQ,'6'
	   write(NCH,'(A,1I2.2,A)')' 2',NEQ,'5 continue'
	   write(NCH,'(A)')
     &	'      CUBS(NA1)=CUBS(NA)+(CUBS(NA)-CUBS(NA-1))*HROA/HRO'
	   write(NCH,'(A,1I2.2,A)')' 2',NEQ,'6 continue'
C	   write(NCH,'(A)')'      if (JIT .lt. 4)	CUBS(j)=0'
	endif
	if (NBEG(LCD).eq.0)	then
			write(NCH,'(A)')'      YWA(J)=CUBS(J)'
	else
			call APPTMP(TMPBUF(NBEG(LCD)),NCH)
			write(NCH,'(A)')'      YWA(J)=CUBS(J)+CD(J)'
	endif
	if (ITYPE.eq.1)	then
	    write(NCH,'(A)')'      YWB(J)=YF*CC(J)*YB/IPOL(J)**2'
	    write(NCH,'(A)')'      YWD(J)=YWA(J)*YD/(IPOL(J)**3*G33(J))'
	endif
	write(NCH,"(' 2',1I2.2,'1 continue')")NEQ
C Current distribution CU is adjusted to have U = Const
C	write(NCH,'(A)')'      YWA(NA1)=0.'
	if (ITYPE.eq.1)		then
C	write(NCH,'(A)')'      YWB(NA1)=-1.'
C	write(NCH,'(A)')'      YWD(NA1)=1.'
	write(NCH,'(A)')
     &		'      FP(NA)=.4*GP*HROA*RTOR/(G22(NA)*IPOL(NA1))'
	write(NCH,'(A)')
     &	'      FP(NA1)=(HROA-.5*HRO)*ROC*CC(NA1)/RTOR/IPOL(NA1)/TAU'
	write(NCH,'(A)')'      FP(NA-1)=1.+FP(NA)*FP(NA1)'
C FP(NA1) - HROA*(d\Psi)/(d\rho)|_{at_the_edge}
	if (NBEG(LMV).eq.0)	then
	  write(NCH,'(A)')'      FP(NA1)=FP(NA)*(IPL+FPO(NA1)*FP(NA1))'
	else
	  write(NCH,'(A)')
     &		'      FP(NA1)=FP(NA)*(IPL+FPO(NA1)*FP(NA1))+YDF'
	  write(NCH,'(A/A)')'      CV(NA1) = CV(NA)+YWD(NA1)+',
     &		'     +((CV(NA)+YWD(NA1))-(CV(NA-1)+YWD(NA)))/HRO*HROA'
	endif
	write(NCH,'(A)')'      FP(NA)=-1.'
	write(NCH,'(A)')
     >	'      call RUNF(G22,YWB,YWC,YWD,FPO,NA,YF*TAU,HRO,HROA,FP,FV)'
Ca	write(NCH,'(A)')'      YCU=1.25/(GP*GP*RTOR)'
	write(NCH,'(A)')'      YFP0=(9.*FP(1)-FP(2))/8.'
	write(NCH,'(A)')'      do	J=1,NA1'
Ca	if (NBEG(LMV).eq.0)	then
Ca	write(NCH,'(A)')'      CU(J)=YCU*G33(J)*IPOL(J)**3*YWC(j)'
Ca	else
Ca	write(NCH,'(A)')'      CU(J)=YCU*G33(J)*IPOL(J)**3*YWC(j)-CV(j)'
Ca	endif
Ca	write(NCH,'(A)')'      MU(J)=YWB(j)/(GP2*BTOR)'
C	write(NCH,'(A)')'      UPL(J)=YWD(J)'
	write(NCH,'(A)')'      UPL(J)=(FP(J)-FPO(J))/TAU/YF'
	write(NCH,'(A)')'      ULON(J)=IPOL(J)*G33(J)*UPL(J)'
	write(NCH,'(A)')'      FP(J)=FP(J)-YFP0'
	write(NCH,'(A)')'      FPO(J)=FP(J)'
	write(NCH,'(A)')'      enddo'	
	write(NCH,'(A)')'      UPL(NA1)=ARRNA1(UPL(NA-3:NA),HROA/HRO)'
	write(NCH,'(A)')'      ULON(NA1)=ULON(NA)'
	write(NCH,'(A)')'      call CUOFP'
	if (NBEG(LMV).gt.0)	then
	write(NCH,'(A)')'      do  J=1,NA1'
	write(NCH,'(A)')'      CU(J)=CU(J)-CV(j)'
	write(NCH,'(A)')'      enddo'
	endif
	if (NEQ.eq.0)	then
	   write(NCH,'(A)')
     >		'      if (BETAJR(0.1*ROC) .gt. 100.) then'
	   write(NCH,'(A)')'      write(*,*)char(7)'
	   write(NCH,'(A)')'      write(*,*)'
	   write(NCH,'(A)')
     >'     > ">>>  ERROR  >>> The initial plasma pressure is too high"'
	   write(NCH,'(A)')'      stop'
	   write(NCH,'(A)')'      endif'
C-----------------------------------------------------------------------
	endif
	write(NCH,"('      if(JCALL.gt.0)  goto	2',1I2.2,'0')")NEQ
	return
	endif		!-------------------  end of setting   U = Const

	if (NEQ.eq.0)	write(NCH,'(A)')'      FP(1)=FV(1)'
	if (NEQ.ne.0)	write(NCH,'(A)')'      FP(1)=FPO(1)+UPL(1)*TAU'
	if (ITYPE.eq. 0)	then
	   write(NCH,"('      do	2',1I2.2,'2	J = 1,NA1')")NEQ
C MU distribution is prescribed:
	   call APPTMP(TMPBUF(NBEG(LMU)),NCH)
	   write(NCH,"(' 2',1I2.2,'2 continue')")NEQ
	   write(NCH,'(A)')'      call	CUOFMU'
	endif

		if (ITYPE.eq.-1)	then
	write(NCH,'(A)')'      YF=GP2*HRO*HRO*BTOR'
	write(NCH,'(A)')'      YC=.4*GP*RTOR/BTOR'
	write(NCH,'(A)')'      YM=0.'
	write(NCH,'(A)')'      YMCD=0.'
	write(NCH,'(A)')'      YHC=HROA/HRO'
	write(NCH,"('      do	2',1I2.2,'2	J = 1,NA1')")NEQ
C Current distribution CU is prescribed:
		call APPTMP(TMPBUF(NBEG(LCU)),NCH)
	write(NCH,"('      if (j.eq.NA1) goto 2',1I2.2,'2')")NEQ
	write(NCH,'(A)')
     &		'      YM  =YM +  CU(J)*RHO(J)/(G33(J)*IPOL(J)**3)'
	write(NCH,'(A)')
     &		'      YMCD=YMCD+YWA(J)*RHO(J)/(G33(J)*IPOL(J)**3)'
	write(NCH,"(' 2',1I2.2,'2 continue')")NEQ
	write(NCH,'(A)')'      YIOH = GP2*YM*HRO*IPOL(NA1)'
	write(NCH,'(A)')'      YICD = GP2*YMCD*HRO*IPOL(NA1)'
	write(NCH,'(A)')'      do j=1,NA1'
C Use the 1st line for prescribing CU profile
C Use the 2nd line for prescribing CUOHM profile
	write(NCH,'(A)')'         CU(J)=CU(J)*(IPL-YICD)/YIOH'
C	write(NCH,'(A)')'         CU(J)=YWA(J)+CU(J)*(IPL-YICD)/YIOH'
	write(NCH,'(A)')'      enddo'
	write(NCH,'(A)')'      YM=0.'
	write(NCH,'(A)')'      do j=1,NA'
	write(NCH,'(A)')
     &		'      YM=YM+YC*CU(J)*RHO(J)/(G33(J)*IPOL(J)**3)'
	write(NCH,'(A)')'      YMF=YM/G22(J)'
	write(NCH,'(A)')'      MU(J)=YMF/J'
	write(NCH,'(A)')'      if (j.eq.NA)	then'
	write(NCH,'(A)')'          YMF=YMF*YHC'
	write(NCH,'(A)')'      endif'
	write(NCH,'(A)')'      FP(J+1)=FP(J)+YF*YMF'
	write(NCH,'(A)')'      enddo'
C	write(NCH,'(A)')'      MU(NA1)=MU(NA)*NA/(NA-0.5+YHC)'
	write(NCH,'(A)')'      MU(NA1)=ARRNA1(MU(NA-3:NA),YHC)'
		endif
C----------------------------------------------------------------------|
C	write(NCH,'(A)')'      CU(NA1)=CU(NA)'
	write(NCH,'(A)')'      YU	=GP2*RTOR'
		if (ITYPE.eq.0)	then
	write(NCH,"('      do	2',1I2.2,'3 J = 1,NA1')")NEQ
		endif
		if (ITYPE.eq.-1)	then
C	write(NCH,'(A)')'      YJ	=5.*BTOR/RTOR*ROC*IPOL(NA1)'
C	write(NCH,'(A)')'      YJ	=IPL/(YJ*MU(NA)*G22(NA))'
	write(NCH,'(A)')
     &'      YJ =(FP(NA1)-FP(NA))/HROA*IPOL(NA1)*G22(NA)/(0.4*GP*RTOR)'
	write(NCH,'(A)')'      YJ = IPL/YJ'
	write(NCH,"('      do	2',1I2.2,'3 J = 1,NA1')")NEQ
	write(NCH,'(A)')'      CU(J)	=YJ*CU(J)'
			if (NBEG(LMV).eq.0)	then
	write(NCH,'(A)')'      MU(J)	=YJ*MU(J)'
	write(NCH,'(A)')'      FP(J)	=YJ*(FP(J)-FP(1))+FP(1)'
			endif
			if (NBEG(LMV).ne.0)	then
C	write(NCH,'(A)')'      CU(J)	=YJ*CU(J)+CV(J)'
	write(NCH,'(A)')'      MU(J)	=YJ*MU(J)+MV(J)'
	write(NCH,'(A)')'      FP(J)	=YJ*(FP(J)-FP(1))+FP(1)+FV(J)'
			endif
		endif
	write(NCH,'(A)')'      if (CC(J).gt.0.0) then'
	write(NCH,'(A)')'         ULON(J)=YU*(CU(J)-YWA(J))/CC(J)'
	write(NCH,'(A)')'      else'
	write(NCH,'(A)')'         ULON(J)=0.0'
	write(NCH,'(A)')'      endif'
	write(NCH,'(A)')'      UPL(J)	=ULON(J)/(IPOL(J)*G33(J))'
	write(NCH,"(' 2',1I2.2,'3 continue')")NEQ
	write(NCH,'(A)')'      ULON(NA1) = ARRNA1(ULON(NA-3:NA),HROA/HRO)'
C	write(NCH,'(A)')'      ULON(NA1) = ULON(NA)'
	write(NCH,'(A)')
     &		'      UPL(NA1)  = ULON(NA1)/(IPOL(NA1)*G33(NA1))'
	end
C CUASN ============= Current iterations ==============================|
