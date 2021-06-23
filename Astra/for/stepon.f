C======================================================================|
	program  astra
C----------------------------------------------------------------------|
	implicit double precision (a-h,o-z)
	include	'for/parameter.inc'
	include	'for/const.inc'	
	include	'for/status.inc'
	include	'for/outcmn.inc'
	include	'for/timeoutput.inc'
	include	'tmp/declar.fml'
	include	'tmp/declar.fnc'
	include	'tmp/declar.usr'
	integer getpid,J,IFKEY,IFIPC,JIT,JEX,JDETV,JCALL,ND,ND1,IPART
	integer	IFSBP(NSBMX)
	double precision YWA(NRD),YWB(NRD),YWC(NRD),YWD(NRD),ARRNA1
C	double precision NE1(NRD),TE1(NRD),TI1(NRD)
	double precision YA,YB,YC,YD,YE,YF,YG,YI,YM,YJ,YU,YYD,YCU
	double precision YMF,YHC,YFP0
	data	JEX/1/ JIT/0/ IPART/0/
C----------------------------------------------------------------------|
	call	fenvex()	!  Enable floating exception handling
C	call	fenvdx()	! Disable floating exception handling
C----------------------- Astra start ----------------------------------|
 1	continue
C-------------------- Initial settings --------------------------------|
	include	'tmp/ininam.tmp'
	call	READAT
	call	SETARX(1)
	include	'tmp/inivar.tmp'
	include	'tmp/setvar.tmp'
	include	'tmp/detvar.tmp'
	include	'tmp/inivar.tmp'
C-------------------- Initial iterations ------------------------------|
	JEX = 20
	do   2	JIT=1,JEX
	   call	DEFARR
	   call	INTVAR
	   include 'tmp/detvar.tmp'
	   call	METRIC
	   call	SETARX(1)
	   include 'tmp/inivar.tmp'
	   include 'tmp/init.inc'
	   j = IFKEY(256)
 2	continue		! End initial iteration loop
	IPART = 1
	call	METRIC
C-------------------- Time step loop ----------------------------------|
 3	continue
C	write(7,*)TIME
C	write(7,'(1P,3(E14.6))')(FP(j),1./MU(j),CU(j),j=1,NA1)
	if (IFKEY(0) .eq. 1)	goto	1
	JIT = 0			! Time step loop
	JEX = max(1,int(ITEREX+0.1))
 4	JIT = JIT+1		! External iteration loop (JIT<=JEX)
	call	DEFARR
	if (JIT .eq. 1) then
	   call	INTVAR		! Set simple exp variables
	   include 'tmp/detvar.tmp'
	   call	METRIC		! 
	   call	SETARX(2)	! Set exp arrays (time interpolation on)
	   call	OLDNEW		! Time advance [setting f(t-tau)=f(t)]
	   call	ADCMP		! 
	endif
	include	'tmp/eqns.inc'	! Sbrs & Eqs call;
	if (JIT .lt. JEX)	goto	4
	goto	3
 100	format(3(2F10.5,2X))
	end
C======================================================================|
	subroutine RADOUT
C	implicit none
	implicit double precision (a-h,o-z)
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	include	'tmp/declar.fml'
	include	'tmp/declar.fnc'
	include	'tmp/declar.usr'
C	call	fenvdx()
	include	'tmp/radout.tmp'
C	call	fenvex()
	end
C======================================================================|
	subroutine TIMOUT
C	implicit none
	implicit double precision (a-h,o-z)
	include	'for/parameter.inc'
	include	'for/const.inc'
	include	'for/status.inc'
	include	'for/outcmn.inc'
	include	'for/timeoutput.inc'
	include	'tmp/declar.fml'
	include	'tmp/declar.fnc'
	include	'tmp/declar.usr'
C	call	fenvdx()
	include	'tmp/timout.tmp'
C	call	fenvex()
	end
C======================================================================|
