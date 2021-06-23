	subroutine TWAVE(ARREXT,AMPLIT,PHASE,DRIVER,YAMPL,YPHASE)
	implicit none
	double precision ARREXT(*),AMPLIT(*),PHASE(*),DRIVER,YAMPL
	include	'for/parameter.inc'
	include  'for/const.inc'
	double precision YA1(NRD),YA2(NRD),YTEMX(NRD),YPHASE,TEMPMN,TST
	double precision FREQW,DRIVE0,TIMEON,TIME0,PERIOD,DF1,DF2
	integer	JCALL,JJ,J
	save	JCALL,JJ,TST
	data	JCALL/1/	JJ/0/
	GOTO (10,8),JCALL

	JJ = JJ+1
	if(JJ.eq.NA1)	JJ = 1
	if( (DRIVE0.lt.0.) .and. (DRIVER.ge.0.) ) then
		TIMEON = TIME0+(TIME-TIME0)*DRIVER/(DRIVER-DRIVE0)
		PERIOD = TIMEON-TST
		FREQW  = GP2/PERIOD
		TST = TIMEON
	endif

	do  1  j = 1,NA1
		DF1	= YA2(j)-YA1(j)
		DF2	= ARREXT(j)-YA2(j)
		if (DF1*DF2.gt.0.)	goto	1
	if (DF1.gt.0.)	then
		YTEMX(j) = YA2(j)
	else
		TEMPMN = YA1(j)
		AMPLIT(j) = YTEMX(j)-YA1(j)
		PHASE(j)  = (TIME-TST)*FREQW
	endif
 1	continue
	YAMPL	= AMPLIT(JJ)
	YPHASE	= PHASE(JJ)

 8	JCALL = 3
	do	9	j = 1,NA1
		YA1(j)	= YA2(j)
		YA2(j)	= ARREXT(j)
 9	continue
	TIME0  = TIME
	DRIVE0 = DRIVER
	return

 10	JCALL = 2
	TST = TIME
	do	11	j = 1,NA1
		PHASE(j) = 0.
		AMPLIT(j)= 0.
		YA2(j)	 = ARREXT(j)
 11	continue
	END
