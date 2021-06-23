C The subroutine is used for smoothing of a radial profile
C 		
C Example: smoothed transfer 'fml/hagb' to 'fml/hagbs'
C
C	if (j.eq.1)	then
C		do  jj	=1,NA
C		j = jj
C		include	'fml/hagb'
C		if (jj .eq. 1)	YY = HAGB
C		call	SMOFML(-1.,JJ,HAGB,'HAGBS')
C		enddo
C		j = 1
C		HAGB = YY
C	endif
C	call	SMOFML(.002,J,HAGBS,'HAGBS')
C
	subroutine	SMOFML(YALPHA,JJ,FMLVAL,FMLNAM)
	implicit none
	include	'for/parameter.inc'
	double precision YALPHA,YAMESH(NRD),YFORIG(NRD),YFML(NRD,10)
	double precision FMLVAL
	integer		jn,jf,jna1,jj,jk,j
	character*6	FMLNUM(10),FMLNAM,FNAM6
	save	JNA1,FMLNUM,YAMESH,YFORIG,YFML,jn
	data	FMLNUM	/10*'      '/

	if (YALPHA.lt.0.0 .and. JJ.gt.1)	goto	14
C assign FMLNAM appended with spaces to FNAM6
	jn	= 0
		do	12	jf	=1,6
		if(jn.eq.1)	goto	11
	if(FMLNAM(jf:jf).eq.' '.or.FMLNAM(jf:jf).eq.'\0')  jn = 1
		if(jn.eq.1)	goto	11
		FNAM6(jf:jf) = FMLNAM(jf:jf)
		goto	12
 11		FNAM6(jf:jf) = ' '
 12		continue
		FNAM6(6:6) = '\0'

C Name identification
	do	jn=1,11
C Arrays are too short:
	if ( jn .eq. 11 )  then
	write(*,*)'>>> SMOFML <<< Warning:',
     .	' Too many formulas for smoothing.\n',
     .		'      Only the following ten are treated:\n',
     .			(FMLNUM(jk),' ',jk=1,10)
		pause 'Increase arrays FMLNUM(10),YFML(...,10)'
			   endif
C  Repeated entries:
	if (FMLNUM(jn).eq.FNAM6)	goto	14

C  First entry with a given name:
	if (FMLNUM(jn)(1:1) .eq. ' ')	then
		FMLNUM(jn) = FNAM6
		goto	14
					endif
	enddo

 14	if(YALPHA .ge. 0.)	GOTO	2
C   First entry.    Original function storage.  Mesh size (NA1) calculation.
		JNA1	=max(0,JJ)+1
		YFORIG(JJ) = FMLVAL
		return

 2	if(JJ .gt. 1)	goto	10
C   First call for smoothing:
	YFORIG(JNA1) = YFORIG(JNA1-1)
C   Radial mesh assignment:
		do	j=1,JNA1-1
		YAMESH(j)=j/(JNA1-0.5)
		enddo
		YAMESH(JNA1)=1.
	call SMOOTH(YALPHA,JNA1,YFORIG,YAMESH,JNA1,YFML(1,jn),YAMESH)

C   Repeated calls.  Smoothed function assignment:
 10	FMLVAL	=YFML(JJ,jn)
	end
