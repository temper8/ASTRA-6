C SDDN [10#19 n/s]:	 Integral {0,R} ( SVD1*NDEUT**2 ) dV
C Neutron rate from thermal D-D reaction
C			(Polevoy 25-JUN-97)
	double precision FUNCTION SHINER(YR)
	implicit none
	double precision Y,YR,Y1
	integer  J,JSRNUM
	include  'for/parameter.inc'
	include  'for/const.inc'
	character*12 SHNAME
	DATA SHNAME/'dat/shth.dat'/
 98	FORMAT(1I1)
 99	FORMAT(1I2)
		JSRNUM	= ABS(CNB1)
		Y	=0.
	DO 	J=1,JSRNUM 
	if(JSRNUM.gt.99) then
	 write(*,*) '>>> Too many NBI sources > 99 (see NBI1TR)'
		return
			endif			
	if(JSRNUM.lt.10) write(SHNAME(12:12),98) JSRNUM
	if(JSRNUM.gt.9.and.JSRNUM.lt.100) write(SHNAME(11:12),99) JSRNUM
	if(JSRNUM.gt.99) then
	 write(*,*) '>>> Too many NBI sources > 99 (see NBI1TR)'
		return
			endif			
			open(38,FILE=SHNAME,ERR=1)
		READ(38,*) Y1
		Y	=Y+Y1
	goto 2
 1	return
 2	close(38)
 	ENDDO
		SHINER	=Y
	END
