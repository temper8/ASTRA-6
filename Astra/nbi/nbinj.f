C SUBROUTINES:	NBI, NBSRSR, NB0, NB1TR, NBION0
C		NBION00, NBSISN, NBSPEC, NBSMTV, NBSMTS 	
C 		NBMENU
C This set of routines together with the basic NBI version (NBIDEF.f) 
C	and time dependent Fokker-Planck solver (NBIBCE.f) provides
C	multisource time dependent multisource calculations
C
C The subroutine NBI is called from ASTRA model. All other are internal
C
C================================================= Version by: 03-SEP-2003
C=========================== corrected for time dependent FP  09-DEC-2007 
C Neutral Beam driven current, momentum and power deposition
C (Default version) including Ripple Losses
C-----------------------------------------------------------------------
C	entry:	AMETR,SHIF,NA1,RTOR,AB,BTOR,NI,HBEAM,
C     		RBMIN,RBMAX,CBMS2,CBMI3,EBEAM,QBEAM,ABEAM,DBM1,DBM2,DBM3,
C	     	MU,TE,TI,NE,NN,NNB,ZEF,AMAIN,PBEAM,SCUBM,CONTR,ELON,TRIAN
C NBI interactive control:
C <VARIABLES>:
C	ABEAM=1,2,3 Hydrogen,Deuterium,Tritium NBI only		(Def=1)
C	CONTR=0-1 fraction of contr-inj. power			(Def=1)
C	DBM1=0-1 ~ full energy component EBEAM power fraction 	(Def=1)
C	DBM2=0-1 ~ half energy component EBEAM/2 power fraction (Def=0)
C	DBM3=0-1 ~third energy component EBEAM/2 power fraction (Def=0) 
C WARNING: DBM1 > 0 must be used; SUM(DBMi).ne.1. could be used
C	Power is renormalized withinNBI: Power(i)~DBMi/SUM(DBMi)	
C	EBEAM [keV] the main component energy			(Def=1)
C	HBEAM [m] NBI footprint center height in respect to Z=0	(Def=0)
C	RBMAX [m] max. major radius of the NBI footprint	(Def=0)
C	RBMIN [m] min. major radius of the NBI footprint	(Def=0)
C	QBEAM [MW] NBI power 					(Def=0)
C
Comments: the footprint position is interpreted as:
C
C	for tangential NBI: (RBMAX+RBMIN)/2 > RTOR - AB 
C	as NBI cross-section by the meridianal plane perpendicular 
C	to the vertical plane of the central NBI pencil beam
C	with tangential radius of (RBMAX+RBMIN)/2 
C
C	for perpendicular NBI: (RBMAX+RBMIN)/2 < RTOR - AB 
C	as NBI cross-section by the verticall plane perpendicular to 
C	the central NBI pencil beam (RBMAX+RBMIN)/2 vertical plane,
C	placed at R = RTOR (the corrections, connected with 
C	Shafranov's shift are calculated into the NBI block
C
C <CONSTANTS>
C	CBM1=1,2,.. number of sorces with different geometry	(Def=1)
C		if CBM1 > 1 source parameters are controlled by 
C		special NBMENU rather than VARIABLE, CONSTANTS 
C		menues for each source separately. NBMENU is called 
C		automatically when you increase CBM1. User can 
C		call NBMENU anytime by changing a sign of CBM1 
C		to negative (CBM1 = 4 to CBM1 = -4, etc.)	
C
C CX heat losses due to NBI control
C	
C	CBM2=1 	heat losses (=0 No losses) from bulk ions due 
C		to NBI charge exchange				(Def=1)
C 		CBM2=1 corresponds to the explicit calculations 
C		of the heat sink from the bulk ion component at
C		the moment when NBI is called.
C		Implicit form is recommended:
C		CBM2=0; PIT=...-.0024*SNNBM
C
CBM3,4 are active for CBMI1=2 only
C
C	CBM3=1 	losses (=0 No losses) from fast ions due to CX
C		with the NBI neutrals itself			(Def=1)
C		The same CX cross-section as for cold neutrals
C		is considered.
C	CBM4=1 	losses (=0 No losses) from fast ions due to CX
C		with the cold neutrals				(Def=1)
C
C NBI Source control:
C
C	CBMH1,2	reserved for the footprint power profile control(Def=0)
C	CBMR1-4	reserved for the footprint power profile control(Def=0)
C
C	CBMS1=0 No averaging over drift orbit
C	CBMS1=1 Averaging over drift oribt
C	CBMS1>1 Zero-width drift orbit
C	CBMS2=1,...,82	number of 'pencils'			(Def=1)
C	CBMS3=tg(angle between central pencil and midplane)	(Def=1)
Comments: CBMS3 is used for perpendicular NBI only		
C	CBMS4=(HBmax-HBmin)/(RBMAX-RBMIN): HBEAM=(HBmax+HBmin)/2(Def=1)
C
C Fast ion solver control
C
C	CBMI1=1(2) steady (time dependent) Fokker-Planck Solver	(Def=1)
C		number of angle mesh points NTET1=50/CBMI1+1
CBMI2-4 are active for CBMI1=2 only
C 
C	CBMI2=TAUNBI/TAU ratio of the NBI time-step to TAU	(Def=1)
C WARNING: TAUNBI does not coinside automatically with the 
C		interval between sequental calls of NBI (that 
C		enables the consumption of calculations at the 
C		NBI-steady state phase) 
C	CBMI3=1,2... NBI X-mesh points number N1=(NA1-1)/CBMI3+1(Def=1)
C		     3(EB,EB/2,EB/3),2(EB,EB/2),1(EB)
C	CBMI4=1,2... free parameter                             (Def=1)
C       NBI V-mesh points number IV1=161	
C
C USER's Functions (NBUSER.f) shold be treated from user's sbr/
C
C	NBFHZ(NBI_No,Z,CBMH1,CBMH2), NBFRY(NBI_No,Y,CBMR1,CBMR2)
C		power distribution in the footprint profile 
C		(see comments in the subroutines)
C
C	RIPRAD(ZUPDWN,J) [m] - Ripple loss cone boundary (major 
C		radius) for each magnetic surface: ZUPDWN [m]- 
C		shift in respect to the midplane of ripple 
C		simmetry, J - surface index
C	RIPRAD -depends on tokamak mag. field coils' and 
C		plasma configurations  (Def: No Ripple losses) (Def=999)
C	exit:	PBEAM,PEBM,PIBM,NIBM,CUFI,CUBM,PBLON,PBPER,
c		SCUBM,SNEBM,SNNBM,NNBM1,2,3	for MAIN
c=======================================================22-MAY-08 Polevoi 
	subroutine  NBINJ(NBINP,
     +			BTOR,RTOR,ABC,AB,ROC,SHIFT,UPDWN,
     +			HRO,TAU,NA1,NB1,
     +			AIM1,AIM2,AIM3,AMJ,ZMJ,NNCL,NNWM,QNBI,
     +			CBM1,CBM2,CBMI3,CBMI1,CBM4,CBMI2,CBM3,CBMI4)
	implicit none
	include	'for/parameter.inc'
C----------------------------------------------------------------------|
C	include 'for/const.inc'
C Quantities defined in const.inc and used below
        double precision	BTOR,RTOR,ABC,AB,ROC,SHIFT,UPDWN,
     +		HRO,TAU,
     +		AIM1,AIM2,AIM3,AMJ,ZMJ,NNCL,NNWM,QNBI,
     +		CBM1,CBM2,CBM3,CBM4,CBMI1,CBMI2,CBMI3,CBMI4
        integer	NA1,NB1
C Quantities used below to be excluded from const.inc
        double precision	CBMH1, 	CBMH2,
     ,		CBMR1, 	CBMR2, 	CBMR3, 	CBMR4,
     ,		CBMS1, 	CBMS2,	CBMS3, 	CBMS4,
     ,		EBEAM,	DBM1,	DBM2,	DBM3,
     ,		ABEAM,	HBEAM,	CONTR,	RBMAX,	RBMIN,	QBEAM
C----------------------------------------------------------------------|
	include 'for/status.inc'
	character 	NBINP*(*),FILNAM*32,STR2*2,STRI*25
	logical	EXI
	integer	KILLBL,N,JSRNUM,J,J2FRST,JN,JN1,JNAX,ERCODE,JWARN,INBMS
        integer	NFIELD,lnbinp
        parameter    (NFIELD=20)
        double precision	ARRAY(NFIELD),YQBEAM,VINT
        double precision	YCBMI3,YCBMI4,YABEAM,YEBEAM,YUD,YHM
c=====Ripple	normalized radius of the ripple boundary
c==== banana with RTCRIT > YRIPLR is lost 
        double precision	PBICX,YRIPLR,RIPRAD,Y,Y1,Y2,YMAX,YQBM
	common	/CRIPL/	YRIPLR(NRD)
	save	J2FRST,INBMS
	data 	J2FRST,INBMS,YCBMI3,YCBMI4 /0,0,1.d0,1.d0/
	data 	ABEAM/1.d0/ DBM1/1.d0/ DBM2/0.d0/ DBM3/0.d0/ 
        data    RBMIN/1.d0/ RBMAX/.5d0/
C----------------------------------------------------------------------|

        lnbinp = len(NBINP)
        FILNAM = NBINP(1:lnbinp)
        if (lnbinp .gt. 32)	goto	94
	inquire(file=FILNAM(1:lnbinp),exist=EXI)
C----------------------------------------------------------------------|
	if(.not.EXI)write(*,*)'File "',FILNAM(1:lnbinp),'" not found'
        J = abs(int(CBM1+0.5d0))
        if (CBM1 .lt. .0)	J = J+1
        if (    J .eq. 0)	goto	99
        if (INBMS .eq. 0)	INBMS = J	! 1st call
	if (INBMS .ne. J)	then		! No. of sources changed
           INBMS = J
           goto	1
        endif
	if (EXI .and. CBM1.gt.0.d0 )	goto	2
 1      continue
C Interactive viewer/editor of the beam source parameter list
	call	NQUERY(2,FILNAM(1:lnbinp),INBMS,ERCODE)
        if (ERCODE .eq. 0)	goto	2
        write(*,*)"NQUERY Errcode =",ERCODE
        goto	96
 2      continue

C CBMI3 to make NBI internal mesh interval > max Larmor raduis
	Y	=.001d0
C*DEC*07
        EBEAM   =.001d0
	open(2,file=FILNAM(1:lnbinp),status='OLD')

	do	j=1,INBMS
C Read a record for one beam source
        call	STREAD(2,20,ARRAY,ERCODE)
C	write(*,*)">>> NBINJ: After STREAD: ercode =",ERCODE
C	write(*,'(1I3,1P,5E12.4)')JN,ARRAY(1),CBM1

c	ABEAM = ARRAY(3)	EBEAM = ARRAY(5)

        goto	(200,96,97,98,95,93),ERCODE+1
 200	        if(j.eq.1) then
 	ABEAM = ARRAY(3)	
        EBEAM = ARRAY(5)
        else
           if(EBEAM.lt.ARRAY(5)) EBEAM=ARRAY(5)
	if(CBMI1.GE.2.d0.and.ABEAM.ne.ARRAY(3))
     &  write(*,*) 'ABEAM must be the same for all NBIs if CBM4 > 1'
        endif
        if(Y.LT.(ARRAY(3)*ARRAY(5))) Y=ARRAY(3)*ARRAY(5)
        
	enddo
	close(2)
c--------------------------------------------
        J = 5.d-3*dsqrt(Y)*(NA1/ABC)/(BTOR*RTOR/(RTOR+SHIFT))
        CBMI3=1.d0
	if(J.gt.1.and.J.lt.(NA1-1)) CBMI3=J
	if(5*J.gt.(NA1-1).and.(NA1-1).gt.5) CBMI3=(NA1-1)/5
	IF(CBMI1.lt.1.d0) CBMI1=1.d0
	if(CBMI1.ne.1.d0.and.J2FRST.eq.0) then
           J2FRST	=1
           YCBMI3	=CBMI3
           YCBMI4	=CBMI4
        endif
	if(CBMI1.ne.1.d0.and.J2FRST.ne.0) then
           if(CBMI3.ne.YCBMI3.or.CBMI4.ne.YCBMI4) 
     .	   write(*,*)'CBMI3,CBMI4 can`t be changed after CBMI1=2'
           CBMI3	=YCBMI3
           CBMI4	=YCBMI4
        endif
        QNBI = 0.d0
        CBM1 = INBMS

c=====Ripple vvvvvvvvvvvvvvvvvvvvvvvvvvvv

        YUD	=UPDWN	
        Y	=(AMETR(NA1)+AMETR(NA1-1))/2.d0
        N	=(NA1-1)/CBMI3
        Y1	=0.
C	write(*,*)N+1,NA1,CBMI3,YUD,CBM1,Y,AMETR(NA1),AMETR(NA1-1)
        do	JN1	=N+1,2,-1
           JN	=JN1-1
           JNAX	=1+CBMI3*(JN-1)
           Y2	=RIPRAD(YUD,JNAX)/Y
           if(Y2.GE.Y1)	YMAX	=Y2
           Y1	=Y2
           YRIPLR(JN)=YMAX
	enddo	
        YRIPLR(N+1)=YRIPLR(N)
c=====Ripple ^^^^^^^^^^^^^^^^^^^^^^^^^^
C...Hot ion' source
	do 	JN	=1,NB1
           PBEAM(JN)=0.
           SCUBM(JN)=0.
           SNNBM(JN)=0.
           SNEBM(JN)=0.
           NNBM1(JN)=0.
           NNBM2(JN)=0.
           NNBM3(JN)=0.
	enddo
C Set plasma composition
        call	NBSPEC(AIM1,AIM2,AIM3,AMJ,ZMJ,NA1,NB1)
C	write(*,*)"After NBSEC"
	if(CBMI1.eq.1.d0)	call	NBION00(NB1)
C	write(*,*)"After NBION00",CBM1,INBMS

C        INBMS	=abs(CBM1)
	if(INBMS.eq.0) 		then
           write(*,*) ' -1 < CBM1 < 1 = no NBI sources '
           call	NBION00(NB1)
           return
        endif
C        j = system('test -d /tmp/ && -w /tmp/')
C        j1 = system('test - /tmp/ && -w /tmp/')
C        write(*,*)j
C        open(NCH+1,status='SCRATCH',file=
        JN	=1
        JSRNUM	=JN

        YEBEAM	=EBEAM		! for double precision
c--------------------------------------------
c	write(*,*) 'Eb ', Ebeam, ' Ab ', abeam, ' c4 ',CBMI4
c      	YEBEAM	=CBMI4

        YABEAM	=ABEAM
	JWARN	=0
        YQBM	=0.d0
	open(2,file=FILNAM(1:lnbinp),status='OLD')
        YEBEAM	=EBEAM
 20     JSRNUM	=JN
C Read a record for one beam source
        call	STREAD(2,20,ARRAY,ERCODE)
C	write(*,*)">>> NBINJ: After STREAD: ercode =",ERCODE
C	write(*,'(1I3,1P,5E12.4)')JN,ARRAY(1),CBM1
        goto	(21,96,97,98,95,93),ERCODE+1
 21	QBEAM = ARRAY(1)
	CONTR = ARRAY(2)
	ABEAM = ARRAY(3)
C	ZBEAM = ARRAY(4)
	EBEAM = ARRAY(5)
	DBM1  = ARRAY(6)
	DBM2  = ARRAY(7)
	DBM3  = ARRAY(8)
	CBMS1 = ARRAY(9)
	CBMS2 = ARRAY(10)
	HBEAM = ARRAY(11)
	RBMAX = ARRAY(12)
	RBMIN = ARRAY(13)
	CBMS3 = ARRAY(14)
	CBMS4 = ARRAY(15)
	CBMH1 = ARRAY(16)
	CBMH2 = ARRAY(17)
	CBMR1 = ARRAY(18)
	CBMR2 = ARRAY(19)
C        write(*,*)JN
C        write(*,'(1P,5(E12.4))')(ARRAY(j),j=1,15)
C Beam footprint drawing:
	write(STRI,'(a,i3,a)')" NBINJ: Beam #",JN,"    "//char(0)
	if (QBEAM .gt. 0.d0) call	DSPMSG(STRI)
C        call	DRFOOT(JN,RBMIN,RBMAX,HBEAM,CBMS4,QBEAM)
	if (YABEAM.ne.ABEAM) JWARN=1
C	write(*,'(1P,5E12.4)')ARRAY
C	write(*,'(1P,4E12.4)')CBMH1, 	CBMH2
C	write(*,'(1P,4E12.4)')CBMR1, 	CBMR2
C	write(*,'(1P,4E12.4)')CBMS1, 	CBMS2,	CBMS3,	CBMS4
C	write(*,'(1P,4E12.4)')EBEAM,	DBM1,	DBM2,	DBM3
C	write(*,'(1P,4E12.4)')ABEAM,	HBEAM,	CONTR
C	write(*,'(1P,4E12.4)')RBMAX,	RBMIN

        QNBI = QNBI+QBEAM
        YQBEAM	=QBEAM
C...Switch off the ion source
        YQBM	=YQBM+YQBEAM
	YHM = (HBEAM-YUD)*100.
	if(CONTR.ne.0.and.CONTR.ne.1.)		then
C...balanced injection
C: coinj. part
           QBEAM	=YQBEAM*(1.-CONTR)
           call	NBSRSR(JSRNUM,  -1.d0,
     ,		NA1,	RTOR,	SHIFT,	AB,	BTOR,
     ,		NNCL,	NNWM,	HRO,	YHM,
     ,		CBMH1, 	CBMH2,	CBMS1, 	CBMS2,	CBMS3, 	CBMS4,
     ,		CBMR1, 	CBMR2, 	CBMI3, 	CBMI1,
     ,		EBEAM,	DBM1,	DBM2,	DBM3,
     ,		ABEAM,	CONTR,	QBEAM,	RBMAX,	RBMIN)
           if(CBMI1.EQ.1.d0)
     >		call	NBION0(NA1,NNCL,NNWM,ABEAM,EBEAM,RTOR,CBMI3)
C: contrinjection part...
           QBEAM	=YQBEAM*CONTR
           call	NBSRSR(JSRNUM,1.d0,
     ,		NA1,	RTOR,	SHIFT,	AB,	BTOR,
     ,		NNCL,	NNWM,	HRO,	YHM,
     ,		CBMH1, 	CBMH2,	CBMS1, 	CBMS2,	CBMS3, 	CBMS4,
     ,		CBMR1, 	CBMR2, 	CBMI3, 	CBMI1,
     ,		EBEAM,	DBM1,	DBM2,	DBM3,
     ,		ABEAM,	CONTR,	QBEAM,	RBMAX,	RBMIN)
           if(CBMI1.EQ.1.d0)
     >		call	NBION0(NA1,NNCL,NNWM,ABEAM,EBEAM,RTOR,CBMI3)
           QBEAM	=YQBEAM
        else
           if(CONTR.EQ.1.d0)	call	NBSRSR(JSRNUM,1.d0,
     ,		NA1,	RTOR,	SHIFT,	AB,	BTOR,
     ,		NNCL,	NNWM,	HRO,	YHM,
     ,		CBMH1, 	CBMH2,	CBMS1, 	CBMS2,	CBMS3, 	CBMS4,
     ,		CBMR1, 	CBMR2, 	CBMI3, 	CBMI1,
     ,		EBEAM,	DBM1,	DBM2,	DBM3,
     ,		ABEAM,	CONTR,	QBEAM,	RBMAX,	RBMIN)
           if(CONTR.EQ.0.d0)	call	NBSRSR(JSRNUM,-1.d0,
     ,		NA1,	RTOR,	SHIFT,	AB,	BTOR,
     ,		NNCL,	NNWM,	HRO,	YHM,
     ,		CBMH1, 	CBMH2,	CBMS1, 	CBMS2,	CBMS3, 	CBMS4,
     ,		CBMR1, 	CBMR2, 	CBMI3, 	CBMI1,
     ,		EBEAM,	DBM1,	DBM2,	DBM3,
     ,		ABEAM,	CONTR,	QBEAM,	RBMAX,	RBMIN)
           if(CBMI1.EQ.1.d0)
     >		call	NBION0(NA1,NNCL,NNWM,ABEAM,EBEAM,RTOR,CBMI3)
	endif
        JN	=JN+1
 	if(JN.le.INBMS) goto 20	
	close(2)
	close(39)
        EBEAM	=YEBEAM
c--------------------------------------------
        QBEAM	=YQBM

C...Fast ion' distribution (3D + t)
	if(CBMI1.GE.2.d0)		then
           if(JWARN.eq.1)	then
              write(*,*) 'ABEAM must be the same for CBMI1#1'
              ABEAM=YABEAM
           endif		
           call	NBIONR(EBEAM,ABEAM,AMJ,RTOR,NA1,TAU,NNCL,NNWM,
     >		CBM1,CBM2,CBM3,CBM4,CBMI1,CBMI2,CBMI3,CBMI4)
	endif
C...Conversion to rough mesh keeping the intagrals

	if(CBMI3.ne.1.d0) 			then
           call	NBSMTV(NIBM,CBMI3,ROC,NA1)
           call	NBSMTV(PIBM,CBMI3,ROC,NA1)
           call	NBSMTV(PEBM,CBMI3,ROC,NA1)
           call	NBSMTV(PBLON,CBMI3,ROC,NA1)
           call	NBSMTV(PBPER,CBMI3,ROC,NA1)	
           call	NBSMTV(PBEAM,CBMI3,ROC,NA1)
           call	NBSMTV(SNNBM,CBMI3,ROC,NA1)
           call	NBSMTV(SNEBM,CBMI3,ROC,NA1)
           call	NBSMTV(SCUBM,CBMI3,ROC,NA1)
           call	NBSMTS(CUFI,CBMI3,ROC,NA1)
           call	NBSMTS(CUBM,CBMI3,ROC,NA1)
        endif

	if(CBMI1.eq.1.d0.and.CBM2.gt.0.d0)		then
	   do j=1,NA1
              include 'fml/pbicx'
              PIBM(J)		=PIBM(J)-PBICX
          
           enddo
	endif
	return
 93     write(*,*)'>>> NBI >>> Wrong NBI configuration file format: '
        write(*,'(A,I2,2A,I2,A)')'      ',JN-1,' beam records ',
     >		'available, ',int(CBM1),' records required.'
        stop
 94     write(*,*)'>>> NBI >>> Input file name is too long'
        stop
 95     write(*,*)'>>> NBI calling STREAD: array out of limits'
        stop
 96     write(*,*)'>>> NBI >>> Error in file "',FILNAM(1:lnbinp),'"'
     >		 , ': unrecognized variable name'
        stop
 97     write(*,*)'>>> NBI >>> File "',FILNAM(1:lnbinp),'" read error'
        stop
 98     write(*,*)'>>> NBI >>> Wrong NBI configuration file format. '
        write(*,*)'            More records expected than available.'
        stop
 99     write(*,*)'>>> NBI >>> Zero number of sources'
        write(*,*)"            Don't know what to do"
        stop
	end
C======================================================================|
	subroutine NQUERY(NCH,FILNAM,NBS,ERCODE)
C----------------------------------------------------------------------|
        implicit none
        integer	 j,j1,jj,i,NBS,NBSMAX,NFIELD,IFNUM,KILLBL,NOBLAN,length
C Input:
C  NBS    - Requested No. of NBI sources
C  NCH    - Logical unit (must be connected to the file NBINP)
C Output:   Data are written into file NBINP
C  ERCODE - Error code (0 for normal exit)
C Internal:
C  NBSMAX - Max. No. of NBI sources (max No. of groupss in *.nbi file)
C  NFIELD - No. of fields in one group of *.nbi file
        parameter    (NBSMAX=16, NFIELD=20)
	character    FILNAM*(*),SFIELD(NFIELD)*12,VALUE*6
	character*80 STR,STRI,ADATA(NBSMAX),BDATA(NBSMAX),CDATA(NBSMAX)
        double precision	YCB
        integer	NCH,ERCODE,lnbinp
        logical	ifexi
C----------------------------------------------------------------------|
C Old NBI file: 25xNBS fields to read
C		CBMH1 	CBMH2	CBMH3 	CBMH4
C		CBMR1 	CBMR2 	CBMR3 	CBMR4
C		CBMS1 	CBMS2	CBMS3 	CBMS4
C		EBEAM	DBM1	DBM2	DBM3
C		ABEAM	HBEAM	CONTR	QBEAMy
C		RBMAX	RBMIN	YCBEAM	YJBEAM
C		SRSNo
C New NBI file: 
C   Parameters
C		QBEAM	CONTR	ABEAM	ZBEAM	EBEAM
C		DBM1	DBM2	DBM3	CBMS1	CBMS2
C     Geometry
C		HBEAM	RBMAX	RBMIN	CBMS3	CBMS4
C		CBMH1	CBMH2	CBMR1	CBMR2   Not_used
C----------------------------------------------------------------------|
        if (NBS .gt. NBSMAX)	then
           write(*,'(2(A,I3,A))')
     >	   '>>> NQUERY: Too many NB sources NBS =',NBS,' requested',
     >	   '            Call ignored,    NBSmax =',NBSMAX,' is allowed'
           return
        endif
        if (NBS .le. 0)	then
           write(*,'(2(A,I3))')
     >	   '>>> NQUERY: Number of sources must be positive, NBS =',NBS
           return
        endif
C Default definition:
        do	jj=1,NBS
           write(ADATA(jj)(1:2),'(1I2)')jj
           write(CDATA(jj)(1:2),'(1I2)')jj
           ADATA(jj)(3:) = "   zrd1     0.     2.     1.   100."
     +			   // "     .7     .2     .1     1.    10."
           CDATA(jj)(3:) = "     .3    1.6     1.     0.     1."
     +			   // "     9.     2.     9.     2.     0."
        enddo
        lnbinp = len(FILNAM)
C        write(*,*)lnbinp
        ifexi = .TRUE.
	open(NCH,file=FILNAM(1:lnbinp),status='OLD',err=1)
        goto	3
 1      ifexi = .FALSE.
        open(NCH,file=FILNAM(1:lnbinp),status='NEW',err=95)
 3      continue
        do	jj=1,NBS
 4         read(NCH,'(A)',end=10)STR ! Skip string
           j  = index(STR,'!')				! Comment
           if (j .ne. 0)	goto	4

C           read(NCH,'(A)',end=10)STR
           read(NCH,'(5A)',err=99)SFIELD
C           read(NCH,'(4(<NFIELD>A12))')SFIELD
C           write(*,'(4(5A12/))')SFIELD
C           write(*,'(4(5I3/))')(IFNUM(SFIELD(j),12),j=1,NFIELD)
           i = -3
           do	j=1,NFIELD
              i = i+7
              call	UPCASE(12,SFIELD(j))
              if (IFNUM(SFIELD(j),12) .ne. 0)	then
                 read(SFIELD(j),*)YCB
                 call	FMTF6(VALUE,YCB)
C                 write(*,*)'"',VALUE,'"',YCB
              else
                 j1 = noblan(SFIELD(j))
                 VALUE = SFIELD(j)(j1:)
C                 j1 = length(VALUE)			! Shift a
C                 if (j1 .eq. 4) VALUE='  '//VALUE(1:4)	! string to
C                 if (j1 .eq. 5) VALUE= ' '//VALUE(1:5)	! the right
              endif
              if (j .le. 10)	then
                 write(ADATA(jj)(i:i+5),'(1A6)')VALUE
              else
                 write(CDATA(jj)(i:i+5),'(1A6)')VALUE
              endif
              if (j .eq. 10)	i = -3
           enddo				! enddo j
        enddo					! enddo jj
 10     continue
        do	jj=1,NBS
           BDATA(jj) = ADATA(jj)
C           write(*,*)'"',BDATA(jj)(1:72),'"'
        enddo
C Template for the 1st string of a table
	STRI = 
C           ---------1---------2---------3---------4---------5---
     +	   "# |QBeam |Contr |ABeam |ZBeam |"// ! QBEAM,CONTR,ABEAM,ZBEAM
     +        "EBeam |DBeam1|DBeam2|DBeam3|"// ! EBEAM,DBM1, DBM2, DBM3
     +	      "Orb_av|Penc.#"//char(0)         ! CBMS1,CBMS2
	STR = "NBI configuration file: "//FILNAM(1:lnbinp)//char(0)
        j = len(ADATA(1))
	call	ASKCOL(STR,STRI,ADATA,j,NBS,0,0)
C        do	jj=1,NBS
C           if (ADATA(jj)(1:72) .ne. BDATA(jj)(1:72))	WRI = .TRUE.
C        enddo
        do	jj=1,NBS
           BDATA(jj) = CDATA(jj)
C           write(*,*)'"',BDATA(jj)(1:72),'"'
        enddo
C        HBEAM  RBMAX  RBMIN  CBMS3  CBMS4
C        CBMH1  CBMH2  CBMR1  CBMR2  
	STRI = 
     +	"# |HBeam |RBmax |RBmin |tg(A) |Aspect|"//
     +  "Cver1 |Cver2 |Chor1 |Chor2 |Unused"//char(0)
C	STR = "NB parameters"//char(0)
	STR="NBI configuration file (cnt.): "//FILNAM(1:lnbinp)//char(0)
        j = len(CDATA(1))
C        call	sleep(1)
C        call	IFKEY(82)
	call	ASKCOL(STR,STRI,CDATA,j,NBS,0,0)
C        do	jj=1,NBS
C           if (CDATA(jj)(1:72) .ne. BDATA(jj)(1:72))	WRI = .TRUE.
C        enddo
C        write(*,'(A)')('"'//ADATA(jj)(1:72)//'"',jj=1,NBS)
C        write(*,'(A)')('"'//CDATA(jj)(1:72)//'"',jj=1,NBS)

C        if (.not.WRI)	goto	77
        rewind(NCH)		! Modify input file

        if (.not.ifexi)	then
           write(NCH,'(1I3)')1
           goto	12
        endif
 11     read(NCH,'(A)',end=99)STR ! Skip string
        j  = index(STR,'!')    ! Comment
        if (j .ne. 0)	goto	11
 12     continue

        do	jj=1,NBS
           i = -3
           do	j=1,NFIELD
              i = i+7
              if (j .le. 10)	then
                 write(VALUE,'(1A6)')ADATA(jj)(i:i+5)
              else
                 write(VALUE,'(1A6)')CDATA(jj)(i:i+5)
              endif
              call	UPCASE(6,VALUE)
              if (IFNUM(VALUE,6) .ne. 0)	then
                 read(VALUE,*)YCB
                 write(SFIELD(j),'(1P,E12.4)')YCB
              else
                 j1 = noblan(VALUE)
                 VALUE = VALUE(j1:)//'     '
C                 SFIELD(j)(1:) ='  '//VALUE//'    '
                 SFIELD(j)(1:) =VALUE//'      '
              endif
              if (j .eq. 10)	i = -3
           enddo
C          write(*,'(1I3)')jj
           if (jj .ne. 1)write(NCH,'(1I3)')jj
C          write(*,'(3(5A12/),5A12)')SFIELD
           write(NCH,'(3(5A12/),5A12)')SFIELD
        enddo
 77     close(NCH)
        ERCODE = 0
        return
 95     write(*,*)'>>> NBI >>> Cannot open file "',FILNAM(1:lnbinp),'"'
        stop
 99     close(NCH)
        ERCODE = 1
        end
C======================================================================|
        subroutine	NBION00(JB1)
c------------------------------------------------------------ 17.09.91
C	All NBI arrays = 0
C---------------------------------------------------------------------
	implicit none
	include  'for/parameter.inc'
	include  'for/status.inc'
        integer	J,JB1
	do 1	J	=1,JB1
		PEBM(J)=0.
		PIBM(J)=0.
		CUBM(J)=0.
		CUFI(J)=0.
		PBPER(J)=0.
		NIBM(J)=0.
 1		PBLON(J)=0.
	RETURN
	end
C<<<<<	NBSMTV,	NBSMTS 		>>>>>>>>>>>>>>>>smoothing of NBI
C================================================================Polevoi
	subroutine	NBSMTV(YFO,YCI3,YROC,JNA1)
C-----------------------------------------------------------------------
C		The subroutine transmits a histogram like function YFO(1:N1)
C		from an N1 grid X(1:N1) to a smooth functon YFO(1:JNA1)
C		keeping the same total volume integral.
C	Input:	YFO(1:N1)
C	Output:	YFO(1:JNA1)
C------------------------------------------------------------ 22-MAY-08
	implicit none
	include  'for/parameter.inc'
	include  'nbi/nbicom.inc'
	double precision	VINT,YFO(*),Y,YOLD,YNEW,ALFA,YCI3,YROC
	integer	j,JNA,JNAC,JNA1,JN,JN1
	do 	j	=1,n1-1
		JNA	=1+YCI3*(j-1)
		JNAC	=JNA-1+YCI3
!		JNAC	=1+YCI3*(j-.5)
		X(j)	=XJ(JNAC)
		DRI(j)	=YFO(JNA)
	enddo
!		DRI(N1)	=YFO(JNA1)
		DRI(N1)	=0.d0
!YFO(JNA1)

                X(N1)=1.d0
                XJ(JNA1)=1.d0

C...Total power normalization
c
		YOLD	=VINT(YFO,YROC)
		ALFA	=0.001d0
c	call transf(n1,dri,x,na1,yfo,xj)
	call SMOOTH(ALFA,n1,dri,x,jna1,yfo,xj)
C*16-MAR-98 vvvvvvvvvvvvvvvvvvvvvvvv
	do	J	=1,JNA1
	if(YFO(j).lt.0.d0)	YFO(J)	=0.d0
	enddo
C*16-MAR-98 ^^^^^^^^^^^^^^^^^^^^^^^^ 
		YNEW	=VINT(YFO,YROC)
	IF(dabs(YNEW).gt.1.d-19) 		then
		Y	=YOLD/YNEW

	do	J	=1,JNA1
		YFO(J)	=Y*YFO(J)
	enddo
						endif
	return
	end
C=======================================================================
	subroutine	NBSMTS(YFO,YCI3,YROC,JNA1)
C-----------------------------------------------------------------------
C	The subroutine transmits a histogram like function YFO(1:N1)
C		from an N1 grid X(1:N1) to a smooth functon YFO(1:JNA1)
C		keeping the same total surface integral.
C	Input:	YFO(1:N1)
C	Output:	YFO(1:JNA1)
C------------------------------------------------------------- 22-MAY-08
	implicit none
	include  'for/parameter.inc'
	include  'nbi/nbicom.inc'
	double precision	IINT,YFO(*),Y,YOLD,YNEW,ALFA,YCI3,YROC
	integer	j,JNA,JNAC,JNA1,JSIGN
	external IINT
        JSIGN=1         ! change of sign? (1= yes, 0= no)
	do 	j	=1,n1-1
		JNA	=1+YCI3*(j-1)
		JNAC	=JNA-1+YCI3
!		JNAC	=1+YCI3*(j-.5)
		X(j)	=XJ(JNAC)
		DRI(J)	=YFO(JNA)
                if(DRI(1)*DRI(J).lt.0.d0) JSIGN=0
	enddo
		dri(n1)	=0.d0
!YFO(JNA1)
                X(N1)=1.d0
                XJ(JNA1)=1.d0
C0.
C...Total current normalization
		YOLD	=IINT(YFO,YROC)
		ALFA	=0.001d0
c	call transf(n1,dri,x,na1,yfo,xj)
	call SMOOTH(ALFA,n1,dri,x,jna1,yfo,xj)
C*22-MAY-22 vvvvvvvvvvvvvvvvvvvvvvvv
        if(JSIGN.eq.1) then
	do	J	=1,JNA1
	if(YOLD.gt.0d0)then
	if(YFO(J).lt.0.d0) YFO(J)	=0.d0 !cut negative
        else
	if(YFO(J).gt.0.d0) YFO(J)	=0.d0 !cut positive
        endif
	enddo
        endif
C*22-MAY-98 ^^^^^^^^^^^^^^^^^^^^^^^^ 
		YNEW	=IINT(YFO,YROC)
	IF(abs(YNEW).gt.1.e-19) 		then
		Y	=YOLD/YNEW

	do	J	=1,JNA1
		YFO(J)	=Y*YFO(J)
	enddo
						endif
	return
	end
C=======================================================================
