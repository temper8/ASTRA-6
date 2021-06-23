C======================================================================|
C Make as
C $> cc -c xpr/neut.c -o xpr/neut.o
C $> cc -c xpr/ipc.c -o xpr/ipc.o
C $> f95f xpr/neut.f xpr/neut.o xpr/ipc.o -o xpr/neut
C----------------------------------------------------------------------|
	program   main
	implicit  none
	integer	  j,iargc,getpid,mampid,mamkey,eignr
	character*132 STRING,eigpath,mampath
C----------------------------------------------------------------------|
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
	call sbp2shm(eigpath, mampath, mampid, mamkey, eignr)
C	do j=0,iargc()
C	   write(*,*)
C	   call	getarg(j,STRING)
C	   write(*,*)j,'"',STRING(1:len_trim(STRING)),'"'
C	   write(*,*)len(STRING),len_trim(STRING)
C	enddo
	end
C----------------------------------------------------------------------|
C======================================================================|
	subroutine A_GLF23P(
     >		IS,
     >		IE,
     >		NA1,
     >		NA1N,
     >		NA1E,
     >		NA1I,
     >		HRO,
     >		HROA,
     >		ROC,
     >		BTOR,
     >		RTOR,
     >		AMJ,
     >		AIM1,
     >		NE,
     >		TE,
     >		NI,
     >		TI,
     >		ZEF,
     >		ZIM1,
     >		AMAIN,
     >		MU,
     >		RHO,
     >		AMETR,
     >		SHIF,
     >		ELON,
     >		VTOR,
     >		ER,
     >		NIBM,
     >		IPOL,
     >		G11,
     >		VRS,
     >		GRADRO,
     >		SHEAR,
C output
     >		CHI,
     >		CHE,
     >		DIF,
     >		VIN,
     >		DPH,
     >		DPL,
     >		DPR,
     >		XTB,
     >		EGM,
     >		GAM,
     >		GM1,
     >		GM2,
     >		OM1,
     >		OM2,
     >		FR1	)
C----------------------------------------------------------------------|
c based on stand-alone driver for the GLF23 model
c       "testglf.f" 18-fev-03 version 1.61
c       written by Jon Kinsey, General Atomics
C----------------------------------------------------------------------|
C WORK(1:NA1,1:13) array is used for output
C                              (when i_delay=0 and egamma_d is not used)
C----------------------------------------------------------------------|
C This subroutine computes gradients and 
C calls GLF23 on the point-by-point basis	(Pereverzev 08-08-2007)
C----------------------------------------------------------------------|
	implicit none
C	include	'for/parameter.inc'
C	include 'for/const.inc'
C	include 'for/status.inc'
C Parameters used:
C In A_vars & A_arrs:
	integer			NA1,IS,IE
	double precision	HRO,ROC,BTOR,RTOR,AMJ,AIM1
	double precision	NE(*),TE(*),NI(*),TI(*),ZEF(*),ZIM1(*)
	double precision	AMAIN(*),ER(*)
	double precision	MU(*),RHO(*),AMETR(*),SHIF(*),ELON(*)
	integer			NA, NA1N, NA1E, NA1I
	double precision	HROA
	double precision VTOR(*),NIBM(*),IPOL(*),G11(*),VRS(*)
	double precision ALMHD,ROTSH,CS,GRADRO(*),SHEAR(*)
	double precision CHI(*),CHE(*),DIF(*),VIN(*),DPH(*),DPL(*)
	double precision DPR(*),XTB(*),EGM(*),GAM(*),GM1(*),GM2(*)
	double precision OM1(*),OM2(*),FR1(*)
C----------------------------------------------------------------------|
      integer jpd,jna
      parameter ( jpd=100 )
      real*8 te_m(0:jpd), ti_m(0:jpd)
     & , ne_m(0:jpd), ni_m(0:jpd), ns_m(0:jpd)
     & , zgte_m(0:jpd), zgti_m(0:jpd), zgne_m(0:jpd), zgni_m(0:jpd)
     & , angrotp_exp(0:jpd), egamma_exp(0:jpd), gamma_p_exp(0:jpd)
     & , vphi_m(0:jpd), vpar_m(0:jpd), vper_m(0:jpd)
     & , zeff_exp(0:jpd), bt_exp, bteff_exp(0:jpd), rh_m(0:jpd),arho_exp
     & , gradrho_exp(0:jpd), gradrhosq_exp(0:jpd)
     & , rmin_exp(0:jpd), rmaj_exp(0:jpd), rmajor_exp
     & , q_exp(0:jpd), shat_exp(0:jpd), alpha_exp(0:jpd)
     & , elong_exp(0:jpd), zimp_exp, amassimp_exp, amassgas_exp
     & , alpha_e, x_alpha
      real*8 y1,y2,zpte_in, zpti_in, zpne_in, zpni_in, drho
      real*8 diffnem, chietem, chiitim
     & , etaphim, etaparm, etaperm, exchm
     & , diff_m(0:jpd), chie_m(0:jpd), chii_m(0:jpd), etaphi_m(0:jpd)
     & , etapar_m(0:jpd), etaper_m(0:jpd), exch_m(0:jpd)
     & , egamma_m(0:jpd), egamma_d(0:jpd,10), gamma_p_m(0:jpd)
     & , anrate_m(0:jpd), anrate2_m(0:jpd)
     & , anfreq_m(0:jpd), anfreq2_m(0:jpd)
c
      integer lprint, nroot, jshoot, jmm, jmaxm, itport_pt(1:5)
     & , igrad, idengrad, i_delay, j, k, leigen, irotstab, bt_flag, iglf
C----------------------------------------------------------------------|
C	write(*,'("--GLF",6I5,1P,7(E12.3))')
C     >		IS,
C     >		IE,
C     >		NA1,
C     >		NA1N,
C     >		NA1E,
C     >		NA1I,
C     >		HRO,
C     >		HROA,
C     >		ROC,
C     >		BTOR,
C     >		RTOR,
C     >		AMJ,
C     >		AIM1

C     >		NE,
C     >		TE,
C     >		NI,
C     >		TI,
C     >		ZEF,
C     >		ZIM1,
C     >		AMAIN,
C     >		MU,
C     >		RHO,
C     >		AMETR,
C     >		SHIF,
C     >		ELON,
C     >		VTOR,
C     >		ER,
C     >		NIBM,
C     >		IPOL,
C     >		G11,
C     >		VRS,
C     >		GRADRO,
C     >		SHEAR,

c git: modern block from ~/glf23_v1.60/testglf.f, some values rechosen
c or adapted to ASTRA needs
	leigen = 1		! 1 for tomsqz, for cgg eigenvalue solver
	leigen = 0		! 1 for tomsqz, for cgg eigenvalue solver
C	nroot  = 8	! n. of roots in eigenvalue solver
	nroot  = 12	! n. of roots in eigenvalue solver (impurity dynamics)
	iglf   = 1		! 1 new, 0 original GLF23 normalization
	jshoot = 0		! for time-dependent code
	igrad  = 1		! 1 input gradients, 0 compute gradients
	jna = max(NA1E,NA1I,NA1N,NA1)
	if (jna .gt. NA1) jna = NA1
	if (jna .eq. 0) then
	   write(*,*)' >>> GLF23: Zero grid encountered. Call ignored'
	   return
	endif
	jmaxm = jna-1		! WARNING! jmaxm cannot exceed jpd: jmaxm < jpd
	NA = NA1-1
	if (jmaxm .ge. jpd)	then
	   write(*,'(A/A,I3,A,I4,A)')
     >		' >>> GLF23: Grid allocation failure. Call ignored',
     >		' >>>        ',jmaxm,' radial points requested,',
     >		jpd,' radial points allowed.'
C----------------------------------------------------------------------|
	   return
	endif
!       idengrad = 2		! simple dilution, 3 
	idengrad = 3		! simple dilution, 3 
	i_delay  = 0
	itport_pt(1) = 1	! 1 particle transport on, 0 off
	itport_pt(2) = 1	! 1 electron heat transport on, 0 off
	itport_pt(3) = 1	! 1 ion heat transport on, 0 off
	itport_pt(4) = 0 ! 1/0/-1 v_phi transport on/off/use egamma_exp
	itport_pt(5) = 0 ! 1/0/-1 v_theta transport on/off/use gamma_p_exp
	irotstab     = 1 ! 1 use internally computed wExB, 0 for prescribed
	bt_exp       = BTOR
	bt_exp       = IPOL(1)*RTOR*BTOR/(RTOR+SHIF(1))
	bt_flag      = 1	! 0 do not use effective B-field
	rmajor_exp   = RTOR+SHIF(1)
	rmajor_exp   = RTOR
	amassgas_exp = AMJ
	zimp_exp     = ZIM1(1)
	if (AIM1 .gt. 1.d0) then
	   amassimp_exp = AIM1
	else
	   amassimp_exp = 1.
	endif
	arho_exp     = RHO(jna)
C       arho_exp     = AMETR(jna)
C       alpha_e      = 1.35	! 1/0 ExB shear stabilization on/off
	alpha_e      = 1.	! 1/0 ExB shear stabilization on/off
C       alpha_e      = 0. 	!     ExB shear stabilization off
	x_alpha      = 1.	! 1/0/-1 alpha stabilization on/off/self-cons

      	do j=1,jmaxm+1
	   rh_m(j-1) = (j-1.d0)/jmaxm 		! GLF grid 0->1
	   te_m(j-1) = TE(j)			! T_e, keV
	   ti_m(j-1) = TI(j)			! T_i, keV
	   ne_m(j-1) = NE(j)		! n_e, electron density, 10^19 m^-3
	   ns_m(j-1) = NIBM(j)		! n_f, fast ion density, 10^19 m^-3
	   zeff_exp(j-1) = ZEF(j)
	   if (zimp_exp .gt. 1.01) then
	      y1 = max(0.d0,(zimp_exp-ZEF(j))/(zimp_exp-1.))
	      ni_m(j-1) = max(0.d0,NE(j)*y1-ns_m(j-1))
	   else
	      ni_m(j-1) = NE(j)-ns_m(j-1)
	   endif
C	   ni_m(j-1) = ni_m(j-1)
	   elong_exp(j-1)= ELON(j)
	   rmin_exp(j-1) = AMETR(j)
	   rmaj_exp(j-1) = RTOR+SHIF(j)+AMETR(j)
	   q_exp(j-1)    = 1/MU(j)
	   include 'fml/almhd'
	   shat_exp(j-1)  = SHEAR(j)
	   alpha_exp(j-1) = ALMHD
	   gradrho_exp(j-1) = GRADRO(j)/arho_exp
!	   gradrho_exp(j-1) = 1./arho_exp
	   bteff_exp(j-1) = BTOR
	   bteff_exp(j-1) = IPOL(1)*RTOR*BTOR/(RTOR+SHIF(1))
	   gradrhosq_exp(j-1) = G11(j)/VRS(j)/arho_exp**2
!	   gradrhosq_exp(j-1) = 1./VRS(j+1)/arho_exp**2
C If (itport_pt(4) == 0  &&  itport_pt(5) == 0)	then
C         vphi_m(j)=rmajor_exp*angrotp_exp(j)
	   angrotp_exp(j-1) = VTOR(j)/RTOR
	   include 'fml/cs'
	   include 'fml/rotsh'
C egamma_exp and gamma_p_exp are used
C if(itport_pt(4) == -1) or (irotstab == 0) 
	   egamma_exp(j-1)  = ROTSH*ROC/(CS+1.d-3)
	   gamma_p_exp(j-1) = 0.0
C	   Y1 = SQEPS(j)
C	   Y2 = sqrt(1.+Y1*Y1)
C	   vpar_m(j-1) = VTOR(j)*Y2 ! calculated if itport_pt(5)=0
C	   vper_m(j-1) = VTOR(j)*Y1 ! calculated if itport_pt(5)=0
C abs(itport_pt(4)) == 1  &&  itport_pt(5) == 0
C this option: vpar is vphi, vexb from neo+vphi
	   vphi_m(j-1) = VTOR(j) ! calculated if itport_pt(4)*itport_pt(5)=0
	   vpar_m(j-1) = 0.0
	   vper_m(j-1) = 0.0
C	 if (j .eq. 14) write(*,'("inGLF",1P,6E12.3)')
C     >			alpha_exp(j-1),CS,ROTSH
	enddo
      	do j=1,jmaxm-1
	   drho=rh_m(j+1)-rh_m(j)+1.d-34
	   zgte_m(j) =-(dlog(te_m(j+1))-dlog(te_m(j)))/drho  ! -dln(T_e)/dr
	   zgti_m(j) =-(dlog(ti_m(j+1))-dlog(ti_m(j)))/drho  ! -dln(T_i)/dr
	   zgne_m(j) =-(dlog(ne_m(j+1))-dlog(ne_m(j)))/drho  ! -dln(n_e)/dr
	   zgni_m(j) =-(dlog(ni_m(j+1))-dlog(ni_m(j)))/drho  ! -dln(n_i)/dr
	   if (abs(zgte_m(j)) .lt. 1.d-6) zgte_m(j) = 1.d-6
	   if (abs(zgti_m(j)) .lt. 1.d-6) zgti_m(j) = 1.d-6
	   if (abs(zgne_m(j)) .lt. 1.d-6) zgne_m(j) = 1.d-6
	   if (abs(zgni_m(j)) .lt. 1.d-6) zgni_m(j) = 1.d-6
	enddo
	diffnem = 0
	chietem = 0
	chiitim = 0
	etaphim = 0
	etaparm = 0
	etaperm = 0
	exchm   = 0
	do j=0,jpd
	   diff_m(j)   = 0.0
	   chie_m(j)   = 0.0
	   chii_m(j)   = 0.0
	   etaphi_m(j) = 0.0
	   etapar_m(j) = 0.0
	   etaper_m(j) = 0.0
	   exch_m(j)   = 0.0
	   egamma_m(j) = 0.0
	   gamma_p_m(j)= 0.0
	   anrate_m(j) = 0.0
	   anrate2_m(j)= 0.0
	   anfreq_m(j) = 0.0
	   anfreq2_m(j)= 0.0
	   do k=1,10
	      egamma_d(j,k) = 0.0
	   enddo
	enddo
C----------------------------------------------------------------------|
C If (jmm == 0)   loop  1<j<jmaxm-1  is used inside CALLGLF2D 
C If (jmm != 0) indices j=jmm,jmm+1 are used inside CALLGLF2D 
c||--------------------------------------------------------------------|
C	write(*,*)jmaxm,jna-1
	IS = max(1,IS)
	IE = min(jmaxm-1,IE)

C	do  j=1,jmaxm-1
C	write(*,'(4(A,I3))')"GLF& loop:",IS," ->",IE,
C     &		';    Total =',jmaxm-1,',  Max =',jpd
	do  j=IS,IE
	   jmm = j
	   call callglf2d( leigen, nroot, iglf
     & ,     jshoot, jmm, jmaxm, itport_pt
     & ,     irotstab, te_m, ti_m, ne_m, ni_m, ns_m
     & ,     igrad, idengrad, zgte_m(j), zgti_m(j), zgne_m(j), zgni_m(j)
     & ,     angrotp_exp, egamma_exp, gamma_p_exp, vphi_m, vpar_m,vper_m
     & ,     zeff_exp, bt_exp, bt_flag, rh_m
     & ,     arho_exp, gradrho_exp, gradrhosq_exp
     & ,     rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp
     & ,     q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp
     & ,     alpha_e, x_alpha, i_delay
     & ,     diffnem, chietem, chiitim, etaphim, etaparm, etaperm
     & ,     exchm, diff_m, chie_m, chii_m, etaphi_m, etapar_m, etaper_m
     & ,     exch_m, egamma_m, egamma_d, gamma_p_m
     & ,     anrate_m, anrate2_m, anfreq_m, anfreq2_m )
C Note! In the local (point) mode callglf2d corrupts arrays as chie_m(j) 
C       or egamma_m(j) because on the entry it sets all arrays to 0
C       but then defines the single value j=jmm only.
	   CHI(j) = chiitim	! \chi_i, m^2/s
	   CHE(j) = chietem	! \chi_e, m^2/s
	   if (diffnem .ge. 0.)	then	! Outward flux -> D > 0, v = 0
	      DIF(j) = diffnem	! D, ion diffusivity, m^2/s
	      VIN(j) = 0.	!
	   else				! Inward flux -> D = 0, v /= 0
	      DIF(j) = 0.	! D, ion diffusivity, m^2/s
	      VIN(j) = diffnem*zgne_m(j)
	   endif
	   DPH(j) = etaphim	! toroidal velocity diffusivity m^2/s
	   DPL(j) = etaparm	! parallel velocity diffusivity m^2/s
	   DPR(j) = etaperm	! perpend. velocity diffusivity m^2/s
	   XTB(j) = exchm	! turbulent e-i equipartition in MW/m^3
	   EGM(j) = egamma_m(j)
	   GAM(j) = gamma_p_m(j) 
	   GM1(j) = anrate_m(j) ! leading mode rate [local csda_m}
	   GM2(j) = anrate2_m(j)! 2nd mode rate [local csda_m}
	   OM1(j) = anfreq_m(j) ! leading mode frequency
	   OM2(j) = anfreq2_m(j)! 2nd mode frequency
	   FR1(j) = 0.d0
	   if  (j .ge. NA1I)	then
	      CHI(j) = 0.
	      DPH(j) = 0.
	   endif
	   if  (j .ge. NA1E)	 CHE(j) = 0.
	   if  (j .ge. NA1N)	 then
	      DIF(j)  = 0.
	      VIN(j) = 0.
	   endif
	enddo			! End of main loop
	if (IE .lt. jmaxm-1)	return

	do j=jna-2,jna-1
	   CHI(j) = CHI(jna-3)
	   CHE(j) = CHE(jna-3)
	   DIF(j) = DIF(jna-3)
	   VIN(j) = VIN(jna-3)
	   DPH(j) = DPH(jna-3)
	   DPL(j) = DPL(jna-3)
	   DPR(j) = DPR(jna-3)
	   XTB(j) = XTB(jna-3)
	   EGM(j) = EGM(jna-3)
	   GAM(j) = GAM(jna-3)
	   GM1(j) = GM1(jna-3)
	   GM2(j) = GM2(jna-3)
	   OM1(j) = OM1(jna-3)
	   OM2(j) = OM2(jna-3)
	   FR1(j) = FR1(jna-3)
	enddo
C	write(*,*)"GLF23P:",jna,jmaxm,IE,DIF(82),DIF(83)
	if (jna-1 .lt. NA)	then
	   do j=jna-1,NA1
	      CHI(j) = 0.
	      CHE(j) = 0.
	      DIF(j) = 0.
	      VIN(j) = 0.
	      DPH(j) = 0.
	      DPL(j) = 0.
	      DPR(j) = 0.
	      XTB(j) = 0.
	      EGM(j) = 0.
	      GAM(j) = 0.
	      GM1(j) = 0.
	      GM2(j) = 0.
	      OM1(j) = 0.
	      OM2(j) = 0.
	      FR1(j) = 0.
	   enddo
	endif
	return
C----------------------------------------------------------------------|
	open(1,file='glf_inp')
 100	format(1P,(5E22.14))
	write(1,*)leigen, nroot, iglf
	write(1,*)jshoot, jmm, jmaxm, itport_pt
	write(1,*)irotstab, igrad, idengrad, i_delay
	write(1,100)(te_m(k),k=1,jmaxm)
	write(1,100)(ti_m(k),k=1,jmaxm)
	write(1,100)(ne_m(k),k=1,jmaxm)
	write(1,100)(ni_m(k),k=1,jmaxm)
	write(1,100)(ns_m(k),k=1,jmaxm)
	write(1,*)"gradients"
C	write(1,'(1I5,1P,1E25.14/)')(k,zgni_m(k),k=1,jmaxm)
	write(1,100)(zgte_m(k),k=1,jmaxm)
	write(1,100)(zgti_m(k),k=1,jmaxm)
	write(1,100)(zgne_m(k),k=1,jmaxm)
	write(1,100)(zgni_m(k),k=1,jmaxm)
	write(1,100)(angrotp_exp(k),k=1,jmaxm)
	write(1,100)(egamma_exp(k),k=1,jmaxm)
	write(1,100)(gamma_p_exp(k),k=1,jmaxm)
	write(1,100)(zeff_exp(k),k=1,jmaxm)
	write(1,*)bt_exp, bt_flag, arho_exp
	write(1,100)(rh_m(k),k=1,jmaxm)
	write(1,100)(gradrho_exp(k),k=1,jmaxm)
	write(1,100)(gradrhosq_exp(k),k=1,jmaxm)
	write(1,100)(rmin_exp(k),k=1,jmaxm)
	write(1,100)(rmaj_exp(k),k=1,jmaxm)
	write(1,*)rmajor_exp, zimp_exp, amassimp_exp
	write(1,100)(elong_exp(k),k=1,jmaxm)
	write(1,100)(q_exp(k),k=1,jmaxm)
	write(1,100)(shat_exp(k),k=1,jmaxm)
	write(1,100)(alpha_exp(k),k=1,jmaxm)
	write(1,*)amassgas_exp, alpha_e, x_alpha
	close(1)
	open(1,file='glf_out')
 200	format(1P,(5E22.10))
C	write(1,200)(vphi_m(k),k=1,jmaxm)
C	write(1,200)(vpar_m(k),k=1,jmaxm)
C	write(1,200)(vper_m(k),k=1,jmaxm)
	write(1,200)(chi(k),k=1,jmaxm-1)
	write(1,200)(che(k),k=1,jmaxm-1)
	write(1,200)(dif(k),k=1,jmaxm-1)
	write(1,200)(chi(k),k=IS,IE)
	write(1,200)(che(k),k=IS,IE)
	write(1,200)(dif(k),k=IS,IE)
	close(1)
 101	format(1P,5E13.3)
	end
C======================================================================|
