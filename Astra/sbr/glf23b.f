C======================================================================|
	subroutine GLF23B
C----------------------------------------------------------------------|
c based on stand-alone driver for the GLF23 model
c       "testglf.f" 18-fev-03 version 1.61
c       written by Jon Kinsey, General Atomics
C----------------------------------------------------------------------|
C WORK(1:NA1,1:13) array is used for output
C WORK1(1:NA1,1:25) array is used internally
C                              (when i_delay=0 and egamma_d is not used)
C----------------------------------------------------------------------|
C Further elaboration of the subroutine GLF161A in order to include
C transfer from transport grid to courser GLF grid
C					(Pereverzev 6-12-2005)
C Note! Simulation result is very sensitive to the parameter "a_sm"
C Namely: if the high gradient in NI propagates from the gradient
C	  pedestal zone inside the core plasma (in other words,
C	  an enhanced gradient grad(NI) appears) then all transport 
C	  coefficients strongly increase throughout the range of
C	  this increased gradient
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
C Parameters used:
C	integer			NA,NA1,NA1N,NA1E,NA1I
C	double precision	HRO,HROA,ROC,BTOR,RTOR,AMJ,AIM1,CBND1
C	double precision	NE(1),TE(1),NI(1),VTOR(1),TI(1)
C	double precision	NIBM(1),ZEF(1),IPOL(1),ZIM1(1)
C	double precision	RHO(1),AMETR(1),SHIF(1),ELON(1)
C	double precision	MU(1),SHEAR(1),G11(1),GRADRO(1),VRS(1)
C	double precision	WORK1(NRD,2*NRD+7),WORK(NRD,2*NRD),
C     1				WORK1D(NRD*(4*NRD+7))
C	common /A_WORKAR/ WORK1D
C	equivalence (WORK1D(1),WORK(1,1)),
C     1		    (WORK1D(NRD*(2*NRD+7)+1),WORK1(1,1))
	double precision ALMHD
C----------------------------------------------------------------------|
      integer jpd,jna,jcall
      data jcall/0/
      save jcall
      parameter ( jpd=100 ) 
      real*8 y1,y2,epsilon,x_glf(1:jpd+1),a_sm,FTAV
      external FTAV
      real*8 te_m(0:jpd), ti_m(0:jpd)
     & , ne_m(0:jpd), ni_m(0:jpd), ns_m(0:jpd)
     & , zpte_m(0:jpd), zpti_m(0:jpd), zpne_m(0:jpd), zpni_m(0:jpd)
     & , angrotp_exp(0:jpd), egamma_exp(0:jpd), gamma_p_exp(0:jpd)
     & , vphi_m(0:jpd), vpar_m(0:jpd), vper_m(0:jpd)
     & , zeff_exp(0:jpd), bt_exp, bteff_exp(0:jpd), rh_m(0:jpd),arho_exp
     & , gradrho_exp(0:jpd), gradrhosq_exp(0:jpd)
     & , rmin_exp(0:jpd), rmaj_exp(0:jpd), rmajor_exp
     & , q_exp(0:jpd), shat_exp(0:jpd), alpha_exp(0:jpd)
     & , elong_exp(0:jpd), zimp_exp, amassimp_exp, amassgas_exp
     & , alpha_e, x_alpha
      real*8 zpte_in, zpti_in, zpne_in, zpni_in, drho
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
c git: following input from ASTRA, not from file,
c take care that all variables are assigned
c      namelist /nlglf/ leigen, lprint, nroot, iglf, jshoot, jmm, jmaxm
c     & , itport_pt, irotstab, te_m, ti_m, ne_m, ni_m, ns_m
c     & , igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in
c     & , angrotp_exp, egamma_exp, gamma_p_exp, vphi_m, vpar_m, vper_m
c     & , zeff_exp, bt_exp, bt_flag, rh_m, arho_exp
c     & , gradrho_exp, gradrhosq_exp
c     & , rmin_exp, rmaj_exp, rmajor_exp, q_exp, shat_exp
c     & , alpha_exp, elong_exp, zimp_exp, amassimp_exp, amassgas_exp
c     & , alpha_e, x_alpha
c
C----------------------------------------------------------------------|
c git: modern block from ~/glf23_v1.60/testglf.f, some values rechosen
c or adapted to ASTRA needs
      epsilon  = 1.d-10
!      leigen   = 1   ! 1 for tomsqz, for cgg eigenvalue solver
      leigen   = 0   ! 1 for tomsqz, for cgg eigenvalue solver
!      nroot    = 8   ! n. of roots in eigenvalue solver (12 impurity dynamics)
      nroot    = 12   ! n. of roots in eigenvalue solver (12 impurity dynamics)
      iglf     = 1   ! 1 new, 0 original GLF23 normalization
      jshoot   = 0   ! for time-dependent code
c      jmm     > 0   ! grid number
      jmm      = 0   ! callglf does full grid
      jna      = max(NA1E,NA1I,NA1N)
      if (jna .eq. 0) jna = NA1
      igrad    = 0   ! 1 input gradients, 0 compute gradients
!      idengrad = 2   ! simple dilution, 3 
      idengrad = 3   ! simple dilution, 3 
      i_delay  = 0
      itport_pt(1) = 1   ! 1 particle transport on, 0 off
      itport_pt(2) = 1   ! 1 electron heat transport on, 0 off
      itport_pt(3) = 1   ! 1 ion heat transport on, 0 off
      itport_pt(4) = 0   ! 1/0/-1 v_phi transport on/off/use egamma_exp
      itport_pt(5) = 0   ! 1/0/- v_theta transport on/off/use gamma_p_exp
      irotstab     = 1   ! 1 use internally computed wExB, 0 for prescribed
      bt_exp       = BTOR
      bt_exp       = IPOL(1)*RTOR*BTOR/(RTOR+SHIF(1))
      bt_flag      = 1    ! 0 do not use effective B-field
c git: next line commented. Isn't bteff_exp an array?
c bteff_exp    = 1.0	! effective B-field (used when bt_flag > 0)
      rmajor_exp   = RTOR+SHIF(1)
      rmajor_exp   = RTOR
      amassgas_exp = AMJ
      zimp_exp     = ZIM1(1)
      amassimp_exp = AIM1
      arho_exp     = RHO(jna)
C      arho_exp     = AMETR(jna)
      alpha_e      = 1.35	! 1/0 ExB shear stabilization on/off
      alpha_e      = 0. 	!     ExB shear stabilization off
      x_alpha      = 1.		! 1/0/-1 alpha stabilization on/off/self-cons
      x_alpha      = 0.		!     alpha stabilization off
      zpte_in      = 0.		! input logaritmic gradients
      zpti_in      = 0.
      zpne_in      = 0.
      zpni_in      = 0.

      jmaxm = jna-1 ! WARNING! jmaxm cannot exceed jpd: jmaxm < jpd
c      jmaxm = 20
      if (jmaxm .ge. jpd)	then
	 write(*,*)' >>> GLF23: Grid allocation failure. Call ignored'
	 return
      endif
      do j=1,jna			! j=1,jna
	 work1(j,1) = RHO(j)/RHO(jna)	! Astra grid
      enddo
      do j=1,jmaxm+1
	 x_glf(j) = (j-1.d0)/jmaxm	! GLF grid
      enddo
C----------------------------------------------------------------------|
C Note! Simulation result is very sensitive to the parameter "a_sm"
C Namely: if the high gradient in NI propagates from the gradient
C	  pedestal zone inside the core plasma (in other words,
C	  an enhanced gradient grad(NI) appears) then all transport 
C	  coefficients strongly increase throughout the range of
C	  this increased gradient
C
	a_sm = 1.d-5
	if (jmaxm .eq. jna-1)	goto	1
C----------------------------------------------------------------------|
	call SMOOTH(a_sm,jna,TE,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   te_m(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,TI,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   ti_m(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,NE,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   ne_m(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,NIBM,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   ns_m(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,ZEF,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   zeff_exp(j-1) = work1(j,2)
	enddo
	if (zimp_exp .gt. 1.01) then
	   do j=1,jna
	      y1 = max(0.d0,(zimp_exp-ZEF(j))/(zimp_exp-1.))
	      work1(j,3) = max(0.d0,NE(j)*y1-NIBM(j))
	   enddo
	else
	   do j=1,jna
	      work1(j,3) = NE(j)-NIBM(j)
	   enddo
	endif
	call SMOOTH(a_sm,jna,work1(1,3),work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   ni_m(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,ELON,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   elong_exp(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,AMETR,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   rmin_exp(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,SHIF,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   rmaj_exp(j-1) = RTOR+work1(j,2)+rmin_exp(j-1)
	enddo
	call SMOOTH(a_sm,jna,MU,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   q_exp(j-1) = 1.d0/work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,SHEAR,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jna
	   include 'fml/almhd'
	   work1(j,4) = ALMHD
	enddo
	do j=1,jmaxm+1
	   shat_exp(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,work1(1,4),work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   alpha_exp(j-1) = work1(j,2)
	enddo
	call SMOOTH(a_sm,jna,GRADRO,work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   gradrho_exp(j-1) = work1(j,2)/arho_exp
	   bteff_exp(j-1) = BTOR
	enddo
	do j=1,jna
	   work1(j,3) = G11(j)/VRS(j)
	   work1(j,4) = VTOR(j)
	enddo
	call SMOOTH(a_sm,jna,work1(1,3),work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   gradrhosq_exp(j-1) = work1(j,2)/arho_exp**2
	enddo
	call SMOOTH(a_sm,jna,work1(1,4),work1,jmaxm+1,work1(1,2),x_glf)
	do j=1,jmaxm+1
	   angrotp_exp(j-1) = work1(j,2)/RTOR
	enddo
	goto	2
 1	continue
	do j=1,jmaxm+1
	   rh_m(j-1) = (j-1.d0)/jmaxm 		! GLF grid
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
	   angrotp_exp(j-1) = VTOR(j)/RTOR
	enddo
C zeffm*nem-nim-ns_m(jm)
C	write(61,*)(zeff_exp(j)*ne_m(j)-ni_m(j)-ns_m(j),j=1,jmaxm)

C	open(unit=61,file="./glf.out",status="unknown")
C	write(61,101)zimp_exp,BTOR,IPOL(1),SHIF(1),arho_exp
C	write(61,*)"rh_m"
C	write(61,101)(rh_m(j),j=1,jmaxm)
C	write(61,*)"te_m"
C	write(61,101)(te_m(j),j=1,jmaxm)
C	write(61,*)"ti_m"
C	write(61,101)(ti_m(j),j=1,jmaxm)
C	write(61,*)"ne_m"
C	write(61,101)(ne_m(j),j=1,jmaxm)
C	write(61,*)"ns_m"
C	write(61,101)(ns_m(j),j=1,jmaxm)
C	write(61,*)"zeff_exp"
C	write(61,101)(zeff_exp(j),j=1,jmaxm)
C	write(61,*)"ni_m"
C	write(61,101)(ni_m(j),j=1,jmaxm)
C	write(61,*)"elong"
C	write(61,101)(elong_exp(j),j=1,jmaxm)
C	write(61,*)"rmin"
C	write(61,101)(rmin_exp(j),j=1,jmaxm)
C	write(61,*)"rmaj"
C	write(61,101)(rmaj_exp(j),j=1,jmaxm)
C	write(61,*)"q"
C	write(61,101)(q_exp(j),j=1,jmaxm)
C	write(61,*)"shat"
C	write(61,101)(shat_exp(j),j=1,jmaxm)
C	write(61,*)"alpha"
C	write(61,101)(alpha_exp(j),j=1,jmaxm)
C	write(61,*)"gradrho"
C	write(61,101)(gradrho_exp(j),j=1,jmaxm)
C	write(61,*)"bteff"
C	write(61,101)(bteff_exp(j),j=1,jmaxm)
C	write(61,*)"gradrhosq"
C	write(61,101)(gradrhosq_exp(j),j=1,jmaxm)
C	write(61,*)"angrotp"
C	write(61,101)(angrotp_exp(j),j=1,jmaxm)
C	close(61)
 2	continue
      do j=1,jmaxm+1
	 egamma_exp(j-1)  = 0.0 ! is used if(itport_pt(4).eq.-1) only 
	 gamma_p_exp(j-1) = 0.0 ! is used if(itport_pt(4).eq.-1) only 
	 vphi_m(j-1)      = 0.0 ! calculated if itport_pt(4)*itport_pt(5)=0
	 vpar_m(j-1)      = 0.0 ! calculated if itport_pt(4)*itport_pt(5)=0
	 vper_m(j-1)      = 0.0 ! calculated if itport_pt(4)*itport_pt(5)=0
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
c||--------------------------------------------------------------------|
c In order to enable calling GLF at each grid point separately
c shift the next two lines inside the do loop. I.e. make
c     do  j=1,jmaxm-1
c        jmm = j
c        call callglf2d(...)
c        work(j,*) = ...
c        ...
c     enddo
c||--------------------------------------------------------------------|
c	jmm = j			! 
C If (jmm == 0)   loop  1<j<jmaxm-1  is used inside CALLGLF2D 
C If (jmm != 0) indices j=jmm,jmm+1 are used inside CALLGLF2D 
C----------------------------------------------------------------------|
        call callglf2d( leigen, nroot, iglf
     & , jshoot, jmm, jmaxm, itport_pt
     & , irotstab, te_m, ti_m, ne_m, ni_m, ns_m
     & , igrad, idengrad, zpte_in, zpti_in, zpne_in, zpni_in
     & , angrotp_exp, egamma_exp, gamma_p_exp, vphi_m, vpar_m, vper_m
     & , zeff_exp, bt_exp, bt_flag, rh_m
     & , arho_exp, gradrho_exp, gradrhosq_exp
     & , rmin_exp, rmaj_exp, rmajor_exp, zimp_exp, amassimp_exp
     & , q_exp, shat_exp, alpha_exp, elong_exp, amassgas_exp
     & , alpha_e, x_alpha, i_delay
     & , diffnem, chietem, chiitim, etaphim, etaparm, etaperm
     & , exchm, diff_m, chie_m, chii_m, etaphi_m, etapar_m, etaper_m
     & , exch_m, egamma_m, egamma_d, gamma_p_m
     & , anrate_m, anrate2_m, anfreq_m, anfreq2_m )

C----------------- Optionally, discard the edge points ----------------|
	chii_m(1)    = chii_m(2)
	chie_m(1)    = chie_m(2)   
	diff_m(1)    = diff_m(2)   
	etaphi_m(1)  = etaphi_m(2) 
	etapar_m(1)  = etapar_m(2) 
	etaper_m(1)  = etaper_m(2) 
	exch_m(1)    = exch_m(2)   
	egamma_m(1)  = egamma_m(2) 
	gamma_p_m(1) = gamma_p_m(2)
	anrate_m(1)  = anrate_m(2) 
	anrate2_m(1) = anrate2_m(2)
	anfreq_m(1)  = anfreq_m(2) 
	anfreq2_m(1) = anfreq2_m(2)
	chii_m(jmaxm)    = chii_m(jmaxm-1)
	chie_m(jmaxm)    = chie_m(jmaxm-1)   
	diff_m(jmaxm)    = diff_m(jmaxm-1)   
	etaphi_m(jmaxm)  = etaphi_m(jmaxm-1) 
	etapar_m(jmaxm)  = etapar_m(jmaxm-1) 
	etaper_m(jmaxm)  = etaper_m(jmaxm-1) 
	exch_m(jmaxm)    = exch_m(jmaxm-1)   
	egamma_m(jmaxm)  = egamma_m(jmaxm-1) 
	gamma_p_m(jmaxm) = gamma_p_m(jmaxm-1)
	anrate_m(jmaxm)  = anrate_m(jmaxm-1) 
	anrate2_m(jmaxm) = anrate2_m(jmaxm-1)
	anfreq_m(jmaxm)  = anfreq_m(jmaxm-1) 
	anfreq2_m(jmaxm) = anfreq2_m(jmaxm-1)
C----------------------------------------------------------------------|
      a_sm = 1.d-3
C      call SMOOTH(a_sm,jmaxm+1,chii_m,x_glf,jna,work(1,1),work1)
C      call SMOOTH(a_sm,jmaxm+1,chie_m,x_glf,jna,work(1,2),work1)
C      call SMOOTH(a_sm,jmaxm+1,diff_m,x_glf,jna,work(1,3),work1)
C      call SMOOTH(a_sm,jmaxm+1,etaphi_m,x_glf,jna,work(1,4),work1)
C      call SMOOTH(a_sm,jmaxm+1,etapar_m,x_glf,jna,work(1,5),work1)
C      call SMOOTH(a_sm,jmaxm+1,etaper_m,x_glf,jna,work(1,6),work1)
C      call SMOOTH(a_sm,jmaxm+1,exch_m,x_glf,jna,work(1,7),work1)
C      call SMOOTH(a_sm,jmaxm+1,egamma_m,x_glf,jna,work(1,8),work1)
C      call SMOOTH(a_sm,jmaxm+1,gamma_p_m,x_glf,jna,work(1,9),work1)
C      call SMOOTH(a_sm,jmaxm+1,anrate_m,x_glf,jna,work(1,10),work1)
C      call SMOOTH(a_sm,jmaxm+1,anrate2_m,x_glf,jna,work(1,11),work1)
C      call SMOOTH(a_sm,jmaxm+1,anfreq_m,x_glf,jna,work(1,12),work1)
C      call SMOOTH(a_sm,jmaxm+1,anfreq2_m,x_glf,jna,work(1,13),work1)
C      do j=jna,2,-1
C	 work(j,1) = work(j-1,1)
C	 work(j,2) = work(j-1,2)
C	 work(j,3) = work(j-1,3)
C	 work(j,4) = work(j-1,4)
C	 work(j,5) = work(j-1,5)
C	 work(j,6) = work(j-1,6)
C	 work(j,7) = work(j-1,7)
C	 work(j,8) = work(j-1,8)
C	 work(j,9) = work(j-1,9)
C	 work(j,10) = work(j-1,10)
C	 work(j,11) = work(j-1,11)
C	 work(j,12) = work(j-1,12)
C	 work(j,13) = work(j-1,13)
C      enddo

C--------------------------------------------------------------------
C      do j=1,jna
C	 work1(j,3) = work(j,3)
C	 if (RHO(j) .gt. CBND1*ROC)	work1(j,3) = 0.d0
C	 if (work1(j,3) .gt.  1.d1)	work1(j,3) = 1.d1
C	 if (work1(j,3) .lt. -1.d1)	work1(j,3) =-1.d1
C      enddo
C      call	FEVEN(32,work1(1,3),work1(1,23))
C      do j=1,jna
C	 work1(j,3) = FTAV(work1(j,23),CV14)
C	 work1(j,4) =-min(0.d0,work(j,3))	! Get negative part
C	 work1(j,5) = max(0.d0,work(j,3))	! Get positive part
C	 work1(j,6) = (NE(j+1)-NE(j))/HRO
C	 work1(j,6) = work1(j,6)/max(1.d-1,NE(j))
C	 work1(j,6) = min(-1.d-1,work1(j,6))
C      enddo
C      call SMOOTH(1.d0,jna,work1(1,4),work1,jna,work1(1,7),work1)
C      call SMOOTH(.5d0,jna,work1(1,6),work1,jna,work1(1,8),work1)
C      do j=1,jna
C	 work1(j,5) = work1(j,5)+CF6*work1(j,7)
C	 work1(j,9) = work1(j,4)+CF6*work1(j,7)
C	 work1(j,9) = max(-RHO(j),work1(j,9)*work1(1,8))
C      enddo
C      do j=1,jna
C	 work(j,33) = max(CV13,work1(j,5))		! -> DN
C	 work(j,30) = work1(j,9)			! -> CN
C      enddo
C--------------------------------------------------------------------
C	write(*,*)jmaxm,jna-1

      do j=1,jna
	 work(j,1) = chii_m(j-1)	! \chi_i, m^2/s
	 work(j,2) = chie_m(j-1)	! \chi_e, m^2/s
	 work(j,3) = diff_m(j-1)	! D, electron diffusivity, m^2/s
	 work(j,4) = etaphi_m(j-1)	!
	 work(j,5) = etapar_m(j-1)	!
	 work(j,6) = etaper_m(j-1)	!
	 work(j,7) = exch_m(j-1)	!
	 work(j,8) = egamma_m(j-1)	!
	 work(j,9) = gamma_p_m(j-1)	!
	 work(j,10) = anrate_m(j-1)	!
	 work(j,11) = anrate2_m(j-1)	!
	 work(j,12) = anfreq_m(j-1)	!
	 work(j,13) = anfreq2_m(j-1)	!
      enddo
c||--------------------------------------------------------------------|
      if (jna .lt. NA1)	then
	 do j=jna+1,NA1
	    work(j,1) = 0.
	    work(j,2) = 0.
	    work(j,3) = 0.
	    work(j,4) = 0.
	    work(j,5) = 0.
	    work(j,6) = 0.
	    work(j,7) = 0.
	    work(j,8) = 0.
	    work(j,9) = 0.
	    work(j,10) = 0.
	    work(j,11) = 0.
	    work(j,12) = 0.
	    work(j,13) = 0.
	 enddo
      endif
C----------------------------------------------------------------------|
	return
C     Added after 18-04-2007 to calculate d(chi)/d(T')
	do j=1,NA
	   y1 = 2.*RTOR*(TI(j+1)-TI(j))/(TI(j+1)+TI(j))/HRO
	   y2 = y1-work(j,32)
	   if (jcall.eq.1)	then
	      work(j,33) = (work(j,1)-work(j,31))/y2	!max(1.d-6,abs(y2))
	   endif
	   work(j,31) = work(j,1)
	   work(j,32) = y1
	enddo
	jcall = 1

	call	SHCOR1(work1(1,4))
	do j=1,NA1
	   if(RHO(j)/ROC .lt. CBND1)	then
	      y1 = work1(j,4)			! CAR30*CAR8
	   else
	      y1 = 0.
	   endif
	   work1(j,1) = max( 0.d0,min(1.d1,work(j,1)))*y1
	   work1(j,2) = max( 0.d0,min(1.d1,work(j,2)))*y1
	   work1(j,3) = max(-1.d1,min(1.d1,work(j,3)))*y1	! CUT=max(-X,min(X,Y))
	enddo
C	call	FEVEN(32,work1(1,1),work(1,21))
C	call	FEVEN(32,work1(1,2),work(1,22))
C	call	FEVEN(32,work1(1,3),work(1,23))
 101	format(1P,5E13.3)
	end
C======================================================================|
