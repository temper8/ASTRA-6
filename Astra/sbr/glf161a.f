C======================================================================|
	subroutine GLF161A
C----------------------------------------------------------------------|
c based on stand-alone driver for the GLF23 model
c       "testglf.f" 18-fev-03 version 1.61
c       written by Jon Kinsey, General Atomics
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c git: renamed rho in rh_m to avoid intereference with ASTRA
C  WORK(1:NA1,1:13) array is used for output
C WORK1(1:NA1,1:25) array is used internally
C                              (when i_delay=0 and egamma_d is not used)
C----------------------------------------------------------------------|
C Implementation of GLF23 model Version 1.61 is made by G.Tardini
C----------------------------------------------------------------------|
C This subroutine is similar to GLF161 
C The difference is in grids alignment: the edge point is ignored
C----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	double precision ALMHD
C----------------------------------------------------------------------|
      integer jpd,jna
      parameter ( jpd=100 )
 
      real*8 epsilon
c
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
c
      real*8 zpte_in, zpti_in, zpne_in, zpni_in, drho
c
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
      epsilon  = 1.e-10
      leigen   = 1   ! 1 for tomsqz, for cgg eigenvalue solver
      nroot    = 8   ! n. of roots in eigenvalue solver (12 impurity dynamics)
      iglf     = 1   ! 1 new, 0 original GLF23 normalization
      jshoot   = 0   ! for time-dependent code
c      jmm     > 0   ! grid number
      jmm      = 0   ! callglf does full grid
      jna      = NA1
      jna      = min(NA1E,NA1I,NA1N)
      if (jna .ne. max(NA1E,NA1I,NA1N)) jna = max(NA1E,NA1I,NA1N)
      if (jna .eq. 0) jna = NA1
      jmaxm    = jna-1 ! WARNING! jmaxm cannot exceed jpd: jmaxm < jpd
      if (jmaxm .ge. jpd)	then
	 write(*,*)' >>> GLF23: Grid allocation failure. Call ignored'
	 return
      endif
      igrad    = 0   ! 1 input gradients, 0 compute gradients
      idengrad = 2   ! simple dilution, 3 
      i_delay  = 0
      itport_pt(1) = 1   ! 1 particle transport on, 0 off
      itport_pt(2) = 1   ! 1 electron heat transport on, 0 off
      itport_pt(3) = 1   ! 1 ion heat transport on, 0 off
      itport_pt(4) = 0   ! 1/0/-1 v_phi transport on/off/use egamma_exp
      itport_pt(5) = 0   ! 1/0/- v_theta transport on/off/use gamma_p_exp
      irotstab     = 1   ! 1 use internally computed wExB, 0 for prescribed
      bt_exp       = BTOR
      bt_flag      = 1    ! 0 do not use effective B-field
c git: next line commented. Isn't bteff_exp an array?
c bteff_exp    = 1.0	! effective B-field (used when bt_flag > 0)
      rmajor_exp   = RTOR+SHIF(1)
      amassgas_exp = AMJ
      zimp_exp     = ZIM1(1)
      amassimp_exp = AIM1
      arho_exp     = RHO(jna)
      alpha_e      = 1.   ! 1/0 ExB shear stabilization on/off
      x_alpha      = 1.   ! 1/0/-1 alpha stabilization on/off/self-cons
      zpte_in      = 0. !input logaritmic gradients
      zpti_in      = 0.
      zpne_in      = 0.
      zpni_in      = 0.
c
      do j=1,jmaxm+1
	 rh_m(j-1)  = RHO(j)/arho_exp ! YRA
	 te_m(j-1) = TE(j)	! T_e, keV
	 ti_m(j-1) = TI(j)	! T_i, keV
	 ne_m(j-1) = NE(j)	! n_e, electron density, 10^19 m^-3
	 ni_m(j-1) = NI(j)	! n_i, ion density, 10^19 m^-3
	 ns_m(j-1) = NIBM(j)	! n_f, fast ion density, 10^19 m^-3
	 zeff_exp(j-1) = ZEF(j)
	 elong_exp(j-1) = ELON(j)
	 rmin_exp(j-1)  = AMETR(j)
	 rmaj_exp(j-1)  = RTOR+SHIF(j)
	 q_exp(j-1)     = 1/MU(j)
	 shat_exp(j-1)   = SHEAR(j)
         include 'fml/almhd'
	 alpha_exp(j-1)  = ALMHD
	 gradrho_exp(j-1) = GRADRO(j)/arho_exp
c git: VR instead of VRS
	 bteff_exp(j-1) = BTOR
	 gradrhosq_exp(j-1) = G11(j)/VR(j)/arho_exp**2
	 angrotp_exp(j-1) = VTOR(j)/RTOR
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
c----------------------------------------------------------------------|
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
      do j=1,jmaxm-1
	 work(j,1) = chii_m(j)		! \chi_i, m^2/s
	 work(j,2) = chie_m(j)		! \chi_e, m^2/s
	 work(j,3) = diff_m(j)		! D, ion diffusivity, m^2/s
	 work(j,4) = etaphi_m(j)	!
	 work(j,5) = etapar_m(j)	!
	 work(j,6) = etaper_m(j)	!
	 work(j,7) = exch_m(j)		!
	 work(j,8) = egamma_m(j)	!
	 work(j,9) = gamma_p_m(j)	!
	 work(j,10) = anrate_m(j)	!
	 work(j,11) = anrate2_m(j)	!
	 work(j,12) = anfreq_m(j)	!
	 work(j,13) = anfreq2_m(j)	!
      enddo
c||--------------------------------------------------------------------|
C       define the edge and the adjacent grid points as identical
C      jmaxm    = jna-1
      do j=jmaxm,jna
	 work(j,1) = work(jmaxm-1,1)
	 work(j,2) = work(jmaxm-1,2)
	 work(j,3) = work(jmaxm-1,3)
	 work(j,4) = work(jmaxm-1,4)
	 work(j,5) = work(jmaxm-1,5)
	 work(j,6) = work(jmaxm-1,6)
	 work(j,7) = work(jmaxm-1,7)
	 work(j,8) = work(jmaxm-1,8)
	 work(j,9) = work(jmaxm-1,9)
	 work(j,10) = work(jmaxm-1,10)
	 work(j,11) = work(jmaxm-1,11)
	 work(j,12) = work(jmaxm-1,12)
	 work(j,13) = work(jmaxm-1,13)
      enddo
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
      return
      end
C======================================================================|
