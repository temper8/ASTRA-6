c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                   WARNING!!!                                        !
c    NEOCL2 uses ASTRA-arrays WORK(j,101)...WORK(j,320)               !
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c=====================================================================|
      subroutine NEOCL4
c
c    Astra_6 interface to Houlberg's code NCLASS
c
c                                     May-2004, A.Zolotukhin
c---------------------------------------------------------------------|
cNCLASS calculates the neoclassical transport properties of a multiple
c  species axisymmetric plasma using k_order parallel and radial force
c  balance equations for each species
cReferences:
c  Houlberg, Shaing, Hirshman, Zarnstorff, Phys Plasmas 4 (1997) 3230
c  Hirshman, Sigmar, Nucl Fusion 21 (1981) 1079
c  W.A.Houlberg 3/99
c---------------------------------------------------------------------|
c	Usage example in an ASTRA model:
c		- for electron density
c		  DN = ... + "WORK(j,102)";
c		  CN = ... + "WORK(j,103)";
c
c		- for electron heat transport
c		  HE = ... + "WORK(j,105)";
c		  CE = ... + "WORK(j,106)";
c
c		- for a density of impurities Z1 
c                                (equation F1 , for example)
c		  DF1 = ... + "WORK(j,242)";
c		  VF1 = ... + "WORK(j,243)";
c
c		- for a bootstrap current
c		    DC=0.;  HC=0.;  XC=0.;
c                   CD=... + "WORK(j,301)";		
c---------------------------------------------------------------------|
c	Below the net fluxes are understood as a flux G_i of the 
c	particular species through the entire flux surface according to
c                 dn_i    d   
c                 ---- = ----[G_i] + Source_i
c                  dt     dV  
c	It coincides with the Astra total fluxes QN, QF1, etc., as
c---------------------------------------------------------------------|
c Output:
c		- electrons
c work(j,101) - particle radial net flux Gamma_electron, [1.e19/s]
c work(j,102) - particle diffusion coefficient Dn_electron, [m**2/s]
c work(j,103) - electron convective velocity Vn_electron, [m/s]
c work(j,104) - radial net heat conduction flux 
c               q_cond_electron, [MW] 
c work(j,105) - heat conductivity Chi_electron, [m**2/s]
c work(j,106) - heat convective velocity V_heat_electron, [m/s]
c work(j,107) - radial energy (conduction+convection) flux 
c               (electrons), [MW] 
c work(j,108) - bootstrap current on (p'/p)_electron, [MA/m**2]
c work(j,109) - bootstrap current on (T'/T)_electron, [MA/m**2]
c work(j,110) - poloidal flow velocity of electrons 
c               on outside midplane, [m/s]
c work(j,111...120) - reserved
c               ------------------------------
c		- main ions (not specified)
c++++++++++++++++  Note!!! ++++++++++++++++++
c This type of species is applied to satisfy quasineutrality condition
c and should be used only (!) in the case when ions H, D, T 
c  or He3 are not specified in the Astra model explicitly.
c  Otherwise use arrays for correspondent species.
c++++++++++++++++++++++++++++++++++++++++++++
c work(j,121) - particle radial net flux Gamma_mainions, [1.e19/s]
c work(j,122) - particle diffusion coefficient Dn_mainions, [m**2/s]
c work(j,123) - particle convective velocity Vn_mainions, [m/s]
c work(j,124) - radial net heat conduction flux 
c               q_cond_mainion, [MW] 
c work(j,125) - heat conductivity Chi_mainions, [m**2/s]
c work(j,126) - heat convective velocity V_heat_mainions, [m/s]
c work(j,127) - radial energy (conduction+convection) flux 
c               (main ions), [MW]
c work(j,128) - bootstrap current on (p'/p)_mainions, [MA/m**2]
c work(j,129) - bootstrap current on (T'/T)_mainions, [MA/m**2]
c work(j,130) - poloidal flow velocity of main ions 
c               on outside midplane, [m/s]
c work(j,131...140) - reserved
c               ------------------------------
c		- protons
c work(j,141) - particle radial net flux Gamma_protons, 
c               [1.e19/s]
c work(j,142) - particle diffusion coefficient Dn_proton, [m**2/s]
c work(j,143) - particle convective velocity Vn_proton, [m/s]
c work(j,144) - radial heat conduction net flux 
c               q_cond_proton, [MW] 
c work(j,145) - heat conductivity Chi_proton, [m**2/s]
c work(j,146) - heat convective velocity V_heat_proton, [m/s]
c work(j,147) - radial energy (conduction+convection) flux 
c               (protons), [MW]
c work(j,148) - bootstrap current on (p'/p)_proton, [MA/m**2]
c work(j,149) - bootstrap current on (T'/T)_proton, [MA/m**2]
c work(j,150) - poloidal flow velocity of protons 
c               on outside midplane, [m/s]
c work(j,151...160) - reserved
c               ------------------------------
c		- deuterons
c work(j,161) - particle radial net flux Gamma_deuterons,
c               [1.e19/s]
c work(j,162) - particle diffusion coefficient Dn_deuteron, [m**2/s]
c work(j,163) - particle convective velocity Vn_deuteron, [m/s]
c work(j,164) - radial heat conduction net flux 
c               q_cond_deuteron, [MW] 
c work(j,165) - heat conductivity Chi_deuteron, [m**2/s]
c work(j,166) - heat convective velocity V_heat_deuteron, [m/s]
c work(j,167) - radial energy (conduction+convection) flux 
c               (deuterons), [MW]
c work(j,168) - bootstrap current on (p'/p)_deuteron, [MA/m**2]
c work(j,169) - bootstrap current on (T'/T)_deuteron, [MA/m**2]
c work(j,170) - poloidal flow velocity of deuterons 
c               on outside midplane, [m/s]
c work(j,171...180) - reserved
c               ------------------------------
c		- tritons
c work(j,181) - particle radial net flux Gamma_tritons,
c               [1.e19/s]
c work(j,182) - particle diffusion coefficient Dn_triton, [m**2/s]
c work(j,183) - particle convective velocity Vn_triton, [m/s]
c work(j,184) - radial heat conduction net flux 
c               q_cond_triton, [MW] 
c work(j,185) - heat conductivity Chi_triton, [m**2/s]
c work(j,186) - heat convective velocity V_heat_triton, [m/s]
c work(j,187) - radial energy (conduction+convection) flux 
c               (tritons), [MW]
c work(j,188) - bootstrap current on (p'/p)_triton, [MA/m**2]
c work(j,189) - bootstrap current on (T'/T)_triton, [MA/m**2]
c work(j,190) - poloidal flow velocity of tritons 
c               on outside midplane, [m/s]
c work(j,191...200) - reserved
c               ------------------------------
c		- He3-particles
c work(j,201) - particle radial net flux Gamma_He3, [1.e19/s]
c work(j,202) - particle diffusion coefficient Dn_He3, [m**2/s]
c work(j,203) - particle convective velocity Vn_He3, [m/s]
c work(j,204) - radial heat conduction net flux 
c               q_cond_He3, [MW] 
c work(j,205) - heat conductivity Chi_He3, [m**2/s]
c work(j,206) - heat convective velocity V_heat_He3, [m/s]
c work(j,207) - radial energy (conduction+convection) flux 
c               (He3's), [MW]
c work(j,208) - bootstrap current on (p'/p)_He3, [MA/m**2]
c work(j,209) - bootstrap current on (T'/T)_He3, [MA/m**2]
c work(j,210) - poloidal flow velocity of He3-particles 
c               on outside midplane, [m/s]
c work(j,211...220) - reserved
c               ------------------------------
c		- alpha-particles
c work(j,221) - particle radial net flux Gamma_alpha, [1.e19/s]
c work(j,222) - particle diffusion coefficient Dn_alpha, [m**2/s]
c work(j,223) - particle convective velocity Vn_alpha, [m/s]
c work(j,224) - radial heat conduction net flux 
c               q_cond_alpha, [MW] 
c work(j,225) - heat conductivity Chi_alpha, [m**2/s]
c work(j,226) - heat convective velocity V_heat_alpha, [m/s]
c work(j,227) - radial energy (conduction+convection) flux 
c               (alphas), [MW]
c work(j,228) - bootstrap current on (p'/p)_alpha, [MA/m**2]
c work(j,229) - bootstrap current on (T'/T)_alpha, [MA/m**2]
c work(j,230) - poloidal flow velocity of alpha-particles 
c               on outside midplane, [m/s]
c work(j,231...240) - reserved
c               ------------------------------
c		- Impurity Z_1
c work(j,241) - particle radial net flux Gamma_impZ1, [1.e19/s]
c work(j,242) - particle diffusion coefficient Dn_impZ1, [m**2/s]
c work(j,243) - particle convective velocity Vn_impZ1, [m/s]
c work(j,244) - radial heat conduction net flux 
c               q_cond_impZ1, [MW]
c work(j,245) - heat conductivity Chi_impZ1, [m**2/s]
c work(j,246) - heat convective velocity V_heat_impZ1, [m/s]
c work(j,247) - radial energy (conduction+convection) flux 
c               (impurity Z1), [MW]
c work(j,248) - bootstrap current on (p'/p)_impZ1, [MA/m**2]
c work(j,249) - bootstrap current on (T'/T)_impZ1, [MA/m**2]
c work(j,250) - poloidal flow velocity of impurity ions Z1 
c               on outside midplane, [m/s]
c work(j,251...260) - reserved
c               ------------------------------
c		- Impurity Z_2
c work(j,261) - particle radial net flux Gamma_impZ2, [1.e19/s]
c work(j,262) - particle diffusion coefficient Dn_impZ2, [m**2/s]
c work(j,263) - particle convective velocity Vn_impZ2, [m/s]
c work(j,264) - radial heat conduction net flux 
c               q_cond_impZ2, [MW]
c work(j,265) - heat conductivity Chi_impZ2, [m**2/s]
c work(j,266) - heat convective velocity V_heat_impZ2, [m/s]
c work(j,267) - radial energy (conduction+convection) flux 
c               (impurity Z2), [MW]
c work(j,268) - bootstrap current on (p'/p)_impZ2, [MA/m**2]
c work(j,269) - bootstrap current on (T'/T)_impZ2, [MA/m**2]
c work(j,270) - poloidal flow velocity of impurity ions Z2 
c               on outside midplane, [m/s]
c work(j,271...280) - reserved
c               ------------------------------
c		- Impurity Z_3
c work(j,281) - particle radial net flux Gamma_impZ3, [1.e19/s]
c work(j,282) - particle diffusion coefficient Dn_impZ3, [m**2/s]
c work(j,283) - particle convective velocity Vn_impZ3, [m/s]
c work(j,284) - radial heat conduction net flux 
c               q_cond_impZ3, [MW]
c work(j,285) - heat conductivity Chi_impZ3, [m**2/s]
c work(j,286) - heat convective velocity V_heat_impZ3, [m/s]
c work(j,287) - radial energy (conduction+convection) flux 
c               (impurity Z3), [MW]
c work(j,288) - bootstrap current on (p'/p)_impZ3, [MA/m**2]
c work(j,289) - bootstrap current on (T'/T)_impZ3, [MA/m**2]
c work(j,290) - poloidal flow velocity of impurity ions Z3 
c               on outside midplane, [m/s]
c work(j,291...300) - reserved
c               ------------------------------
c               - Miscellanious parameters
c work(j,301) - bootstrap current, [MA/m**2]
c work(j,302) - external current, [MA/m**2]
c work(j,303) - current conductivity, [MS/m=1/(microOhm*m)]
c work(j,304) - density of main ions as they are defined above
c work(j,305...320) - reserved
c---------------------------------------------------------------------|
      implicit none
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      include 'sbr/pamx_mi.inc'
      include 'sbr/pamx_ms.inc'
      include 'sbr/pamx_mz.inc'
      double precision  YGRRdB2(1000),YNGRTHETA(1000),YFM(3,1000)
      equivalence	(WORK1(1,1),   YGRRdB2(1)),
     1			(WORK1(1,10),YNGRTHETA(1)),
     2			(WORK1(1,20),  YFM(1,1))
      double precision          y_grrho2,   grti,  y_den,  yh
      real          dq_s(mx_ms),   vq_s(mx_ms)
      real          p_eps
      integer k_electron,	k_mainion,	k_proton,
     #        k_deuteron,	k_triton,	k_he3,
     #        k_alpha,		k_impZ1,	k_impZ2,
     #        k_impZ3,		narray
      real        ybbmax, ybbmax2, yftupper, yftlower
      real        z_coulomb,               z_electronmass,
     #               z_j7kv,                  z_mu0,
     #               z_pi,                    z_protonmass
      real           rdum(8),                 RARRAY_SUM
cDeclaration of input to NCLASS
      integer        k_order,                 k_potato
      integer        m_i,                     m_z
      real           c_den,                   c_potb,
     #               c_potl
      real           p_b2,                    p_bm2,
     #               p_eb,                    p_fhat,
     #               p_fm(3),                 p_ft,
     #               p_grbm2,                 p_grphi,
     #               p_gr2phi,                p_ngrth
      real           amu_i(mx_mi),            grt_i(mx_mi),
     #               temp_i(mx_mi)
      real           den_iz(mx_mi,mx_mz),     fex_iz(3,mx_mi,mx_mz),
     #               grp_iz(mx_mi,mx_mz)
      real           YNMAIN,                  YNMAIN1
cDeclaration of output from NCLASS
      integer        iflag,                   m_s
      integer        jm_s(mx_ms),             jz_s(mx_ms)
      real           p_bsjb,                  p_etap,
     #               p_exjb
      real           calm_i(3,3,mx_mi)
      real           caln_ii(3,3,mx_mi,mx_mi),
     #               capm_ii(3,3,mx_mi,mx_mi),
     #               capn_ii(3,3,mx_mi,mx_mi)
      real           bsjbp_s(mx_ms),          bsjbt_s(mx_ms),
     #               dn_s(mx_ms),             gfl_s(5,mx_ms),
     #               qfl_s(5,mx_ms),          sqz_s(mx_ms),
     #               upar_s(3,3,mx_ms),       utheta_s(3,3,mx_ms),
     #               vn_s(mx_ms),             veb_s(mx_ms),
     #               qeb_s(mx_ms),            xi_s(mx_ms),
     #               ymu_s(3,3,mx_ms)
      real	     chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms),
     #               dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
cDeclaration of local variables
      character      label*120
      integer        i,		iza,	j,    jj,	j1,	k
      integer        idum(8),   nout
C----------------------------------------------------------------------|
C used for output only
 	integer	im
C----------------------------------------------------------------------|
      if (mx_mi .lt. 6)	then
         write(*,*)"NEOCL warning: call ignored. mx_mi < 6"
         return
      endif
c  Physical and conversion constants
      z_coulomb=1.6022e-19
      z_electronmass=9.1095e-31
      z_j7kv=1.6022e-16
      z_mu0=1.2566e-06
      z_pi=ACOS(-1.0)
      z_protonmass=1.6726e-27
C Set control and radially independent data
c  k_order-order of v moments to be solved [-]
c         =2 u and p_q
c         =3 u, p_q, and u2
c         =else error
c  k_potato-option to include potato orbits [-]
c          =0 off
c          =else on
      k_order    = 2
      k_potato   = 1
c  c_den-density cutoff below which species is ignored (/m**3)
c  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
c  c_potl-q(0)*R(0) (m)
      c_den      = 1.0e10
      y_den      = 1.e-19*c_den                  ! c_den in ASTRA units
      c_potb     = -0.5*ELON(1)*BTOR*MU(1)**2
C Parameters:
C  mx_mi=9  - max number of isotopes
C  mx_mz=18 - max charge
C  mx_ms=40 - max number of species
c  m_i-number of isotopes (1<mi<mx_mi+1)
c  m_z-highest charge state of all species (0<mz<mx_mz+1)
c  grt_i(i)-temperature gradient of i (keV/rho)
c  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)

c========  Initialization of all output arrays  =======================|
	do j=1,NA1
		do narray=101,320
		 work(j,narray) = 0.
		enddo
	enddo
c======================================================================|
	CALL ZBFAUX(yGRRdB2,yNGRTheta,YFM)
c----------------------------------------------------------------------|
	do j=1,NA
C Set radially dependent data
c  p_eps-inverse aspect ratio [-]
c  p_grphi-radial electric field Phi' (V/rho)
c  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
c  p_q-safety factor [-] (-1./MU(j))
c  p_eb-<E.B> (V*T/m)
c  p_b2-<B**2> (T**2)
c  p_bm2-<1/B**2> (/T**2)
c  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
c  p_ft-trapped fraction [-]
c  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
c  p_ngrth-<n.grad(Theta)> (1/m)
c  c_potl   = -(RTOR+SHIF(j))/MU(1)

	k_electron = 0
	k_mainion  = 0
	k_proton   = 0
	k_deuteron = 0
	k_triton   = 0
	k_he3      = 0
	k_alpha    = 0
	k_impZ1    = 0
	k_impZ2    = 0
	k_impZ3    = 0

	 if (j .lt. NA)	then
            YH = HRO
         else
            YH = HROA
         endif
         c_potl   = -(RTOR+SHIF(j))/MU(1)

 	  y_grrho2 = G11(j)/VRS(j)

         p_eps    = HRO*j/ROC*ABC/(RTOR+SHIF(j))
         p_grphi  = ER(j)
         p_gr2phi = MU(j)*j*
     >              (ER(j+1)/MU(j+1)/((j+1)*YH) - ER(j)/MU(j)/(j*YH))
C Warning: The NCLASS version 1.2 returns NaN resistivity
C          if p_eb = 0.
C         p_eb     = max(.00001,-ULON(j)*BTOR/(GP2*RTOR))
         p_eb     = -ULON(j)*BTOR/(GP2*RTOR)
         if (abs(p_eb) .lt. 1.d-3)	p_eb = sign(0.001,p_eb)

          p_b2     = BTOR**2*BDB02(j)
          p_bm2    = B0DB2(j)/BTOR**2
          p_fhat	 = -IPOL(j)*RTOR/(MU(j)*j*YH)
c Trapped particle fraction
          ybbmax   = BTOR/BMAXT(j)*BDB0(j)      !->  <h>=<B/B_max>
 	 ybbmax2  = (BTOR/BMAXT(j))**2*BDB02(j)!->  <h2>=<(B/B_max)2>
 	 yftupper = 1.- ybbmax2/ybbmax**2*
     >                    (1.- sqrt(1. - ybbmax)*(1. + 0.5*ybbmax))
 	 yftlower = 1. - ybbmax2*(BMAXT(j)/BTOR)**2*FOFB(j)
 	 p_ft     = 0.75*yftupper + 0.25*yftlower

!         p_grbm2 = <grad(rho)**2/B**2> ~

	 p_grbm2 =  yGRRdB2(j)

!         p_ngrth - <n.grad(Theta)>     ~
	 p_ngrth =  yNGRTheta(j)

c  p_fm(3)-poloidal moments of geometric factor for PS viscosity [-]
         do i=1,3
            p_fm(i)=0.0
         enddo
         do i=1,3
c            p_fm(i)=i*((1.0-SQRT(1.0-p_eps**2))/p_eps)**(2.0*i)
c     #           *(1.0+i*SQRT(1.0-p_eps**2))/((1.0-p_eps**2)**1.5
c     #           )*p_ngrth**2
	    p_fm(i) = YFM(i,j)
        enddo


         do j1=1,mx_mi
            do jj=1,mx_mz
               den_iz(j1,jj) = 0.
               grp_iz(j1,jj) = 0.
               do i=1,3
                  fex_iz(i,j1,jj) = 0.
               enddo
            enddo
         enddo

 	 m_i = 1                ! Reserved for electrons
	 k_electron = 1
         m_z = 1
         amu_i(m_i)  = 5.4463e-4
         temp_i(m_i) = TE(j)
         grt_i(m_i)  = (TE(j+1)-TE(j))/YH
         den_iz(m_i,m_z) = 1.e19*NE(j)
         grp_iz(m_i,m_z) = 1.e19/YH*(TE(j+1)*NE(j+1)-TE(j)*NE(j)) 

         grti = (TI(j+1)-TI(j))/YH

	if (NHYDR(j) .gt. y_den)	then
	    m_i = m_i+1
	    k_proton = m_i
            m_z = 1
            amu_i(m_i) = 1.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            den_iz(m_i,m_z)=1.e19*NHYDR(j)
            grp_iz(m_i,m_z) = 1.e19/YH*
     >            (TI(j+1)*NHYDR(j+1)-TI(j)*NHYDR(j))
         endif


	 if (NDEUT(j) .gt. y_den)	then
            m_i = m_i+1
	    k_deuteron = m_i
            amu_i(m_i) = 2.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            den_iz(m_i,m_z) = 1.e19*NDEUT(j)
            grp_iz(m_i,m_z) = 1.e19/YH*
     >                        (TI(j+1)*NDEUT(j+1)-TI(j)*NDEUT(j)) 
         endif
         if (NTRIT(j) .gt. y_den)	then
            m_i = m_i+1
	    k_triton = m_i
            amu_i(m_i) = 3.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            den_iz(m_i,m_z) = 1.e19*NTRIT(j)
            grp_iz(m_i,m_z) = 1.e19/YH*
     >                    (TI(j+1)*NTRIT(j+1)-TI(j)*NTRIT(j))
         endif

         if (m_i .eq. 1 .and. mx_mz .lt. 2 .and. ZMJ .ge. 2 )	then
            write(*,*)">>> NEOCL: ion composition error:"
            write(*,*)" mx_mz < 2 and no hydrogen. Call ignored."
            return
         endif
c It is assumed here that the ion species are not specified explicitly
c in the model
c or/and
c is used to satisfy the quasineutrality condition
	    YNMAIN = 1./ZMJ * ( NE(j) 
     >                           - NHYDR(j) - NDEUT(j) - NTRIT(j)
     >                           - 2.*NHE3(j) - 2.*NALF(j)
     >                           - ZIM1(j)*NIZ1(j)
     >                           - ZIM2(j)*NIZ2(j)
     >                           - ZIM3(j)*NIZ3(j) )
            work(j,304) = YNMAIN
	    YNMAIN1 = 1./ZMJ * ( NE(j+1) 
     >                           - NHYDR(j+1) - NDEUT(j+1) - NTRIT(j+1)
     >                           - 2.*NHE3(j+1) - 2.*NALF(j+1)
     >                           - ZIM1(j+1)*NIZ1(j+1)
     >                           - ZIM2(j+1)*NIZ2(j+1)
     >                           - ZIM3(j+1)*NIZ3(j+1) )
          if (YNMAIN .gt. y_den)        then
	    m_i = m_i + 1
	    k_mainion = m_i
	    m_z = ZMJ
            amu_i(m_i) = AMJ
            temp_i(m_i) = TI(j)
            grt_i(m_i) = grti
            den_iz(m_i,m_z) = 1.e19*YNMAIN
	    grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*YNMAIN1-TI(j)*YNMAIN) 
	  endif
         if (NHE3(j) .gt. y_den)	then
            m_i = m_i+1
	    k_he3 = m_i
	    amu_i(m_i) = 3.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            m_z = 2
            den_iz(m_i,m_z) = 1.e19*NHE3(j)
            grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*NHE3(j+1)-TI(j)*NHE3(j)) 
         endif
         
	 if (m_i .eq. mx_mi)	goto	5

         if (NALF(j) .gt. y_den)	then
            m_i = m_i+1
	    k_alpha = m_i
            amu_i(m_i) = 4.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            m_z = 2
            den_iz(m_i,m_z) = 1.e19*NALF(j)
            grp_iz(m_i,m_z) = 1.e19/YH*(TI(j+1)*NALF(j+1)-TI(j)*NALF(j)) 
         endif

c consider up to three impurity species
         jj = ZIM1(j)+0.5
         if (jj .gt. mx_mz)	goto	3
         if (NIZ1(j) .gt. y_den)	then
            m_i = m_i+1
	    k_impZ1 = m_i
            amu_i(m_i) = AIM1
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj .gt. m_z)	m_z = jj
            den_iz(m_i,jj) = 1.e19*NIZ1(j)
            grp_iz(m_i,jj) = 1.e19/YH*(TI(j+1)*NIZ1(j+1)-TI(j)*NIZ1(j)) 
         endif
         if (m_i .eq. mx_mi)	goto	5
 3       jj = ZIM2(j)+0.5
         if (jj .gt. mx_mz)	goto	4
         if (NIZ2(j) .gt. y_den)	then
            m_i = m_i+1
	    k_impZ2 = m_i
            amu_i(m_i) = AIM2
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj .gt. m_z)	m_z = jj
            den_iz(m_i,jj) = 1.e19*NIZ2(j)
            grp_iz(m_i,jj) = 1.e19/YH*(TI(j+1)*NIZ2(j+1)-TI(j)*NIZ2(j)) 
         endif
         if (m_i .eq. mx_mi)	goto	5
 4       jj = ZIM3(j)+0.5
         if (jj .gt. mx_mz)	goto	5
         if (NIZ3(j) .gt. y_den)	then
            m_i = m_i+1
	    k_impZ3 = m_i
            amu_i(m_i) = AIM3
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj .gt. m_z)	m_z = jj
            den_iz(m_i,jj) = 1.e19*NIZ3(j)
            grp_iz(m_i,jj) = 1.e19/YH*(TI(j+1)*NIZ3(j+1)-TI(j)*NIZ3(j)) 
         endif
	        
 5       continue

         CALL NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,
     #        p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi,
     #        p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,
c output:
     #        m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,capm_ii,
     #        capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s,
     #        utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,chit_ss,
     #        dp_ss,dt_ss,iflag)


         IF (iflag.gt.0) THEN
            nout=6
            OPEN(unit=nout,
     >           file='user')	    
c     >          ,file='out_nclass_pt.dat'
c     >           ,status='unknown'
c     >            )
            IF(iflag.eq.1) THEN
               label='ERROR:NCLASS-k_order must be 2 or 3, k_order='
               idum(1)=k_order
               CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
            ELSEIF(iflag.eq.2) THEN
               label='ERROR:NCLASS-require 1<m_i<mx_mi, m_i='
               idum(1)=m_i
               CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
            ELSEIF(iflag.eq.3) THEN
               label='ERROR:NCLASS-require 0<m_z<mx_mz, m_z='
               idum(1)=m_z
               CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
            ELSEIF(iflag.eq.4) THEN
               label='ERROR:NCLASS-require 0<m_s<mx_ms, m_s='
               idum(1)=m_s
               CALL WRITE_LINE_IR(nout,label,1,idum,0,rdum,0)
            ELSEIF(iflag.eq.5) THEN
               label='ERROR:NCLASS-inversion of flow matrix failed'
               CALL WRITE_LINE(nout,label,0,0)
            ENDIF
            GOTO 1000
         ENDIF
c=====================================================================|
c================   Output  ==========================================|
c=====================================================================|
c=================  Electrons  =======================================|
	 	IF (k_electron .gt. 0) then
c--------------> Total radial particle flux (electrons)
		  CALL RARRAY_COPY(5,gfl_s(1,k_electron),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,101) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (electrons)
		  im  = jm_s(k_electron)
		  iza = IABS(jz_s(k_electron))
		  work(j,102) = dn_s(k_electron)/y_grrho2
		  work(j,103) = (vn_s(k_electron)  +
     >                           veb_s(k_electron) +
     >                           gfl_s(5,k_electron)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (electrons)
		  CALL RARRAY_COPY(5,qfl_s(1,k_electron),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,104) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (electrons)
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_electron),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_electron)=chit_ss(k_electron,k_electron) +
     >                              chip_ss(k_electron,k_electron)
                   vq_s(k_electron)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                           dq_s(k_electron)*grt_i(im)/temp_i(im)
                   work(j,105) = dq_s(k_electron)/y_grrho2
                   work(j,106) = vq_s(k_electron)/y_grrho2
c--------------> Total radial energy flux (electrons)
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_electron) +
     >                     2.5*gfl_s(k,k_electron)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,107) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (electrons)
                  rdum(1)=bsjbp_s(k_electron)
                  rdum(2)=bsjbt_s(k_electron)
	          work(j,108) = -1.e-6*rdum(1)/BTOR
	          work(j,109) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (electrons)
                  work(j,110) = ( utheta_s(1,1,k_electron) + 
     >                            utheta_s(1,2,k_electron) +
     >                            utheta_s(1,3,k_electron)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Main Ions  =======================================|
	 	IF (k_mainion .gt. 0) then
c--------------> Total radial particle flux (main ions)
		  CALL RARRAY_COPY(5,gfl_s(1,k_mainion),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,121) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (main ions)
		  im  = jm_s(k_mainion)
		  iza = IABS(jz_s(k_mainion))
		  work(j,122) = dn_s(k_mainion)/y_grrho2
		  work(j,123) = (vn_s(k_mainion)  +
     >                           veb_s(k_mainion) +
     >                           gfl_s(5,k_mainion)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (main ions)
		  CALL RARRAY_COPY(5,qfl_s(1,k_mainion),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,124) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (main ions)
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_mainion),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_mainion)=chit_ss(k_mainion,k_mainion) +
     >                              chip_ss(k_mainion,k_mainion)
                   vq_s(k_mainion)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                            dq_s(k_mainion)*grt_i(im)/temp_i(im)
                   work(j,125) = dq_s(k_mainion)/y_grrho2
                   work(j,126) = vq_s(k_mainion)/y_grrho2
c--------------> Total radial energy flux (main ions)
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_mainion) +
     >                     2.5*gfl_s(k,k_mainion)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,127) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (main ions)
                  rdum(1)=bsjbp_s(k_mainion)
                  rdum(2)=bsjbt_s(k_mainion)
	          work(j,128) = -1.e-6*rdum(1)/BTOR
	          work(j,129) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (main ions)
                  work(j,130) = ( utheta_s(1,1,k_mainion) + 
     >                            utheta_s(1,2,k_mainion) +
     >                            utheta_s(1,3,k_mainion)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Protons    =======================================|
	 	IF (k_proton .gt. 0) then
c--------------> Total radial particle flux (protons  )
		  CALL RARRAY_COPY(5,gfl_s(1,k_proton),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,141) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (protons  )
		  im  = jm_s(k_proton)
		  iza = IABS(jz_s(k_proton))
		  work(j,142) = dn_s(k_proton)/y_grrho2
		  work(j,143) = (vn_s(k_proton)  +
     >                           veb_s(k_proton) +
     >                           gfl_s(5,k_proton)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (protons  )
		  CALL RARRAY_COPY(5,qfl_s(1,k_proton),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,144) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (protons  )
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_proton),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_proton)=chit_ss(k_proton,k_proton) +
     >                              chip_ss(k_proton,k_proton)
                   vq_s(k_proton)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                             dq_s(k_proton)*grt_i(im)/temp_i(im)
                   work(j,145) = dq_s(k_proton)/y_grrho2
                   work(j,146) = vq_s(k_proton)/y_grrho2
c--------------> Total radial energy flux (protons  )
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_proton) +
     >                     2.5*gfl_s(k,k_proton)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,147) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (protons  )
                  rdum(1)=bsjbp_s(k_proton)
                  rdum(2)=bsjbt_s(k_proton)
	          work(j,148) = -1.e-6*rdum(1)/BTOR
	          work(j,149) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (protons)
                  work(j,150) = ( utheta_s(1,1,k_proton) + 
     >                            utheta_s(1,2,k_proton) +
     >                            utheta_s(1,3,k_proton)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Deuterons  =======================================|
	 	IF (k_deuteron .gt. 0) then
c--------------> Total radial particle flux (deuterons  )
		  CALL RARRAY_COPY(5,gfl_s(1,k_deuteron),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,161) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (deuterons  )
		  im  = jm_s(k_deuteron)
		  iza = IABS(jz_s(k_deuteron))
		  work(j,162) = dn_s(k_deuteron)/y_grrho2
		  work(j,163) = (vn_s(k_deuteron)  +
     >                           veb_s(k_deuteron) +
     >                           gfl_s(5,k_deuteron)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (deuterons  )
		  CALL RARRAY_COPY(5,qfl_s(1,k_deuteron),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,164) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (deuterons  )
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_deuteron),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_deuteron)=chit_ss(k_deuteron,k_deuteron) +
     >                              chip_ss(k_deuteron,k_deuteron)
                   vq_s(k_deuteron)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                          dq_s(k_deuteron)*grt_i(im)/temp_i(im)
                   work(j,165) = dq_s(k_deuteron)/y_grrho2
                   work(j,166) = vq_s(k_deuteron)/y_grrho2
c--------------> Total radial energy flux (deuterons  )
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_deuteron) +
     >                     2.5*gfl_s(k,k_deuteron)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,167) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (deuterons  )
                  rdum(1)=bsjbp_s(k_deuteron)
                  rdum(2)=bsjbt_s(k_deuteron)
	          work(j,168) = -1.e-6*rdum(1)/BTOR
	          work(j,169) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (deuterons)
                  work(j,170) = ( utheta_s(1,1,k_deuteron) + 
     >                            utheta_s(1,2,k_deuteron) +
     >                            utheta_s(1,3,k_deuteron)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Tritons    =======================================|
	 	IF (k_triton .gt. 0) then
c--------------> Total radial particle flux (tritons  )
		  CALL RARRAY_COPY(5,gfl_s(1,k_triton),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,181) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (tritons  )
		  im  = jm_s(k_triton)
		  iza = IABS(jz_s(k_triton))
		  work(j,182) = dn_s(k_triton)/y_grrho2
		  work(j,183) = (vn_s(k_triton)  +
     >                           veb_s(k_triton) +
     >                           gfl_s(5,k_triton)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (tritons  )
		  CALL RARRAY_COPY(5,qfl_s(1,k_triton),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,184) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (tritons  )
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_triton),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_triton)=chit_ss(k_triton,k_triton) +
     >                              chip_ss(k_triton,k_triton)
                   vq_s(k_triton)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                          dq_s(k_triton)*grt_i(im)/temp_i(im)
                   work(j,185) = dq_s(k_triton)/y_grrho2
                   work(j,186) = vq_s(k_triton)/y_grrho2
c--------------> Total radial energy flux (tritons  )
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_triton) +
     >                     2.5*gfl_s(k,k_triton)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,187) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (tritons  )
                  rdum(1)=bsjbp_s(k_triton)
                  rdum(2)=bsjbt_s(k_triton)
	          work(j,188) = -1.e-6*rdum(1)/BTOR
	          work(j,189) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (tritons)
                  work(j,190) = ( utheta_s(1,1,k_triton) + 
     >                            utheta_s(1,2,k_triton) +
     >                            utheta_s(1,3,k_triton)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  He3 Particles  ===================================|
	 	IF (k_he3 .gt. 0) then
c--------------> Total radial particle flux (He3 particles)
		  CALL RARRAY_COPY(5,gfl_s(1,k_he3),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,101) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (He3 particles)
		  im  = jm_s(k_he3)
		  iza = IABS(jz_s(k_he3))
		  work(j,202) = dn_s(k_he3)/y_grrho2
		  work(j,203) = (vn_s(k_he3)  +
     >                           veb_s(k_he3) +
     >                           gfl_s(5,k_he3)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (He3 particles)
		  CALL RARRAY_COPY(5,qfl_s(1,k_he3),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,204) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (He3 particles)
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_he3),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_he3)=chit_ss(k_he3,k_he3) +
     >                              chip_ss(k_he3,k_he3)
                   vq_s(k_he3)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                          dq_s(k_he3)*grt_i(im)/temp_i(im)
                   work(j,205) = dq_s(k_he3)/y_grrho2
                   work(j,206) = vq_s(k_he3)/y_grrho2
c--------------> Total radial energy flux (He3 particles)
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_he3) +
     >                     2.5*gfl_s(k,k_he3)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,207) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (He3 particles)
                  rdum(1)=bsjbp_s(k_he3)
                  rdum(2)=bsjbt_s(k_he3)
	          work(j,208) = -1.e-6*rdum(1)/BTOR
	          work(j,209) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (He3 particles)
                  work(j,210) = ( utheta_s(1,1,k_he3) + 
     >                            utheta_s(1,2,k_he3) +
     >                            utheta_s(1,3,k_he3)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Alpha Particles  =================================|
	 	IF (k_alpha .gt. 0) then
c--------------> Total radial particle flux (alpha particles)
		  CALL RARRAY_COPY(5,gfl_s(1,k_alpha),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,221) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (alpha particles)
		  im  = jm_s(k_alpha)
		  iza = IABS(jz_s(k_alpha))
		  work(j,222) = dn_s(k_alpha)/y_grrho2
		  work(j,223) = (vn_s(k_alpha)  +
     >                           veb_s(k_alpha) +
     >                           gfl_s(5,k_alpha)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (alpha particles)
		  CALL RARRAY_COPY(5,qfl_s(1,k_alpha),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,224) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (alpha particles)
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_alpha),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_alpha)=chit_ss(k_alpha,k_alpha) +
     >                              chip_ss(k_alpha,k_alpha)
                   vq_s(k_alpha)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                          dq_s(k_alpha)*grt_i(im)/temp_i(im)
                   work(j,225) = dq_s(k_alpha)/y_grrho2
                   work(j,226) = vq_s(k_alpha)/y_grrho2
c--------------> Total radial energy flux (alpha particles)
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_alpha) +
     >                     2.5*gfl_s(k,k_alpha)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,227) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (alpha particles)
                  rdum(1)=bsjbp_s(k_alpha)
                  rdum(2)=bsjbt_s(k_alpha)
	          work(j,228) = -1.e-6*rdum(1)/BTOR
	          work(j,229) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (alpha particles)
                  work(j,230) = ( utheta_s(1,1,k_alpha) + 
     >                            utheta_s(1,2,k_alpha) +
     >                            utheta_s(1,3,k_alpha)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Impurity Z1      =================================|
	 	IF (k_impZ1 .gt. 0) then
c--------------> Total radial particle flux (impZ1 particles)
		  CALL RARRAY_COPY(5,gfl_s(1,k_impZ1),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,241) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (impZ1 particles)
		  im  = jm_s(k_impZ1)
		  iza = IABS(jz_s(k_impZ1))
		  work(j,242) = dn_s(k_impZ1)/y_grrho2
		  work(j,243) = (vn_s(k_impZ1)  +
     >                           veb_s(k_impZ1) +
     >                           gfl_s(5,k_impZ1)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (impZ1 particles)
		  CALL RARRAY_COPY(5,qfl_s(1,k_impZ1),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,244) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (impZ1 particles)
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_impZ1),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_impZ1)=chit_ss(k_impZ1,k_impZ1) +
     >                              chip_ss(k_impZ1,k_impZ1)
                   vq_s(k_impZ1)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                          dq_s(k_impZ1)*grt_i(im)/temp_i(im)
                   work(j,245) = dq_s(k_impZ1)/y_grrho2
                   work(j,246) = vq_s(k_impZ1)/y_grrho2
c--------------> Total radial energy flux (impZ1 particles)
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_impZ1) +
     >                     2.5*gfl_s(k,k_impZ1)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,247) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (impZ1 particles)
                  rdum(1)=bsjbp_s(k_impZ1)
                  rdum(2)=bsjbt_s(k_impZ1)
	          work(j,248) = -1.e-6*rdum(1)/BTOR
	          work(j,249) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (impZ1 particles)
                  work(j,250) = ( utheta_s(1,1,k_impZ1) + 
     >                            utheta_s(1,2,k_impZ1) +
     >                            utheta_s(1,3,k_impZ1)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Impurity Z2      =================================|
	 	IF (k_impZ2 .gt. 0) then
c--------------> Total radial particle flux (impZ2 particles)
		  CALL RARRAY_COPY(5,gfl_s(1,k_impZ2),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,261) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (impZ2 particles)
		  im  = jm_s(k_impZ2)
		  iza = IABS(jz_s(k_impZ2))
		  work(j,262) = dn_s(k_impZ2)/y_grrho2
		  work(j,263) = (vn_s(k_impZ2)  +
     >                           veb_s(k_impZ2) +
     >                           gfl_s(5,k_impZ2)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (impZ2 particles)
		  CALL RARRAY_COPY(5,qfl_s(1,k_impZ2),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,264) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (impZ2 particles)
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_impZ2),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_impZ2)=chit_ss(k_impZ2,k_impZ2) +
     >                              chip_ss(k_impZ2,k_impZ2)
                   vq_s(k_impZ2)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                          dq_s(k_impZ2)*grt_i(im)/temp_i(im)
                   work(j,265) = dq_s(k_impZ2)/y_grrho2
                   work(j,266) = vq_s(k_impZ2)/y_grrho2
c--------------> Total radial energy flux (impZ2 particles)
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_impZ2) +
     >                     2.5*gfl_s(k,k_impZ2)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,267) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (impZ2 particles)
                  rdum(1)=bsjbp_s(k_impZ2)
                  rdum(2)=bsjbt_s(k_impZ2)
	          work(j,268) = -1.e-6*rdum(1)/BTOR
	          work(j,269) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (impZ2 particles)
                  work(j,270) = ( utheta_s(1,1,k_impZ2) + 
     >                            utheta_s(1,2,k_impZ2) +
     >                            utheta_s(1,3,k_impZ2)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c=================  Impurity Z3      =================================|
	 	IF (k_impZ3 .gt. 0) then
c--------------> Total radial particle flux (impZ3 particles)
		  CALL RARRAY_COPY(5,gfl_s(1,k_impZ3),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,281) = rdum(6)*VRS(j)*1.e-19
c--------------> Particle diffusion and velocity (impZ3 particles)
		  im  = jm_s(k_impZ3)
		  iza = IABS(jz_s(k_impZ3))
		  work(j,282) = dn_s(k_impZ3)/y_grrho2
		  work(j,283) = (vn_s(k_impZ3)  +
     >                           veb_s(k_impZ3) +
     >                           gfl_s(5,k_impZ3)/den_iz(im,iza))
     >                           /y_grrho2
c--------------> Radial conduction flux (impZ3 particles)
		  CALL RARRAY_COPY(5,qfl_s(1,k_impZ3),1,rdum,1)
		  rdum(6) = RARRAY_SUM(5,rdum,1)
		  work(j,284) = rdum(6)*VRS(j)*1.e-6
c--------------> Heat conduction and velocity (impZ3 particles)
c                  Conduction total is sum of components
                   rdum(1)=RARRAY_SUM(5,qfl_s(1,k_impZ3),1)
c                  Diagonal conductivity plus convective velocity
                   dq_s(k_impZ3)=chit_ss(k_impZ3,k_impZ3) +
     >                              chip_ss(k_impZ3,k_impZ3)
                   vq_s(k_impZ3)=rdum(1)/(den_iz(im,iza)*
     >                                    z_j7kv*temp_i(im)) +
     >                          dq_s(k_impZ3)*grt_i(im)/temp_i(im)
                   work(j,285) = dq_s(k_impZ3)/y_grrho2
                   work(j,286) = vq_s(k_impZ3)/y_grrho2
c--------------> Total radial energy flux (impZ3 particles)
                  DO k=1,5
                   rdum(k)=qfl_s(k,k_impZ3) +
     >                     2.5*gfl_s(k,k_impZ3)*temp_i(im)*z_j7kv
                  ENDDO
                  rdum(6)=RARRAY_SUM(5,rdum,1)
	          work(j,287) = rdum(6)*VRS(j)*1.e-6
c--------------> Bootstrap current on p'/p (impZ3 particles)
                  rdum(1)=bsjbp_s(k_impZ3)
                  rdum(2)=bsjbt_s(k_impZ3)
	          work(j,288) = -1.e-6*rdum(1)/BTOR
	          work(j,289) = -1.e-6*rdum(2)/BTOR
c--------------> Poloidal flow velociry on outside midplane (impZ3 particles)
                  work(j,290) = ( utheta_s(1,1,k_impZ3) + 
     >                            utheta_s(1,2,k_impZ3) +
     >                            utheta_s(1,3,k_impZ3)  ) *
     >                          BTOR/(1.0 + p_eps)/p_fhat
		ENDIF
c================ Miscellaneous Parameters ===========================|
c--------------> Bootstrap current
                  rdum(1)=p_bsjb
                  work(j,301)  = -1.e-6*rdum(1)/BTOR
c--------------> External current
                  rdum(2)=p_exjb
	          work(j,302)  = -1.e-6*rdum(2)/BTOR
c--------------> Current conductivity
                  work(j,303)  = 1.e-6/p_etap


	ENDDO   !<<<<--------- End of j-loop
c=====================================================================|
		do narray=101,320
		 work(NA1,narray) = work(NA,narray)
c		 work(NA1,narray) = work(NA,narray)*(ROC/HRO-(NA-1))+
c     >                              work(NA-1,narray)*(NA-ROC/HRO)
		enddo
c		do narray=101,320
c		 work(NA-1,narray) = work(NA-2,narray)
c		 work(NA,narray) = work(NA-2,narray)
c		 work(NA1,narray) = work(NA-2,narray)
c		enddo
 1000	end

c#####################################################################|
c################    END OF THE DRIVER   !!!!    #####################|
c#####################################################################|
	SUBROUTINE ZBFAUX(GRRdB2,NGRTheta,YFM)
c Returns <(grad(rho)/B)2>
c         <n.grad(theta)>
c         and F_m - poloidal moments of geometric factor for PS viscosity
c                                  A.V.Zolotukhin 09-09-2003
c                  corrected by I.Yu. Senichenkov (IYS) February 2008
c-----------------------------------------------------------------------

	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'

	integer	j_theta, j_theta1, ntheta, j_rho, j_mvisc
	double precision YH,yametr,ydroda,yrho,I_str, THETA_cap, 
     #                   y_gamma, yr_theta, theta_s
	double precision yLambda, yLambdap, yDelta, yDeltap, 
     #               ytriang, ytriangp
	double precision COS_th,COS2_th,yd,y_sqg,y_sqg1,grr2dB2,yBm1
	double precision Ynablarho2,GRRdB2(1000),NGRTheta(1000),
     #                   YFM(3,1000)
	double precision yBmod,yBtheta,yBlnBCOS(3),yBlnBSIN(3),
     #                   yB2COS(3),yB2SIN(3)
	double precision THETA_gen,theta1_s
	double precision COS_th1, COS2_th1, yr_theta1, yd1, yBmod1, 
     #                   yBtheta1,y_sqg11
	double precision yVRR(1000),y_gamma1,yBm(1000)
	parameter( NTHETA=32 )


	DO j_rho=1, NA                                   ! j_rho loop begins
!!**!!	  if (j_rho .lt. NA)  then
            YH = AMETR(j_rho+1) - AMETR(j_rho)
!!**!!          else
!!**!!            YH = HROA
!!**!!          endif

	  yametr    = AMETR(j_rho)   !!**!!  0.5*(AMETR(j_rho)+AMETR(j_rho+1))
	  yrho      = RHO(j_rho)     !!**!!
	  ydroda    = DRODA(j_rho)   !!**!!
	  I_str     = IPOL(j_rho)*RTOR*BTOR
	  THETA_cap = BTOR*yrho*ydroda*MU(j_rho)/I_str   !!**!!

c--------> Magnetic surface ellipticity and its derivative
	  yLambda  = ELON(j_rho)
	  yLambdap = (ELON(j_rho+1)-ELON(j_rho))/YH

c--------> Magnetic axis shift and its derivative
	  yDelta   = SHIF(j_rho)
	  yDeltap  = (SHIF(j_rho+1)-SHIF(j_rho))/YH

c--------> Magnetic surface triangularity and its derivative
	  ytriang  = TRIA(j_rho)                     !!**!!
	  ytriangp = (TRIA(j_rho+1)-TRIA(j_rho))/YH  !!**!!


	  GRRdB2(j_rho)   = 0.
	  NGRTheta(j_rho) = 0.
	  yVRR(j_rho)     = 0.
	  yBm(j_rho)      = 0.
		do j_mvisc=1,3
		   YFM(j_mvisc,j_rho) = 0.
		   yBlnBCOS(j_mvisc) = 0.
		   yBlnBSIN(j_mvisc) = 0.
		   yB2COS(j_mvisc)   = 0.
		   yB2SIN(j_mvisc)   = 0.
		enddo

	  DO j_theta=NTHETA,1,-1                                 ! j_theta loop begins
	    theta_s = GP2*j_theta/NTHETA
	    COS_th  = COS(theta_s)
	    COS2_th = COS_th*COS_th

!!**!!	    yr_theta = RTOR + yDelta + yametr*COS_th
!------- changed by IYS Feb 2008 ---------
	    yr_theta = RTOR + yDelta + 
     #         yametr*(COS_th - ytriang * (1.-COS2_th))        ! corrected by IYS Feb 2008

!------- changed by IYS Feb 2008 ---------
!!**!!	    yd = yLambda + yLambda*yDeltap*COS_th +
!!**!!     >                      yametr*yLambdap*(1.-COS2_th)
!------- changed by IYS Feb 2008 ---------
          yd = yLambda * (1.0 + yDeltap * Cos_th) + (1.-COS2_th) *
     #       (yametr * yLambdap * (1.0 + 2.0 * ytriang * Cos_th) + 
     #       yLambda * (ytriang - yametr * ytriangp) * Cos_th)
!------- changed by IYS Feb 2008 ---------

c----------> Square root of the determinant g
	    y_sqg = yr_theta*yametr*yd

c----------> sqrt(g)/rho
	    y_sqg1 = yr_theta*yd


c----------> B - magnetic field
!!**!!	    yBmod = I_str/(yr_theta*yd)*
!!**!!     >	      sqrt(yd**2 + THETA_cap**2*(1. +(yLambda**2-1.)*COS2_th))
!------- changed by IYS Feb 2008 ---------
	    yBmod = I_str/(yr_theta*yd)*sqrt(yd**2 + THETA_cap**2 *
     >	       ((1.0 + 2.0 * ytriang * Cos_th)**2 * (1.-COS2_th) +
     #            yLambda**2 * Cos2_th))
!------- changed by IYS Feb 2008 ---------

c----------> 1/B
	    yBm1 = 1./yBmod

c----------> B_theta - contravariant component of the magnetic field
!!**!!		yBtheta = -BTOR*MU(j_rho)/y_sqg1
		yBtheta = BTOR*MU(j_rho)*yrho*ydroda/y_sqg

c----------> nabla(rho)2
!!**!!		ynablarho2 = (1+(yLambda**2-1.)*COS2_th)/yd**2
c----------> nabla(rho)2 -- corrected by IYS Feb 2008
		ynablarho2 = ydroda**2 * ((1.0 + 2.0 * ytriang * Cos_th)**2 * 
     #            (1.-COS2_th) + yLambda**2 * Cos2_th) / yd**2


c----------> nabla(rho)2/B2
	    grr2dB2 = ynablarho2/(yBmod*yBmod)

c----------> Integral[sqrt(g)/rho, 0, 2pi]
	    yVRR(j_rho) = yVRR(j_rho) + y_sqg1*GP2/NTHETA


c----------> Integral[(nabla(rho)/B)2*SQRT(g)/rho, 0, 2pi]
	    GRRdB2(j_rho)   = GRRdB2(j_rho)+ grr2dB2*y_sqg1*GP2/NTHETA

c----------> Integral[1/B, 0, 2pi]
!!**!!	    NGRTheta(j_rho) = NGRTheta(j_rho) + yBm1*GP2/NTHETA

	    THETA_gen = 0.
	      DO j_theta1=NTHETA,1,-1                            ! j_theta1 loop begins
	        theta1_s   = theta_s*j_theta1/NTHETA
		COS_th1    = COS(theta1_s)
		COS2_th1   = COS_th1*COS_th1

!!**!!		yr_theta1  = RTOR + yDelta + yametr*COS_th1
	    yr_theta1 = RTOR + yDelta + 
     #         yametr*(COS_th1 - ytriang * (1.-COS2_th1))        ! corrected by IYS Feb 2008

!!**!!	        yd1        = yLambda + yLambda*yDeltap*COS_th1   +
!!**!!     >                                 yametr*yLambdap*(1.-COS2_th1)
 !------- changed by IYS Feb 2008 ---------
          yd1 = yLambda * (1.0 + yDeltap * Cos_th1) + (1.-COS2_th1) *
     #       (yametr * yLambdap * (1.0 + 2.0 * ytriang * Cos_th) + 
     #       yLambda * (ytriang - yametr * ytriangp) * Cos_th1)
!------- changed by IYS Feb 2008 ---------

               y_sqg11    = yr_theta1*yd1

c----------> B1 - magnetic field
!!**!!	        yBmod1 = I_str/(yr_theta1*yd1)*
!!**!!     >	                  sqrt(yd1**2 +
!!**!!     >                     THETA_cap**2*(1. +(yLambda**2-1.)*COS2_th1))
!------- changed by IYS Feb 2008 ---------
	       yBmod1 = I_str/(yr_theta1*yd1)*sqrt(yd1**2 + THETA_cap**2 *
     >	       ((1.0 + 2.0 * ytriang * Cos_th1)**2 * (1.-COS2_th1) +
     #            yLambda**2 * Cos2_th1))
!------- changed by IYS Feb 2008 ---------

c----------> B_theta1 - contravariant component of the magnetic field
!!**!!		yBtheta1 = -BTOR*MU(j_rho)/y_sqg11
		yBtheta1 = BTOR*MU(j_rho)*yrho*ydroda/y_sqg11/yametr

                THETA_gen = THETA_gen + yBmod1/yBtheta1*theta_s/NTHETA
		if (j_theta.eq.NTHETA) y_gamma1 = THETA_gen
              ENDDO                                            ! j_theta1 loop ends

		if (j_theta.eq.NTHETA) y_gamma = GP2/y_gamma1

	      THETA_gen = y_gamma*THETA_gen


	      DO j_mvisc=1,3                                   ! j_mvisc loop begins
		yBlnBCOS(j_mvisc) = yBlnBCOS(j_mvisc) +
     >                     yBmod*LOG(yBmod)*COS(j_mvisc*THETA_gen)
     >                                           *y_sqg1*GP2/NTHETA
		yBlnBSIN(j_mvisc) = yBlnBSIN(j_mvisc) +
     >                     yBmod*LOG(yBmod)*SIN(j_mvisc*THETA_gen)
     >                                           *y_sqg1*GP2/NTHETA
		yB2COS(j_mvisc)   = yB2COS(j_mvisc) +
     >                     yBmod*yBmod*COS(j_mvisc*THETA_gen)
     >                                           *y_sqg1*GP2/NTHETA
		yB2SIN(j_mvisc)   = yB2SIN(j_mvisc) +
     >                     yBmod*yBmod*SIN(j_mvisc*THETA_gen)
     >                                           *y_sqg1*GP2/NTHETA


	      ENDDO                                            ! j_mvisc loop ends
	      yBm(j_rho) = yBm(j_rho) + yBmod*y_sqg1*GP2/NTHETA
	  ENDDO                                                ! j_theta loop ends

	   GRRdB2(j_rho)   = GRRdB2(j_rho)/yVRR(j_rho)

!!**!!	   NGRTheta(j_rho) = -NGRTheta(j_rho)/
!!**!!     >                       yVRR(j_rho)*BTOR*MU(j_rho)
         NGRTheta(j_rho) = y_gamma 

           yBm(j_rho) = yBm(j_rho)/yVRR(j_rho)

	   DO j_mvisc=1,3
		   yBlnBCOS(j_mvisc) = yBlnBCOS(j_mvisc)/yVRR(j_rho)
		   yBlnBSIN(j_mvisc) = yBlnBSIN(j_mvisc)/yVRR(j_rho)
		   yB2COS(j_mvisc)   = yB2COS(j_mvisc)/yVRR(j_rho)
		   yB2SIN(j_mvisc)   = yB2SIN(j_mvisc)/yVRR(j_rho)


		   YFM(j_mvisc,j_rho) = 2.*(j_mvisc*y_gamma)**2/
     >                             (BTOR**3*BDB02(j_rho)*BDB0(j_rho))*
     >                       (  yBlnBCOS(j_mvisc)*yB2COS(j_mvisc)  +
     >                          yBlnBSIN(j_mvisc)*yB2SIN(j_mvisc)     )
           ENDDO


	 ENDDO                                                 ! j_rho loop ends

	GRRdB2(NA1) = GRRdB2(NA)*(ROC/HRO-(NA-1))+
     >                         GRRdB2(NA-1)*(NA-ROC/HRO)
	NGRTheta(NA1) = NGRTheta(NA)*(ROC/HRO-(NA-1))+
     >                         NGRTheta(NA-1)*(NA-ROC/HRO)
	DO j_mvisc=1,3
	  YFM(j_mvisc,NA1) = YFM(j_mvisc,NA)*(ROC/HRO-(NA-1))+
     >                        YFM(j_mvisc,NA-1)*(NA-ROC/HRO)
        ENDDO
         RETURN
	 END
c#####################################################################|
c#####################################################################|
