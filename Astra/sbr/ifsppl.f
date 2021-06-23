C======================================================================|
C First-Principles chi_i & chi_e bazed on 
C  M. Kotschenreuther, W. Dorland, M.A. Beer, and G.W. Hammett,
C 	"Quantitative Predictions of Tokamak Energy Confinement 
C  	 from First-Principles Simulations with Kinetic Effects"
C  	     Physics of Plasmas, Vol. 2, p. 2381, (1995). 
C
C  subroutine ip_chi2 by W. Dorland
C		bdorland@zonker.ph.utexas.edu or bdorland@pppl.gov
C
C======================================================================|
	subroutine IFSPPL(YOUT1,YOUT2,YOUT3,YOUT4,CHI_I,CHI_E)
C Input: No input parameter in the subroutine head
C
C Output: (orig_name)
C  YOUT1 (RLTcrit):  R/L_Tcrit for ITG mode
C  YOUT2 (RLTcritz): R/L_Tcrit for carbon branch
C  YOUT3 (g):        L_Tc/L_T, where L_Tc is the critical temperature gradient
C                        scale length for the deuterium branch of the ITG mode
C  YOUT4 (gamma):    linear growth rate (see comment at the end)
C  CHI_I (chi_i):  Anomalous ion thermal diffusivity from toroidal ITG mode
C  CHI_E (chi_e):  Anomalous electron thermal diffusivity from toroidal ITG mode
C
C More: (from letter of 28-Jan-99)
c The inputs to the subroutine include rho_i, v_ti, and rmajor.  These inputs
c determine the units of the output, as follows.
c 
c If you pass values of rho_i, v_ti, and rmajor to the subroutine in SI
c units, then the units of the chi's that are returned will be m**2/sec.  If
c you pass rho_i and rmajor in units of cm, and v_ti in units of cm/sec, the
c chi's will be reported in units of cm**2/sec.
c 
c In our notation, v_ti == sqrt(T_i/m_i); i.e., there is a factor of sqrt(2)
c missing compared to some other conventions.  Similarly,
c rho_i == v_ti / Omega, using the definition of v_ti that does not include
c the factor of sqrt(2).
c 
c One common source of confusion has been the 'gnu' input parameter, and I
c would double-check that.  However, it sounds like you are doing that right,
c since the agreement is much poorer when there is an error in this input.
c For example, the critical gradient was negative for a few people because of
c errors in this input.
c 
c There is also a factor of 3/2 that is confusing in some of the
c documentation.  If your transport equation is of the form
c 
c 3/2 d(nT)/dt + ...
c 
c then the chi that is reported has this factor included correctly.  However,
c if you are solving an equation that looks like
c 
c d(nT)/dt + ...
c 
c then there is a factor of 2/3 missing.
c----------------------------------------------------------------------|
	implicit none
	include	'for/parameter.inc'
	include 'for/const.inc'
	include 'for/status.inc'
	integer	itest_chi,j
	double precision RLTCR,SHAT,VTI,  YTAUB,YSHAT,YRLN,YNU,GAITG
	double precision YOUT1(*),YOUT2(*),YOUT3(*),YOUT4(*),CHI_I(*)
	double precision LTI,LNI,LNE,RLI,YF,ROTSH,RLTKD,CHI_E(*)
	real*8	zRLT, zRLN, zRLNe, zq, zkappa,
     .          zshat, zth, znbeam, ztau, zeps, zgnu,
     .          g_perp, zrmajor, zrho_i, zv_ti, zRLTcrit,
     .          zRLTcrc, zchi_0, zg, zgamma, zchi_i, zchi_e
	itest_chi = 0
	do	j=1,NA
	include	'fml/lti'
	    zRLT = RTOR/LTI
	include	'fml/lni'
	    zRLN = RTOR/LNI
	include	'fml/lne'
	    zRLNe = RTOR/LNE
	    zq = 1./MU(j)
	    zkappa = ELON(j)
	include	'fml/shat'
	    zshat = SHAT
c	    zshat =-AMETR(j)/MU(j)*(MU(j+1)-MU(j))/(AMETR(j+1)-AMETR(j))
c          zth == (n_i + 36 n_C)/(n_e - n_beam)
	    zth = ZEF(j)		! *NE(j)/(NE(j)-NIBM(j))
	    znbeam = NIBM(j)/NE(j)
	    ztau = TI(j)/TE(j)
	    zeps = AMETR(j)/RTOR
	    zgnu = 2.5*RTOR*NE(j)/sqrt(TI(j)*TE(j)**3)
	include	'fml/rotsh'
	    g_perp = 0.
	    g_perp = ROTSH
	    zrmajor = RTOR
	    RLI=0.003235*sqrt(AMAIN(j)*TI(j))/BTOR
	    zrho_i = RLI
	    zv_ti = 3.1E5*sqrt(TI(j)/AMAIN(j))
	    call ip_chi2(itest_chi, zRLT, zRLN, zRLNe, zq, zkappa,
     .          	 zshat, zth, znbeam, ztau, zeps, zgnu,
     .          	 g_perp, zrmajor, zrho_i, zv_ti, 
C Output:
     .          zRLTcrit, zRLTcrc, zchi_0, zg, zgamma, zchi_i, zchi_e)
C	     YF = VR(j)*HRO/(AMETR(J+1)-AMETR(J))*GRADRO(j)/G11(j)
	     YF = VR(j)*DRODA(j)*GRADRO(j)/G11(j)
	     YOUT1(j) = zRLTcrit
	     YOUT2(j) = zRLTcrc
	     YOUT3(j) = zg
	     YOUT4(j) = zgamma
	     CHI_I(j) = YF*zchi_i
	     CHI_E(j) = YF*zchi_e
	enddo
	YOUT1(NA) = YOUT1(NA-1)
	YOUT2(NA) = YOUT2(NA-1)
	YOUT3(NA) = YOUT3(NA-1)
	YOUT4(NA) = YOUT4(NA-1)
	CHI_I(NA) = CHI_I(NA-1)
	CHI_E(NA) = CHI_E(NA-1)
	YOUT1(NA1) = YOUT1(NA)
	YOUT2(NA1) = YOUT2(NA)
	YOUT3(NA1) = YOUT3(NA)
	YOUT4(NA1) = YOUT4(NA)
	CHI_I(NA1) = CHI_I(NA)
	CHI_E(NA1) = CHI_E(NA)
	end
C======================================================================|
      subroutine ip_chi2 (itest_chi, RLT, RLN, RLNe, q, kappa,
     .                    shat, zth, nbeam, tau, eps, gnu,
     .                    g_perp, rmajor, rho_i, v_ti, RLTcrit,
     .                    RLTcritz, chi_0, g, gamma, chi_i, chi_e)
c
c  modified to use REAL*8  HSJ 2/11/9
c  modified to use tshat everywhere instead of shat 2/20/97 HSJ
c
c Note: Watch apostrophes in comments
c
c The formulas embodied in this subroutine are documented in the Physics
c of Plasmas article entitled Quantitative Predictions of Tokamak Energy
c Confinement from First-Principles Simulations with Kinetic Effects,
c by M. Kotschenreuther, W. Dorland, M.A. Beer, and G.W. Hammett,
c Vol. 2, p. 2381, (1995). Extensions to non-circular cross-sections are
c described below.
c
c There is a significant typographical error in that paper.  R/Ln* is
c defined to be max(6,R/Ln); it should be min(6,R/Ln) as defined in
c this subroutine.
c
c Also, note that in deriving these formulas, we assumed that the
c density gradient scale lengths for the different species were equal.
c This is an approximation that needs to be relaxed.
c
c As emphasized in the paper, these formulas were derived numerically and
c are therefore not trustworthy outside a particular region of parameter
c space.  For example, we did not parameterize the heat flux in the weak
c magnetic shear limit; thus, one should not use the model in this limit.
c I have attempted to reduce related strange numerical behaviors by
c limiting some inputs to be roughly within their range of validity.
c
c Questions, problems, errors, etc. should be reported to
c bdorland@zonker.ph.utexas.edu or bdorland@pppl.gov.
c
c I will reply as quickly as possible.
c
c Stiffness:
c **********
c For many cases that we have simulated, the transport equations
c tend to be very stiff.  That is, the plasma temperature gradient scale
c length tends to adjust itself to be close to the critical gradient scale
c length over some region of the plasma, because chi becomes very large
c very fast for shorter temperature gradient scale lengths.  Typically,
c we have had to be very careful with the numerical algorithm used in the
c transport equation solver with this experience in mind.  The details
c of our implementation are available to anyone that is interested.
c
c Geometry:
c ********
c
c The nonlinear simulations that were done to obtain these formulas were
c mostly done in a simplified geometry, using a shifted circle, low beta,
c high aspect ratio expansion.  Some modifications due to more sophisticated
c geometrical models have been calculated and have been included here,
c but should be considered preliminary.  There are two important issues
c that must be noted.  First, we derived our formulas using a different
c radial coordinate.  Second, since we are actually calculating the
c transport coefficients in general geometry, we require less assumptions
c for the form of the transport equation to be solved.
c
c Let me describe the <|grad rho|> issue first:
c
c     The database standard modeling assumptions that were agreed upon for
c     this exercise include the assumption that the anomalous fluxes for
c     non-circular cross-sections are simply related to the anomalous
c     fluxes for related circular cross-section plasmas.  That is, in order
c     to get the factor of < |grad rho|**2 > that appears as a coefficient
c     of chi in the energy transport equations, one assumes that
c
c     chi_anom_general = chi_anom_circular * (grad rho).
c
c     One need not make this assumption; one can just calculate the quantity
c     chi_anom_general directly.  One would then have a transport equation
c     of the form
c
c   (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|> n chi d/drho(T)] + ...
c
c     in which (grad rho) appears to the first power, rather than the second,
c     and chi is the thermal diffusivity from a general geometry theory.
c
c     This is arguably the better way to proceed, since
c
c         Vprime <|grad rho|> = A
c
c     where A is the surface area.  In this form, the quantity
c
c          -n chi dT/drho
c
c     can be identified as the heat flux per unit area, a natural
c     quantity from a theoretical/simulation point of view.  This chi is
c     the quantity returned by this subroutine.
c
c     If you are solving the transport equations in the form that
c     the ITER Expert Group agreed upon, e.g.,
c
c (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|**2> n chi d/drho(T)] + ...
c
c     then you need to multiply the chi_i and chi_e reported by this
c     subroutine by the factor <|grad rho|>/<|grad rho|**2>.  This should
c     result in only small corrections to the predicted profiles.
c
c The choice of radial coordinate is more difficult to resolve:
c
c     We did not use the sqrt(toroidal flux) radial coordinate in our
c     non-circular cross-section simulations.  Instead, we used "rho_d",
c     where rho_d is defined to be the average horizontal minor radius
c     at the elevation of the magnetic axis, normalized to the value of
c     this quantity at the LCFS.
c
c     In other words, denote by "d" the horizontal diameter of a given flux
c     surface measured at the elevation of the magnetic axis.  Denote by
c     "D" the horizontal diameter of the last closed flux surface at the
c     elevation of the magnetic axis.  Then rho_d = d/D.  I believe this
c     is variable number 67 (RMINOR) in the ITER Profile Database Standard
c     List.
c
c     It is not difficult to allow for an arbitrary radial coordinate
c     in a transport code.  One must obtain all of the radial
c     quantities as functions of rho_d rather than rho via interpolation.
c
c     However, I do not expect everyone to go to this length to test our
c     model, since you agreed to use the sqrt(toroidal flux) definition
c     of the radial coordinate.  Thus, I suggest the following alternative:
c     Simply use the rho_d coordinate to define the scale lengths that
c     appear in the formulas below.  For most quantities (such as R/LT),
c     this simply amounts to including an additional factor d rho/d rho_d
c     in the expressions passed to this subroutine.  While not completely
c     correct, this workaround captures the dominant effect, related to
c     evaluating the flux near the critical gradient.
c
c ****** Summary of comments:  ***********
c
c (1) The general geometry extensions to the IFS/PPPL model were derived
c     using rho = d/D = RMINOR as the radial coordinate.  To be
c     most accurate, the transport equation should be solved using d/D
c     as the radial coordinate.
c
c     If you use rho proportional to sqrt(toroidal flux) instead of rho=d/D
c     as your radial coordinate, you should at least carefully define the
c     scale lengths as indicated below (using rho=d/D).
c
c (2) This routine should be used to return the thermal transport
c     coefficients (chi_i, chi_e) for energy transport equations of the form
c
c (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|> n chi d/drho(T)] + ...
c
c     Note that <|grad rho|> only appears to the first power according
c     to this definition of chi.  If your code is hardwired to solve an
c     equation of the form
c
c (3/2) d(n T)/dt = (1/V') d/drho V' <|grad rho|**2> n chi d/drho(T)] + ...
c
c     then multiply the chi_i and chi_e obtained from this routine by
c     the factor <|grad rho|>/<|grad rho|**2>.
c
c *****************************************************************
c
c     RLT  R/L_Ti, where R = the major radius and
c                  1/L_Ti = -1/T_i dT_i/drho_d
c     RLN  R/L_ni, where 1/L_ni = -1/n_i dn_i/drho_d
c     RLNe R/L_ne, where 1/L_ne = -1/n_e dn_e/drho_d
c     q    The safety factor.
c     kappa  The elongation, defined here to be the
c          kappa = maximum height/maximum diameter
c     shat == rho_d/q dq/drho_d
c     zth  Thermal Z_eff.  The simulations that were carried out to
c          generate the formulae in this subroutine assumed the plasma
c          was composed of a thermal hydrogenic species, thermal carbon,
c          a hydrogenic beam species, and electrons.  We found that low-Z
c          impurities primarily act to dilute the main ion concentration,
c          and can accounted for to first order by modifying the
c          definition of Z_eff.  Some of the more important effects of
c          the fast ions in the plasma are also partially accounted
c          for by this parameter, which is:
c          zth == (n_i + 36 n_C)/(n_e - n_beam)
c
c     nbeam == local fast ion (beam) density normalized to the electron
c          density.
c     tau  == T_i/T_e.  Note that this is opposite to a widely
c                       used convention.
c     eps  == rho_d/R, the local minor radius normalized to the major radius.
c     gnu   Dimensionless collisionality parameter.
c          gnu == 2.5e-7 * n_e / (T_e**1.5 T_i**0.5) * rmajor
c          where n_e is in units of cm**-3, T_e and T_i are in eV, and
c          rmajor is in units of m.  For an R = 2.4 m, 100 eV, 1e13 plasma,
c          gnu=600.
c     beta  Local total beta; presently ignored.
c     g_perp == velocity shear parameter.  Use g_perp=0 (actual value
c          discussed in the Waltz paper cited below).
c     rmajor == major radius of the plasma
c     rho_i == local thermal gyroradius of thermal hydrogenic species.
c     v_ti == sqrt(T_i/m_i) where T_i and m_i are the local thermal
c          hydrogenic temperature and average thermal hydrogenic mass.
c
c     rho_e,v_te,etg,rlte: dummy parameters at present. Ignore.
c
c     Units: The only dimensional parameters in the inputs are the major
c     radius, rho_i, and v_t.  Their units should be consistent; the chis
c     that are returned will be in units of rho_i**2 v_ti / rmajor.
c
c OUTPUT:
c     RLTcrit: R/L_Tcrit for ITG mode
c     RLTcritz: R/L_Tcrit for carbon branch
c     chi_0: normalized chi (ignore)
c     g: L_Tc/L_T, where L_Tc is the critical temperature gradient
c        scale length for the deuterium branch of the ITG mode.
c     chi_i: Anomalous ion thermal diffusivity from toroidal ITG mode.
c     chi_e: Anomalous electron thermal diffusivity from toroidal ITG mode.
c
c     This parameterization of chi is not complete.  There are significant
c     neglected physical processes that are known to be important in
c     many operational regimes.
c
c     The most significant problems are:
c
c     (1) Trapped ion/long wavelength ITG modes.  These modes are known
c     to be unstable for typical edge tokamak parameters.  However, until
c     we have nonlinear estimates of the associated thermal diffusivity,
c     these modes are ignored, leading to overly optimistic predictions of
c     edge thermal confinement.
c     (2) Trapped electron modes, which can alter the stability boundary
c     significantly for low collisionality.  At high collisionality these
c     modes are generally stable and thus largely irrelevant.  When present,
c     they are associated most strongly with particle transport, although
c     there is also an associated heat transport.
c     (3) Minority ion density gradients, which can strongly change chi
c     and LT_crit.
c     (4) Sheared flows, which are stabilizing.  This includes diamagnetic
c     and ExB shear flows.
c     (5) Finite beta effects, generally stabilizing.
c
c jek 1/21/97 added gamma to argument list and compute rot shear
c     outside of routine (gperp = 0)
c
      implicit none
c
      integer itest_chi
c
c     replaced REAL with REAL*8 in the following HSJ 2/11/96
c
      real*8 RLT,RLN,RLNe,shat,zth,tau,rmajor,chi_i,q,
     .       nbeam,taub,rho_i,v_ti,eps,chi_e,nu,chi_0,g,g_perp,gnu
      real*8 gamma,rot
      real*8 RLTcrit,RLTcritz,chi0,a_0,b_0
      real*8 f_0,f_z,chiz,g_facz,c_0
      real*8 c1,trln,trlne,tshat,kappa
      real*8 chie1,chie2
      real*8 g_fac1,gfac
      real*8 zero, sixth, fourth, half, one, two, six  ! for portability
c
      data a_0/ 0.0/
      data b_0/ 0.0/
      data c_0/ 1.0/
****  data a_0/-1.0/
****  data b_0/-1.0/
c
      save a_0, b_0, c_0
c
      zero   = 0.0
      sixth  = 1.0 / 6.0
      fourth = 0.25
      half   = 0.5
      one    = 1.0
      two    = 2.0
      six    = 6.0
c
      if (c_0 .eq. -1) then
        write (*, *) 'C_0 multiplier'
        read  (*, *)  c_0
      end if
c
      if (b_0 .eq. -1) then
        write (*, *) 'beta=?'
        read  (*, *)  b_0
      end if
      taub  = tau / (1.0-nbeam)
      nu    = gnu * 0.84
c
      tRLN  = MIN (ABS (RLN ), six) * SIGN (one, RLN)
      tRLNe = MIN (ABS (RLNe), six)
c
c Formula is not applicable for shat<0.5; doesn't matter most places.
c
      tshat = MAX (shat, half)
c
c     critical ion temperature gradient:
c
      RLTcrit = 2.46*(1.+2.78/q**2)**0.26*(zth/2.)**0.7*taub**0.52
     .          *( (0.671+0.570*tshat-0.189*tRLN)**2
     .          +0.335*tRLN+0.392-0.779*tshat+0.210*tshat**2)
     .          *( 1.-0.942*(2.95*eps**1.257/nu**0.235-0.2126)
     .          *zth**0.516 / ABS (tshat)**0.671)
c
      RLTcritz = 0.75 * (1.0+taub) * (1.0+tshat)
     .                * MAX (one, 3.0 - 2.0 * tRLNe / 3.0)
     .                * (1.0 + 6.0 * MAX (zero, 2.9-zth))
c
      c1 = 1.0
      if (zth .gt. 3.0)  c1 = (3.0/zth)**1.8
c
      f_0 = 11.8*c1*q**1.13/(1.+tshat**0.84)/taub**1.07
     .      *(1.+6.72*eps/nu**0.26/q**0.96)
     .      /(1.+((kappa-1.)*q/3.6)**2)
c
      f_z = 7.88/(1.+tshat) * MAX (fourth, zth-3) / taub**0.8
     .      /(1.+((kappa-1.)*q/3.6)**2)
c
      chi0 = f_0 * rho_i**2*v_ti/rmajor
      chiz = f_z * rho_i**2*v_ti/rmajor
c
Cgrp      if (RLT-RLTcrit .gt. 0.0) then
Cgrp         g_fac1 = MIN ((RLT-RLTcrit)**0.5, (RLT-RLTcrit))
Cgrp      else
Cgrp         g_fac1 = 0.0
Cgrp      end if
c
Cgrp      if (RLT-RLTcritz.gt.0.) then
Cgrp         g_facz = MIN ((RLT-RLTcritz)**0.5, (RLT-RLTcritz))
Cgrp      else
Cgrp         g_facz = 0.0
Cgrp      end if
         g_fac1 = GFAC(RLT-RLTcrit,half)
         g_facz = GFAC(RLT-RLTcritz,half)
c
      chi_i = MAX (chi0*g_fac1, chiz*g_facz)
c
      g     = RLT / RLTcrit
      chi_0 = chi0 * SQRT (ABS (RLTcrit))
c
      chie1 = chi0*g_fac1*1.44*tau**0.4*(q/tshat)**0.3*nu**0.14
     .      * MAX (sixth, eps)
      chie1 = 0.5 * chie1*(1.+tRLNe/3.0)
c
      chie2 = 0.5 * MAX (two, (1.0+RLNe/3.0))
      chie2 = chie2*0.526*tau*nu**0.22
      chie2 = chie2*chiz*g_facz
c
c Correction for n_i/n_e and ratio of heat fluxes rather than chi_s:
c
      chi_e = MAX (chie1, chie2) * (7.0-zth)/6.0
c
c     **Preliminary** model of rotational stabilization.  Based on
c     work described in Waltz, et al., Phys. of Plasmas, Vol. 1,
c     p. 3138 (1992).  Here, the ExB shearing rate is denoted
c     as g_perp and the linear growth rate gamma_l is parameterized as:
c
      gamma = 0.25/(1+0.5*tshat**2)/tau
     .      * (1 + 3.0 * MAX (zero, eps-0.16667))
     .      / (1 + MAX (zero, q-3)/15.0)
     .      * (RLT-RLTcrit)*v_ti/rmajor
c
      rot   = MIN (one, MAX (zero, (1.0 - ABS (g_perp)/gamma)))
      if (gamma .le. 0) rot = 0.
c
      chi_i = MAX (zero, c_0*chi_i*rot)
      chi_e = MAX (zero, c_0*chi_e*rot)
c
c test
c
      if (itest_chi .eq. 1) then
       chi_i = 1.0
       chi_e = 1.0
      end if
      return
c
      end
	real*8	function GFAC(x,r)
	real*8	x,r,s,c,xl
	s = 0.382683432
	c = 0.707106781
	xl = -r*s
	xr = -xl*c
	if (x .le. xl)	then
	   GFAC = 0
	   return
	endif
	if (x .le. xr)	then
	   GFAC = r-sqrt(r*r-(x-xl)**2)
	   return
	endif
	if (x .le. 1.)	then
	   GFAC = x
	   return
	endif
	GFAC = sqrt(x)
	end
