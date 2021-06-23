!======================================================================|
      subroutine mmm4a
!		Interface MMM95 for Astra	Pereverzev 14/05/2002
!     		MMM95 was downloaded from the NTCC web page in May 2002
! To define:
!   ywexb
!
! Questions:
!  1) What to do if no impurity species is present?
!  2) Can dimension of gamma() or omega() exceed matdim=jmxdim=8?
!
! Input:	See commentaries by G.Bateman below
! 
! Output:	Astra work array
!		WORK(1:NRD,16:29+(4+jmxdim)*jmxdim=221 < 2*nrd)
!		currently jmxdim=12 must be <19
!----------------------------------------------------------------------|
      implicit	 none
      include	'for/parameter.inc'
      include	'for/const.inc'
      include	'for/status.inc'
      double precision     rotsh
      real	 yrminor,  yrmajor,  yelong
     &, ydense,  ydensh,   ydensimp, ydensfe
     &, yzeff,   ytekev,   ytikev,   yq,      ybtor
     &, avezimp, amassimp, amasshyd, aimass,  ywexbs
     &, ygrdne,  ygrdni,   ygrdnh,   ygrdnz,  ygrdte,  ygrdti,  ygrdq
     &, ythiig,  ythdig,   ytheig,   ythzig
     &, ythirb,  ythdrb,   ytherb,   ythzrb
     &, ythikb,  ythdkb,   ythekb,   ythzkb
! jmxdim  <-> matdim;      jpoints <-> npoints
! jprint <-> lprint;       jreset <-> lreset;          jprout <-> nprout
      integer    jmxdim,   jpoints,  jprout,  jmdim,   j,  j1,  j2, jna
      parameter( jmxdim=12, jpoints=1 )
      real       y_den,    yrda,     ygamma(jmxdim),   yomega(jmxdim)
     &, yvelthi(jmxdim),   yvflux(jmxdim),    ydifthi(jmxdim,jmxdim)
c
c    Note, if the mmm95 module is used in another code, the minimum 
c    dimension required for cswitch is 23, the minimum dimension for 
c    lswitch is 5 and fig, frb, and fkb must all be dimensioned to a 
c    minimum of 4.
      integer	 jprint,  jerr , jsuper,  jreset,   jswitch(5)
      real       yswitch(23),    yfig(4), yfrb(4),  yfkb(4)
!----------------------------------------------------------------------|
c
c  Input integers:
c  ---------------
c
c  matdim  = first and second dimension of transport matricies
c            difthi(j1,j2,jz) and the first dimension of 
c            velthi(j1,jz), vflux(j1,jz), gamma(j1,jz), and omega(j1,jz).
c            matdim must be at least 5
c
c  npoints = number of values of jz in all of the above arrays
c
c  nprout  = output unit number for long printout
c
      jprout = 6		! std output ! Output to /dev/tty
c
c  Input switches
c  --------------
c
c  lprint      controls the amount of printout (0 => no printout)
c              higher values yield more diagnostic output
c
      jprint = 0
c  lsuper   = 0 for simulations of all other discharges
c           > 0 for supershot simulations; substantially reduces 
c               contribution from kinetic ballooning mode
c
      jsuper = 0
c  lreset  = 0 to use only internal settings for lswitch, cswitch
c              and for the coefficients fig, frb, and fkb that control
c              the contributions form the various instability modes
c
      jreset = 0
c    Note that when lreset = 0, the values of the switches and
c    coefficients in the argument list are ignored and all the 
c    switches and coefficients are set internally.
c
c    WARNING:  use lreset > 0 only if you want to pass all the switches
c              lswitch, cswitch, fig, frb, and fkb through the 
c              argument list.
c
c    WARNING:  NTCC users desiring to use anything other than lreset = 0
c              should consult with the mmm95 code authors first.
c
      jmdim = jmxdim
      y_den = 1.e-5
!----------------------------------------------------------------------|
      jna = max(NA1N,NA1E,NA1I)
      if (jna .eq. 0)	return
      jna = min(NA1,jna)-1
      do   j = 1,na1
         if (j .gt. jna)	goto	10
c
c  Input arrays:
c  -------------
c
c      All the following 1-D arrays are assumed to be defined on flux
c      surfaces called zone boundaries where the transport fluxes are
c      to be computed.  The number of flux surfaces is given by npoints
c      (see below).  For example, if you want to compute the transport
c      on only one flux surface, set npoints = 1.
c  rminor(jz)   = minor radius (half-width) of zone boundary [m]
         yrminor = ametr(j)
c  rmajor(jz)    = major radius to geometric center of zone bndry [m]
         yrmajor = rtor+shif(j)
c  elong(jz)    = local elongation of zone boundary
         yelong = elon(j)
c  dense(jz)    = electron density [m^-3]
         ydense = 1.e19*ne(j)
c  densh(jz)    = sum over thermal hydrogenic ion densities [m^-3]
         if ( nhydr(j)+ndeut(j)+ntrit(j) .gt. y_den )	then
            ydensh = nhydr(j)+ndeut(j)+ntrit(j)
         else
            ydensh = ni(j)
         endif
         if (ydensh .lt. y_den)	ydensh = y_den
c  densimp(jz)  = sum over impurity ion densities [m^-3]
         if ( nhe3(j)+nalf(j)+niz1(j)+niz2(j)+niz3(j) .gt. y_den )  then
            ydensimp = nhe3(j)+nalf(j)+niz1(j)+niz2(j)+niz3(j)
         else
            ydensimp = y_den
         endif
         if (ydensimp .lt. y_den)	ydensimp = y_den
c  densfe(jz)   = electron density from fast (non-thermal) ions [m^-3]
         ydensfe= 1.e19*nibm(j)
c  xzeff(jz)    = Z_eff
         yzeff  = zef(j)
c  tekev(jz)    = T_e (electron temperature) [keV] 
         ytekev = te(j)
c  tikev(jz)    = T_i (temperature of thermal ions) [keV]
         ytikev = ti(j)
c  q(jz)        = magnetic q-value
         yq     = 1./mu(j)
c  btor(jz)     = ( R B_tor ) / rmajor(jz)  [tesla]
         ybtor  = rtor*btor/yrmajor
c  avezimp(jz)  = average density weighted charge of impurities
c               = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
c                 sum_imp = sum over impurity ions with charge state Z_imp
         avezimp= (zim1(j)*niz1(j)+zim2(j)*niz2(j)+zim3(j)*niz3(j)
     &		  +2.*(nhe3(j)+nalf(j))
     &            )/ydensimp
c  amassimp(jz) = average density weighted atomic mass of impurities
c               = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where 
c                 sum_imp = sum over impurity ions, each with mass M_imp
         amassimp=(aim1*niz1(j)+aim2*niz2(j)+aim3*niz3(j)
     &		  +3.*nhe3(j)+4.*nalf(j)
     &            )/ydensimp
c  amasshyd(jz) = average density weighted atomic mass of hydrogen ions
c               = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
c                 sum_hyd = sum over hydrogenic ions, each with mass M_hyd
         amasshyd=(nhydr(j)+2.*ndeut(j)+3.*ntrit(j))/ydensh
c  aimass(jz)   = mean atomic mass of thermal ions [AMU]
c               = ( sum_i n_i M_i ) / ( sum_i n_i ) where
c                 sum_i = sum over all ions, each with mass M_i
         aimass = (nhydr(j)+2.*ndeut(j)+3.*(ntrit(j)+nhe3(j))+4.*nalf(j)
     &            +aim1*niz1(j)+aim2*niz2(j)+aim3*niz3(j)
     &            )/ydensh
c  wexbs(jz)    = ExB shearing rate in [rad/s].  See  K.H. Burrell,
c                 "Effects of {ExB} velocity shear and magnetic shear 
c                 on turbulence and transport in magnetic confinement 
c                 devices", Phys. of Plasmas, 4, 1499 (1997).
         include 'fml/rotsh'
         ywexbs = rotsh
c
c    All of the following normalized gradients are at zone boundaries.
c    r   = half-width, R = major radius to center of flux surface
c    n_i = thermal ion density (sum over hydrogenic and impurity)
c    n_h = thermal hydrogenic density (sum over hydrogenic species)
c    n_Z = thermal impurity density sumed over all impurities
c    Z   = average impurity charge  sumed over all impurities
c
         yrda = -2.*(rtor+shif(j))/(ametr(j+1)-ametr(j))
c  grdne(jz) = -R ( d n_e / d r ) / n_e
         ygrdne = yrda*(ne(j+1)-ne(j))/(ne(j+1)+ne(j))
c  grdni(jz) = -R ( d n_i / d r ) / n_i
         ygrdni = yrda*(ni(j+1)-ni(j))/(ni(j+1)+ni(j))
c  grdnh(jz) = -R ( d n_h / d r ) / n_h
         ygrdnh = nhydr(j+1)+ndeut(j+1)+ntrit(j+1)
         ygrdnh = yrda*(ygrdnh-ydensh)/(ygrdnh+ydensh)
c  grdnz(jz) = -R ( d Z n_Z / d r ) / ( Z n_Z )
         ygrdte = 2.*(nhe3(j+1)+nalf(j+1))+zim1(j+1)*niz1(j+1)
     &                +zim2(j+1)*niz2(j+1)+zim3(j+1)*niz3(j+1)
         ygrdnz = 2.*(nhe3(j)+nalf(j))+zim1(j)*niz1(j)
     &		      +zim2(j)*niz2(j)+zim3(j)*niz3(j)
         ygrdnz = yrda*(ygrdte-ygrdnz)/(ygrdte+ygrdnz)
c  grdte(jz) = -R ( d T_e / d r ) / T_e
         ygrdte = yrda*(te(j+1)-te(j))/(te(j+1)+te(j))
c  grdti(jz) = -R ( d T_i / d r ) / T_i
         ygrdti = yrda*(ti(j+1)-ti(j))/(ti(j+1)+ti(j))
c  grdq (jz) =  R ( d q   / d r ) / q    related to magnetic shear
c            = -R ( d mu  / d r ) / mu
         ygrdq  = yrda*(mu(j+1)-mu(j))/(mu(j+1)+mu(j))
         ydensh = 1.e19*ydensh
         ydensimp = 1.e19*ydensimp
C         write(*,'(I3,1P,5E12.3)')j
C     & , yrminor,  yrmajor,  yelong
!----------------------------------------------------------------------|
         call mmm95 (
     &   yrminor,  yrmajor,  yelong
     & , ydense,  ydensh,   ydensimp, ydensfe
     & , yzeff,   ytekev,   ytikev,   yq,      ybtor
     & , avezimp, amassimp, amasshyd, aimass,  ywexbs
     & , ygrdne,  ygrdni,   ygrdnh,   ygrdnz,  ygrdte,  ygrdti,  ygrdq
c output:
     & , ythiig,  ythdig,   ytheig,   ythzig
     & , ythirb,  ythdrb,   ytherb,   ythzrb
     & , ythikb,  ythdkb,   ythekb,   ythzkb
     & , ygamma,  yomega,   ydifthi,  yvelthi, yvflux
c control:
     & , jmdim,   jpoints,  jprout,   jprint,  jerr
     & , jsuper,  jreset,   jswitch,  yswitch, yfig,    yfrb,    yfkb)
!----------------------------------------------------------------------|
 10      continue
c  Output:
c  -------
c
c    The following effective diffusivities represent contributions
c    to the total diffusivity matrix (difthi and velthi given below)
c    from each of the models that contribute to the Multi-Mode model.
c    Generally, these arrays are used for diagnostic output only.
c
c  thiig(jz) = ion thermal diffusivity from the Weiland model
c  thdig(jz) = hydrogenic ion diffusivity from the Weiland model
c  theig(jz) = elelctron thermal diffusivity from the Weiland model
c  thzig(jz) = impurity ion diffusivity from the Weiland model
         work(j,20) = ythiig
         work(j,21) = ythdig
         work(j,22) = ytheig
         work(j,23) = ythzig
c	    
c  thirb(jz) = ion thermal diffusivity from resistive ballooning modes
c  thdrb(jz) = hydrogenic ion diffusivity from resistive ballooning modes
c  therb(jz) = elelctron thermal diffusivity from resistive ballooning modes
c  thzrb(jz) = impurity ion diffusivity from resistive ballooning modes
         work(j,24) = ythirb
         work(j,25) = ythdrb
         work(j,26) = ytherb
         work(j,27) = ythzrb
c	    
c  thikb(jz) = ion thermal diffusivity from kinetic ballooning modes
c  thdkb(jz) = hydrogenic ion diffusivity from kinetic ballooning modes
c  thekb(jz) = elelctron thermal diffusivity from kinetic ballooning modes
c  thzkb(jz) = impurity ion diffusivity from kinetic ballooning modes
         work(j,16) = ythikb
         work(j,17) = ythdkb
         work(j,18) = ythekb
         work(j,19) = ythzkb
c
c    The following are growth rates and mode frequencies from the
c    Weiland model for drift modes such as ITG and TEM.
c    These arrays are intended for diagnostic output.
c
c  gamma(jm,jz) = growth rate for mode jm at point jz ( 1/sec )
c  omega(jm,jz) = frequency for mode jm at point jz ( radians/sec )
c
c    All of the transport coefficients are given in the following two
c    matricies for diffusion difthi and convection velthi in MKS units.
c    See the LaTeX documentation for difthi and velthi just below.
c
c    NOTE:  difthi and velthi include all of the anomalous transport.
c    There are no additional contributions to the heat fluxs from
c    charged particle convection.
c
c  difthi(j1,j2,jz) = full matrix of anomalous transport diffusivities
c  velthi(j1,jz)    = convective velocities
c  vflux(j1,jz)     = flux matrix
         do  j1	= 1,jmxdim
! maximum work(1:nrd,29+(4+jmxdim)*jmxdim<=2*nrd) ! currently jmxdim<19
            work(j,29+j1)       = ygamma(j1)
            work(j,29+j1+jmxdim) = yomega(j1)
            work(j,29+j1+jmxdim*2) = yvflux(j1)
            work(j,29+j1+jmxdim*3) = yvelthi(j1)
            do  j2	= 1,jmxdim
               work(j,29+j1+jmxdim*(3+j2)) = ydifthi(j1,j2)
            enddo
         enddo
! Truncated output for most unstable mode:
         work(j,28) = 0.
         work(j,29) = 0.
         do  j1	= 1,jmxdim
            if (work(j,28) .lt. ygamma(j1))	then
                work(j,28) = ygamma(j1)
                work(j,29) = yomega(j1)
            endif
         enddo

c  Output Integer
c  --------------
c
c  nerr        status code returned; 0 = OK; .ne. 0 indicates error
c
         if (jmdim.ne.jmxdim) write(*,*)">>> MMM4A: ?????"
         if (jerr .ne. 0)	then
            write(*,*)'>>> MMM4A: Error at j =',j
         endif
      enddo
c
c  Internal control variables:
c  ---------------------------
c
c  lswitch(j), j=1,8   integer control variables: 
c
c  cswitch(j), j=1,25   general control variables:
c
c  lswitch(1)  controls which version of the Weiland model is used
c                  Default lswitch(1) = 10
c             = 2  2 eqn  Weiland model Hydrogen \eta_i mode only
c             = 4  4 eqn  Weiland model with Hydrogen and trapped electrons
c             = 5  5 eqn  Weiland model with trapped electrons, FLR effects, 
c                         and parallel ion motion
c             = 6  6 eqn  Weiland model Hydrogen, trapped electrons,
c                    and one impurity species
c             = 7  7 eqn   Weiland model Hydrogen, trapped electrons,
c                  one impurity species, and collisions
c             = 8  8 eqn  Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, and parallel
c                  ion (hydrogenic) motion
c             = 9  9 eqn  Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, and finite beta
c             = 10 10 eqn Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, parallel
c                  ion (hydrogenic) motion, and finite beta
c             = 11 11 eqn Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, parallel
c                  ion (hydrogenic, impurity) motion, and finite beta
c
c  lswitch(2) = 0  full matrix representation for difthi and velthi
c                  Default lswitch(2) = 2
c             = 1  set diagonal matrix elements difthi and velthi
c             = 2  set diagonal matrix elements = effective diffusivities
c
c  lswitch(3)  controls \kappa scaling
c                  Default lswitch(3) = 0
c             = 0  use \kappa scaling raised to
c                  exponents (cswitch(3) - cswitch(5))
c             = 1  use (1+\kappa^2)/2 instead of \kappa scaling
c               
c  lswitch(4) > 0  to replace negative diffusivity with velocity
c                  Default lswitch(4) = 1
c
c  lswitch(5) = 1  to limit magnitude of all normalized gradients
c                  to ( major radius ) / ( ion Larmor radius )
c                  Default lswitch(5) = 1
c
c  cswitch(1)   0.5  minimum value of shear
c  cswitch(2)   3.5  coeff in exponential (fbeta-th) of kinetic ballooning model
c  cswitch(3)  -4.0  exponent of local elongation multiplying drift waves
c  cswitch(4)  -4.0  exponent of local elongation multiplying resistive
c                     ballooning modes
c  cswitch(5)  -4.0  exponent of local elongation multiplying
c                     kinetic balllooning modes
c  cswitch(6)   0.0  k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
c  cswitch(8)   1.0  coeff of beta_prime_1 in kinetic ballooning mode
c  cswitch(9)  0.15  alpha in diamagnetic stabil. in kinetic ballooning model
c  cswitch(10)  0.0  rel fract of ion thermal diffusivity given to convection 
c  cswitch(11)  0.0  rel fract of hydrogen particle diffusivity given to convection 
c  cswitch(12)  0.0  rel fract of el thermal diffusivity given to convection 
c  cswitch(13)  0.0  rel fract of impurity particle diffusivity given to convection 
c  cswitch(14)  1.0  coef of finite beta effect in weiland14 = cetain(20) 
c  cswitch(15)  0.0  min value of impurity charge state zimpz
c  cswitch(16)  0.0  coef of fast particle fraction (superthermal ions) 
c                    in weiland model -- coef of densfe
c  cswitch(17)  1.0  coeff of k_\parallel (parallel ion motion) in weiland14 
c                    = cetain(10)
c  cswitch(18)  0.0  coeff of nuhat (effect of collisions) in weiland14 
c                    = cetain(15)
c  cswitch(19)  0.0  coeff for including v_parallel in strong ballooning limit
c                    = cetain(12); cswitch(19) = 1 for inclusion of v_par effect
c  cswitch(20)  0.0  trapping fraction used in weiland14 (when > 0.0)
c                    multiplies electron trapping fraction when < 0.0
c                    no effect when cswitch(20) = 0.0
c  cswitch(21)  1.0  multiplier for wexbs (flow shear rate) in Weiland model
c  cswitch(22)  0.0  ranges from 0.0 to 1.0 adds impurity heat flow to total 
c                    ionic heat flow for the weiland model
c  cswitch(23)  0.0  controls finite diff to construct the zgm matrix 
c                    = cetain(30)
c
c     contributions to vfluxes and interchanges: 
c
c  fig(1)   hydrogen particle transport from ITG (eta_i) mode
c  fig(2)   electron thermal  transport from ITG (eta_i) mode
c  fig(3)   ion      thermal  transport from ITG (eta_i) mode
c  fig(4)   impurity particle transport from ITG (eta_i) mode
c
c  frb(1)   hydrogen particle transport from resistive ballooning mode
c  frb(2)   electron thermal  transport from resistive ballooning mode
c  frb(3)   ion      thermal  transport from resistive ballooning mode
c  frb(4)   impurity particle transport from resistive ballooning mode
c
c  fkb(1)   hydrogen particle transport from kinetic ballooning mode
c  fkb(2)   electron thermal  transport from kinetic ballooning mode
c  fkb(3)   ion      thermal  transport from kinetic ballooning mode
c  fkb(4)   impurity particle transport from kinetic ballooning mode
c
      end
!======================================================================|
!| %
!| %  mmm95.f is the source code for the MMM95 Multi-Mode transport model
!| %  to produce LaTeX documentation, type:
!| %  s2tex.py mmm95.f
!| 
!| %
!| % The following lines control how the LaTeX document is typeset
!|  
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt  \topmargin 0pt  \oddsidemargin 0pt
!| \textheight 9.0in \textwidth 6.5in
!| \begin{document}           % End of preamble and beginning of text.
!| 
!| \begin{center}
!| \Large {\bf Multi-Mode Transport Model MMM95 {\large Version 1.1}} \\
!| \vspace{1pc} \normalsize
!| Glenn Bateman, Arnold H.~Kritz \\
!| Lehigh University Physics Department \\
!| 16 Memorial Drive East, Bethlehem PA 18015 \\
!| bateman@plasma.physics.lehigh.edu \\
!| kritz@plasma.physics.lehigh.edu
!| \end{center}
!| 
!| This file documents a subroutine called {\tt mmm95}, which computes
!| plasma transport coefficients using the Multi-Mode transport model
!| which has been held fixed since 1995.  A complete derivation of the
!| MMM95 model is given in reference \cite{bate98a}.  Note, the mmm95 
!| subroutine has been tested and used when compiled with flags such that
!| double precision is used.  Note also, if this module is used in 
!| another code, the minimum dimension required for cswitch is 23, 
!| the minimum dimension for lswitch is 5 and fig, frb, and fkb 
!| must all be dimensioned to a minimum of 4.
!| 
!| When sbrtn mmm95 is used in the BALDUR transport code to compute
!| particle and thermal fluxes, smoothing of the gradients is normally
!| needed for numerical stability.  Note that the lower bounds of the
!| gradient scale lengths are limited by the poloidal ion Larmor radius
!| (zlarpo) in sbrtn mmm95.
!| 
c
c  Revision History
c  ----------------
c       date            Description
c
c   18-Dec-2001         Resistive ballooning diffusion = abs ( ... )
c
c   13-Aug-2001         Fixed normalization of gamma and omega
c
c   26-Mar-1999         Added statement at beginning of do 300 loop
c                       to avoid computation of fluxes at 
c                       rminor(jz) .lt. 1.e-4 * rmajor(jz)
c
c   09-Mar-1999         Revamped comments and included changes as 
c                 	suggested by D. McCune, module reviewer for
c                       NTCC.
c

c@mmm95.tex
c--------1---------2---------3---------4---------5---------6---------7-c
c
      subroutine mmm95 (
     &   rminor,  rmajor,   elong
     & , dense,   densh,    densimp,  densfe
     & , xzeff,   tekev,    tikev,    q,       btor
     & , avezimp, amassimp, amasshyd, aimass,  wexbs
     & , grdne,   grdni,    grdnh,    grdnz,   grdte,   grdti,  grdq
     & , thiig,   thdig,    theig,    thzig
     & , thirb,   thdrb,    therb,    thzrb
     & , thikb,   thdkb,    thekb,    thzkb
     & , gamma,   omega,    difthi,   velthi,  vflux
     & , matdim,  npoints,  nprout,   lprint,  nerr
     & , lsuper,  lreset,   lswitch,  cswitch, fig,    frb,     fkb)
c
c
c    All the following 1-D arrays are assumed to be defined on flux
c    surfaces called zone boundaries where the transport fluxes are
c    to be computed.  The number of flux surfaces is given by npoints
c    (see below).  For example, if you want to compute the transport
c    on only one flux surface, set npoints = 1.
c
c    Note, if the mmm95 module is used in another code, the minimum 
c    dimension required for cswitch is 23, the minimum dimension for 
c    lswitch is 5 and fig, frb, and fkb must all be dimensioned to a 
c    minimum of 4.

c  Input arrays:
c  -------------
c
c  rminor(jz)   = minor radius (half-width) of zone boundary [m]
c  rmajor(jz)   = major radius to geometric center of zone bndry [m]
c  elong(jz)    = local elongation of zone boundary
c
c  dense(jz)    = electron density [m^-3]
c  densh(jz)    = sum over thermal hydrogenic ion densities [m^-3]
c  densimp(jz)  = sum over impurity ion densities [m^-3]
c  densfe(jz)   = electron density from fast (non-thermal) ions [m^-3]
c
c  xzeff(jz)    = Z_eff
c  tekev(jz)    = T_e (electron temperature) [keV] 
c  tikev(jz)    = T_i (temperature of thermal ions) [keV]
c  q(jz)        = magnetic q-value
c  btor(jz)     = ( R B_tor ) / rmajor(jz)  [tesla]
c
c  avezimp(jz)  = average density weighted charge of impurities
c               = ( sum_imp n_imp Z_imp ) / ( sum_imp n_imp ) where
c                 sum_imp = sum over impurity ions with charge state Z_imp
c
c  amassimp(jz) = average density weighted atomic mass of impurities
c               = ( sum_imp n_imp M_imp ) / ( sum_imp n_imp ) where 
c                 sum_imp = sum over impurity ions, each with mass M_imp
c
c  amasshyd(jz) = average density weighted atomic mass of hydrogen ions
c               = ( sum_hyd n_hyd M_hyd ) / ( sum_hyd n_hyd ) where
c                 sum_hyd = sum over hydrogenic ions, each with mass M_hyd
c
c  aimass(jz)   = mean atomic mass of thermal ions [AMU]
c               = ( sum_i n_i M_i ) / ( sum_i n_i ) where
c                 sum_i = sum over all ions, each with mass M_i
c
c  wexbs(jz)    = ExB shearing rate in [rad/s].  See  K.H. Burrell,
c                 "Effects of {ExB} velocity shear and magnetic shear 
c                 on turbulence and transport in magnetic confinement 
c                 devices", Phys. of Plasmas, 4, 1499 (1997).
c
c    All of the following normalized gradients are at zone boundaries.
c    r = half-width, R = major radius to center of flux surface
c
c  grdne(jz) = -R ( d n_e / d r ) / n_e
c  grdni(jz) = -R ( d n_i / d r ) / n_i
c  grdnh(jz) = -R ( d n_h / d r ) / n_h
c  grdnz(jz) = -R ( d Z n_Z / d r ) / ( Z n_Z )
c  grdte(jz) = -R ( d T_e / d r ) / T_e
c  grdti(jz) = -R ( d T_i / d r ) / T_i
c  grdq (jz) =  R ( d q   / d r ) / q    related to magnetic shear
c
c  where:
c    n_i     = thermal ion density (sum over hydrogenic and impurity)
c    n_h     = thermal hydrogenic density (sum over hydrogenic species)
c    n_Z     = thermal impurity density,  Z = average impurity charge
c                      sumed over all impurities
c
c  Output:
c  -------
c
c    The following effective diffusivities represent contributions
c    to the total diffusivity matrix (difthi and velthi given below)
c    from each of the models that contribute to the Multi-Mode model.
c    Generally, these arrays are used for diagnostic output only.
c
c  thiig(jz) = ion thermal diffusivity from the Weiland model
c  thdig(jz) = hydrogenic ion diffusivity from the Weiland model
c  theig(jz) = elelctron thermal diffusivity from the Weiland model
c  thzig(jz) = impurity ion diffusivity from the Weiland model
c	    
c  thirb(jz) = ion thermal diffusivity from resistive ballooning modes
c  thdrb(jz) = hydrogenic ion diffusivity from resistive ballooning modes
c  therb(jz) = elelctron thermal diffusivity from resistive ballooning modes
c  thzrb(jz) = impurity ion diffusivity from resistive ballooning modes
c	    
c  thikb(jz) = ion thermal diffusivity from kinetic ballooning modes
c  thdkb(jz) = hydrogenic ion diffusivity from kinetic ballooning modes
c  thekb(jz) = elelctron thermal diffusivity from kinetic ballooning modes
c  thzkb(jz) = impurity ion diffusivity from kinetic ballooning modes
c
c    The following are growth rates and mode frequencies from the
c    Weiland model for drift modes such as ITG and TEM.
c    These arrays are intended for diagnostic output.
c
c  gamma(jm,jz) = growth rate for mode jm at point jz ( 1/sec )
c  omega(jm,jz) = frequency for mode jm at point jz ( radians/sec )
c
c    All of the transport coefficients are given in the following two
c    matricies for diffusion difthi and convection velthi in MKS units.
c    See the LaTeX documentation for difthi and velthi just below.
c
c    NOTE:  difthi and velthi include all of the anomalous transport.
c    There are no additional contributions to the heat fluxs from
c    charged particle convection.
c
c  difthi(j1,j2,jz) = full matrix of anomalous transport diffusivities
c  velthi(j1,jz)    = convective velocities
c  vflux(j1,jz)     = flux matrix
!| 
!| The full matrix form of anomalous transport has the form
!| $$ \frac{\partial}{\partial t}
!|  \left( \begin{array}{c} n_H T_H  \\ n_H \\ n_e T_e \\
!|     n_Z \\ n_Z T_Z \\ \vdots
!|     \end{array} \right)
!|  = - \nabla \cdot
!| \left( \begin{array}{l} {\rm vFlux}_1 \; n_H T_H \\
!|  {\rm vFlux}_2 \; n_H \\
!|  {\rm vFlux}_3 \; n_e T_e \\
!|  {\rm vFlux}_4 \; n_Z \\
!|  {\rm vFlux}_5 \; n_Z T_Z \\
!|  \vdots \end{array} \right) 
!|  + \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  = \nabla \cdot
!| \left( \begin{array}{llll}
!| D_{1,1} n_H & D_{1,2} T_H & D_{1,3} n_H T_H / T_e \\
!| D_{2,1} n_H / T_H & D_{2,2} & D_{2,3} n_H / T_e \\
!| D_{3,1} n_e T_e / T_H & D_{3,2} n_e T_e / n_H & D_{3,3} n_e & \vdots \\
!| D_{4,1} n_Z / T_H & D_{4,2} n_Z / n_H & D_{4,3} n_Z / T_e \\
!| D_{5,1} n_Z T_Z / T_H & D_{5,2} n_Z T_Z / n_H &
!|         D_{5,3} n_Z T_Z / T_e \\
!|  & \ldots & & \ddots
!| \end{array} \right)
!|  \nabla
!|  \left( \begin{array}{c}  T_H \\ n_H \\  T_e \\
!|    n_Z \\  T_Z \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  + \nabla \cdot
!| \left( \begin{array}{l} {\bf v}_1 \; n_H T_H \\ {\bf v}_2 \; n_H \\
!|    {\bf v}_3 \; n_e T_e \\
!|    {\bf v}_4 \; n_Z \\ {\bf v}_5 \; n_Z T_Z \\
!|     \vdots \end{array} \right) +
!|  \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right) $$
!| Note that all the diffusivities are in units of m$^2$/sec while the
!| convective velocities and vfluxes are in units of m/sec.
!| 
!| WARNING:  Do not add separate convective transport terms to this
!| anomalous transport model.  All the anomalous transport 
!| predicted by this Multi-Mode model is contained
!| in the diffusion coefficients {\tt difthi} and {\tt velthi} given
!| above.
!| 
c
c  Input integers:
c  ---------------
c
c  matdim  = first and second dimension of transport matricies
c            difthi(j1,j2,jz) and the first dimension of 
c            velthi(j1,jz), vflux(j1,jz), gamma(j1,jz), and omega(j1,jz).
c            matdim must be at least 5
c
c  npoints = number of values of jz in all of the above arrays
c
c  nprout  = output unit number for long printout
c
c
c  Input switches
c  --------------
c
c  lprint      controls the amount of printout (0 => no printout)
c              higher values yield more diagnostic output
c
c  lsuper   = 0 for simulations of all other discharges
c           > 0 for supershot simulations; substantially reduces 
c               contribution from kinetic ballooning mode
c
c  lreset  = 0 to use only internal settings for lswitch, cswitch
c              and for the coefficients fig, frb, and fkb that control
c              the contributions form the various instability modes
c
c    Note that when lreset = 0, the values of the switches and
c    coefficients in the argument list are ignored and all the 
c    switches and coefficients are set internally.
c
c    WARNING:  use lreset > 0 only if you want to pass all the switches
c              lswitch, cswitch, fig, frb, and fkb through the 
c              argument list.
c
c    WARNING:  NTCC users desiring to use anything other than lreset = 0
c              should consult with the mmm95 code authors first.
c
c
c  Output Integer
c  --------------
c
c  nerr        status code returned; 0 = OK; .ne. 0 indicates error
c
c
c  Internal control variables:
c  ---------------------------
c
c  lswitch(j), j=1,8   integer control variables: 
c
c  cswitch(j), j=1,25   general control variables:
c
c  lswitch(1)  controls which version of the Weiland model is used
c                  Default lswitch(1) = 10
c             = 2  2 eqn  Weiland model Hydrogen \eta_i mode only
c             = 4  4 eqn  Weiland model with Hydrogen and trapped electrons
c             = 5  5 eqn  Weiland model with trapped electrons, FLR effects, 
c                         and parallel ion motion
c             = 6  6 eqn  Weiland model Hydrogen, trapped electrons,
c                    and one impurity species
c             = 7  7 eqn   Weiland model Hydrogen, trapped electrons,
c                  one impurity species, and collisions
c             = 8  8 eqn  Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, and parallel
c                  ion (hydrogenic) motion
c             = 9  9 eqn  Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, and finite beta
c             = 10 10 eqn Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, parallel
c                  ion (hydrogenic) motion, and finite beta
c             = 11 11 eqn Weiland model Hydrogen, trapped electrons,
c                  one impurity species, collisions, parallel
c                  ion (hydrogenic, impurity) motion, and finite beta
c
c  lswitch(2) = 0  full matrix representation for difthi and velthi
c                  Default lswitch(2) = 2
c             = 1  set diagonal matrix elements difthi and velthi
c             = 2  set diagonal matrix elements = effective diffusivities
c
c  lswitch(3)  controls \kappa scaling
c                  Default lswitch(3) = 0
c             = 0  use \kappa scaling raised to
c                  exponents (cswitch(3) - cswitch(5))
c             = 1  use (1+\kappa^2)/2 instead of \kappa scaling
c               
c  lswitch(4) > 0  to replace negative diffusivity with velocity
c                  Default lswitch(4) = 1
c
c  lswitch(5) = 1  to limit magnitude of all normalized gradients
c                  to ( major radius ) / ( ion Larmor radius )
c                  Default lswitch(5) = 1
c
c  cswitch(1)   0.5  minimum value of shear
c  cswitch(2)   3.5  coeff in exponential (fbeta-th) of kinetic ballooning model
c  cswitch(3)  -4.0  exponent of local elongation multiplying drift waves
c  cswitch(4)  -4.0  exponent of local elongation multiplying resistive
c                     ballooning modes
c  cswitch(5)  -4.0  exponent of local elongation multiplying
c                     kinetic balllooning modes
c  cswitch(6)   0.0  k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
c  cswitch(8)   1.0  coeff of beta_prime_1 in kinetic ballooning mode
c  cswitch(9)  0.15  alpha in diamagnetic stabil. in kinetic ballooning model
c  cswitch(10)  0.0  rel fract of ion thermal diffusivity given to convection 
c  cswitch(11)  0.0  rel fract of hydrogen particle diffusivity given to convection 
c  cswitch(12)  0.0  rel fract of el thermal diffusivity given to convection 
c  cswitch(13)  0.0  rel fract of impurity particle diffusivity given to convection 
c  cswitch(14)  1.0  coef of finite beta effect in weiland14 = cetain(20) 
c  cswitch(15)  0.0  min value of impurity charge state zimpz
c  cswitch(16)  0.0  coef of fast particle fraction (superthermal ions) 
c                    in weiland model -- coef of densfe
c  cswitch(17)  1.0  coeff of k_\parallel (parallel ion motion) in weiland14 
c                    = cetain(10)
c  cswitch(18)  0.0  coeff of nuhat (effect of collisions) in weiland14 
c                    = cetain(15)
c  cswitch(19)  0.0  coeff for including v_parallel in strong ballooning limit
c                    = cetain(12); cswitch(19) = 1 for inclusion of v_par effect
c  cswitch(20)  0.0  trapping fraction used in weiland14 (when > 0.0)
c                    multiplies electron trapping fraction when < 0.0
c                    no effect when cswitch(20) = 0.0
c  cswitch(21)  1.0  multiplier for wexbs (flow shear rate) in Weiland model
c  cswitch(22)  0.0  ranges from 0.0 to 1.0 adds impurity heat flow to total 
c                    ionic heat flow for the weiland model
c  cswitch(23)  0.0  controls finite diff to construct the zgm matrix 
c                    = cetain(30)
c
c     contributions to vfluxes and interchanges: 
c
c  fig(1)   hydrogen particle transport from ITG (eta_i) mode
c  fig(2)   electron thermal  transport from ITG (eta_i) mode
c  fig(3)   ion      thermal  transport from ITG (eta_i) mode
c  fig(4)   impurity particle transport from ITG (eta_i) mode
c
c  frb(1)   hydrogen particle transport from resistive ballooning mode
c  frb(2)   electron thermal  transport from resistive ballooning mode
c  frb(3)   ion      thermal  transport from resistive ballooning mode
c  frb(4)   impurity particle transport from resistive ballooning mode
c
c  fkb(1)   hydrogen particle transport from kinetic ballooning mode
c  fkb(2)   electron thermal  transport from kinetic ballooning mode
c  fkb(3)   ion      thermal  transport from kinetic ballooning mode
c  fkb(4)   impurity particle transport from kinetic ballooning mode
c
c
c***********************************************************************
c
c-----------------------------------------------------------------------
c
c  Compile this routine and routines that it calls with a compiler 
c  option, such as -r8, to convert real to double precision when used on 
c  workstations.
c
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: MMM95 calls the following routines
c
c  WEILAND14       - Computes diffusion matrix and convect velocities
c                        for the Weiland transport model
c    WEILAND14FLUX - Calculates fluxes and effective diffusivities
c      TOMSQZ      - Wrapper for QZ algorithm solving Ax = lambda Bx
c         CQZHES   - First step in QZ algorithm 
c         CQZVAL   - Second and third step in QZ algorithm
c         CQZVEC   - Fourth step in QZ algorithm
c
c-----------------------------------------------------------------------

      implicit none
c
      integer km, klswitch, kcswitch
c
      integer  matdim,  npoints, nprout,  lprint,   nerr
c
      integer  lsuper,  lreset,  lswitch(*)
c
      parameter ( km = 12, klswitch = 8, kcswitch = 25 )
c
      real
     &   rminor(*),  rmajor(*),   elong(*)
     & , dense(*),   densh(*),    densimp(*),  densfe(*)
     & , xzeff(*),   tekev(*),    tikev(*),    q(*),       btor(*)
     & , avezimp(*), amassimp(*), amasshyd(*), aimass(*),  wexbs(*)
     & , grdne(*),   grdni(*),    grdnh(*),    grdnz(*)
     & , grdte(*),   grdti(*),    grdq(*)
c
      real  
     &   thiig(*),   thdig(*),    theig(*),    thzig(*)
     & , thirb(*),   thdrb(*),    therb(*),    thzrb(*)
     & , thikb(*),   thdkb(*),    thekb(*),    thzkb(*)
     & , omega(matdim,*),         gamma(matdim,*)
     & , difthi(matdim,matdim,*), velthi(matdim,*)
     & , vflux(matdim,*)
c
      real     cswitch(*)
c
      real     fig(*),  fkb(*),  frb(*)
c
c..physical constants
c
      real zpi,  zcc,  zcmu0,  zceps0,  zckb,  zcme,  zcmp,  zce
c
c  zpi     = pi
c  zcc     = speed of light                  [m/sec]
c  zcmu0   = vacuum magnetic permeability    [henrys/m]
c  zceps0  = vacuum electrical permittivity  [farads/m]
c  zckb    = energy conversion factor        [Joule/keV]
c  zcme    = electron mass                   [kg]
c  zcmp    = proton mass                     [kg]
c  zce     = electron charge                 [Coulomb]
c
c..computer constants
c
      real  zepslon, zlgeps
c
c  zepslon = machine epsilon [smallest number so that 1.0+zepslon>1.0]
c  zlgeps  = ln ( zepslon )
c
c
c..local variables
c
      integer  jz, j1, j2, jm

      real  zelong,  zelonf,  zai,     zne,     zni,    zte,    zti
     & ,    zq,      zeff,    zgne,    zgni,    zgnh,   zgnz,   zgte
     & ,    zgth,    zgtz,    zshear,  zrmin,   zrmaj,  zbtor,  zep
     & ,    zgyrfi,  zbeta,   zvthe,   zvthi,   zsound, zlog,   zcf
     & ,    znuei,   znueff,  zlari,   zlarpo,  zrhos,  zwn,    znude
     & ,    znuhat,  zgpr,    zscyl,   zsmin,   zshat,  zgmax
c
c.. variables for Weiland model
c
c  iletai(j1) and cetain(j1) are control variables
c
      integer        iletai(32)
c
      real  cetain(32), zomega(km), zgamma(km), zchieff(km)
c
      real           zdfthi(km,km),    zvlthi(km),     zflux(km)
c
      integer        idim,    ieq,     imodes
c
      real  zthte,   zbetae,  znz,     zmass,  zimpz,  ztzte
     & ,    zfnzne,  zmzmh,   zfnsne,  zftrap, zkyrho, zomegde
     & ,    zwexb,   znormd,  znormv
c
c  zexb    = local copy of ExB shearing rate
c  znormd  = factor to convert normalized diffusivities
c  znormv  = factor to convert normalized convective velocities
c
c..local variables for resistive ballooning modes
c
      real  zgyrfe, zlare,   zgddia, zgdp
c
c..local variables for kinetic ballooning modes
c
      real  zbprim, zbcoef1, zbc1,   zelfkb,  zfbthn, zdk
c
c-----------------------------------------------------------------------
c
c..initialize imodes
	imodes  = 0
c..physical constants
c
        zpi     = atan2 ( 0.0, -1.0 )
        zcc     = 2.997925e+8
        zcmu0   = 4.0e-7 * zpi
        zceps0  = 1.0 / ( zcc**2 * zcmu0 )
        zckb    = 1.60210e-16
        zcme    = 9.1091e-31
        zcmp    = 1.67252e-27
        zce     = 1.60210e-19
c
c..computer constants
c
        zepslon = 1.0e-34
        zlgeps  = log ( zepslon )
c
c
c..initialize arrays
c
      do jz = 1, npoints
        thiig(jz)  = 0.
        thdig(jz)  = 0.
        theig(jz)  = 0.
        thzig(jz)  = 0.
        therb(jz)  = 0.
        thirb(jz)  = 0.
        thdkb(jz)  = 0.
        thekb(jz)  = 0.
        thzkb(jz)  = 0.
        thikb(jz)  = 0.
        thdkb(jz)  = 0.
        thekb(jz)  = 0.
        thzkb(jz)  = 0.
      enddo
c
      do jz = 1, npoints
        do j1 = 1, matdim
          velthi(j1,jz) = 0.0
          vflux(j1,jz) = 0.0
          gamma(j1,jz) = 0.0
          omega(j1,jz) = 0.0
          do j2 = 1, matdim
            difthi(j1,j2,jz) = 0.0
          enddo
        enddo
      enddo
c
      nerr = 0
c
c..if lreset < 1, use internal settings for switches and coefficients
c  otherwise, use values passed through the argument list above
c
      if ( lreset .lt. 1 ) then
c
c..initialize switches
c
      do j1=1,kcswitch
        cswitch(j1) = 0.0
      enddo
c
      do j1=1,klswitch
        lswitch(j1) = 0
      enddo
c
c
c  Multi Mode Model in sbrtn THEORY version MMM95
c  for use in the BALDUR transport code
c
      lswitch(1) = 10 ! Weiland ITG model weiland14 (10 eqns, no collisions)
      lswitch(2) = 2  ! use effective diffusivities
      lswitch(3) = 0  ! use kappa instead of (1+\kappa^2)/2
      lswitch(4) = 1  ! replace -ve diffusivity with convective velocity
      lswitch(5) = 1  ! limit gradients by major radius / ion Larmor radius
c
c  misc. parameters for subroutine mmm95
c
      cswitch(1)  =  0.5  ! minimum value of shear
      cswitch(2)  =  3.5  ! coef in exponential (fbeta-th) in kinetic ballooning
      cswitch(3)  = -4.0  ! elongation scaling for drift wave mode
      cswitch(4)  = -4.0  ! elongation scaling for resistive ballooning mode
      cswitch(5)  = -4.0  ! elongation scaling for kinetic ballooning mode
      cswitch(6)  =  0.0  ! k_y \rho_s (= 0.316 if abs(cswitch(6)) < zepslon)
      cswitch(8)  =  1.0  ! coeff of beta_prime_1 in kinetic ballooning mode
      cswitch(9)  = 0.15  ! alpha in diamagnetic stabil. in kinetic balloon model 
      cswitch(10) =  0.0  ! relative fraction of ion thermal diffusivity 
                          ! given to convection
      cswitch(11) =  0.0  ! relative fract of hydrogen particle diffusivity 
                          ! given to convection
      cswitch(12) =  0.0  ! relative fraction of el thermal diffusivity 
                          ! given to convection
      cswitch(13) =  0.0  ! relative fract of impurity particle diffusivity 
                          ! given to convection
      cswitch(14) =  1.0  ! coef of finite beta effect in weiland14 = cetain(20)
      cswitch(15) =  0.0  ! min value of impurity charge state zimpz
      cswitch(16) =  1.0  ! coef of fast particle fraction (superthermal ions) 
                          ! in weiland model -- coef of densfe
      cswitch(17) =  1.0  ! coeff of k_\parallel (parallel ion motion) in 
                          ! weiland14 = cetain(10)
      cswitch(18) =  0.0  ! coeff of nuhat in weiland14 = cetain(15)
      cswitch(19) =  0.0  ! coeff for including v_parallel in strong ballooning
                          ! limit = cetain(12); cswitch(19) = 1 for inclusion 
                          ! of v_par effect
      cswitch(20) =  0.0  ! trapping fraction used in weiland14 (when > 0.0)
                          ! multiplies electron trapping fraction when < 0.0
                          ! no effect when cswitch(20) = 0.0
      cswitch(21) =  1.0  ! multiplier for wexbs (flow shear rate) 
                          ! in Weiland model
      cswitch(22) =  0.0  ! multiplier to impurity heat flux
      cswitch(23) =  0.0  ! controls finite diff to construct the 
                          ! zgm matrix = cetain(30)

c  contributions to hydrogenic particle, elec-energy, ion-energy,
c    and impurity ion fluxes

        fig(1) = 0.80
        fig(2) = 0.80
        fig(3) = 0.80
        fig(4) = 0.80
c
        fkb(1) = 1.00
       	fkb(2) = 0.65
        fkb(3) = 0.65
        fkb(4) = 1.00
c
        if ( lsuper .gt. 0 ) then
          fkb(1) = 0.045
          fkb(2) = 0.010
          fkb(3) = 0.010
          fkb(4) = 0.045
        endif
c
        frb(1) = 1.00
        frb(2) = 1.00
        frb(3) = 1.00
        frb(4) = 1.00
c
      endif

!| 
!| We then enter a loop over the spatial zones, and establish local
!| variables for all the input arrays for more compact notation.
!| 
c
c.. start the main do-loop over the radial index "jz"..........
c
c
      do 300 jz = 1, npoints

c avoid computation of fluxes at 
c rminor(jz) .lt. 1.e-4 * rmajor(jz)
c
      if (rminor(jz) .lt. (1.e-4 * rmajor(jz))) go to 300
c
c  transfer common to local variables to compact the notation
c
      zelong = max (zepslon,elong(jz))
      if ( lswitch(3) .eq. 1 ) then
        zelonf = ( 1. + zelong**2 ) / 2.
      else
        zelonf = zelong
      endif
c
      zai    = aimass(jz)
      zne    = dense(jz)
      zni    = densh(jz) + densimp(jz)
      znz    = densimp(jz)
      zte    = tekev(jz)
      zti    = tikev(jz)
      zq     = q(jz)
      zeff   = xzeff(jz)
c
c  normalized gradients
c
      zgne   = grdne(jz)
      zgni   = grdni(jz)
      zgnh   = grdnh(jz)
      zgnz   = grdnz(jz)
      zgte   = grdte(jz)
      zgth   = grdti(jz)
      zgtz   = grdti(jz)

      zrmin  = max( rminor(jz), zepslon )
      zrmaj  = rmajor(jz)
      zshear = grdq(jz) * zrmin / zrmaj
      zbtor  = btor(jz)
c
c  compute inverse aspect ratio
c
      zep    = max( zrmin/zrmaj, zepslon )
c
!| 
!| To complete the rest of the calculation we then compute various
!| quantities needed for the transport flux formulas (as in Table 1 of
!| the Comments paper, from which
!| $\omega_{ce}$ was inadvertantly omitted).
!| To begin with, we compute only quantities
!| which do not involve scale heights.
!| In the order in which they are computed, algebraic notation for
!| these quantities is:
!| $$ \omega_{ci}=eB_{o}/(m_{p}A_{i}) \eqno{\tt zgyrfi} $$
!| $$ \beta=(2\mu_{o}k_{b}/B_{o}^{2})(n_{e}T_{e}+n_{i}T_{i})
!|  \eqno{\tt zbeta} $$
!| $$ v_{e}=(2k_{b}T_{e}/m_{e})^{1/2} \eqno{\tt zvthe} $$
!| $$ v_{i}=(2k_{b}T_{i}/m_{p}A_{i})^{1/2} \eqno{\tt zvthi} $$
!| $$ c_{s}=[k_{b}T_{e}/(m_{p}A_{i})]^{1/2} \eqno{\tt zsound} $$
!| $$ \ln (\lambda)=37.8 - \ln (n_{e}^{1/2}T_{e}^{-1}) \eqno{\tt zlog} $$
!| $$ \nu_{ei}=4(2\pi)^{1/2}n_{e}(\ln \lambda)e^{4}Z_{eff}
!|                /[3(4\pi \epsilon_{o})^{2}m_{e}^{1/2}(k_{b}T_{e})^{3/2}]
!|  \eqno{\tt znuei} $$
!| $$ \eta=\nu_{ei}/(2\epsilon_{o}\omega_{pe}^{2}) \eqno{\tt zresis} $$
!| $$ \nu_{eff}=\nu_{ei}/\epsilon \eqno{\tt znueff} $$
!| $$ \nu_{e}^{*}=\nu_{ei}qR_{o}/(\epsilon^{3/2}v_{e}) \eqno{\tt thnust} $$
!| $$ \hat{\nu}=\nu_{eff}/\omega_{De} \eqno{\tt znuhat} $$
!| $$ \rho_{\theta i}=\rho_{i}q/\epsilon \eqno{\tt zlari} $$
!| $$ \rho_{i}=v_{i}/\omega_{ci} \eqno{\tt zlarpo} $$
!| $$ \rho_{s}=c_{s}/\omega_{ci} \eqno{\tt zrhos} $$
!| $$ k_{\perp}=0.3/\rho_{s} \eqno{\tt zwn} $$
!| 
!| The corresponding coding is:
!| 
c
      zgyrfi = zce * zbtor / (zcmp * zai)
      zbeta  = (2. * zcmu0 * zckb / zbtor**2) * (zne * zte + zni * zti)
      zvthe  = sqrt(2. * zckb * zte / zcme)
      zvthi  = sqrt(2. * zckb * zti / (zcmp * zai))
      zsound = sqrt(zckb * zte / (zcmp * zai))
      zlog   = 37.8-log(sqrt(zne) / zte)
      zcf    = (4. * sqrt(zpi) / 3.)
      zcf    = zcf * (zce / (4. * zpi * zceps0))**2
      zcf    = zcf * (zce / zckb) * sqrt( (zce/zcme) * (zce/zckb) )
      znuei  = zcf * sqrt(2.) * zne * zlog * zeff / (zte * sqrt(zte))
c
      znueff = znuei / zep
      zlari  = zvthi / zgyrfi
      zlarpo = max(zlari * zq / zep, zepslon)
      zrhos  = zsound / zgyrfi
      zwn    = 0.3 / zrhos
      znude  = 2 * zwn * zrhos * zsound / zrmaj
      znuhat = znueff / znude
c
c..if lswitch(5) = 1, limit magnitude of normalized gradients
c                    to ( major radius ) / ( ion Larmor radius )
c
      zgmax = zrmaj / zlarpo
c
      if ( lswitch(5) .eq. 1 ) then
c
        zgne = sign ( min ( abs ( zgne ), zgmax ), zgne )
        zgni = sign ( min ( abs ( zgni ), zgmax ), zgni )
        zgnh = sign ( min ( abs ( zgnh ), zgmax ), zgnh )
        zgnz = sign ( min ( abs ( zgnz ), zgmax ), zgnz )
        zgte = sign ( min ( abs ( zgte ), zgmax ), zgte )
        zgth = sign ( min ( abs ( zgth ), zgmax ), zgth )
        zgtz = sign ( min ( abs ( zgtz ), zgmax ), zgtz )
c
      endif
c
c  zgpr = -R ( d p   / d r ) / p    for thermal pressure
c
c  Compute the pressure scale length using smoothed and bounded
c  density and temperature
c
      zgpr = ( zne * zte * ( zgne + zgte )
     &         + zni * zti * ( zgni + zgth ) )
     &         / ( zne * zte + zni * zti )
c
      if ( lswitch(5) .eq. 1 )
     &  zgpr = sign ( min ( abs ( zgpr ), zgmax ), zgpr )
c
c
!| 
!| Our formulas for the shear begin with
!| $$ {\hat s}_{cyl}=|(r/q)(\partial q/\partial r)| \eqno{\tt zscyl} $$
!| computed earlier in this subroutine.
!| The minimum prescribed shear is 
!| $$ {\hat s}_{min}=max(c_{1},0) \eqno{\tt zsmin} $$
!| where $c_1=$ = {\tt cswitch(1)} so that shear is then given by
!| $$ {\hat s}=max({\hat s}_{min},{\hat s}_{cyl}) \eqno{\tt zshat} $$
!| 
!| The relevant coding for the calculations just described is:
!| 
c
      zscyl=max(abs(zshear),zepslon)
      zsmin=max(cswitch(1),zepslon)
      zshat=max(zsmin,zscyl)
c
!| 
!| %**********************************************************************c
!| 
!| \section{Transport Models}
!| 
!| The computation of the anomalous transport coefficients is
!| now described.  Please note that all the heat flux is included in the
!| thermal diffusion and velocity coefficients.  There are no additional
!| ``convective velocities''.
!| The mode abbreviations used here are
!| \begin{center}
!| \begin{tabular}{llll}
!|     &             &                                         &        \\
!|     & {\tt ig}    & $\eta_i$-mode and drift wave modes      &        \\
!|     & {\tt rb}    & resistive ballooning                    &        \\
!|     & {\tt kb}    & kinetic ballooning                      &        \\
!|     &             &                                         &
!| \end{tabular}
!| \end{center}
!| 
!| %**********************************************************************c
!| 
!| \subsection{$\eta_i$ Modes}
!| %%%%%
!| 
!| The $\eta_i$ and trapped electron mode model 
!| by Weiland et al\cite{nord90a} is implemented when
!| ${\tt lswitch(1)}$ is set greater than 1.
!| When $ {\tt lswitch(1)} = 2 $, only the hydrogen equations are used
!| (with no trapped electrons or impurities) to compute only the 
!| $ \eta_i $ mode.
!| When $ {\tt lswitch(1)} = 4 $, trapped electrons are included,
!| but not impurities.
!| When $ {\tt lswitch(1)} = 6 $, a single species of impurity ions is
!| included as well as trapped electrons.
!| When $ {\tt lswitch(1)} = 7 $, the effect of collisions is included.
!| When $ {\tt lswitch(1)} = 8 $, parallel ion (hydrogenic) motion and 
!| the effect of collisions are included.
!| When $ {\tt lswitch(1)} = 9 $, finite beta effects and collisions are
!| included.
!| When $ {\tt lswitch(1)} = 10 $, parallel ion (hydrogenic) motion, 
!| finite beta effects, and the effect of collisions are included.
!| When $ {\tt lswitch(1)} = 11 $, parallel ion (hydrogenic and impurity) motion, 
!| finite beta effects, and the effect of collisions are included.
!| Finite Larmor radius corrections are included in all cases.
!| Values of {\tt lswitch(1)} greater than 11 are reserved for extensions
!| of this Weiland model.
!| 
!| The mode growth rate, frequency, and effective diffusivities are
!| computed in subroutine {\tt weiland14}.
!| Frequencies are normalized by $\omega_{De}$ and diffusivities are
!| normalized by $ \omega_{De} / k_y^2 $.
!| The order of the diffusivity equations is 
!| $ T_H $, $ n_H $, $ T_e $, $ n_Z $, $ T_Z $, \ldots
!| Note that the effective diffusivities can be negative.
!| 
!| The diffusivity matrix $ D = {\tt difthi(j1,j2)}$ 
!| is given above.
!| 
!| 
!| The impurity density gradient scale length is defined as 
!| $$g_{nz}=-R{{d\ }\over {dr}}\left(Zn_z\right)/(Zn_z)$$
!| The electron density gradient scale length is defined as
!| $$g_{ne}=(1-Zf_z-f_s)g_{nH}+Zf_zg_{nz}+f_sg_{ns}$$
!| where $ f \equiv n_Z / n_e $ and $ n_e = n_H + Z n_Z +n_s$.
!| For this purpose, all the impurity species are lumped together as 
!| one effective impurity species and all the hydrogen isotopes are lumped 
!| together as one effective hydrogen isotope.
!| 
!| 
c
        do j1=1,32
          iletai(j1) = 0
          cetain(j1) = 0.0
        enddo
c
        thiig(jz) = 0.0
        theig(jz) = 0.0
        thdig(jz) = 0.0
        thzig(jz) = 0.0
c
c..set the number of equations to use in the Weiland model
c
        if ( (lswitch(1) .lt. 2) .or. (lswitch(1) .gt. 11 )) then
          nerr = -10
          return
        elseif (lswitch(1) .eq. 3) then
	  nerr = -10
          return
        else
          ieq = lswitch(1)
        endif
c
        cetain(11) = 1.0
c
c.. coefficient of k_parallel for parallel ion motion
c.. cswitch(19) for v_parallel in strong ballooning limit
c.. in 9 eqn model
c
        cetain(10) = cswitch(17)
        cetain(12) = cswitch(19)
        cetain(15) = cswitch(18)
        cetain(20) = cswitch(14)
c
        iletai(10) = 0
c
        idim   = km
c
c  Hydrogen species
c
        zthte  = zti / zte

        zbetae = 2. * zcmu0 * zckb * zne * zte / zbtor**2
c
c  Impurity species (use only impurity species 1 for now)
c  assume T_Z = T_H throughout the plasma here
c
        znz    = densimp(jz)
        zmass  = amassimp(jz)
        zimpz  = avezimp(jz)
        zimpz  = max ( zimpz, cswitch(15) )
c
        ztzte  = zti / zte
        zfnzne = znz / zne
        zmzmh  = zmass / amasshyd(jz)
c
c  superthermal ions
c
c  zfnsne = ratio of superthermal ions to electrons
c  L_ns   = gradient length of superthermal ions
c
        zfnsne = max ( cswitch(16) * densfe(jz) / dense(jz) , 0.0 )
c
        zftrap = sqrt ( 2. * zrmin / ( zrmaj * ( 1. + zrmin / zrmaj )))
        if ( cswitch(20) .gt. zepslon ) zftrap = cswitch(20)
        if ( cswitch(20) .lt. -zepslon )
     &       zftrap = abs(cswitch(20))*zftrap
c
        if ( abs(cswitch(6)) .lt. zepslon ) then
          zkyrho = 0.316
        else
          zkyrho = cswitch(6)
        endif
c
c
c...Define a local copy of normalized ExB shearing rate : pis
c
        zomegde = 2.0 * zkyrho * zsound / zrmaj 
c
        zwexb = cswitch(21) * wexbs(jz) / zomegde 
c
c
        cetain(30) = cswitch(23)
        iletai(6)  = 0
        if ( lswitch(2) .lt. 1 ) iletai(7) = 1
c
c  if lswitch(2) .lt. 1, compute only the effective diffusivities
c
        iletai(9) = lswitch(2)
c
        call weiland14 ( 
     &   iletai,   cetain,   lprint,   ieq,      nprout,   zgne
     & , zgnh,     zgnz,     zgte,     zgth,     zgtz,     zthte
     & , ztzte,    zfnzne,   zimpz,    zmzmh,    zfnsne,   zbetae
     & , zftrap,   znuhat,   zq,       zshat,    zkyrho,   zwexb
     & , idim,     zomega,   zgamma,   zdfthi,   zvlthi,   zchieff
     & , zflux,    imodes,   nerr )
c
c  If nerr not equal to 0 an error has occured
c
	if (nerr .ne. 0) return
c
c
c  Growth rates for diagnostic output
c    Note that all frequencies are normalized by \omega_{De}
c      consequently, trapped electron modes rotate in the positive
c      direction (zomega > 0) while eta_i modes have zomega < 0.
c
        jm = 0
        do j1=1,imodes
          if ( zgamma(j1) .gt. zepslon ) then
            jm = jm + 1
            gamma(jm,jz) = zgamma(j1) * zomegde
            omega(jm,jz) = zomega(j1) * zomegde
          endif
        enddo
c
c  conversion factors for diffusion and convective velocity
c
        znormd = zelonf**cswitch(3) *
     &    2.0 * zsound * zrhos**2 / ( zrmaj * zkyrho )
        znormv = zelonf**cswitch(3) *
     &    2.0 * zsound * zrhos**2 / ( zrmaj**2 * zkyrho )
c
c  compute effective diffusivites for diagnostic purposes only
c
        thdig(jz) = fig(1) * znormd * zchieff(2)
        theig(jz) = fig(2) * znormd * zchieff(3)
        thiig(jz) = fig(3) * znormd * zchieff(1)
     &  + fig(3) * znormd * zchieff(5) * cswitch(22) * znz / zni
        thzig(jz) = fig(4) * znormd * zchieff(4)
c
c  start computing the fluxes
c
        vflux(1,jz) = vflux(1,jz) + thiig(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdig(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + theig(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzig(jz) * zgnz / zrmaj
c
c  compute diffusivity matrix
c
        do j1=1,matdim
          velthi(j1,jz) = 0.0
          vflux(j1,jz) = 0.0
          do j2=1,matdim
            difthi(j1,j2,jz) = 0.0
          enddo
        enddo
c
c..set difthi and velthi
c
        if ( lswitch(2) .gt. 1 ) then
c
c  diagonal elements of matrix = effective diffusivities
c
          difthi(1,1,jz) = difthi(1,1,jz) + thiig(jz)
          difthi(2,2,jz) = difthi(2,2,jz) + thdig(jz)
          difthi(3,3,jz) = difthi(3,3,jz) + theig(jz)
          difthi(4,4,jz) = difthi(4,4,jz) + thzig(jz)
c
        else
c
c..full matrix form of model
c
          if ( ieq .eq. 2 ) then
            difthi(1,1,jz) = difthi(1,1,jz) +
     &        fig(3) * znormd * zdfthi(1,1)
            velthi(1,jz)   = velthi(1,jz) +
     &        fig(3) * znormv * zvlthi(1)
          elseif ( ieq .eq. 4 ) then
            do j2=1,3
              difthi(1,j2,jz) = difthi(1,j2,jz) +
     &          fig(3) * znormd * zdfthi(1,j2)
              difthi(2,j2,jz) = difthi(2,j2,jz) +
     &          fig(1) * znormd * zdfthi(2,j2)
              difthi(3,j2,jz) = difthi(3,j2,jz) + 
     &          fig(2) * znormd * zdfthi(3,j2)
              difthi(4,j2,jz) = difthi(4,j2,jz) + 
     &          fig(4) * znormd * zdfthi(2,j2)
            enddo
              velthi(1,jz)    = velthi(1,jz) +
     &          fig(3) * znormv * zvlthi(1)
              velthi(2,jz)    = velthi(2,jz) +
     &          fig(1) * znormv * zvlthi(2)
              velthi(3,jz)    = velthi(3,jz) +
     &          fig(2) * znormv * zvlthi(3)
              velthi(4,jz)    = velthi(4,jz) +
     &          fig(4) * znormv * zvlthi(2)
          else
            do j2=1,4
              difthi(1,j2,jz) = difthi(1,j2,jz) +
     &          fig(3) * znormd * zdfthi(1,j2)
              difthi(2,j2,jz) = difthi(2,j2,jz) +
     &          fig(1) * znormd * zdfthi(2,j2)
              difthi(3,j2,jz) = difthi(3,j2,jz) +
     &          fig(2) * znormd * zdfthi(3,j2)
              difthi(4,j2,jz) = difthi(4,j2,jz) +
     &          fig(4) * znormd * zdfthi(4,j2)
            enddo
              velthi(1,jz)    = velthi(1,jz) +
     &          fig(3) * znormv * zvlthi(1)
              velthi(2,jz)    = velthi(2,jz) +
     &          fig(1) * znormv * zvlthi(2)
              velthi(3,jz)    = velthi(3,jz) +
     &          fig(2) * znormv * zvlthi(3)
              velthi(4,jz)    = velthi(4,jz) +
     &          fig(4) * znormv * zvlthi(4)
          endif
c
        endif
c
c
c..transfer from diffusivity to convective velocity
c
        if ( lswitch(4) .gt. 0 ) then
c
          if ( thiig(jz) .lt. 0.0 ) then
            velthi(1,jz) = velthi(1,jz) - thiig(jz) * zgth / zrmaj
            thiig(jz) = 0.0
            do j2=1,4
              difthi(1,j2,jz) = 0.0
            enddo
          endif
c
          if ( thdig(jz) .lt. 0.0 ) then
            velthi(2,jz) = velthi(2,jz) - thdig(jz) * zgnh / zrmaj
            thdig(jz) = 0.0
            do j2=1,4
              difthi(2,j2,jz) = 0.0
            enddo
          endif
c
          if ( theig(jz) .lt. 0.0 ) then
            velthi(3,jz) = velthi(3,jz) - theig(jz) * zgte / zrmaj
            theig(jz) = 0.0
            do j2=1,4
              difthi(3,j2,jz) = 0.0
            enddo
          endif
c
          if ( thzig(jz) .lt. 0.0 ) then
            velthi(4,jz) = velthi(4,jz) - thzig(jz) * zgnz / zrmaj
            thzig(jz) = 0.0
            do j2=1,4
              difthi(4,j2,jz) = 0.0
            enddo
          endif
c
        else
c
c..shift from diffusion to convective velocity
c
          if ( abs(cswitch(10)) + abs(cswitch(11)) + abs(cswitch(12))
     &       + abs(cswitch(13)) .gt. zepslon ) then
c
            velthi(1,jz) = velthi(1,jz)
     &       + cswitch(10) * thiig(jz) * zgth / zrmaj
            velthi(2,jz) = velthi(2,jz)
     &       + cswitch(11) * thdig(jz) * zgnh / zrmaj
            velthi(3,jz) = velthi(3,jz)
     &       + cswitch(12) * theig(jz) * zgte / zrmaj
            velthi(4,jz) = velthi(4,jz)
     &       + cswitch(13) * thzig(jz) * zgnz / zrmaj
c
c..alter the effective diffusivities 
c  if they are used for more than diagnostic purposes
c
            thiig(jz) = ( 1.0 - cswitch(10) ) * thiig(jz)
            thdig(jz) = ( 1.0 - cswitch(11) ) * thdig(jz)
            theig(jz) = ( 1.0 - cswitch(12) ) * theig(jz)
            thzig(jz) = ( 1.0 - cswitch(13) ) * thzig(jz)
c
            do j2=1,4
              difthi(1,j2,jz) = ( 1.0 - cswitch(10) ) * difthi(1,j2,jz)
              difthi(2,j2,jz) = ( 1.0 - cswitch(11) ) * difthi(2,j2,jz)
              difthi(3,j2,jz) = ( 1.0 - cswitch(12) ) * difthi(3,j2,jz)
              difthi(4,j2,jz) = ( 1.0 - cswitch(13) ) * difthi(4,j2,jz)
            enddo
c
          endif
c
        endif
c
c..end of Weiland model
c
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
c
!| 
!| %**********************************************************************c
!| %%%%%
!| 
!| \subsection{Resistive Ballooning}
!| 
!| \subsubsection{Guzdar-Drake Drift-Resistive Ballooning Model}
!| 
!| The 1993 $\bf E \times B$ drift-resistive ballooning mode model by
!| Guzdar and Drake \cite{drake93} is selected
!| $$
!|   D^{RB} = F_1^{RB}
!|   \abs \left( 2\pi q_{a}^2 \right) \rho_e^2 \nu_{ei}\left( {R\over p}
!|   \right)\left(- {{d p }\over {d r}}\right)$$
!| where $F_1^{RB}=$ {\tt frb(1)}.  Here the normalize presure gradient
!| scale length has been substituted for the density gradient scale 
!| given in their paper following a comment made by Drake at the
!| 1995 TTF Workshop\cite{drakecom2}.  Including diamagnetic and
!| elongation effects, the particle diffusivity is
!| $$
!|   D^{RB} =  F_1^{RB} \abs \left( 2\pi q_{a}^2 \right) \rho_e^2 \nu_{ei} 
!|   \left( {R\over p}\right)\left(- {{d p }\over {d r}}\right)
!|   c_9 \kappa^{c_{4}}  \eqno{\tt zgddb}
!| $$
!| where $\rho_e=v_e/\omega_{ce}$, $c_9 =$ {\tt cswitch(9)} and 
!| and $c_4 =$ {\tt cswitch(4)}.
!| 
!| The electron  and ion thermal diffusivities are taken equal taken to be
!| an adjustable fraction of the particle diffusivity.
!| $$
!|  \chi_e^{RB} =  F_2^{RB}
!|     \abs \left( 2\pi q_{a}^2 \right) \rho_e^2 \nu_{ei} 
!|   \left( {R\over p}\right)\left(- {{d p }\over {d r}}\right)
!|   c_9 \kappa^{c_{4}}  \eqno{\tt therb}
!| $$
!| $$
!|  \chi_i^{RB} =  F_3^{RB}
!|     \abs \left( 2\pi q_{a}^2 \right) \rho_e^2 \nu_{ei} 
!|   \left( {R\over p}\right)\left(- {{d p }\over {d r}}\right)
!|   c_9 \kappa^{c_{4}}   \eqno{\tt thirb}
!| $$  where $F_2^{RB}=$ {\tt frb(2)} and $F_3^{RB}=$ {\tt frb(3)} 
c
c..Guzdar-Drake theory (Phys Fluids B 5 (1993) 3712
c..L_p used instead of L_n
c
        zgyrfe = zce * zbtor / zcme  ! electron plasma frequency
        zlare = zvthe / zgyrfe    ! electron Larmor radius
c
c..   Diamagnetic stabilization
c
          zgddia = cswitch(9)
c
c..   Diffusivities
c
        zgdp = abs ( 2. * zpi * ((zq * zlare)**2.) * znuei
     &    * zgpr * 100. * zgddia )

        thdrb(jz) = frb(1) * zgdp * zelonf**cswitch(4)
        therb(jz) = frb(2) * zgdp * zelonf**cswitch(4)
        thirb(jz) = frb(3) * zgdp * zelonf**cswitch(4)
        thzrb(jz) = frb(4) * zgdp * zelonf**cswitch(4)
c
c  add to the fluxes
c
        vflux(1,jz) = vflux(1,jz) + thirb(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdrb(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + therb(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzrb(jz) * zgnz / zrmaj
c
        difthi(1,1,jz) = difthi(1,1,jz) + thirb(jz)
        difthi(2,2,jz) = difthi(2,2,jz) + thdrb(jz)
        difthi(3,3,jz) = difthi(3,3,jz) + therb(jz)
        difthi(4,4,jz) = difthi(4,4,jz) + thzrb(jz)
c

!| %**********************************************************************c
!| %%%%%%
!| 
!| \subsection{Kinetic Ballooning}
!| 
!| For transport due to the kinetic ballooning mode, we compute $D^{KB}$ and
!| the thermal diffusivities in terms of the pressure gradient $\beta'$.
!| $$  D^{KB} = \frac{c_s \rho_s^2}{p}\left(-{{dp}\over{dr}}\right) f_{\beta th} \quad 
!| {\rm where}\qquad  f_{\beta th} = \exp \left[c_{2}\left( \frac{\beta '}
!| {\beta_{c1}'} -1 \right) \right] \eqno{\tt zdk}$$
!| and where $\beta_{cl}'$ is the ideal pressure gradient threshold for 
!| the onset of the ideal ballooning mode in $s-\alpha$ geometry, 
!| $$  \beta_{c1}' = c_{8} \hat{s}/(1.7 q^{2}R_{o}) \eqno{\tt zbc1} $$ \\[-2mm]
!| Here, $c_8=$ {\tt cswitch(8)=1}  by default, but is included for flexibility.
!| The coefficent $c_2$ in the expression for $f_{\beta th}$ is set equal to 
!| {\tt cswitch(2)}.  
!| 
!| The diffusivities are then given as:
!| $$ D_{a}^{KB}=D^{KB}F^{KB}_{1}  \kappa^{c_{5}} \eqno{\tt thdkb} $$
!| $$ Q_{e}^{KB}\frac{L_{Te}}{n_{e}T_{e}}=D^{KB}F^{KB}_{2} \kappa^{c_{5}}  \eqno{\tt thekb} $$
!| $$ Q_{i}^{KB}\frac{L_{Ti}}{n_{i}T_{i}}=D^{KB}F^{KB}_{3} \kappa^{c_{5}}
!| \eqno{\tt thikb} $$  
!| where $F_1^{KB}=$ {\tt fkb(1)}, $F_2^{KB}=$ {\tt fkb(2)},
!| $F_3^{KB}=$ {\tt fkb(3)}, and $c_5=$ {\tt cswitch(5)}. 
!|  Note that the new version does not include
!| the (5/2) factor in the thermal diffusivities.
!| 
!| \noindent
!| The relevant coding is:
!| 
c ..................................
c .  the kinetic ballooning model  .
c ..................................
c
c       zbprim and zbc1 computed above under drift model
c
      if (  abs(cswitch(2)) .gt. zepslon
     &   .and.  zgpr .gt. 0.0  ) then
c
      zbprim = abs( zbeta * zgpr / zrmaj )
      zbcoef1 = 1.0
      if ( abs(cswitch(8)) .gt. zepslon ) zbcoef1 = cswitch(8)
      zbc1   = zbcoef1 * abs(zshat)/(1.7*zq**2*zrmaj)
      zelfkb = zelonf**cswitch(5)
c
        zfbthn = exp( min(abs(zlgeps),
     &     max(-abs(zlgeps),cswitch(2)*(zbprim/zbc1 - 1.))) )
c
        zdk = abs( zsound * zrhos**2 * zfbthn * zgpr / zrmaj )
c
        thdkb(jz) = zdk*fkb(1)*zelfkb
        thekb(jz) = zdk*fkb(2)*zelfkb
        thikb(jz) = zdk*fkb(3)*zelfkb
        thzkb(jz) = zdk*fkb(4)*zelfkb
c
c  add to the fluxes
c
        vflux(1,jz) = vflux(1,jz) + thikb(jz) * zgth / zrmaj
        vflux(2,jz) = vflux(2,jz) + thdkb(jz) * zgnh / zrmaj
        vflux(3,jz) = vflux(3,jz) + thekb(jz) * zgte / zrmaj
        vflux(4,jz) = vflux(4,jz) + thzkb(jz) * zgnz / zrmaj
c
        difthi(1,1,jz) = difthi(1,1,jz) + thikb(jz)
        difthi(2,2,jz) = difthi(2,2,jz) + thdkb(jz)
        difthi(3,3,jz) = difthi(3,3,jz) + thekb(jz)
        difthi(4,4,jz) = difthi(4,4,jz) + thzkb(jz)
c
      endif
c
 300  continue
c
c
c   end of the main do-loop over the radial index, "jz"----------
c
      return
      end
!======================================================================|
!|  
!| %**********************************************************************c
!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{bate98a}
!| Glenn Bateman, Arnold~H. Kritz, Jon~E. Kinsey, Aaron~J. Redd, and Jan Weiland,
!| ``Predicting temperature and density profiles in tokamaks,''
!| {\em Physics of Plasmas,} {\bf 5} (1998) 1793--1799.
!| 
!| \bibitem{drake93}
!| P.~N. Guzdar, J.~F. Drake, D. McCarthy, A.~B. Hassam, and C.~S. Liu,
!| ``Three-dimensional Fluid Simulations of the Nonlinear Drift-resistive
!| Ballooning Modes in Tokamak Edge Plasmas,'' Physics of Fluids B {\bf 5}
!| (1993) 3712.
!| 
!| \bibitem{guzcomm} P.N. Guzdar, University of Maryland, personal communication
!|  (11/94).
!| 
!| %\bibitem{drakecom} J.F. Drake, University of Maryland, personal
!| %communication (6/94).
!| 
!| \bibitem{drakecom2} J.F. Drake, University of Maryland, personal
!| communication (3/95).
!| 
!| %\bibitem{Comments} C. E. Singer, ``Theoretical Particle and Energy
!| %Flux Formulas for Tokamaks,'' Comments on Plasma Physics and Controlled
!| %Fusion {\bf 11} (1988) 165.
!| 
!| \bibitem{nord90a} H. Nordman, J. Weiland, and A. Jarmen, 
!| ``Simulation of toroidal drift mode turbulence driven by 
!| temperature gradients and electron trapping,'' 
!| Nucl. Fusion {\bf 30} (1990) 983--996.
!| 
!| %\bibitem{Singer} C.E.Singer, G.Bateman, and D.D.Stotler,
!| %``Boundary Conditions for OH, L, and H-mode Simulations,''
!| %Princeton University Plasma Physics Report PPPL-2527 (1988).
!| 
!| \end{thebibliography}
!| 
!| %**********************************************************************c
!| \end{document}             % End of document.
!| %---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
!| \documentstyle{article}
!| \headheight 0pt \headsep 0pt          
!| \topmargin 0pt  \textheight 9.0in
!| \oddsidemargin 0pt \textwidth 6.5in
!| 
!| \newcommand{\Partial}[2]{\frac{\partial #1}{\partial #2}}
!| \newcommand{\jacobian}{{\cal J}}
!| 
!| \begin{document}
!| 
!| \begin{center} 
!| {\bf 
!| Weiland Model for Transport Driven by
!| Toroidal Ion Temperature Gradient and \\
!| Trapped Electron Modes  \\
!| {\tt weiland14.tex} \\
!| \vspace{1pc}
!| Glenn Bateman and Arnold Kritz \\
!| Lehigh University, Physics Department \\
!| 16 Memorial Drive East, Bethlehem, PA 18015 \\
!| \vspace{1pc}
!| Jan Weiland, Hans Nordman and P{\"a}r Strand\\
!| Department of Electromagnetics \\
!| Chalmers University of Technology \\
!| S-412 96 G\"{o}teborg, Sweden \\
!| \vspace{1pc}
!| Jon Kinsey \\
!| General Atomics \\
!| P.O. Box 85608, San Diego, CA 92186} \\ 
!| \vspace{1pc}
!| \end{center}
!| This subroutine evaluates the transport matrix for $\eta_i$ and trapped 
!| electron modes derived by Jan Weiland, H. Nordman and their group in 
!| G\"{o}teborg Sweden \cite{bate98a} -- \cite{jarm87a}.
!| The equations in this routine include fast Hydrogenic ions, impurities, 
!| trapped electron and finite Larmor radius effects. New options include 
!| parallel ion motion, finite beta and collisional effects together with an 
!| approximate treatment of ${\bf E}\times {\bf B}$ flow shear reduction of 
!| transport.
!| 
!| The dimensionless diffusivity matrix $ D_{j_1,j_2} = {\tt difthi(j1,j2)}$
!| and convective velocity array $ v_{j_1} = {\tt velthi(j1)} $
!| are given as:
!| $$ \frac{\partial}{\partial t}
!|  \left( \begin{array}{c} n_H T_H  \\ n_H \\ n_e T_e \\
!|     n_Z \\ n_Z T_Z \\ \vdots
!|     \end{array} \right)
!|  = - \nabla \cdot
!| \left( \begin{array}{l} {\rm vFlux}_1 n_H T_H \\
!|  {\rm vFlux}_2 n_H \\
!|  {\rm vFlux}_3 n_e T_e \\
!|  {\rm vFlux}_4 n_Z \\
!|  {\rm vFlux}_5 n_Z T_Z \\
!|  \vdots \end{array} \right) 
!|  + \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  = \nabla \cdot
!| \left( \begin{array}{llll}
!| D_{1,1} n_H & D_{1,2} T_H & D_{1,3} n_H T_H / T_e \\
!| D_{2,1} n_H / T_H & D_{2,2} & D_{2,3} n_H / T_e \\
!| D_{3,1} n_e T_e / T_H & D_{3,2} n_e T_e / n_H & D_{3,3} n_e & \vdots \\
!| D_{4,1} n_Z / T_H & D_{4,2} n_Z / n_H & D_{4,3} n_Z / T_e \\
!| D_{5,1} n_Z T_Z / T_H & D_{5,2} n_Z T_Z / n_H &
!|         D_{5,3} n_Z T_Z / T_e \\
!|  & \ldots & & \ddots
!| \end{array} \right)
!|  \nabla
!|  \left( \begin{array}{c}  T_H \\ n_H \\  T_e \\
!|    n_Z \\  T_Z \\ \vdots
!|     \end{array} \right)
!| $$
!| $$
!|  + \nabla \cdot
!| \left( \begin{array}{l} {\bf v}_1 n_H T_H \\ {\bf v}_2 n_H \\
!|    {\bf v}_3 n_e T_e \\
!|    {\bf v}_4 n_Z \\ {\bf v}_5 n_Z T_Z \\ \vdots \end{array} \right) +
!|  \left( \begin{array}{c} S_{T_H} \\ S_{n_H} \\ S_{T_e} \\
!|     S_{n_Z} \\ S_{T_Z} \\ \vdots
!|     \end{array} \right) $$
!| Note that all the diffusivities in this routine are normalized by
!| $ \omega_{De} / k_y^2 $,
!| convective velocities are normalized by $ \omega_{De} / k_y $,
!| and all the frequencies are normalized by $ \omega_{De} $.
!| \newpage
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine weiland14 (
     &   letain,   cetain,   lprint,   neq,      nout,     gnein
     & , gnhin,    gnzin,    gtein,    gthin,    gtzin,    tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    difthi,   velthi,   chieff
     & , vflux,    nmodes,   nerr )
c
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Names of variables in the argument list:}
!| \begin{tabular}{lllp{3.0in}}
!| variable & status & symbol & meaning \\
!| {\tt letain(j)} & input & & Integer control variables
!|                             (see table below).\\
!| {\tt cetain(j)} & input & & Real-valued control variables
!|                             (see table below).\\
!| {\tt lprint}    & input & & Controls printout. 
!|  Higher values produce more. \\
!| {\tt neq} & input & & number of equations \\
!| {\tt nout} & input & & output unit for error messages \\
!| {\tt gnein} & input & $  $ & $ - R\hat{r} \cdot \nabla n_e / n_e $ \\
!| {\tt gnhin} & input & $  $ & $ - R\hat{r} \cdot \nabla n_H / n_H $ \\
!| {\tt gnzin} & input & $  $ & $ - R\hat{r} \cdot \nabla n_Z / n_Z $ \\
!| {\tt gtein} & input & $  $ & $ - R\hat{r} \cdot \nabla T_e / T_e $ \\
!| {\tt gthin} & input & $  $ & $ - R\hat{r} \cdot \nabla T_H / T_H $ \\
!| {\tt gtzin} & input & $  $ & $ - R\hat{r} \cdot \nabla T_Z / T_Z $ \\
!| {\tt tauhin} & input & $\tau_H$ & $ \tau_H = T_H / T_e $ \\
!| {\tt tauzin} & input & $\tau_Z$ & $ \tau_Z = T_Z / T_e $ \\
!| {\tt fnzin} & input & $ f_{nZ} $ & $f_{nZ} = n_Z / n_e $ \\
!| {\tt czin}  & input & $ Z $ & impurity charge number \\
!| {\tt azin}  & input & $ m_Z / m_H $
!|      & impurity mass to hydrogen isotope mass. \\
!|  & & & Note that ``hydrogen'' may include a deuterium or tritium mix. \\
!| {\tt fnsin}  & input & $ f_s $
!|   & $ f_s = n_s / n_e $ fraction of superthermal hydrogenic ions \\
!| {\tt betaein} & input & $\beta_e$ &
!|      $ = n_e T_e / ( B^2 / 2 \mu_0 ) $ \\
!| {\tt ftrapein}  & input & $f_{trap} $ &
!|     fraction of trapped electrons \\
!| {\tt vef}       & input & $ \nu_{th} / \omega_{De} $ &
!|      thermal collision frequency, normalized \\
!| {\tt q}         & input & $ q $ & magnetic q-value \\
!| {\tt shear}     & input & $ s $ & $ d \ln q / d \ln r $ \\
!| {\tt ekyrhoin}  & input & $ k_y \rho_s $ & normalized poloidal
!|                     wave number \\
!| {\tt wexb}      & input & $\omega_{E\times B}$& $E\times B$ shearing rate 
!| (normalized with $\omega_{D_e}$) \\
!| {\tt ndim} & input & & first dimension of the 2-D array difthi
!|                and the maximum number of unstable modes allowed \\
!| {\tt omega(j)}  & output & $\omega / \omega_{De} $ &
!|      real part of the frequencies normalized by $ \omega_{De} $ \\
!| {\tt gamma(j)}  & output & $\gamma / \omega_{De} $ &
!|      growth rates normalized by $ \omega_{De} $ \\
!| {\tt difthi(i,j)}      & output & $ D \omega_{De} / k_y^2 $
!|       & diffusivity matrix normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt velthi(j)}      & output & $ v \omega_{De} / k_y $
!|       & convective velocities normalized by $ k_y / \omega_{De} $ \\
!| {\tt chieff(j)} & output & $ \chi_{\rm eff} \omega_{De} / k_y^2 $ 
!|       & effective total diffusivities
!|         for $ n_H T_H $, $ n_H $, $ n_e T_e $, 
!|         $ n_Z $, $ n_Z T_Z $, \ldots
!|         normalized by $ k_y^2 / \omega_{De} $ \\
!| {\tt vflux(j)}   & output & $ ({F}/{nT})({R k_y^2}/{\omega_{De}}) $
!|    & array containing computed fluxes \\
!| {\tt nmodes} & output & & number of unstable modes \\
!| {\tt nerr}       & output & & nerr $\neq 0 \rightarrow$ error\\  
!| \end{tabular}
!| \end{center}
!| \newpage
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Integer control variables in the argument list:}
!| \begin{tabular}{lp{6.0in}}
!| {\bf variable} & {\bf meaning} \\
!| {\tt letain(2)} & = number of elements computed for transport matrix
!|                     (only when $> 0$) \\
!| {\tt letain(7)} &$ > 0 \rightarrow$ rescale transport matrix with velthi(j) = 0 \\
!| {\tt letain(9)} &$ > 1 \rightarrow$ do not produce transport matrix, only 
!|                        effective diffusivities needed \\
!|                 & $ = 1 \rightarrow$ only diagonal elements of diffusion
!|                           matrix \\
!|                 & $ < 1 \rightarrow$ full diffusion matrix \\
!| {\tt letain(29)} & $ > 0 $ to print frequencies and fluxes mode by mode \\
!| \end{tabular}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Real-valued control variables in the argument list:}
!| \begin{tabular}{lp{6.0in}}
!| {\bf variable} & {\bf meaning} \\
!| {\tt cetain(10)} & = coefficient of $ k_\parallel $ [default to 1.0 in mmm95 model cswitch(17)] \\
!| {\tt cetain(11)}  & = 1./kpc, used to normalize $A_\parallel$ and $K$ [default to 1.0 in mmm95 model] \\
!| {\tt cetain(12)} & = coefficient in expression for $H$ [default to 0.0 in mmm95 model cswitch(19)] \\
!| {\tt cetain(15)} & = coefficient of $ \hat{\nu} $ [default to 0.0 in mmm95 model cswitch(18)] \\
!| {\tt cetain(20)} & = coefficient of $ \beta_{e,h,z} $ [default to 1.0 in mmm95 model cswitch(14)] \\
!| {\tt cetain(29)} & = radius used in printouts in weiland14flux\\
!| {\tt cetain(30)} & = finite difference used to construct
!|                    transport matrix \\
!| \end{tabular}
!| \end{center}
!| 
!| \renewcommand{\arraystretch}{1.0}
!| \begin{center}
!| {\bf Effects included with different numbers of equations:}
!| \begin{tabular}{lp{5.0in}}
!| {\bf neq} & {\bf effects} \\
!| 2 & Hydrogenic Ion Temperature Gradient mode only \\
!| 4 & Hydrogenic ITG + trapped electron modes \\
!| 5 & Hydrogenic ITG with parallel ion motion + TEM \\
!| 6 & Hydrogenic + impurity ITG + TEM \\
!| 7 & Hydrogenic + impurity ITG + TEM + collisions \\
!| 8 & Hydrogenic ITG with parallel hydrogen ion motion  
!|   + impurity ITG (without parallel impurity ion motion)
!|   + TEM + collisions\\
!| 9 & Hydrogenic ITG + ipurity ITG (without any parallel ion motion)
!|   + TEM + electromagnetic (finite $\beta$) effects \\
!| 10 & Hydrogenic ITG with parallel ion motion
!|   + impurity ITG without parallel ion motion
!|   + TEM + electromagnetic (finite $\beta$) effects + collisions \\
!| 11 & Hydrogenic ITG with parallel ion motion
!|   + impurity ITG with parallel ion motion
!|   + TEM + electromagnetic (finite $\beta$) effects + collisions \\
!| 
!| \end{tabular}
!| \end{center}
!| 
!| This routine {\tt weiland14} calls the routine {\tt weiland14flux}
!| to compute transport fluxes {\tt vflux} and the effective diffusivities
!| {\tt chieff}.
!| If a full transport matrix is requested ($ {\tt letain(9)} < 1 $),
!| derivatives of the transport fluxes are taken with respect to the
!| normalized gradients to compute {\tt difthi(i,j)}.
!| These diffusive fluxes are then subtracted from the transport fluxes
!| to compute the convective velocities {\tt velthi(j)}.
!| 
c-----------------------------------------------------------------------
c
c  External dependencies:
c
c  Call tree: WEILAND14 calls the following routines
c
c    WEILAND14FLUX - Calculates fluxes and effective diffusivities
c      TOMSQZ      - Wrapper for QZ algorithm solving Ax = lambda Bx
c        CQZHES    - First step in QZ algorithm 
c        CQZVAL    - Second and third step in QZ algorithm
c        CQZVEC    - Fourth step in QZ algorithm
c
c-----------------------------------------------------------------------

      implicit none
c
      integer idp
      parameter (idp = 15)

      integer 
     &   letain,   lprint, neq,    ndim,    nmodes,   nerr,   nout
     & , imatrx,   j1,     j2,     j
c
c ndim  = first dimension of the 2-D array difthi
c         and the maximum number of unstable modes allowed
c nmodes = number of unstable modes
c imatrx= the number of elements computed along each row and
c          column of the transport matrix

      logical zdiffinit
      data zdiffinit /.true./
c
      dimension 
     &   letain(32),        cetain(32),     omega(*),         gamma(*)
     & , chieff(*),         vflux(*),       difthi(ndim,*),   velthi(*)
c
      real 
     &   cetain,   gnein,  gnhin,  gnzin,   gtein,    gthin,  gtzin
     & , tauhin,   tauzin, fnzin,  czin,    azin,     fnsin,  betaein
     & , ftrapein, vef,    q,      shear,   ekyrhoin, wexb,   omega
     & , gamma,    chieff, vflux,  difthi,  velthi
c
c
      real  
     &   zgne,     zgnh,   zgnz,    zgns,    zgte,      zgth, zgtz
     & , zone,     zdg,    ztemp,   zepsqrt, zepsmach,  zdgflux(idp)
     & , ztwo,     zhalf,  chidum(idp),      zfluxmax

      save 
     &   zone,     ztwo,   zhalf,   zepsqrt, zepsmach, zdiffinit

      if (zdiffinit) then 
        zone  = 1.0
        ztwo  = 2.0
        zhalf = zone/ztwo
        zepsmach = zhalf
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
        zepsqrt = sqrt(zepsmach)
        zdiffinit = .false.
      endif

c
c..initialize arrays 
c
      do j1=1,ndim
        omega(j1)  = 0.0
        gamma(j1)  = 0.0
        chieff(j1) = 0.0
        vflux(j1)  = 0.0
        velthi(j1) = 0.0
        do j2=1,ndim
          difthi(j1,j2) = 0.0
        enddo
      enddo
c
      do j1=1,idp
        zdgflux(j1) = 0.0
        chidum(j1)  = 0.0
      enddo
c
c...Define Transport matrix sizes
c
      imatrx = min ( neq - 1, 4 )
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, letain(2)-1 )

c
c...copy gradients to local variables
c

      zgne = gnein
      zgnh = gnhin
      zgnz = gnzin
      zgth = gthin
      zgte = gtein
      zgtz = gtzin
      zgns = zgne-zgnh*( zone-czin*fnzin-fnsin ) - zgnz*czin*fnzin
c
c... Establish fluxes
c

      call weiland14flux (
     &   letain,   cetain,   lprint,   neq,      nout,     zgne
     & , gnhin,    gnzin,    gtein,    gthin,    gtzin,    tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    chieff,   vflux,    nmodes
     & , nerr )
c
c  return on error condition
c
      if ( nerr .ne. 0 ) then 
         nerr = -7
         return
      endif
c
c  return if only chieff array is needed
c
      if (letain(9) .gt. 1) return
c
c  return if there is no flux
c
      zfluxmax = 0.0
      do j = 1, imatrx + 1
        zfluxmax = zfluxmax + abs( vflux(j) )
      enddo
      if ( zfluxmax .lt. zepsqrt ) return
c
c... Define forward differencing step
c
      zdg = cetain(30)
      if ( abs(zdg) .lt. zepsqrt ) zdg = 0.01
c
c... Take the derivative w.r.t gth
c

      zgth = zgth + sign(zdg,zgth)
      call weiland14flux (
     &   letain,   cetain,   lprint,   neq,      nout,     zgne 
     & , gnhin,    gnzin,    gtein,    zgth,     gtzin,    tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes
     & , nerr )
c
c  return on error condition
c
      if ( nerr .ne. 0 ) then 
         nerr = -7
         return
      endif

c
      do j = 1, imatrx +1 
         difthi(j,1) = (zdgflux(j) - vflux(j))/sign(zdg,zgth)
      enddo
      zgth = gthin

c
c... Take derivatives w.r.t gnh
c
      zgnh = zgnh + sign(zdg,zgnh)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns

      call weiland14flux (
     &   letain,   cetain,   lprint,   neq,      nout,     zgne
     & , zgnh,     gnzin,    gtein,    gthin,    gtzin,    tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes
     & , nerr )
c
c  return on error condition
c
      if ( nerr .ne. 0 ) then 
         nerr = -7
         return
      endif
c
      do j = 1, imatrx +1
         difthi(j,2) = (zdgflux(j) - vflux(j))/sign(zdg,zgnh)
      enddo
      zgnh = gnhin
c
c... Take derivatives w.r.t gte
c
      zgte = zgte + sign(zdg,zgte)
      zgne = zgnh*( zone-czin*fnzin-fnsin )+gnzin*czin*fnzin + zgns

      call weiland14flux (
     &   letain,   cetain,   lprint,   neq,      nout,     zgne
     & , gnhin,    gnzin,    zgte,     gthin,    gtzin,    tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes
     & , nerr )
c
c  return on error condition
c
      if ( nerr .ne. 0 ) then 
         nerr = -7
         return
      endif
c
      do j = 1, imatrx +1
         difthi(j,3) = (zdgflux(j) - vflux(j))/sign(zdg,zgte)
      enddo
      zgte = gtein
c
c... Take derivatives w.r.t gnz
c

      if (neq .gt. 4) then

        zgnz = zgnz + sign(zdg,zgnz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns

        call weiland14flux (
     &   letain,   cetain,   lprint,   neq,      nout,     zgne
     & , gnhin,    zgnz,     gtein,    gthin,    gtzin,    tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes
     & , nerr )
c
c  return on error condition
c
      if ( nerr .ne. 0 ) then 
         nerr = -7
         return
      endif
c
        do j = 1, imatrx + 1
           difthi(j,4) = (zdgflux(j) - vflux(j))/sign(zdg,zgnz)
        enddo
        zgnz = gnzin

c
c... Take derivatives w.r.t gtz
c
        zgtz = zgtz + sign(zdg,zgtz)
        zgne = zgnh*( zone-czin*fnzin-fnsin ) + zgnz*czin*fnzin + zgns

        call weiland14flux (
     &   letain,   cetain,   lprint,   neq,      nout,     zgne
     & , gnhin,    gnzin,    gtein,    gthin,    zgtz,     tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    chidum,   zdgflux,  nmodes
     & , nerr )
c
c  return on error condition
c
      if ( nerr .ne. 0 ) then 
         nerr = -7
         return
      endif
c
        do j = 1, imatrx + 1
          difthi(j,5) = (zdgflux(j) - vflux(j))/sign(zdg,zgtz)
        enddo
        zgtz = gtzin
      else 
        do j = 1,imatrx+1
          difthi(j,4) = 0.0
          difthi(j,5) = 0.0
        enddo
      endif
           

c
c... Convective velocities
c
      velthi(1) = vflux(1) - difthi(1,1) * gthin - difthi(1,2) * gnhin
     &                     - difthi(1,3) * gtein - difthi(1,4) * gnzin
c
      velthi(2) = vflux(2) - difthi(2,1) * gthin - difthi(2,2) * gnhin
     &                     - difthi(2,3) * gtein - difthi(2,4) * gnzin
c
      velthi(3) = vflux(3) - difthi(3,1) * gthin - difthi(3,2) * gnhin
     &                     - difthi(3,3) * gtein - difthi(3,4) * gnzin
c
      if (neq .gt. 4) then 
c
      velthi(4) = vflux(4) - difthi(4,1) * gthin - difthi(4,2) * gnhin
     &                     - difthi(4,3) * gtein - difthi(4,4) * gnzin
c
      velthi(5) = vflux(5) - difthi(5,1) * gthin - difthi(5,2) * gnhin
     &                     - difthi(5,3) * gtein - difthi(5,4) * gnzin
      else
         velthi(4) = 0.0
         velthi(5) = 0.0
      endif
c
c... Temporary lines Needed for correspondance with etaw14a and earlier
c
      do j = 1, imatrx+1
         difthi(j,5) = 0.0
         difthi(5,j) = 0.0
      enddo
      velthi(5) = 0.0
c
c
c... Rescale diffusivity matrix and set velthi = 0 if that is requested
c    through letain(7) > 0
c      
      if (letain(7) .gt. 0) then
         do j = 1, imatrx +1
            ztemp = vflux(j) - velthi(j)
            do j1 = 1, imatrx +1
                difthi(j,j1) = difthi(j,j1) * vflux(j) / ztemp
            enddo
            velthi(j) = 0.0
         enddo
      endif
c
      return
c
      end
c@weiland14flux
c---:----1----:----2----:----3----:----4----:----5----:----6----:----7-c
      subroutine weiland14flux (
     &   letain,   cetain,   lprint,   neq,      nout,     gnein
     & , gnhin,    gnzin,    gtein,    gthin,    gtzin,    tauhin
     & , tauzin,   fnzin,    czin,     azin,     fnsin,    betaein
     & , ftrapein, vef,      q,        shear,    ekyrhoin, wexb
     & , ndim,     omega,    gamma,    chieff,   vfluxout, nmodes
     & , nerr )
c
c..See description of the argument list at the beginning of this file
c
!| 
!| This subroutine evaluates the eigenvalue problem for $\eta_i$
!| and trapped electron modes derived by 
!| Jan Weiland, H. Nordman and their group in G\"{o}teborg Sweden
!| \cite{bate98a} -- \cite{jarm87a}.
!| Transport fluxes are then computed using a quasi-linear estimate.
!| The equations in this routine include fast Hydrogenic ions,
!| impurities, trapped electron and finite Larmor radius effects.
!| New options include parallel ion motion, finite beta and collisional effects
!| together with an approximate treatment of ${\bf E}\times {\bf B}$ flow 
!| shear reduction of transport.
!| 
!| The essential idea is to linearize the fluid equations 
!| with magnetic drifts for a given Fourier harmonic
!| \[ u = u_0(x) + \tilde{u} \exp\left[ -i (\omega t
!|    + \vec{k} \cdot \vec{x} ) \right] \]
!| Compute the eigenvalues and eigenvectors from these equations and then
!| compute the quasi-linear heat and particle fluxes.
!| 
!| 
!| The fundamental equations used in this model are the fluid equations.
!| For ions (hydrogen or impurity), the equation of continuity is
!| \[ \Partial{n_i}{t} + \nabla \cdot ( n_i v_i ) = 0 \]
!| where
!| \[ v_i = v_E + v_{*i} + v_{Pi} + v_{\pi i}
!|    + \hat{B} v_{\parallel i} \]
!| are the fluid flows where
!| \[ v_E = E \times \hat{B} / B = \hat{B} \times \nabla \phi / B
!|  = - i k_y \phi \hat{x} / B \]
!| is the E cross B drift;
!| \[ v_{*i} = \frac{ \hat{B} \times \nabla ( n_i T_i )}{Z_i e n_i B}
!|  = i ( 1 + \eta_i ) \omega_{*i} / k_y \qquad{\rm where}\qquad
!|  \omega_{*i} = \frac{ - k_y T_i g_{ni} }{ Z_i e B R } \]
!| is the diamagnetic drift ($ Z_i = 1 $ for hydrogen isotopes and
!| $ Z_i = -1 $ for electrons);
!| \[ g_{ni} \equiv - R \hat{x} \cdot \nabla n_i / n_i \]
!| is the density gradient scale length
!| \[ \eta_i \equiv g_{Ti} / g_{ni}
!|  = \frac{n_i \hat{x} \cdot \nabla T_i}{T_i \hat{x} \cdot \nabla n_i}, \]
!| \[ v_{Pi} = \frac{ d E }{ dt } / ( B \Omega_i ) \]
!| is the polarization drift with
!| \[ \Omega_i = \frac{ Z_i e B }{ m_i}; \]
!| and
!| \[ v_{\pi i} = \frac{ \hat{B} \times \nabla \cdot \pi_i}{
!|   Z_i e n_i B } \]
!| is the drift due to the off-diagonal elements of the stress tensor, $\pi_i$.
!| 
!| The parallel ion motion $ v_{\parallel i} $ is determined by the parallel
!| ion momentum equation driven by electromagnetic forces and ion pressure
!| gradient along the field lines\cite{weil92a}:
!| \begin{equation}
!|   m_i n_i \Partial{v_{\parallel i}}{t}
!|  =  - e Z_i n_i \left[ \nabla_{\parallel} \phi + \frac{1}{c} \left(
!|     \Partial{A_\parallel}{t}
!|      - ( {\mathbf v}_{*i} \times {\mathbf B}_{\perp} ) \cdot \hat{{\mathbf B}}
!|          \right) \right]  - \nabla_{\parallel} ( n_i T_i ) 
!| \end{equation}
!| where $A_\parallel$ is the parallel component of the vector potential and
!| $ {\mathbf B}_{\perp} = {\mathbf \nabla} \times
!|     ( A_\parallel \hat{{\mathbf B}} )
!|    \approx {\mathbf \nabla} A_\parallel \times  \hat{{\mathbf B}} $.
!| 
!| \noindent
!| The ion energy balance equation is
!| \[ \frac{3}{2} n_i \left( \Partial{ }{t} + \vec{v}_i \cdot \nabla \right) T_i
!|   + n_i T_i \nabla \cdot \vec{v}_i = - \nabla \cdot \vec{q}_{*i}
!|   = \frac{5}{2} n_i ( \vec{v}_{*i} - \vec{v}_{Di} ) \cdot \nabla T_i \]
!| where $ \vec{q}_{*i} $ is the diamagnetic ion heat flow and
!| \[ \vec{v}_{Di} = \frac{T_i}{m_i \Omega_i} \hat{B} \times 
!|   \left( \frac{\nabla B}{B} + \vec{\kappa} \right) \]
!| is the drift due to $\nabla |B|$ and magnetic curvature
!| $ \vec{\kappa} = \hat{B} \cdot \nabla \hat{B} $.  Note
!| \[ \vec{k} \cdot \vec{v}_{Di} = \omega_{Di} = \frac{-2 k_y T_i}{Z_i e B R} \]
!| is the diamagnetic drift frequency.
!| 
!| The equations for trapped electron continuity and energy flow have
!| the same form except that the polarization drift and the drift
!| due to the stress tensor
!| can be neglected and the electron fluid velocity can be written
!| \[ v_e = v_E +  v_{*e}\]
!| where
!| \[ v_{*e} = -\frac{ \hat{B} \times \nabla ( n_e T_e )}{ e n_e B}
!|  = i ( 1 + \eta_e ) \omega_{*e} / k_y \qquad{\rm where}\qquad
!|  \omega_{*e} = \frac{  k_y T_e g_{ne} }{  e B R } \]
!| 
!| Let $f_t$ be the fraction of trapped electrons, $ f_t = n_{et} / n_e $,
!| where $n_{et}$ is the density of trapped electrons and
!| $n_{ef}$ and the density of free (circulating) electrons.  Then
!| \[ n_e = n_{et} + n_{ef}. \]
!| The perturbed density can be written
!| \[ \frac{\tilde{n}_e}{n_e}
!|  =  \frac{\tilde{n}_{et}}{n_e} + \frac{\tilde{n}_{ef}}{n_e}
!|  =  f_t \frac{\tilde{n}_{et}}{n_{et}}
!|        + ( 1 - f_t ) \frac{\tilde{n}_{ef}}{n_{ef}} \]
!| %  =  f_t \frac{\tilde{n}_{et}}{n_{et}}
!| %        + ( 1 - f_t ) \frac{ e \tilde{\phi} }{T_e}. \]
!| 
!| Let $ f_Z $ be the fraction of impurity ions with charge $Z$
!| relative to the electron density $ f_Z = n_Z / n_e $.
!| Charge neutrality gives
!| \[ n_e = n_H + Z n_Z + n_s \]
!| where $n_H$ is the density of the hydrogen isotopes and
!| $n_Z$ is the density of the impurity species and $n_s$ is the density of
!| superthermal hydrogenic ions.  Then
!| \[ \frac{\tilde{n}_e}{n_e}
!|  =  \frac{\tilde{n}_H}{n_e} + \frac{Z \tilde{n}_Z}{n_e}
!|  =  ( 1 - Z f_Z - f_s ) \frac{\tilde{n}_H}{n_H}
!|        + Z f_Z \frac{\tilde{n}_Z}{n_Z}  \]
!| where $ f_s \equiv n_s / n_e $.
!| For now, it is assumed that the superthermal ions 
!| do not take part in the perturbation $ \tilde{n}_s = 0.$
!| 
!| A relation between the perturbed free (circulating) electron density 
!| $ n_{ef} $ and the perturbed electric and magnetic potentials 
!| can be obtained from the momentum equation for free electrons parallel 
!| to the unperturbed magnetic field
!| \begin{equation}
!|  m_e \frac{\partial v_{\parallel e}}{\partial t}  = 
!|    e \left( \hat{{\mathbf B}} \cdot {\mathbf \nabla}\phi
!|    + \frac{1}{c} \frac{\partial A_{\parallel}}{\partial t} \right)
!|    - \frac{e}{c} \left( {\mathbf v}_{*ep} \times \delta {\mathbf B}_{\perp}
!|       \right) \cdot \hat{{\mathbf B}}
!|    - \frac{1}{n} \hat{{\mathbf e}}_{\parallel} \cdot {\mathbf \nabla} p 
!| \end{equation}
!| where 
!| the unperturbed perpendicular free electron velocity is the 
!| diamagnetic velocity (including the full electron pressure gradient)
!| $ {\mathbf v}_{*ep} = - \hat{{\mathbf B}} \times {\mathbf \nabla}
!|    ( n_{ef} T_e ) / ( e n_{ef} B ) $.
!| The inertial term on the left of this momentum equation can be neglected
!| if $ k_\parallel v_{the} \gg \omega $.
!| 
!| Assuming the electron temperature is nearly uniform along
!| the perturbed magnetic field 
!| $ {\mathbf B} \cdot {\mathbf \nabla} T_e \approx 0 $,
!| if follows that perturbed electron temperature is approximately
!| \begin{equation}
!|  \hat{T}_e = 
!|    \frac{ \omega_{*ep} }{ k_{\parallel} c }
!|    \frac{ g_{Te} }{  g_{ne} }
!|       \hat{A}_{\parallel}
!| \qquad{\rm where}\qquad  \hat{A}_\parallel \equiv e A_\parallel / T_e .
!| \end{equation}
!| 
!| The perturbed free electron density then follows by neglecting electron
!| inertia in Eq. (2)
!| \begin{equation}
!| {{\tilde{n}_{ef}}\over n_e}\equiv \hat{n}_{ef} =
!|    \hat{\phi} - ( \omega - \omega_{*e} ) \hat{A}_\parallel
!|     / ( k_\parallel c ).
!| \end{equation}
!| 
!| In order to relate $\phi$ and $A_{\parallel}$, use the 
!| electron continuity equation for the free electrons,
!| \begin{equation}
!| \frac{\partial n_{ef}}{\partial t}  + 
!|   {\mathbf v}_{E}\cdot{\mathbf \nabla} n_{ef} 
!|   + n_{ef}{\mathbf \nabla}\cdot{\mathbf v}_{E} 
!|    +  {\mathbf \nabla}\cdot \left( n_{ef} {\mathbf v}_{*e} \right) 
!|   - \frac{1}{e}{\mathbf \nabla}\cdot
!|    \left( J_{\parallel ef}\hat{{\mathbf B}} \right)
!|    = 0
!| \end{equation}
!| Only the free electrons contribute to the parallel electron current.
!| When parallel ion motion is included in the current,
!| Ampere's law implies,
!| \begin{equation}
!| J_{\parallel} = J_{\parallel ef} + J_{\parallel i} = - \frac{c}{4 \pi}
!|    \Delta A_{\parallel}
!| \end{equation}
!| Combining the last four equations,
!| eliminating $\tilde{n}_{ef}$, using the curvature relations
!| ${\mathbf \nabla} \cdot {\mathbf v}_E$ and
!|  ${\mathbf \nabla} \cdot \left( n_e {\mathbf v}_{*e} \right)$,
!| and assuming the same unperturbed temperature for the free and trapped
!| electrons, we obtain
!| \begin{equation}
!| \hat{A}_{\parallel}
!|  = \frac{k_{\parallel} c \left( \omega_{*e} - \omega \right) \hat{\phi}}
!| { \omega \left( \omega_{*e} - \omega \right)
!|    + \left( k_{\perp} \rho k_{\parallel} v_{A} \right)^2
!|    + \omega_{De} \left( \omega - \omega_{*eT} \right) }
!| \end{equation}
!| where $\omega_{*eT} = \omega_{De} ( g_{Te} + g_{ne} ) / 2 $.
!| 
!| Justification for the use of all of the above equations 
!| and more complete derivations
!| can be found in the references.\cite{weil92a,nord90a}
!| 
!| Note that 
!| \[ n_i \nabla \cdot \vec{v}_i = - \Partial{n_i}{t} - \vec{v}_i \cdot \vec{n}_i \]
!| from the equation of continuity can be used to eliminate
!| $ \nabla \cdot \vec{v}_i $ in the energy equation (divided by 3/2):
!| \[ n_i \left( \Partial{}{t} + \vec{v}_i \cdot \nabla \right) T_i
!| -\frac{2}{3} T_i \Partial{n_i}{t}
!|  - \frac{2}{3} T_i \vec{v}_i \cdot \nabla n_i
!|  = \frac{5}{3} n_i ( \vec{v}_{*i} - \vec{v}_{Di} ) \cdot \nabla T_i \]
!| and a similar equation for the trapped electrons.
!| The $\vec{v}_{*i}$ (and corresponding $\vec{v}_{*e}$)
!| contributions drop out because
!| \[ n_i \vec{v}_{*i} \cdot \nabla T_i
!|  - \frac{2}{3} T_i \vec{v}_{*i} \cdot \nabla n_i
!|  - \frac{5}{3} n_i \vec{v}_{*i} \cdot \nabla T_i
!|  = \frac{2}{3} \vec{v}_{*i} \cdot \nabla ( n_i T_i ) = 0 \]
!| since the diamagnetic drift
!| $ \vec{v}_{*i} = \hat{B} \times \nabla ( n_i T_i ) / ZeB n_i $
!| is perpendicular to the gradient of the pressure
!| $ \nabla ( n_i T_i ) $.
!| 
!| See the derivations in the notes by Weiland\cite{weil92a} for
!| the relations
!| \[ \nabla \cdot \delta ( {n_i v_{*i}} ) 
!|  = \vec{v}_{Di} \cdot \nabla \delta ( {n_i T_i} ) / T_i \]
!| (Ref. [2], page 132 Eq. (4.20)),
!| \[ \nabla \cdot \tilde{v}_E
!|    = \frac{Z e}{T_i} \vec{v}_{Di} \cdot \nabla \tilde{\phi} \]
!| (Ref. [2], page 132 Eq. (4.21)), and
!| \[ \nabla \cdot [ n_i ( \vec{v}_{pi} + \vec{v}_{\pi i} ) ]
!|  \approx - i n_i k_y^2 \rho_{si}^2
!|  [ \omega - \omega_{*i} ( 1 + \eta_i ) ] \frac{e \tilde{\phi}}{T_e} \]
!| (Ref. [2], page 25 Eq. (1.28)) where
!| \[ \rho_{si}^2 = \frac{T_e}{m_i \Omega_i^2} \]
!| with $ v_{thi}^2 = 2 T_i / m_i $ and $ \Omega_i = Z e B / m_i $.
!| 
!| Using these relations in the continuity, momentum, and energy equations, 
!| we obtain the following perturbed relations,
!| \[ (-\omega + \omega_{Di}) \hat{n}_i + \omega_{Di} \hat{T}_i
!|  + [ ( \omega_{Di} - \omega_{*i} ) Z_i T_e/T_i
!|  - k_y^2 \rho_{si}^2 ( \omega - \omega_{*i} ( 1 + \eta_i ) )
!|  ] \hat{\phi} + k_\parallel v_{\parallel i} = 0 \]
!| \[ - \omega m_i v_{\parallel i} + k_\parallel Z_i T_e \hat{\phi}
!|  + k_\parallel T_i ( \hat{n}_i + \hat{T}_i ) = 0 \]
!| \[ (-\omega + \frac{5}{3} \omega_{Di} ) \hat{T}_i
!|    + \frac{2}{3} \omega \hat{n}_i
!|    + \omega_{*e} \frac{g_{ni}}{g_{ne}}
!|    ( \eta_i - \frac{2}{3} ) \hat{\phi} = 0 \]
!| with corresponding equations for trapped electrons ($i=et$)
!| and impurities ($i=Z$).
!| Here $\hat{n} \equiv \tilde{n} / n$, 
!| $\hat{T} \equiv \tilde{T} / T$, and
!| $\hat{\phi} \equiv e \tilde{\phi} / T_e$ are dimensionless forms
!| of the perturbation.
!| 
!| Normalizing all frequencies by 
!| $\omega_{De} = \vec{k} \cdot \vec{v}_{De} = 2 k_y T_e /  e B R $,
!| and normalizing the parallel ion velocity by the speed of sound 
!| in that ion,
!| $ \hat{v}_{\parallel i} \equiv v_{\parallel i} / c_{si} $ 
!| where $ c_{si} \equiv \sqrt{T_e / m_i} $,
!| we obtain
!| \[ ( \hat{\omega} + \frac{T_i}{Z_i T_e} ) \hat{n}_i
!|   + \frac{T_i}{Z_i T_e} \hat{T}_i
!|   + \left[ 1 - \frac{g_{ni}}{2}
!|   + k_y^2 \rho_{si}^2 \left( \hat{\omega}
!|     + \frac{T_i}{Z_i T_e} \frac{g_{ni}+g_{Ti}}{2}
!|     \right) \right] \hat{\phi}
!|   - \frac{k_\parallel c_{si}}{\omega_{De}} \hat{v}_{\parallel i}
!|   = 0 \]
!| \[ \hat{\omega} \hat{v}_{\parallel i}
!|  - \frac{k_\parallel c_{si}}{\omega_{De}} Z_i \left[ \hat{\phi}
!|    + \frac{T_i}{Z_i T_e} ( \hat{n}_i + \hat{T}_i ) \right] = 0 \]
!| \[ \left( - \hat{\omega} - \frac{5}{3} \frac{T_i}{Z_i T_e}
!|    \right) \hat{T}_i + \frac{2}{3} \hat{\omega} \hat{n}_i
!|    + \frac{1}{2} \left( g_{Ti} - \frac{2}{3} g_{ni}
!|      \right) \hat{\phi} = 0 \]
!| where $ \hat{\omega} \equiv \omega / \omega_{De} $
!| and
!| \[ g_{ni} \equiv - R \hat{x} \cdot \nabla n_i / n_i \]
!| is the normalized gradient.
!| 
!| Now use
!| \[ f_t \hat{n}_{et} = ( 1 - Z f_Z - f_s ) \hat{n}_H
!|    + Z f_Z \hat{n}_Z - ( 1 - f_t ) \hat{\phi} \]
!| to eliminate the perturbed trapped electron density $\hat{n}_{et}$
!| in favor of the perturbed ion densities and perturbed potential.
!| The perturbed hydrogen and impurity equations remain as above
!| while the trapped electron density and energy equations become:
!| \[ (1-\hat{\omega}) ( 1 - Z f_Z - f_s ) \hat{n}_H
!|  + (1-\hat{\omega}) Z f_Z \hat{n}_Z  + f_t \hat{T}_{et}
!|  - [ 1 - f_t g_{ne} / 2 - ( 1 - f_t ) \hat{\omega} ] \hat{\phi}
!|  = 0 \]
!| \[ \frac{2}{3} \hat{\omega} ( 1 - Z f_Z - f_s ) \hat{n}_H
!|  + \frac{2}{3} \hat{\omega} Z f_Z \hat{n}_H
!|  + \left( \frac{5}{3} - \hat{\omega} \right) f_t \hat{T}_{et}
!|  + \left[ f_t \frac{g_{Te} - (2/3) g_{ne}}{2}
!|    - \frac{2}{3} \hat{\omega} ( 1 - f_t ) \right] \hat{\phi}
!|  = 0  \]
!| 
!| The parallel wavenumber $ k_\parallel = \sqrt{\alpha} / ( q R ) $
!| is estimated in reference \cite{weil92a}
!| in terms of the ballooning parameter $\alpha$ given by
!| \begin{equation}
!| \alpha = \left( k_1 + \sqrt{k_1 + S^2 k_2} \right) / 2
!| \end{equation}
!| with
!| \begin{equation}
!| k_1 = \frac{q^2}{4} \left( k_{\perp}\rho \right)^2
!|    \sqrt{ \frac{T_i}{T_e} \frac{ ( g_{Ti} + g_{ni} ) }{ 2 ( 1 - f_t ) } }
!| \qquad{\rm and}\qquad
!| k_2 = q^2 \left( k_{\perp}\rho \right)^4
!|     \frac{T_i}{T_e} \frac{ ( g_{Ti} + g_{ni} ) }{ 2 ( 1 - f_t ) } 
!| \end{equation}
!| where $ S = d \ln q / d \ln r $ is the magnetic shear.
!| 
!| In all of these expressions we use the notation
!| $$ \omega_{*e} = k_y T_e g_{ne} / ( e B R ) $$
!| $$ \omega_{*i} = - k_y T_i g_{ni} / ( Z_i e B R )
!|     = - T_i \omega_{*e} / ( Z_i T_e ) $$
!| $$ \omega_{De} = 2 k_y T_e / e B R = 2 \omega_{*e} / ( g_{ne} ) $$
!| $$ \omega_{Di} = 2 k_y T_i / Z_i e B R = 2 \omega_{*i} / ( g_{ni} )
!|    = - 2 T_i \omega_{*e} / (  g_{ni} T_e ) $$
!| $$ g_{ni} = - R \hat{x} \cdot \nabla n_i / n_i $$
!| $$ \eta_i = g_{Ti} / g_{ni}  $$
!| $$ \eta_e = g_{Te} / g_{ne} $$
!| and
!| $$ f_t \approx \sqrt{ \frac{ 2 r/R }{ 1 + r/R } } $$
!| is the trapped electron fraction.
!| Note that here I define 
!| \[ \epsilon_{ni} \equiv 1 / ( g_{ni} ) \]
!| while $ \epsilon_{ni} $ in the Weiland papers is defined to be
!| twice this value.
!| 
!| The basic technique is to set up the generalized eigenvalue problem
!| \[ A v = \lambda B v \]
!| where $ \lambda = \hat{\omega} + i \hat{\gamma} $ is the eigenvalue
!| and $ v $ is the corresponding eigenvector.
!| Hence the eigenvalues give the frequency and growth rates of the
!| modes while the eigenvectors give the phase of the perturbed
!| variables relative to one another.
!| 
!| In order to approximate the reduction of transport with $E\times B$ 
!| velocity shear the $E\times B $ shearing rate is subtracted from the 
!| linear growth rates obtained as above. Currently this is implemented as
!| follows: The generalized eigenvalue problem is redefined using  
!|  \[ (A - {\it i} \omega_{E\times B} B) v' = {\lambda}' B v' \]
!| and fluxes are calculated rom the new growth rates   ${\lambda}'$ 
!| and eigenvectors $v'$. This method do give the same results as a 
!| direct reduction of the growth rates only would give but is more 
!| integrated with the current framework of the model.
!| 
!|   
!| In this routine, the perturbed variables are always computed in
!| the following order:
!| \[ \hat{\phi}, \hat{T}_H, \hat{n}_H, \hat{T}_{et},
!|     \hat{n}_Z, \hat{T}_Z, \ldots \]
!| 
c
c.. the table of argument input and output variable names and 
c   definitions is givin in tex documentation
c
c  Note that external subroutines need to be provided:
c
c  tomsqz   Code wrapper for ACM/TOMS routine 535 implementing the QZ algorithm
c           complex generailized eigenvalue problem (eigenvalues and 
c           eigenvectors) 
c           The algorithm itself consists of the three routines 
c           cqzhes, cqzval and cqzvec
c
c
      implicit none
c
      integer idp
      parameter ( idp = 15 )
c
      integer 
     &   letain,   lprint,  neq,     ndim,   nmodes, nerr,      nout
     & , ieq,      idim,    imatrx,  j1,     j2,     j
c
c ieq    = number of equations
c
c imatrx = the number of elements computed along each row and
c          column of the transport matrix
c
      logical inital
      data inital /.true./
c
      dimension   
     &   letain(32),        cetain(32),      omega(*)
     & ,  gamma(*),         vfluxout(*),     chieff(*)
c
      real 
     &   cetain,   gnein,   gnhin,   gnzin,  gtein,  gthin,     gtzin
     & , tauhin,   tauzin,  fnzin,   czin,   azin,   fnsin,     betaein
     & , ftrapein, vef,     q,       shear,  h,      ekyrhoin,  wexb
     & , omega,    gamma,   chieff,  vfluxout 
c
c  Note:  Normally the transport matrix is  ieq-1 by ieq-1
c    however, if there are 6 equations, compute only a 4 by 4 matrix
c    since the impurity temperature equation is not used by most
c    transport codes including BALDUR
c
      real 
     &   zamr(idp,idp),       zami(idp,idp),     zbmr(idp,idp)
     & , zbmi(idp,idp),       zamrt(idp,idp),    zamit(idp,idp)
     & , zbmrt(idp,idp),      zbmit(idp,idp),    zvr(idp,idp)
     & , zvi(idp,idp),        zomega(idp),       zgamma(idp)
     & , zalfr(idp),          zalfi(idp),        zbeta(idp)
c
      integer ifail
c
c ( zamr(i,j), zami(i,j) ) = matrix A
c ( zbmr(i,j), zbmi(i,j) ) = matrix B
c
c Note that the eigenvalues are
c
c   zomega(j) = zalfr(j) / zbeta(j)
c   zgamma(j) = zalfi(j) / zbeta(j)
c
c and the eigenvectors are
c
c   zevec(j) = cmplx( zvr(j), zvi(j) )
c
c Here, zbeta(j) will be 0.0 in the case of an infinite eigenvalue
c
      complex  zevec(idp,idp)
c
      real  
     &   zepsmach, zepsqrt,   zone,  ztwo,   zthree,   zfour,   zfive
     & , zhalf,    zquarter,  ztvr,  zftr,   ztwohlf,  zgne,    zgnh
     & , zgnz,     zgns,      zgte,  zgth,   zgtz,     ztauh,   ztauz
     & , zft,      zimp,      zfnz,  zmass,  zfns,     zflh,    zflz
     & , zgamax,   zetae,     zetah
c
c zepsmach = machine epsilon
c zepsqrt  = sqrt ( machine epsilon )
c zone     = 1.0
c ztwo     = 2.0
c zthree   = 3.0
c zfour    = 4.0
c zfive    = 5.0
c zhalf    = 0.5
c zquarter = 0.25
c ztwohlf  = 2.5
c ztvr     = 2. / 3.
c zftr     = 5. / 3.
c
      real 
     &   zflxph,   zflxnh,   zflxpe, zflxnz, zflxpz,   zphsph,  zphsnh
     & , zphspe,   zphsnz,   zphspz, ztemp1, zreal,    zimag   
     & , ztempa(idp),        ztempb(idp),    zerreal,  zerimag, zerrmax
c
c  These are the thermal and particle fluxes and effective diffusivities
c
      real zflxm(idp), zchim(idp)
c
c  Normalized transport fluxes         Eff. Diffusivities
c       zflxm(1)   n_H T_H             zchim(1)    chi_i                
c       zflxm(2)   n_H                 zchim(2)    D_i
c       zflxm(3)   n_e T_e             zchim(3)    chi_e
c       zflxm(4)   n_Z                 zchim(4)    D_q
c       zflxm(5)   n_Z T_Z             zchim(5)    chi_q
c
      real 
     &   bt,        bt1,      zeni,  k1,      k2,      zkpsh,   zkpsz
     & , zalp,      zalf,     zanorm,         zbetae
c
!| 
!| Definitions of some of the internal variables:
!| 
!| \renewcommand{\arraystretch}{1.4}
!| \begin{center}
!| \begin{tabular}{lp{4.0in}}
!| variable & meaning \\
!| \\
!| {\tt zgth} & $ - \frac{R}{T_H} \Partial{T_H}{r} $ \\
!| {\tt zgte} & $ - \frac{R}{T_e} \Partial{T_e}{r} $ \\
!| {\tt zgtz} & $ - \frac{R}{T_Z} \Partial{T_Z}{r} $ \\
!| {\tt zgnh} & $ - \frac{R}{n_H} \Partial{n_H}{r} $ \\
!| {\tt zgnz} & $ - \frac{R}{n_Z} \Partial{n_Z}{r} $ \\
!| \end{tabular}  \end{center}
!| 
c
      save 
     &   idim,  zepsmach, zepsqrt,  zone,   ztwo,   zthree,   zfour
     & , zfive, zhalf,    zquarter, ztvr,   zftr,   ztwohlf, inital
c
c..initialize variables
c
c
      ieq = max ( 2, neq )
c
      vef = cetain(15) * vef
      imatrx = min ( ieq - 1, 4 )
      if ( letain(2) .gt. 0 ) imatrx = min ( imatrx, letain(2)-1 )
c
      if ( inital ) then
c
        idim = idp
c
        zone   = 1.0
        ztwo   = 2.0
        zthree = 3.0
        zfour  = 4.0
        zfive  = 5.0
c
        zhalf    = zone  / ztwo
        ztwohlf  = ztwo  + zhalf
        zquarter = zone  / zfour
        ztvr     = ztwo  / zthree
        zftr     = zfive / zthree
c
        zepsmach = zhalf
  2     if ( zhalf * zepsmach + zone .gt. zone ) then
          zepsmach = zhalf * zepsmach
          go to 2
        endif
c
        zepsqrt = sqrt ( zepsmach )
c
        inital = .false.
c
      endif
c
c..end of initialization
c
      nerr = 0
c
c..print header
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*)
     & 'Weiland-Nordman eigenvalue equations, subroutine weiland14flux'
        write (nout,*) '(all frequencies normalized by omega_{De})'
        write (nout,*) '(all diffusivities normalized by '
     &    ,'omega_{De} / k_y^2'
c
        write (nout,108) 'gnein', gnein
        write (nout,108) 'gnhin', gnhin
        write (nout,108) 'gnzin', gnzin
        write (nout,108) 'gtein', gtein
        write (nout,108) 'gthin', gthin
        write (nout,108) 'gtzin', gtzin
        write (nout,108) 'tauhin', tauhin
        write (nout,108) 'tauzin', tauzin
        write (nout,108) 'ftrapein', ftrapein
        write (nout,108) 'fnzin', fnzin
        write (nout,108) 'czin', czin
        write (nout,108) 'azin', azin
        write (nout,108) 'fnsin', fnsin
        write (nout,108) 'betaein', betaein
        write (nout,108) 'vefin', vef
        write (nout,108) 'qin', q
        write (nout,108) 'shearin', shear
        write (nout,108) 'ekyrhoin', ekyrhoin
        write (nout,108) 'wexb', wexb
 108    format (1x,a8,' = ',1pe14.6,',')
c
      endif
c
c...copy to local variables
c
      zgne = gnein
      zgnh = gnhin
      zgnz = gnzin
      zgth = gthin
      zgte = gtein
      zgtz = gtzin
c
c..check validity of input data
c
      if ( ndim .gt. idim ) then
         nerr = -6
         return
      endif
c
c..initialize arrays
c
      do j1=1,ndim
        omega(j1)    = 0.0
        gamma(j1)    = 0.0
        vfluxout(j1) = 0.0
        chieff(j1)   = 0.0
      enddo
c
      do j1=1,idp
        zomega(j1)   = 0.0
        zgamma(j1)   = 0.0
        zalfr (j1)   = 0.0
        zalfi (j1)   = 0.0
        zbeta (j1)   = 0.0
        zchim(j1)    = 0.0
        zflxm(j1)    = 0.0
        do j2=1,idp
          zevec(j1,j2)  = ( 0.0, 0.0 )
          zamr (j1,j2)  = 0.0
          zami (j1,j2)  = 0.0
          zbmr (j1,j2)  = 0.0
          zbmi (j1,j2)  = 0.0
          zamrt(j1,j2)  = 0.0
          zamit(j1,j2)  = 0.0
          zbmrt(j1,j2)  = 0.0
          zbmit(j1,j2)  = 0.0
          zvr  (j1,j2)  = 0.0
          zvi  (j1,j2)  = 0.0
        enddo
      enddo
c
c..set up initial gradients
c
      zimp   = czin
      zmass  = azin
      zfnz   = fnzin * zimp
      zfns   = fnsin
c
c..Save zgns = - R ( d n_s / d r ) / n_e
c    where n_s is the fast ion density
c
      zgns = zgne - zgnh * ( zone - zfnz - zfns ) - zgnz * zfnz
c
c..compute the rest of the dimensionless variables needed
c
c     ztauz = T_Z / ( Z T_e ),     zfnz  = Z n_Z / n_e
c     zimp  = Z,                   zflh  = k_y^2 \rho_{sH}^2
c     zmass = m_Z / m_H,           zflz  = k_y^2 \rho_{sZ}^2
c

      ztauh  = tauhin
      zeni = zhalf * zgne
      zbetae = max ( cetain(20) * betaein, zepsqrt )
      zimp   = max ( czin, zone )
      zmass  = max ( azin, zone )
      zfnz   = fnzin * zimp
      zft    = ftrapein
      zflh   = ekyrhoin**2
c
      if ( neq .gt. 4 ) then
        ztauz  = tauzin / czin
        zflz   = zmass * zflh / zimp**2
      else
        ztauz  = zone
        zflz   = 0.0
      endif
      zetae = zgte / zgne
      zetah = zgth / zgne
c
c..diagnostic output
c
      if ( lprint .gt. 2 ) then
        write (nout,*)
        write (nout,*) '--------------------------------------'
        write (nout,*)
        write (nout,*)
        write (nout,*) zgnh,' = zgnh'
        write (nout,*) zgne,' = zgne'
        write (nout,*) zgnz,' = zgnz'
        write (nout,*) zgns,' = zgns'
        write (nout,*) zgth,' = zgth'
        write (nout,*) zgte,' = zgte'
        write (nout,*) zgtz,' = zgtz'
        write (nout,*) zetae,' = zetae'
        write (nout,*) zetah,' = zetah'
        write (nout,*)
c
        write (nout,*) ztauh,' = ztauh'
        write (nout,*) ztauz,' = ztauz'
        write (nout,*)
        write (nout,*) zft,' = zft'
        write (nout,*) zimp,' = zimp'
        write (nout,*) zmass,' = zmass'
        write (nout,*) fnzin,' = fnz'
        write (nout,*) zfnz,' = zfnz'
        write (nout,*) zfns,' = zfns'
        write (nout,*) zflh,' = zflh'
        write (nout,*) zflz,' = zflz'
        write (nout,*)
        write (nout,*) zepsqrt,' = zepsqrt'
        write (nout,*) zepsmach,' = zepsmach'
        write (nout,109) 'gnsin', zgns
 109  format (1x,a8,' = ',1pe14.6,', computed from quasi-neutrality')

      endif

      if ( ieq .gt. 4 ) then
c
c..constants from Nilsson and Weiland, NF 34 (1994) 803  section 2
c  and from Weiland and Hirose NF 32 (1992) 151
c
c  ** Note:  zanorm = 1./kpc is used to normalize A_\parallel and K
c
        zanorm = cetain(11)
c
        bt  = 1.5
        bt1 = bt - ztwohlf
c
c  Expressions from Weiland and Hirose, NF 32 (1992) 151  Eq. (6)
c  Note that k1 rather than k1**2 appears under sqrt in zalp.
c  This is the result of a numerical fit rather than analytic solution.
c  Also, only free electrons are used from zbetae in zalf
c  Hence, zbetae -> zbetae * (zone-zft).
c
        if ( zgnh + zgth .lt. zepsqrt ) then
c
          zalp  = 0.0
          zalf  = 0.0
          zkpsh = 0.0
          zkpsz = 0.0
c
        else
c
          k1 = zquarter * q * q * zflh *
     &         sqrt( zhalf * (zgnh + zgth) * ztauh / (zone-zft) )
          k2 = q * q * zflh * zflh * zhalf * (zgnh + zgth) *
     &         ztauh / ( zone-zft)

          zalp = zhalf * ( k1 + sqrt( k1 + shear * shear * k2) )
          zalf = zalp / (ztwo * zflh * q * q * zbetae * (zone-zft))

          zkpsh = cetain(10) * zhalf * sqrt( zalp / zflh ) / q

          if ( zmass .gt. zepsqrt ) then
            zkpsz = zkpsh / sqrt ( zmass )
          else
            zkpsz = zkpsh
          endif
c
        endif
c
        if (lprint .gt. 5 ) then
          write (nout,*)
          write (nout,*) k1,'  = k1'
          write (nout,*) k2,'  = k2'
          write (nout,*) zalp,' = zalp'
          write (nout,*) zalf,' = zalf'
          write (nout,*) zkpsh,' = zkpsh'
          write (nout,*) zkpsz,' = zkpsz'
          write (nout,*)
          write (nout,*) bt1,' = bt1'
          write (nout,*) zanorm,' = zanorm'
        endif
c
      endif
c
c..NOTE! Below only the section, appropriate for the number of equations 
c        in the version of the Weiland model that is chosen, is used.
c        In each section the A and the B matrices are constructed.
c
c  For the MMM95 transport model used by the Lehigh Group, ieq = 10, and
c  collisions are turned off.
c
      if (ieq .eq. 11) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Eleven equations with impurities, trapped electrons, FLR,
c  collisional effects, parallel hydrogenic and impurity ion motion,
c  electromagnetic (finite beta) effects
c
c  equations for
c  e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vph, Av, K, Vpz
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &  'Eleven eqns for '
     &  ,'e phi / T_e, T_H, n_H, T_et, n_Z, T_Z, F, Vph, Av, K, Vpz'
      endif
c
c
c  hydrogenic density
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  hydrogenic energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  trapped electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - 2.5 * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c      impurity density
c
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
      zamr(5,11) = zkpsz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c      impurity energy
c
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
c
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
c
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
c
c  Parallel hydrogenic ion motion Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
c
      zbmr(8,8) = zone
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = ( zhalf * (zgne + zgte) - zalf * zflh ) * zanorm
      zamr(9,10) = - zeni * zanorm
      zamr(9,11) = zkpsz * zfnz / ( zone - zft )
c
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
c
c  time derivative of Av
c
      zamr(10,10) = zone
c
      zbmr(10,9) = zone
c
c..Impurity parallel ion motion
c
      zamr(11,1) = zimp * zkpsz
      zamr(11,5) = zimp * zkpsz * ztauz
      zamr(11,6) = zimp * zkpsz * ztauz
c
      zbmr(11,11) = zone
c
c
      elseif (ieq .eq. 10) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Ten equations with impurities, trapped electrons, FLR,
c  collisional effects, parallel ion motion, and electromagnetic effects
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp, Av, K
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &'Ten eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
c
c
c  hydrogenic density
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zamr(1,2) = - ztauh
      zamr(1,3) = - ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c  hydrogenic energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c  trapped electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,9) =  - ( zone - zft ) * zeni * zanorm
      zami(3,9) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,10) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,10) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,10) = ( zone - zft ) * zanorm
c
c  trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - 2.5 * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c      impurity density
c
      zamr(5,1) = - zone + zhalf*(zgnz - zflz * ztauz * (zgnz + zgtz))
      zamr(5,5) = - ztauz
      zamr(5,6) = - ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c      impurity energy
c
      zamr(6,1) = zhalf * (zgtz - ztvr * zgnz)
      zamr(6,6) = - ztauz * zftr
c
      zbmr(6,5) = - ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf * zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = - vef
c
      zbmr(7,1) = - zone
      zbmr(7,7) = zone
c
c  Parallel ion motion Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh * ztauh
      zamr(8,3) = zkpsh * ztauh
c
      zbmr(8,8) = zone
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(9,1) = zeni
      zamr(9,8) = zkpsh * ( zone - zfnz - zfns ) / ( zone - zft )
      zamr(9,9) = ( zhalf * (zgne + zgte) - zalf * zflh ) * zanorm
      zamr(9,10) = - zeni * zanorm
c
      zbmr(9,1) = zone
      zbmr(9,9) = zanorm
      zbmr(9,10) = - zanorm
c
c  time derivative of Av
c
      zamr(10,10) = zone
c
      zbmr(10,9) = zone
c
c
      elseif (ieq .eq. 9) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Nine equations with impurities, trapped electrons, FLR,
c  collisional effects, and electromagnetic effects
c  and parallel ion motion in strong ballooning approximation
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Av, K
c
      if (lprint .gt. 5 ) then
      write(nout, *)
     &'Nine eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Av, K'
      endif
c
c
      H = cetain(12) * 0.5 * abs( shear ) / max ( q, zepsqrt )
c
c     ion continuity
c
      zamr(1,1) = - zone + zhalf * (zgnh - zflh * ztauh * (zgnh + zgth))
      zami(1,1) = - H
      zamr(1,2) = - ztauh
      zami(1,2) = - ztauh*H
      zamr(1,3) = - ztauh
      zami(1,3) = - ztauh*H
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c      ion energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
      zamr(2,2) = - ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = - ztvr
c
c total electron density expressed through quasineutrality
c
      zamr(3,1) = zft * zeni - zone
      zami(3,1) = vef * (zone-zft)
      zamr(3,3) = zone - zfnz - zfns
      zami(3,3) = - vef * (zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = - vef * zfnz
      zami(3,7) = vef * zft
      zamr(3,8) = - ( zone - zft ) * zeni * zanorm
      zami(3,8) = vef * ( zone - zft ) * zeni * zanorm
      zamr(3,9) = ( zone - zft ) * (zone + zeni) * zanorm
      zami(3,9) = - vef * ( zone - zft ) * zanorm
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
      zbmr(3,9) = ( zone - zft ) * zanorm
c
c      trapped electron energy
c
      zamr(4,1) = zft * zhalf * ( zgte - ztvr * zgne )
      zami(4,1) = vef * ztvr * (bt - 2.5 * (zone-zft))
      zami(4,3) = - vef * ztvr * bt1 * (zone - zfnz - zfns)
      zamr(4,4) = zft * zftr
      zami(4,5) = - vef * ztvr * bt1 * zfnz
      zami(4,7) = - zftr * vef * zft
c
      zbmr(4,1) = (zone - zft) * ztvr
      zbmr(4,3) = - (zone - zfnz - zfns) * ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = - zfnz * ztvr
c
c      impurity density
c
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c      impurity energy
c
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
c
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
c
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
c
c  electromagnetic parallel vector potential Av = e A_par /Te
c
      zamr(8,1) = zeni
      zamr(8,8) = ( zhalf*(zgne+zgte)-zalf*zflh ) * zanorm
      zamr(8,9) = - zeni * zanorm
c
      zbmr(8,1) = zone
      zbmr(8,8) = zanorm
      zbmr(8,9) = - zanorm
c
c   time derivative of Av
c
      zamr(9,9) = zone
c
      zbmr(9,8) = zone
c
c
      elseif (ieq .eq. 8) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Eight equations with impurities, trapped electrons, FLR,
c  collisional effects, and parallel ion motion
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z, F, Vp
c
      if (lprint .gt. 5 ) then
      write(nout,*)
     &'Eight eqns for e phi/T_e, Th, n_h, Tet, n_z, T_z, F, Vp'
      endif
c
c
c     ion continuity
c
      zamr(1,1) = -zone + zhalf*(zgnh-zflh*ztauh*(zgnh+zgth))
      zamr(1,2) = -ztauh
      zamr(1,3) = -ztauh
      zamr(1,8) = zkpsh
c
      zbmr(1,1) = zflh
      zbmr(1,3) = zone
c
c      ion energy
c
      zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
      zamr(2,2) = -ztauh * zftr
c
      zbmr(2,2) = zone
      zbmr(2,3) = -ztvr
c
c total electron density expressed through quasineutrality
c
      zamr(3,1) = zft*zeni - zone
      zami(3,1) = vef*(zone-zft)
      zamr(3,3) = zone - zfnz- zfns
      zami(3,3) = -vef*(zone - zfnz - zfns)
      zamr(3,4) = zft
      zamr(3,5) = zfnz
      zami(3,5) = -vef*zfnz
      zami(3,7) = vef*zft
c
      zbmr(3,1) = zft - zone
      zbmr(3,3) = zone - zfnz - zfns
      zbmr(3,5) = zfnz
c
c      trapped electron energy
c
      zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
      zami(4,1) = vef*ztvr*(bt-2.5*(zone-zft))
      zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
      zamr(4,4) = zft*zftr
      zami(4,5) = -vef*ztvr*bt1*zfnz
      zami(4,7) = -zftr*vef*zft
c
      zbmr(4,1) = (zone-zft)*ztvr
      zbmr(4,3) = - (zone - zfnz - zfns)*ztvr
      zbmr(4,4) = zft
      zbmr(4,5) = -zfnz*ztvr
c
c      impurity density
c
      zamr(5,1) = -zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
      zamr(5,5) = -ztauz
      zamr(5,6) = -ztauz
c
      zbmr(5,1) = zflz
      zbmr(5,5) = zone
c
c      impurity energy
c
      zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
      zamr(6,6) = -ztauz*zftr
c
      zbmr(6,5) = -ztvr
      zbmr(6,6) = zone
c
c  variable F
c
      zamr(7,1) = zhalf*zgte - zone
      zami(7,1) = vef
      zamr(7,7) = zone
      zami(7,7) = -vef
c
      zbmr(7,1) = -zone
      zbmr(7,7) = zone
c
c  Parallel ion motion Vpi/Cs
c
      zamr(8,1) = zkpsh
      zamr(8,2) = zkpsh*ztauh
      zamr(8,3) = zkpsh*ztauh
c
      zbmr(8,8) = zone
c
c
      elseif ( ieq .eq. 7 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Seven equations with impurities, trapped electrons, FLR, 
c  and collisions
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, T_Z and F
c  Here, F is defined as F = GM*e phi/T_e
c  where GM=1+etae/(epsn*(omega-1+i*vef))
c
      if ( lprint .gt. 5 ) then
      write (nout,*)
      write (nout,*)
     & 'Seven eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z and F'
      endif
c
c
c  hydrogen density
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  hydrogen energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
c
c  trapped electron density
c
        zamr(3,1) = - zone + zft * zeni 
        zami(3,1) = vef*(zone-zft)
        zamr(3,3) = zone - zfnz - zfns
        zami(3,3) = -vef*(zone - zfnz - zfns)
        zamr(3,4) = zft
        zamr(3,5) = zfnz
        zami(3,5) = -vef*zfnz
        zami(3,7) = vef*zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft*zhalf*(zgte-ztvr*zgne)
        zami(4,1) = vef*ztvr*(bt-2.5*(zone-zft))
        zami(4,3) = -vef*ztvr*bt1*(zone - zfnz - zfns)
        zamr(4,4) = zft * zftr
        zami(4,5) = -vef*ztvr*bt1*zfnz
        zami(4,7) = -zftr*vef*zft
c
        zbmr(4,1) = ( zone - zft ) *ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) *ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
        zamr(5,1) = - zone+zhalf*(zgnz-zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
c
c  impurity energy
c
        zamr(6,1) = zhalf*(zgtz-ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
c
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
c
c  variable F
c
        zamr(7,1) = zhalf*zgte - zone
        zami(7,1) = vef
        zamr(7,7) = zone 
        zami(7,7) = -vef
c
        zbmr(7,1) = -zone
        zbmr(7,7) = zone
c
c
      elseif ( ieq .eq. 6 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..Six equations with impurities, trapped electrons, and FLR
c
c  equations for e \phi / T_e, T_H, n_H, T_{et}, n_Z, and T_Z
c
      if ( lprint .gt. 1 ) then
        write (nout,*)
        write (nout,*) 
     &   'Six eqns for e phi/T_e, T_H, n_H, T_{et}, n_Z, and T_Z'
      endif
c
c
c  hydrogen density
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  hydrogen energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) = zone
        zbmr(2,3) = - ztvr
c
c  trapped electron density
c
        zamr(3,1) = zft*zeni - zone
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
        zamr(3,5) = zfnz
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
        zbmr(3,5) = zfnz
c
c  trapped electron energy
c
        zamr(4,1) = zft * zhalf * (zgte - ztvr*zgne)
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - ( zone - zfnz - zfns ) * ztvr
        zbmr(4,4) = zft
        zbmr(4,5) = - zfnz * ztvr
c
c  impurity density
c
        zamr(5,1) = -zone + zhalf * (zgnz - zflz*ztauz*(zgnz+zgtz))
        zamr(5,5) = - ztauz
        zamr(5,6) = - ztauz
c
        zbmr(5,1) = zflz
        zbmr(5,5) = zone
c
c  impurity energy
c
        zamr(6,1) = zhalf * (zgtz - ztvr*zgnz)
        zamr(6,6) = - ztauz * zftr
c
        zbmr(6,5) = - ztvr
        zbmr(6,6) = zone
c
c
      elseif ( ieq .eq. 5 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..5 equations with trapped electrons, FLR effects, and parallel ion motion
c
c  equations for e phi/T_e, T_H, n_i, T_e, and Vp
c
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*)
     &    ' Five eqns for e phi/T_e, T_H, n_H, T_e, and Vp'
       endif
c
c
c  ion continuity
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
        zamr(1,5) = zkpsh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  ion energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   by the ion density perturbation.
c   The dilution factor 1-zfnz has now been added.
c
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
c
c
c  Parallel ion motion Vpi/Cs
c
        zamr(5,1) = zkpsh
        zamr(5,2) = zkpsh*ztauh
        zamr(5,3) = zkpsh*ztauh
c
        zbmr(5,5) = zone
c
      elseif ( ieq .eq. 4 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..4 equations with trapped electrons and FLR effects
c
c  equations for e phi/T_e, T_H, n_i, and T_e
c
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*) ' Four eqns for e phi/T_e, T_H, n_H, and T_e'
       endif
c
c
c  ion continuity
c
        zamr(1,1) = -zone + zhalf * (zgnh - zflh * ztauh * (zgnh+zgth))
        zamr(1,2) = - ztauh
        zamr(1,3) = - ztauh
c
        zbmr(1,1) = zflh
        zbmr(1,3) = zone
c
c  ion energy
c
        zamr(2,1) = zhalf * ( zgth - ztvr*zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(2,2) =   zone
        zbmr(2,3) = - ztvr
c
c  trapped electron continuity
c
c   Calculates the total electron density perturbation and replaces it
c   The dilution factor 1-zfnz has now been added.
c
        zamr(3,1) = - zone + zft * zhalf * zgne
        zamr(3,3) = zone - zfnz - zfns
        zamr(3,4) = zft
c
        zbmr(3,1) = zft - zone
        zbmr(3,3) = zone - zfnz - zfns
c
c  trapped electron energy
c
        zamr(4,1) = zft * ( zgte - ztvr * zgne ) * zhalf
        zamr(4,4) = zft * zftr
c
        zbmr(4,1) = ( zone - zft ) * ztvr
        zbmr(4,3) = - zone * ztvr
        zbmr(4,4) = zft
c
      elseif ( ieq .eq. 2 ) then
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c..two equations when trapped particles and FLR effects omitted
c
c
c  equations for e phi/T_e and T_H
c
       if ( lprint .gt. 1 ) then
         write (nout,*)
         write (nout,*) ' Two eqns for e phi/T_e and T_H'
       endif
c
        zamr(1,1) = zhalf * zgnh - ztauh - zone
        zamr(1,2) = - ztauh
        zamr(2,1) = zhalf * ( zgth - ztvr * zgnh )
        zamr(2,2) = - ztauh * zftr
c
        zbmr(1,1) = zone
        zbmr(1,2) = 0.0
        zbmr(2,1) = - ztvr
        zbmr(2,2) = zone
c
      endif
c
c
c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c
c..find the eigenvalues and eigenvectors using ACM/TOMS routine 535
c
      ifail = -1
c
c..Update A matrix w.r.t wexb and save local copies of the matrices
c
      do j1=1,ieq
       do j2=1,ieq
        zamr(j1,j2)  = zamr(j1,j2) + abs(wexb)*zbmi(j1,j2)
        zami(j1,j2)  = zami(j1,j2) - abs(wexb)*zbmr(j1,j2)
        zamrt(j1,j2) = zamr(j1,j2) 
        zamit(j1,j2) = zami(j1,j2) 
        zbmrt(j1,j2) = zbmr(j1,j2) 
        zbmit(j1,j2) = zbmi(j1,j2) 
       enddo
      enddo

c
c..diagnostic output for defining matrices
c
      if ( lprint .gt. 6 ) then
c
        write (nout,*)
        write (nout,*) ieq
c
        write (nout,*)
        write (nout,*) ' zamr(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zamr(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zami(j1,j2)  j2 ->'
        do j1=1,ieq
          write (nout,192) (zami(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zbmr(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmr(j1,j2),j2=1,ieq)
        enddo
c
        write (nout,*)
        write (nout,*) ' zbmi(j1,j2)  j2->'
        do j1=1,ieq
          write (nout,192) (zbmi(j1,j2),j2=1,ieq)
        enddo
 192  format (1p10e12.4)
      endif
c
c... Solution of the generalized eigenvalue problem
c
      call tomsqz(idim,ieq,zamr,zami,zbmr,zbmi,ZALFR,ZALFI,ZBETA, 
     &            ZVR,ZVI,IFAIL)
c
c... Storing the results -  eigenvalues and eigenvectors
c
      zgamax = 0.0
      do j=1,ieq
        ztemp1 = zbeta(j)
        if ( abs(zbeta(j)) .lt. zepsqrt ) ztemp1 = zepsqrt
        zomega(j) = zalfr(j)  / ztemp1
        zgamma(j) = zalfi(j)  / ztemp1
        omega(j) = zomega(j)
        gamma(j) = zgamma(j)
        zgamax = max ( zgamax, zgamma(j) )
        do j1=1,ieq 
          zevec(j1,j) = cmplx ( zvr(j1,j), zvi(j1,j))
        end do
      enddo

      nerr = ifail

c
c... If no unstable roots are found set fluxes to zero and leave
c        
      if ( zgamax .lt. zepsqrt ) go to 90
c
c..check the eigenfunctions
c
      if ( lprint .gt. 12 ) then
        write (nout,*)
        write (nout,*) ' Checking eigenfunctions'
      endif
c
c  Real and imaginary parts
c
      zerrmax = 0.0
      nmodes = 0
c
c Loop over number of possible modes
c
      do j=1,ieq
c
c Start an error check on the eigensolution comparin RHS w LHS
c
        do j1=1,ieq
          ztempa(j1) = 0.0
          ztempb(j1) = 0.0
          do j2=1,ieq
            zerreal =
     &            zamrt(j1,j2) * real (zevec(j2,j))
     &          - zamit(j1,j2) * aimag(zevec(j2,j))
     &          - zbmrt(j1,j2) * real (zevec(j2,j)) * zomega(j)
     &          + zbmrt(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * real (zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zomega(j)

            zerimag =
     &            zamrt(j1,j2) * aimag(zevec(j2,j))
     &          + zamit(j1,j2) * real (zevec(j2,j))
     &          - zbmrt(j1,j2) * aimag(zevec(j2,j)) * zomega(j)
     &          - zbmrt(j1,j2) * real (zevec(j2,j)) * zgamma(j)
     &          + zbmit(j1,j2) * aimag(zevec(j2,j)) * zgamma(j)
     &          - zbmit(j1,j2) * real (zevec(j2,j)) * zomega(j)

            ztempa(j1) = ztempa(j1) + zerreal
            ztempb(j1) = ztempb(j1) + zerimag

          enddo
          zerrmax = max ( zerrmax, abs(ztempa(j1)), abs(ztempb(j1)) )
        enddo

        if ( lprint .gt. 12 ) then
          write (nout,*)
          write (nout,*) ' LHS - RHS for j =  ',j
          do j1=1,ieq
            write (nout,142) ztempa(j1), ztempb(j1)
          enddo
        endif
c
!| 
!| The effective diffusivities can be computed directly from the 
!| eigenvectors in the following way:
!| Consider, for example, the flux of hydrogen particles produced by
!| the perturbed $ E \times B $ motion of the plasma
!| \[ \Gamma_H = \tilde{n}_H \tilde{v}_E^* + c.c. 
!|  = 2 ( {\rm Re} \tilde{n}_H {\rm Im} \tilde{\phi}
!|      - {\rm Im} \tilde{n}_H {\rm Re} \tilde{\phi} ) k_y / B. \]
!| Using the quasi-linear assumption, the saturation level is 
!| approximated by
!| \[ \frac{e \tilde{\phi}}{T_e}
!|  \approx \frac{1}{k_x \rho_{sH}} \frac{ \gamma }{ k_y c_{sH} }
!|  = \frac{2}{R k_x} \frac{\gamma}{\omega_{De}}. \]
!| Hence, the hydrogen particle flux is given by
!| \[ \frac{ F_H }{ n_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \hat{\gamma}^2
!|   ( {\rm Re} \hat{n}_H {\rm Im} \hat{\phi}
!|   - {\rm Im} \hat{n}_H {\rm Re} \hat{\phi} ) / | \hat{\phi} |^2.
!|    \]
!| Correspondingly, the hydrogen heat (here, $ n_H T_H $) flux is given by
!| \[ \frac{ F_{p_H} }{ n_H T_H } \frac{ R k_x^2 }{ \omega_{De} }
!|   = 2 \hat{\gamma}^2
!|   ( {\rm Re} ( \hat{n}_H + \hat{T}_H ) {\rm Im} \hat{\phi}
!|   - {\rm Im} ( \hat{n}_H + \hat{T}_H ) {\rm Re} \hat{\phi} )
!|     / | \hat{\phi} |^2.
!|    \]
!| 
!| The case of the electron heat flux (here, the flux of $ n_e T_e $)
!| is a little more complicated.  Note that
!| 
c
c..compute effective diffusivities directly from eigenvalues
c  assume eigenvectors are arranged in the order of
c  e\phi/T_e, T_H, n_H, T_{et}, n_Z, T_Z
c
c
c... Calculate norm of solution in terms of (e\phi/T_e) **2
c
        ztemp1 =  real(zevec(1,j)) *  real(zevec(1,j))
     &         + aimag(zevec(1,j)) * aimag(zevec(1,j))
c
c... Check if current mode j is unstable
c
        if ( zgamma(j).gt.zepsqrt .and. ztemp1.gt.zepsqrt ) then

           nmodes = nmodes + 1
c
c Define fluxes : Thermal hydrogen flux
c
          zreal =  real(zevec(2,j)) +  real(zevec(3,j))
          zimag = aimag(zevec(2,j)) + aimag(zevec(3,j))
c
c...phase difference
c
          zphsph = - ( zimag * real(zevec(1,j))
     &             -   zreal * aimag(zevec(1,j)) ) / ztemp1
c
c...flux from j:th mode
c
          zflxph = 2.0 * zphsph * zgamma(j) * zgamma(j)

c
c...flux summed over all unstable modes
c
          zflxm(1) = zflxm(1) + zflxph
c
c Define hydrogen density flux - phase shift
c
          zphsnh = - ( aimag(zevec(3,j)) * real(zevec(1,j))
     &             - real(zevec(3,j)) * aimag(zevec(1,j)) ) / ztemp1
c
c...Flux from j:th mode
c
          zflxnh = 2.0 * zphsnh * zgamma(j) * zgamma(j)
c
c...Flux summed over all unstable modes
c
          zflxm(2) = zflxm(2) + zflxnh
c
c Define thermal electron flux
c
          if ( ieq .gt. 3 ) then
            zreal =  real(zevec(4,j))
     &            + ( zone - zfnz - zfns ) * real(zevec(2,j))
     &            + zfnz * real(zevec(5,j))
            zimag = aimag(zevec(4,j))
     &            + ( zone - zfnz - zfns ) * aimag(zevec(2,j))
     &            + zfnz * aimag(zevec(5,j))
c
c  Note, the electron heat flux is reduced by the fraction of
c  trapped electrons - phase shift
c
            zphspe = - zft * ( zimag * real(zevec(1,j))
     &                       - zreal * aimag(zevec(1,j)) ) / ztemp1

c
c... Flux from the j:th mode
c
            zflxpe = 2.0 * zphspe * zgamma(j)* zgamma(j)
c
c... Flux summed over all unstable modes
c
            zflxm(3) = zflxm(3) + zflxpe
c
          endif
c
c... Impurity density flux
c
          if ( ieq .gt. 4 ) then
c
c... phase difference
c
            zphsnz = - (aimag(zevec(5,j))*real (zevec(1,j))
     &               -  real (zevec(5,j))*aimag(zevec(1,j)) ) / ztemp1

c
c... Flux from the j:th mode
c
            zflxnz =  2.0 * zphsnz * zgamma(j) * zgamma(j)
c
c... Flux summed over all unstable modes
c
            zflxm(4) = zflxm(4) + zflxnz
c
          endif
c
c Define the impurity heat flux
c
          if ( ieq .gt. 5 ) then
c
            zreal =  real(zevec(5,j)) +  real(zevec(6,j))
            zimag = aimag(zevec(5,j)) + aimag(zevec(6,j))
c
c... phase difference
c
            zphspz = - ( zimag * real (zevec(1,j))
     &               -   zreal * aimag(zevec(1,j)) ) / ztemp1
c
c... Flux from the j:th mode
c
            zflxpz = 2.0 * zphspz * zgamma(j) * zgamma(j) 
c
c... Flux summed over all unstable modes
c
            zflxm(5) = zflxm(5) + zflxpz
c
          endif
c			
c..header for diagnostic printout of frequencies, fluxes, and phases
c
        if ( letain(29) .gt. 0 .and. lprint .gt. 0) then
          write (nout,134)
c
 134  format (//' Diagnostic printout of frequencies, fluxes and phases'
     &  /' Note: fluxph = flux of hydrogenic thermal energy'
     &  ,' = 2.0 * gamma^2 * phaseph'
     &  /9x,' (phaseph is related to the phases of the perturbations)'
     &  /7x,'fluxnh = flux of hydrogenic ions = 2.0 * gamma^2 * phasenh'
     &  /7x,'fluxpe = flux of electron thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepe'
     &  /7x,'fluxnz = flux of impurity ions = 2.0 * gamma^2 * phasenz'
     &  /7x,'fluxpz = flux of impurity thermal energy'
     &  ,' = 2.0 * gamma^2 * phasepz'
     &  //1x,'radius',t10,'omega',t20,'gamma'
     &  ,t30,'fluxph',t40,'phaseph',t50,'fluxnh',t60,'phasenh'
     &  ,t70,'fluxpe',t80,'phasepe',t90,'fluxnz',t100,'phasenz'
     &  ,t110,'fluxpz',t120,'phasepz  #m')
c
c..diagnostic printout of frequencies and fluxes mode by mode
c
        write (nout,135) cetain(29), zomega(j), zgamma(j)
     &    , zflxph, zphsph, zflxnh, zphsnh, zflxpe, zphspe
     &    , zflxnz, zphsnz, zflxpz, zphspz
        endif
c
 135  format (f7.3,1p12e10.2,' #m')
c
        endif
c
c...End of flux and diffusivity definitions loop 
c
      enddo

c
c... Error diagnostics
c
      if ( lprint .gt. 0 ) then
c        if (abs(zerrmax) .gt. 100*zepsqrt) nerr = max (nerr , 1)
        if (abs(zerrmax) .gt. zepsqrt) nerr = max (nerr , 1)
        write (nout,*) 
        write (nout,*) zerrmax,' = zerrmax'
        write (nout,*) nerr,
     &       ' = nerr, error in eigenvalue in weiland14flux'
      endif

 142  format (1p10e12.4)
c
c..compute effective total diffusivities
c
      zchim(1) = zflxm(1) / sign( max( abs( zgth ), zepsqrt ),  zgth )
      zchim(2) = zflxm(2) / sign( max( abs( zgnh ), zepsqrt ),  zgnh )
      zchim(3) = zflxm(3) / sign( max( abs( zgte ), zepsqrt ),  zgte )
      zchim(4) = zflxm(4) / sign( max( abs( zgnz ), zepsqrt ),  zgnz )
      zchim(5) = zflxm(5) / sign( max( abs( zgtz ), zepsqrt ),  zgtz )
c
c..save effective diffusivities and fluxes
c
      do j1=1,min(ieq-1,5)
        chieff(j1)  = zchim(j1)
        vfluxout(j1) = zflxm(j1)
      enddo

c
      if ( lprint .gt. 2 ) then
c
c..print eigenvalues and eigenfunctions
c
        write (nout,121)
        do j=1,ieq
          write (nout,122) zomega(j), zgamma(j)
        enddo

 121    format (/' Solution of the eigenvalue equations'
     &   /t4,'zomega',t18,'zgamma')
 122    format (1p2e14.5,i5)
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (zchim(j1),j1=1,ieq-1)
c
      endif
c
      if ( lprint .gt. 99 ) then
c
        write (nout,*)
        write (nout,*) ' Eigenvectors zevec(j1,j2) j2->'
        do j1=1,ieq
          write (nout,124) (zevec(j1,j2),j2=1,ieq)
 124      format (2(1p12e11.3))
        enddo
c
        write (nout,*)
      endif
c
      if ( ifail .gt. 0 ) then
        nerr = 5
        return
      endif
c
c
c
c
      if ( imatrx .lt. 1 ) go to 90
c
c..diagnostic printout
c
      if ( lprint .gt. 6 ) then
c
        write (nout,*)
        write (nout,*) ' vector zflxm(j1)'
        do j1=1,ieq-1
          write (nout,132) zflxm(j1)
        enddo
c
      endif
c
  90  continue
c
      if ( lprint .gt. 0 ) then
c
        write (nout,*)
        write (nout,*) '--------------------------------------'
     &    ,'  Final results from sbrtn weiland14flux'
c
        write (nout,*)
        write (nout,*) ' Effective diffusivities from eigenvectors'
     &    ,' normalized by omega_{De} / k_y^2'
c
        write (nout,130)
        write (nout,132) (chieff(j1),j1=1,imatrx)
c
      endif
c
c
      return
c
 130    format (/t3,'chi_i',t16,'dif_h',t29,'chi_e'
     &    ,t42,'dif_Z',t55,'chi_Z')
 132    format (1p11e12.4)
c
      end
!| 
!| Note that
!| $$ \omega_{De} = 2 k_\perp T_e / e B R = 2 k_\perp \rho_s c_s / R $$
!| $$ D \propto \gamma / k^2
!|   \propto \rho_s^2 \omega_{De}
!|     \frac{ \gamma / \omega_{De} }{ k^2 \rho_s^2 }. $$
!| Note, $ \omega_{*e} = k_\perp \rho_s c_s g_{ne} $,
!| $ \rho_s = c_s / \omega_{ci} $, $ c_s = \sqrt{ 2 T_e / m_i} $,
!| and $ \omega_{ci} = Z_i e B / m_i $, 
!| $ \omega_{De} R = k_y \rho_s c_s $,
!| $v_i = \sqrt{T_i / m_i}$, and $ g_n = - R \Partial{n}{r} / n $
!| in the above normalizations.
!| Note that all the diffusivities in this routine are normalized by
!| $ \omega_{De} / k_y^2 $, 
!| convective velocities are normalized by $ \omega_{De} / k_y $, 
!| and all the frequencies are normalized by $ \omega_{De} $.
!| 
!| 
!| \begin{thebibliography}{99}
!| 
!| \bibitem{bate98a}
!| Glenn Bateman, Arnold~H. Kritz, Jon~E. Kinsey, Aaron~J. Redd, and Jan Weiland,
!| ``Predicting temperature and density profiles in tokamaks,''
!| {\em Physics of Plasmas,} {\bf 5} (1998) 1793--1799.
!| 
!| \bibitem{weil92a} J. Weiland,
!| ``Low Frequency Modes Associated with Drift Motions in
!| Inhomogeneous Plasmas,''
!| report CTH--IEFT/PP-1992-17,
!| Institute for Electromagnetic Field Theory and Plasma Physics,
!| Chalmers University of Technology,
!| G\"{o}teborg, Sweden.
!| 
!| \bibitem{froj92a} M. Fr\"{o}jdh, M. Liljestr\"{o}m, H. Nordman,
!| ``Impurity effects on $\eta_i$ mode stability and transport,''
!| Nuclear Fusion {\bf 32} (1992) 419--428.
!| 
!| \bibitem{nord92a} H. Nordman and J. Weiland, ``Comments on
!| `Ion-temperature-gradient-driven 
!| transport in a density modification experiment on the tokamak fusion test
!| reactor [Phys. Fluids {\bf B4} (1992) 953]' ''.
!| 
!| \bibitem{weil92b} J. Weiland, 
!| ``Nonlinear effects in velocity space and drift wave
!| transport in tokamaks,'' Phys. Fluids {\bf B4} (1992) 1388--1390.
!| 
!| \bibitem{weil92c} J. Weiland and H. Nordman, ``Drift wave model for inward
!| energy transport in tokamak plasmas,'' Institute for Electromagnetic Field
!| Theory and Plasma Physics, Gothenburg, Sweden, (1992) CTH-IEFT/PP-1992-13 ISSN.
!| 
!| \bibitem{weil92d} J.Weiland and A. Hirose, ``Electromagnetic and kinetic
!| effects on the ion temperature gradient mode,'' Nucl. Fusion {\bf 32} (1992)
!| 151--155.
!| 
!| \bibitem{nord91a} H. Nordman and J. Weiland, ``The concept of marginal
!| stability and recent experimental results from the TFTR tokamak,'' Institute
!| for Electromagnetic Field Theory and Plasma Physics, Gothenburg, Sweden, (1991)
!| CTH-IEFT/PP-1991-26 ISSN.
!| 
!| \bibitem{weil91a} J. Weiland and H. Nordman, ``Enhanced confinement regimes in
!| transport code simulations of toroidal drift wave transport,'' Nucl. Fusion
!| {\bf 31} (1991) 390--394.
!| 
!| \bibitem{nord90a} H. Nordman, J. Weiland, and A. Jarmen, ``Simulation of
!| toroidal drift mode turbulence driven by temperature gradients and electron
!| trapping,'' Nucl. Fusion {\bf 30} (1990) 983--996.
!| 
!| \bibitem{weil89a} J. Weiland, A.B. Jarm\'{e}n, and H. Nordman, ``Diffusive
!| particle and heat pinch effects in toroidal plasmas,'' Nucl. Fusion {\bf 29}
!| (1989) 1810--1814.
!| 
!| \bibitem{nord89a} H. Nordman and J. Weiland, ``Transport due to toroidal
!| $\eta_i$ mode turbulence in tokamaks,'' Nucl. Fusion {\bf 29} (1989) 251--263.
!| 
!| \bibitem{ande88a} P. Andersson and J. Weiland, 
!| ``A fully toroidal fluid analysis
!| of the magnetohydrodynamic ballooning mode branch in tokamaks,'' 
!| Phys. Fluids {\bf 31} (1988) 359--365.
!| 
!| \bibitem{jarm87a} A. Jarm\'{e}n, P. Andersson, and J. Weiland, 
!| ``Fully toroidal ion temperature gradient driven drift modes,'' 
!| Nucl. Fusion {\bf 27} (1987) 941--949.
!| \end{thebibliography}
!| 
!| \end{document}
!======================================================================|
      SUBROUTINE TOMSQZ(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI,IFAIL)
C-----------------------------------------------------------------------
C TOMSQZ  written by P. Strand 27-apr-98,       elfps@elmagn.chalmers.se
C----------------------------------------------------------------------- 
C CQZHES, CQZVEC, and CQZVAL: Fortran subroutines implementing the QZ 
C algorithm for solving the generalized eigenvalue problem for complex 
C matrices. (See B.S.C Garbow, ACM TOMS 4 (1978) pp. 404-410.).
C-----------------------------------------------------------------------
C 
C ON INPUT THE GENERALIZED EIGENVALUE PROBLEM IS DEFINED THROUGH THE 
C COMPLEX MATRICES
C
C       A = cmplx (AR, AI)  AND   B = cmplx (BR, BI)
C
C WHERE LEADING DIMENSION N IS AS DEFINED IN THE CALLING ROUTINE AND 
C WHERE NA IS THE ROW  RANGE IN THE CURRENT PROBLEM. THE EIGENVALUE 
C PROBLEM IS THEN DEFINED THROUGH 
C
C       A x = w B x
C
C WHERE  THE COMPLEX EIGENVECTORS 
C
C       x = cmplx (ZVR, ZVI) 
C
C TOGETHER WITH THE COMPLEX EIGENVALUE
C
C        w = cmplx(alfr, alfi)/beta
C
C IS OUTPUT FROM THE ROUTINE
C
C IFAIL WILL BE NONZERO IF CONVERGENCE HAS NOT BEEN REACH WITHIN 50 
C ITERATIONS
C-----------------------------------------------------------------------
C DECLARATIONS FOR INPUT VARIABLES
C-----------------------------------------------------------------------

      INTEGER N, NA
      REAL AR(N,NA),AI(N,NA),BR(N,NA),BI(N,NA)

C-----------------------------------------------------------------------
C DECALRATIONS FOR OUTPUT VARIABLES
C-----------------------------------------------------------------------

      REAL ALFR(N),ALFI(N),BETA(N)
      REAL ZVR(N,NA), ZVI(N,NA)

C-----------------------------------------------------------------------
C LOCAL VARIABLES
C-----------------------------------------------------------------------

      LOGICAL WANTX
      REAL EPS1

C-----------------------------------------------------------------------
C START OF ACTUAL CODING
C-----------------------------------------------------------------------

      WANTX = .TRUE.
      EPS1  = -0.0

      CALL CQZHES(N,NA,AR,AI,BR,BI,WANTX,ZVR,ZVI)
      CALL CQZVAL(N,NA,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,WANTX,
     &            ZVR,ZVI,IFAIL)
      CALL CQZVEC(N,NA,AR,AI,BR,BI,ALFR,ALFI,BETA,ZVR,ZVI)
      RETURN
      END
                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE CQZHES(NM,N,AR,AI,BR,BI,MATZ,ZR,ZI)                    
C
      INTEGER I,J,K,L,N,K1,LB,L1,NM,NK1,NM1
      REAL AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ZR(NM,N),ZI(NM,N)
      REAL R,S,T,TI,U1,U2,XI,XR,YI,YR,RHO,U1I
      REAL SQRT,CABS,ABS
      LOGICAL MATZ
      COMPLEX CMPLX
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FIRST STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX GENERAL MATRICES AND
C     REDUCES ONE OF THEM TO UPPER HESSENBERG FORM WITH REAL (AND NON-
C     NEGATIVE) SUBDIAGONAL ELEMENTS AND THE OTHER TO UPPER TRIANGULAR
C     FORM USING UNITARY TRANSFORMATIONS.  IT IS USUALLY FOLLOWED BY
C     CQZVAL  AND POSSIBLY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX GENERAL MATRIX,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER HESSENBERG FORM.  THE ELEMENTS
C          BELOW THE FIRST SUBDIAGONAL HAVE BEEN SET TO ZERO, AND THE
C          SUBDIAGONAL ELEMENTS HAVE BEEN MADE REAL (AND NON-NEGATIVE),
C
C        B HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        Z=(ZR,ZI) CONTAINS THE PRODUCT OF THE RIGHT HAND
C          TRANSFORMATIONS IF MATZ HAS BEEN SET TO .TRUE.
C          OTHERWISE, Z IS NOT REFERENCED.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
C     ********** INITIALIZE Z **********

      ZERO = 0.0
      ZONE = 1.0
      IF (.NOT. MATZ) GO TO 10
C
      DO 3 I = 1, N
C
         DO 2 J = 1, N
            ZR(I,J) = 0.0
            ZI(I,J) = 0.0
    2    CONTINUE
C
         ZR(I,I) = 1.0
    3 CONTINUE
C     ********** REDUCE B TO UPPER TRIANGULAR FORM WITH
C                TEMPORARILY REAL DIAGONAL ELEMENTS **********
   10 IF (N .LE. 1) GO TO 170
      NM1 = N - 1
C
      DO 100 L = 1, NM1
         L1 = L + 1
         S = 0.0
C
         DO 20 I = L, N
            S = S + ABS(BR(I,L)) + ABS(BI(I,L))
   20    CONTINUE
C
         IF (S .EQ. ZERO) GO TO 100
         RHO = 0.0
C
         DO 25 I = L, N
            BR(I,L) = BR(I,L) / S
            BI(I,L) = BI(I,L) / S
            RHO = RHO + BR(I,L)**2 + BI(I,L)**2
   25    CONTINUE
C
         R = SQRT(RHO)
         XR = CABS(CMPLX(BR(L,L),BI(L,L)))
         IF (XR .EQ. ZERO) GO TO 27
         RHO = RHO + XR * R
         U1 = -BR(L,L) / XR
         U1I = -BI(L,L) / XR
         YR = R / XR + 1.0
         BR(L,L) = YR * BR(L,L)
         BI(L,L) = YR * BI(L,L)
         GO TO 28
C
   27    BR(L,L) = R
         U1 = -1.0
         U1I = 0.0
C
   28    DO 50 J = L1, N
            T = 0.0
            TI = 0.0
C
            DO 30 I = L, N
               T = T + BR(I,L) * BR(I,J) + BI(I,L) * BI(I,J)
               TI = TI + BR(I,L) * BI(I,J) - BI(I,L) * BR(I,J)
   30       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 40 I = L, N
               BR(I,J) = BR(I,J) - T * BR(I,L) + TI * BI(I,L)
               BI(I,J) = BI(I,J) - T * BI(I,L) - TI * BR(I,L)
   40       CONTINUE
C
            XI = U1 * BI(L,J) - U1I * BR(L,J)
            BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
            BI(L,J) = XI
   50    CONTINUE
C
         DO 80 J = 1, N
            T = 0.0
            TI = 0.0
C
            DO 60 I = L, N
               T = T + BR(I,L) * AR(I,J) + BI(I,L) * AI(I,J)
               TI = TI + BR(I,L) * AI(I,J) - BI(I,L) * AR(I,J)
   60       CONTINUE
C
            T = T / RHO
            TI = TI / RHO
C
            DO 70 I = L, N
               AR(I,J) = AR(I,J) - T * BR(I,L) + TI * BI(I,L)
               AI(I,J) = AI(I,J) - T * BI(I,L) - TI * BR(I,L)
   70       CONTINUE
C
            XI = U1 * AI(L,J) - U1I * AR(L,J)
            AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
            AI(L,J) = XI
   80    CONTINUE
C
         BR(L,L) = R * S
         BI(L,L) = 0.0
C
         DO 90 I = L1, N
            BR(I,L) = 0.0
            BI(I,L) = 0.0
   90    CONTINUE
C
  100 CONTINUE
C     ********** REDUCE A TO UPPER HESSENBERG FORM WITH REAL SUBDIAGONAL
C                ELEMENTS, WHILE KEEPING B TRIANGULAR **********
      DO 160 K = 1, NM1
         K1 = K + 1
C     ********** SET BOTTOM ELEMENT IN K-TH COLUMN OF A REAL **********
         IF (AI(N,K) .EQ. ZERO) GO TO 105
         R = CABS(CMPLX(AR(N,K),AI(N,K)))
         U1 = AR(N,K) / R
         U1I = AI(N,K) / R
         AR(N,K) = R
         AI(N,K) = 0.0
C
         DO 103 J = K1, N
            XI = U1 * AI(N,J) - U1I * AR(N,J)
            AR(N,J) = U1 * AR(N,J) + U1I * AI(N,J)
            AI(N,J) = XI
  103    CONTINUE
C
         XI = U1 * BI(N,N) - U1I * BR(N,N)
         BR(N,N) = U1 * BR(N,N) + U1I * BI(N,N)
         BI(N,N) = XI
  105    IF (K .EQ. NM1) GO TO 170
         NK1 = NM1 - K
C     ********** FOR L=N-1 STEP -1 UNTIL K+1 DO -- **********
         DO 150 LB = 1, NK1
            L = N - LB
            L1 = L + 1
C     ********** ZERO A(L+1,K) **********
            S = ABS(AR(L,K)) + ABS(AI(L,K)) + AR(L1,K)
            IF (S .EQ. ZERO) GO TO 150
            U1 = AR(L,K) / S
            U1I = AI(L,K) / S
            U2 = AR(L1,K) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            AR(L,K) = R * S
            AI(L,K) = 0.0
            AR(L1,K) = 0.0
C
            DO 110 J = K1, N
               XR = AR(L,J)
               XI = AI(L,J)
               YR = AR(L1,J)
               YI = AI(L1,J)
               AR(L,J) = U1 * XR + U1I * XI + U2 * YR
               AI(L,J) = U1 * XI - U1I * XR + U2 * YI
               AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110       CONTINUE
C
            XR = BR(L,L)
            BR(L,L) = U1 * XR
            BI(L,L) = -U1I * XR
            BR(L1,L) = -U2 * XR
C
            DO 120 J = L1, N
               XR = BR(L,J)
               XI = BI(L,J)
               YR = BR(L1,J)
               YI = BI(L1,J)
               BR(L,J) = U1 * XR + U1I * XI + U2 * YR
               BI(L,J) = U1 * XI - U1I * XR + U2 * YI
               BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
               BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  120       CONTINUE
C     ********** ZERO B(L+1,L) **********
            S = ABS(BR(L1,L1)) + ABS(BI(L1,L1)) + ABS(BR(L1,L))
            IF (S .EQ. ZERO) GO TO 150
            U1 = BR(L1,L1) / S
            U1I = BI(L1,L1) / S
            U2 = BR(L1,L) / S
            R = SQRT(U1*U1+U1I*U1I+U2*U2)
            U1 = U1 / R
            U1I = U1I / R
            U2 = U2 / R
            BR(L1,L1) = R * S
            BI(L1,L1) = 0.0
            BR(L1,L) = 0.0
C
            DO 130 I = 1, L
               XR = BR(I,L1)
               XI = BI(I,L1)
               YR = BR(I,L)
               YI = BI(I,L)
               BR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               BI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               BR(I,L) = U1 * YR - U1I * YI - U2 * XR
               BI(I,L) = U1 * YI + U1I * YR - U2 * XI
  130       CONTINUE
C
            DO 140 I = 1, N
               XR = AR(I,L1)
               XI = AI(I,L1)
               YR = AR(I,L)
               YI = AI(I,L)
               AR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               AI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               AR(I,L) = U1 * YR - U1I * YI - U2 * XR
               AI(I,L) = U1 * YI + U1I * YR - U2 * XI
  140       CONTINUE
C
            IF (.NOT. MATZ) GO TO 150
C
            DO 145 I = 1, N
               XR = ZR(I,L1)
               XI = ZI(I,L1)
               YR = ZR(I,L)
               YI = ZI(I,L)
               ZR(I,L1) = U1 * XR + U1I * XI + U2 * YR
               ZI(I,L1) = U1 * XI - U1I * XR + U2 * YI
               ZR(I,L) = U1 * YR - U1I * YI - U2 * XR
               ZI(I,L) = U1 * YI + U1I * YR - U2 * XI
  145       CONTINUE
C
  150    CONTINUE
C
  160 CONTINUE
C
  170 RETURN
C     ********** LAST CARD OF CQZHES **********
      END
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE CQZVAL(NM,N,AR,AI,BR,BI,EPS1,ALFR,ALFI,BETA,           
     X                                       MATZ,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,N,EN,K1,K2,LL,L1,NA,NM,ITS,KM1,LM1,
     X        ENM2,IERR,LOR1,ENORN
      REAL AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL R,S,A1,A2,EP,SH,U1,U2,XI,XR,YI,YR,ANI,A1I,A33,A34,A43,A44,
     X       BNI,B11,B33,B44,SHI,U1I,A33I,A34I,A43I,A44I,B33I,B44I,
     X       EPSA,EPSB,EPS1,ANORM,BNORM,B3344,B3344I
      REAL SQRT,CABS,ABS
      INTEGER MAX0
      LOGICAL MATZ
      COMPLEX Z3
      COMPLEX CSQRT,CMPLX
      REAL REAL,AIMAG
C
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF STEPS 2 AND 3 OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART,
C     AS MODIFIED IN TECHNICAL NOTE NASA TN E-7305(1973) BY WARD.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES, ONE OF THEM
C     IN UPPER HESSENBERG FORM AND THE OTHER IN UPPER TRIANGULAR FORM,
C     THE HESSENBERG MATRIX MUST FURTHER HAVE REAL SUBDIAGONAL ELEMENTS.
C     IT REDUCES THE HESSENBERG MATRIX TO TRIANGULAR FORM USING
C     UNITARY TRANSFORMATIONS WHILE MAINTAINING THE TRIANGULAR FORM
C     OF THE OTHER MATRIX AND FURTHER MAKING ITS DIAGONAL ELEMENTS
C     REAL AND NON-NEGATIVE.  IT THEN RETURNS QUANTITIES WHOSE RATIOS
C     GIVE THE GENERALIZED EIGENVALUES.  IT IS USUALLY PRECEDED BY
C     CQZHES  AND POSSIBLY FOLLOWED BY  CQZVEC.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER HESSENBERG MATRIX
C          WITH REAL SUBDIAGONAL ELEMENTS,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        EPS1 IS A TOLERANCE USED TO DETERMINE NEGLIGIBLE ELEMENTS.
C          EPS1 = 0.0 (OR NEGATIVE) MAY BE INPUT, IN WHICH CASE AN
C          ELEMENT WILL BE NEGLECTED ONLY IF IT IS LESS THAN ROUNDOFF
C          ERROR TIMES THE NORM OF ITS MATRIX.  IF THE INPUT EPS1 IS
C          POSITIVE, THEN AN ELEMENT WILL BE CONSIDERED NEGLIGIBLE
C          IF IT IS LESS THAN EPS1 TIMES THE NORM OF ITS MATRIX.  A
C          POSITIVE VALUE OF EPS1 MAY RESULT IN FASTER EXECUTION,
C          BUT LESS ACCURATE RESULTS,
C
C        MATZ SHOULD BE SET TO .TRUE. IF THE RIGHT HAND TRANSFORMATIONS
C          ARE TO BE ACCUMULATED FOR LATER USE IN COMPUTING
C          EIGENVECTORS, AND TO .FALSE. OTHERWISE,
C
C        Z=(ZR,ZI) CONTAINS, IF MATZ HAS BEEN SET TO .TRUE., THE
C          TRANSFORMATION MATRIX PRODUCED IN THE REDUCTION
C          BY  CQZHES, IF PERFORMED, OR ELSE THE IDENTITY MATRIX.
C          IF MATZ HAS BEEN SET TO .FALSE., Z IS NOT REFERENCED.
C
C     ON OUTPUT-
C
C        A HAS BEEN REDUCED TO UPPER TRIANGULAR FORM.  THE ELEMENTS
C          BELOW THE MAIN DIAGONAL HAVE BEEN SET TO ZERO,
C
C        B IS STILL IN UPPER TRIANGULAR FORM, ALTHOUGH ITS ELEMENTS
C          HAVE BEEN ALTERED.  IN PARTICULAR, ITS DIAGONAL HAS BEEN SET
C          REAL AND NON-NEGATIVE.  THE LOCATION BR(N,1) IS USED TO
C          STORE EPS1 TIMES THE NORM OF B FOR LATER USE BY  CQZVEC,
C
C        ALFR AND ALFI CONTAIN THE REAL AND IMAGINARY PARTS OF THE
C          DIAGONAL ELEMENTS OF THE TRIANGULARIZED A MATRIX,
C
C        BETA CONTAINS THE REAL NON-NEGATIVE DIAGONAL ELEMENTS OF THE
C          CORRESPONDING B.  THE GENERALIZED EIGENVALUES ARE THEN
C          THE RATIOS ((ALFR+I*ALFI)/BETA),
C
C        Z CONTAINS THE PRODUCT OF THE RIGHT HAND TRANSFORMATIONS
C          (FOR BOTH STEPS) IF MATZ HAS BEEN SET TO .TRUE.,
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF AR(J,J-1) HAS NOT BECOME
C                     ZERO AFTER 50 ITERATIONS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      ZTWO = 2.0
      ZONE = 1.0
      ZERO = 0.0

      IERR = 0
C     ********** COMPUTE EPSA,EPSB **********
      ANORM = 0.0
      BNORM = 0.0
C
      DO 30 I = 1, N
         ANI = 0.0
         IF (I .NE. 1) ANI = ABS(AR(I,I-1))
         BNI = 0.0
C
         DO 20 J = I, N
            ANI = ANI + ABS(AR(I,J)) + ABS(AI(I,J))
            BNI = BNI + ABS(BR(I,J)) + ABS(BI(I,J))
   20    CONTINUE
C
         IF (ANI .GT. ANORM) ANORM = ANI
         IF (BNI .GT. BNORM) BNORM = BNI
   30 CONTINUE
C
      IF (ANORM .EQ. ZERO) ANORM = 1.0
      IF (BNORM .EQ. ZERO) BNORM = 1.0
      EP = EPS1
      IF (EP .GT. ZERO) GO TO 50
C     ********** COMPUTE ROUNDOFF LEVEL IF EPS1 IS ZERO **********
      EP = 1.0
   40 EP = EP / ZTWO
      IF (ZONE + EP .GT. ZONE) GO TO 40
   50 EPSA = EP * ANORM
      EPSB = EP * BNORM
C     ********** REDUCE A TO TRIANGULAR FORM, WHILE
C                KEEPING B TRIANGULAR **********
      LOR1 = 1
      ENORN = N
      EN = N
C     ********** BEGIN QZ STEP **********
   60 IF (EN .EQ. 0) GO TO 1001
      IF (.NOT. MATZ) ENORN = EN
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     ********** CHECK FOR CONVERGENCE OR REDUCIBILITY.
C                FOR L=EN STEP -1 UNTIL 1 DO -- **********
   70 DO 80 LL = 1, EN
         LM1 = EN - LL
         L = LM1 + 1
         IF (L .EQ. 1) GO TO 95
         IF (ABS(AR(L,LM1)) .LE. EPSA) GO TO 90
   80 CONTINUE
C
   90 AR(L,LM1) = 0.0
C     ********** SET DIAGONAL ELEMENT AT TOP OF B REAL **********
   95 B11 = CABS(CMPLX(BR(L,L),BI(L,L)))
      IF (B11     .EQ. ZERO) GO TO 98
      U1 = BR(L,L) / B11
      U1I = BI(L,L) / B11
C
      DO 97 J = L, ENORN
         XI = U1 * AI(L,J) - U1I * AR(L,J)
         AR(L,J) = U1 * AR(L,J) + U1I * AI(L,J)
         AI(L,J) = XI
         XI = U1 * BI(L,J) - U1I * BR(L,J)
         BR(L,J) = U1 * BR(L,J) + U1I * BI(L,J)
         BI(L,J) = XI
   97 CONTINUE
C
      BI(L,L) = 0.0
   98 IF (L .NE. EN) GO TO 100
C     ********** 1-BY-1 BLOCK ISOLATED **********
      ALFR(EN) = AR(EN,EN)
      ALFI(EN) = AI(EN,EN)
      BETA(EN) = B11
      EN = NA
      GO TO 60
C     ********** CHECK FOR SMALL TOP OF B **********
  100 L1 = L + 1
      IF (B11 .GT. EPSB) GO TO 120
      BR(L,L) = 0.0
      S = ABS(AR(L,L)) + ABS(AI(L,L)) + ABS(AR(L1,L))
      U1 = AR(L,L) / S
      U1I = AI(L,L) / S
      U2 = AR(L1,L) / S
      R = SQRT(U1*U1+U1I*U1I+U2*U2)
      U1 = U1 / R
      U1I = U1I / R
      U2 = U2 / R
      AR(L,L) = R * S
      AI(L,L) = 0.0
C
      DO 110 J = L1, ENORN
         XR = AR(L,J)
         XI = AI(L,J)
         YR = AR(L1,J)
         YI = AI(L1,J)
         AR(L,J) = U1 * XR + U1I * XI + U2 * YR
         AI(L,J) = U1 * XI - U1I * XR + U2 * YI
         AR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         AI(L1,J) = U1 * YI + U1I * YR - U2 * XI
         XR = BR(L,J)
         XI = BI(L,J)
         YR = BR(L1,J)
         YI = BI(L1,J)
         BR(L1,J) = U1 * YR - U1I * YI - U2 * XR
         BR(L,J) = U1 * XR + U1I * XI + U2 * YR
         BI(L,J) = U1 * XI - U1I * XR + U2 * YI
         BI(L1,J) = U1 * YI + U1I * YR - U2 * XI
  110 CONTINUE
C
      LM1 = L
      L = L1
      GO TO 90
C     ********** ITERATION STRATEGY **********
  120 IF (ITS .EQ. 50) GO TO 1000
      IF (ITS .EQ. 10) GO TO 135
C     ********** DETERMINE SHIFT **********
      B33 = BR(NA,NA)
      B33I = BI(NA,NA)
      IF (CABS(CMPLX(B33,B33I)) .GE. EPSB) GO TO 122
      B33 = EPSB
      B33I = 0.0
  122 B44 = BR(EN,EN)
      B44I = BI(EN,EN)
      IF (CABS(CMPLX(B44,B44I)) .GE. EPSB) GO TO 124
      B44 = EPSB
      B44I = 0.0
  124 B3344 = B33 * B44 - B33I * B44I
      B3344I = B33 * B44I + B33I * B44
      A33 = AR(NA,NA) * B44 - AI(NA,NA) * B44I
      A33I = AR(NA,NA) * B44I + AI(NA,NA) * B44
      A34 = AR(NA,EN) * B33 - AI(NA,EN) * B33I
     X    - AR(NA,NA) * BR(NA,EN) + AI(NA,NA) * BI(NA,EN)
      A34I = AR(NA,EN) * B33I + AI(NA,EN) * B33
     X     - AR(NA,NA) * BI(NA,EN) - AI(NA,NA) * BR(NA,EN)
      A43 = AR(EN,NA) * B44
      A43I = AR(EN,NA) * B44I
      A44 = AR(EN,EN) * B33 - AI(EN,EN) * B33I - AR(EN,NA) * BR(NA,EN)
      A44I = AR(EN,EN) * B33I + AI(EN,EN) * B33 - AR(EN,NA) * BI(NA,EN)
      SH = A44
      SHI = A44I
      XR = A34 * A43 - A34I * A43I
      XI = A34 * A43I + A34I * A43
      IF (XR .EQ. ZERO .AND. XI .EQ. ZERO) GO TO 140
      YR = (A33 - SH) / 2.0
      YI = (A33I - SHI) / 2.0
      Z3 = CSQRT(CMPLX(YR**2-YI**2+XR,2.0*YR*YI+XI))
      U1 = REAL(Z3)
      U1I = AIMAG(Z3)
      IF (YR * U1 + YI * U1I .GE. ZERO) GO TO 125
      U1 = -U1
      U1I = -U1I
  125 Z3 = (CMPLX(SH,SHI) - CMPLX(XR,XI) / CMPLX(YR+U1,YI+U1I))
     X   / CMPLX(B3344,B3344I)
      SH = REAL(Z3)
      SHI = AIMAG(Z3)
      GO TO 140
C     ********** AD HOC SHIFT **********
  135 SH = AR(EN,NA) + AR(NA,ENM2)
      SHI = 0.0
C     ********** DETERMINE ZEROTH COLUMN OF A **********
  140 A1 = AR(L,L) / B11 - SH
      A1I = AI(L,L) / B11 - SHI
      A2 = AR(L1,L) / B11
      ITS = ITS + 1
      IF (.NOT. MATZ) LOR1 = L
C     ********** MAIN LOOP **********
      DO 260 K = L, NA
         K1 = K + 1
         K2 = K + 2
         KM1 = MAX0(K-1,L)
C     ********** ZERO A(K+1,K-1) **********
         IF (K .EQ. L) GO TO 170
         A1 = AR(K,KM1)
         A1I = AI(K,KM1)
         A2 = AR(K1,KM1)
  170    S = ABS(A1) + ABS(A1I) + ABS(A2)
         U1 = A1 / S
         U1I = A1I / S
         U2 = A2 / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
C
         DO 180 J = KM1, ENORN
            XR = AR(K,J)
            XI = AI(K,J)
            YR = AR(K1,J)
            YI = AI(K1,J)
            AR(K,J) = U1 * XR + U1I * XI + U2 * YR
            AI(K,J) = U1 * XI - U1I * XR + U2 * YI
            AR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            AI(K1,J) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(K,J)
            XI = BI(K,J)
            YR = BR(K1,J)
            YI = BI(K1,J)
            BR(K,J) = U1 * XR + U1I * XI + U2 * YR
            BI(K,J) = U1 * XI - U1I * XR + U2 * YI
            BR(K1,J) = U1 * YR - U1I * YI - U2 * XR
            BI(K1,J) = U1 * YI + U1I * YR - U2 * XI
  180    CONTINUE
C
         IF (K .EQ. L) GO TO 240
         AI(K,KM1) = 0.0
         AR(K1,KM1) = 0.0
         AI(K1,KM1) = 0.0
C     ********** ZERO B(K+1,K) **********
  240    S = ABS(BR(K1,K1)) + ABS(BI(K1,K1)) + ABS(BR(K1,K))
         U1 = BR(K1,K1) / S
         U1I = BI(K1,K1) / S
         U2 = BR(K1,K) / S
         R = SQRT(U1*U1+U1I*U1I+U2*U2)
         U1 = U1 / R
         U1I = U1I / R
         U2 = U2 / R
         IF (K .EQ. NA) GO TO 245
         XR = AR(K2,K1)
         AR(K2,K1) = U1 * XR
         AI(K2,K1) = -U1I * XR
         AR(K2,K) = -U2 * XR
C
  245    DO 250 I = LOR1, K1
            XR = AR(I,K1)
            XI = AI(I,K1)
            YR = AR(I,K)
            YI = AI(I,K)
            AR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            AI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            AR(I,K) = U1 * YR - U1I * YI - U2 * XR
            AI(I,K) = U1 * YI + U1I * YR - U2 * XI
            XR = BR(I,K1)
            XI = BI(I,K1)
            YR = BR(I,K)
            YI = BI(I,K)
            BR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            BI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            BR(I,K) = U1 * YR - U1I * YI - U2 * XR
            BI(I,K) = U1 * YI + U1I * YR - U2 * XI
  250    CONTINUE
C
         BI(K1,K1) = 0.0
         BR(K1,K) = 0.0
         BI(K1,K) = 0.0
         IF (.NOT. MATZ) GO TO 260
C
         DO 255 I = 1, N
            XR = ZR(I,K1)
            XI = ZI(I,K1)
            YR = ZR(I,K)
            YI = ZI(I,K)
            ZR(I,K1) = U1 * XR + U1I * XI + U2 * YR
            ZI(I,K1) = U1 * XI - U1I * XR + U2 * YI
            ZR(I,K) = U1 * YR - U1I * YI - U2 * XR
            ZI(I,K) = U1 * YI + U1I * YR - U2 * XI
  255    CONTINUE
C
  260 CONTINUE
C     ********** SET LAST A SUBDIAGONAL REAL AND END QZ STEP **********
      IF (AI(EN,NA) .EQ. ZERO) GO TO 70
      R = CABS(CMPLX(AR(EN,NA),AI(EN,NA)))
      U1 = AR(EN,NA) / R
      U1I = AI(EN,NA) / R
      AR(EN,NA) = R
      AI(EN,NA) = 0.0
C
      DO 270 J = EN, ENORN
         XI = U1 * AI(EN,J) - U1I * AR(EN,J)
         AR(EN,J) = U1 * AR(EN,J) + U1I * AI(EN,J)
         AI(EN,J) = XI
         XI = U1 * BI(EN,J) - U1I * BR(EN,J)
         BR(EN,J) = U1 * BR(EN,J) + U1I * BI(EN,J)
         BI(EN,J) = XI
  270 CONTINUE
C
      GO TO 70
C     ********** SET ERROR -- BOTTOM SUBDIAGONAL ELEMENT HAS NOT
C                BECOME NEGLIGIBLE AFTER 50 ITERATIONS **********
 1000 IERR = EN
C     ********** SAVE EPSB FOR USE BY CQZVEC **********
 1001 IF (N .GT. 1) BR(N,1) = EPSB
      RETURN
C     ********** LAST CARD OF CQZVAL **********
      END
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE CQZVEC(NM,N,AR,AI,BR,BI,ALFR,ALFI,BETA,ZR,ZI)          
C
      INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN
      REAL AR(NM,N),AI(NM,N),BR(NM,N),BI(NM,N),ALFR(N),ALFI(N),
     X       BETA(N),ZR(NM,N),ZI(NM,N)
      REAL R,T,RI,TI,XI,ALMI,ALMR,BETM,EPSB
      REAL CABS
      COMPLEX Z3
      COMPLEX CMPLX
      REAL REAL,AIMAG
C
C
C
C
C
C     THIS SUBROUTINE IS A COMPLEX ANALOGUE OF THE FOURTH STEP OF THE
C     QZ ALGORITHM FOR SOLVING GENERALIZED MATRIX EIGENVALUE PROBLEMS,
C     SIAM J. NUMER. ANAL. 10, 241-256(1973) BY MOLER AND STEWART.
C
C     THIS SUBROUTINE ACCEPTS A PAIR OF COMPLEX MATRICES IN UPPER
C     TRIANGULAR FORM, WHERE ONE OF THEM FURTHER MUST HAVE REAL DIAGONAL
C     ELEMENTS.  IT COMPUTES THE EIGENVECTORS OF THE TRIANGULAR PROBLEM
C     AND TRANSFORMS THE RESULTS BACK TO THE ORIGINAL COORDINATE SYSTEM.
C     IT IS USUALLY PRECEDED BY  CQZHES  AND  CQZVAL.
C
C     ON INPUT-
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT,
C
C        N IS THE ORDER OF THE MATRICES,
C
C        A=(AR,AI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX,
C
C        B=(BR,BI) CONTAINS A COMPLEX UPPER TRIANGULAR MATRIX WITH REAL
C          DIAGONAL ELEMENTS.  IN ADDITION, LOCATION BR(N,1) CONTAINS
C          THE TOLERANCE QUANTITY (EPSB) COMPUTED AND SAVED IN  CQZVAL,
C
C        ALFR, ALFI, AND BETA ARE VECTORS WITH COMPONENTS WHOSE
C          RATIOS ((ALFR+I*ALFI)/BETA) ARE THE GENERALIZED
C          EIGENVALUES.  THEY ARE USUALLY OBTAINED FROM  CQZVAL,
C
C        Z=(ZR,ZI) CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTIONS BY  CQZHES  AND  CQZVAL, IF PERFORMED.
C          IF THE EIGENVECTORS OF THE TRIANGULAR PROBLEM ARE
C          DESIRED, Z MUST CONTAIN THE IDENTITY MATRIX.
C
C     ON OUTPUT-
C
C        A IS UNALTERED,
C
C        B HAS BEEN DESTROYED,
C
C        ALFR, ALFI, AND BETA ARE UNALTERED,
C
C        Z CONTAINS THE EIGENVECTORS.  EACH EIGENVECTOR IS NORMALIZED
C          SO THAT THE MODULUS OF ITS LARGEST COMPONENT IS 1.0 .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
C
C     ------------------------------------------------------------------
C
      ZERO = 0.0
      IF (N .LE. 1) GO TO 1001
      EPSB = BR(N,1)
C     ********** FOR EN=N STEP -1 UNTIL 2 DO -- **********
      DO 800 NN = 2, N
         EN = N + 2 - NN
         NA = EN - 1
         ALMR = ALFR(EN)
         ALMI = ALFI(EN)
         BETM = BETA(EN)
C     ********** FOR I=EN-1 STEP -1 UNTIL 1 DO -- **********
         DO 700 II = 1, NA
            I = EN - II
            R = 0.0
            RI = 0.0
            M = I + 1
C
            DO 610 J = M, EN
               T = BETM * AR(I,J) - ALMR * BR(I,J) + ALMI * BI(I,J)
               TI = BETM * AI(I,J) - ALMR * BI(I,J) - ALMI * BR(I,J)
               IF (J .EQ. EN) GO TO 605
               XI = T * BI(J,EN) + TI * BR(J,EN)
               T = T * BR(J,EN) - TI * BI(J,EN)
               TI = XI
  605          R = R + T
               RI = RI + TI
  610       CONTINUE
C
            T = ALMR * BETA(I) - BETM * ALFR(I)
            TI = ALMI * BETA(I) - BETM * ALFI(I)
            IF (T .EQ. ZERO .AND. TI .EQ. ZERO) T = EPSB
            Z3 = CMPLX(R,RI) / CMPLX(T,TI)
            BR(I,EN) = REAL(Z3)
            BI(I,EN) = AIMAG(Z3)
  700    CONTINUE
C
  800 CONTINUE
C     ********** END BACK SUBSTITUTION.
C                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
C                FOR J=N STEP -1 UNTIL 2 DO -- **********
      DO 880 JJ = 2, N
         J = N + 2 - JJ
         M = J - 1
C
         DO 880 I = 1, N
C
            DO 860 K = 1, M
               ZR(I,J) = ZR(I,J) + ZR(I,K) * BR(K,J) - ZI(I,K) * BI(K,J)
               ZI(I,J) = ZI(I,J) + ZR(I,K) * BI(K,J) + ZI(I,K) * BR(K,J)
  860       CONTINUE
C
  880 CONTINUE
C     ********** NORMALIZE SO THAT MODULUS OF LARGEST
C                COMPONENT OF EACH VECTOR IS 1 **********
      DO 950 J = 1, N
         T = 0.0
C
         DO 930 I = 1, N
            R = CABS(CMPLX(ZR(I,J),ZI(I,J)))
            IF (R .GT. T) T = R
  930    CONTINUE
C
         DO 940 I = 1, N
            ZR(I,J) = ZR(I,J) / T
            ZI(I,J) = ZI(I,J) / T
  940    CONTINUE
C
  950 CONTINUE
C
 1001 RETURN
C     ********** LAST CARD OF CQZVEC **********
      END



















!======================================================================|
