C======================================================================|
      subroutine NEOCL1(YBS,YCC,YHE,YXI)
C----------------------------------------------------------------------|
      implicit none
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      INCLUDE 'sbr/pamx_mi.inc'
      INCLUDE 'sbr/pamx_ms.inc'
      INCLUDE 'sbr/pamx_mz.inc'
      double precision YBS(NA1),YCC(NA1),YHE(NA1),YXI(NA1),ddum(4)
      double precision GRAD, TPF, yth2, yh, grti, y_den, z_coulomb
      double precision z_j7kv, ybtot, ybtor, ybpol, yuthai, yppr, p_eps
      REAL           rdum
!Declaration of input to NCLASS
      INTEGER        k_order,                 k_potato
      INTEGER        m_i,                     m_z
      REAL           c_den,                   c_potb,
     #               c_potl
      REAL           p_b2,                    p_bm2,
     #               p_eb,                    p_fhat,
     #               p_fm(3),                 p_ft,
     #               p_grbm2,                 p_grphi,
     #               p_gr2phi,                p_ngrth
      REAL           amu_i(mx_mi),            grt_i(mx_mi),
     #               temp_i(mx_mi)
      REAL           den_iz(mx_mi,mx_mz),     fex_iz(3,mx_mi,mx_mz),
     #               grp_iz(mx_mi,mx_mz)
!Declaration of output from NCLASS
      INTEGER        iflag,                   m_s
      INTEGER        jm_s(mx_ms),             jz_s(mx_ms)
      REAL           p_bsjb,                  p_etap,
     #               p_exjb
      REAL           calm_i(3,3,mx_mi)
      REAL           caln_ii(3,3,mx_mi,mx_mi),capm_ii(3,3,mx_mi,mx_mi),
     #               capn_ii(3,3,mx_mi,mx_mi)
      REAL           bsjbp_s(mx_ms),          bsjbt_s(mx_ms),
     #               dn_s(mx_ms),             gfl_s(5,mx_ms),
     #               qfl_s(5,mx_ms),          sqz_s(mx_ms),
     #               upar_s(3,3,mx_ms),       utheta_s(3,3,mx_ms),
     #               vn_s(mx_ms),             veb_s(mx_ms),
     #               qeb_s(mx_ms),            xi_s(mx_ms),
     #               ymu_s(3,3,mx_ms)
      REAL           chip_ss(mx_ms,mx_ms),    chit_ss(mx_ms,mx_ms),
     #               dp_ss(mx_ms,mx_ms),      dt_ss(mx_ms,mx_ms)
!Declaration of local variables
      CHARACTER      label*120
      INTEGER        i,                       j,
     #               jj,                      j1,
     #               k_out
      INTEGER        idum(8),                 nout
C----------------------------------------------------------------------|
C used for output only
	integer	im, iz, iza
C----------------------------------------------------------------------|

      if (mx_mi .lt. 6)	then
         write(*,*)"NEOCL warning: call ignored. mx_mi < 6"
         return
      endif
C Set control and radially independent data
!  k_out-option for output to nout [-]
!       =1 errors only
!       =2 errors and results 
!       =else no output
!  k_order-order of v moments to be solved [-]
!         =2 u and p_q
!         =3 u, p_q, and u2
!         =else error
!  k_potato-option to include potato orbits [-]
!          =0 off
!          =else on
      k_out      = 1
      k_order    = 2
      k_potato   = 1
!  c_den-density cutoff below which species is ignored (/m**3)
!  c_potb-kappa(0)*Bt(0)/[2*q(0)**2] (T)
!  c_potl-q(0)*R(0) (m)
      c_den      = 1.0e10
      y_den      = 1.d-19*c_den                  ! c_den in ASTRA units
      c_potb     = -0.5*ELON(1)*BTOR*MU(1)**2
C Parameters:
C  mx_mi=9  - max number of isotopes
C  mx_mz=18 - max charge
C  mx_ms=40 - max number of species
!  m_i-number of isotopes (1<mi<mx_mi+1)
!  m_z-highest charge state of all species (0<mz<mx_mz+1)
!  grt_i(i)-temperature gradient of i (keV/rho)
!  grp_iz(i,z)-pressure gradient of i,z (keV/m**3/rho)
C----------------------------------------------------------------------|
      do	j=1,NA
C Set radially dependent data
!  p_eps-inverse aspect ratio [-]
!  p_grphi-radial electric field Phi' (V/rho)
!  p_gr2phi-radial electric field gradient Psi'(Phi'/Psi')' (V/rho**2)
!  p_q-safety factor [-] (-1./MU(j))
!  p_eb-<E.B> (V*T/m)
!  p_b2-<B**2> (T**2)
!  p_bm2-<1/B**2> (/T**2)
!  p_fhat-mu_0*F/(dPsi/dr) (rho/m)
!  p_ft-trapped fraction [-]
!  p_grbm2-<grad(rho)**2/B**2> (rho**2/m**2/T**2)
!  p_ngrth-<n.grad(Theta)> (1/m)
         c_potl   = -(RTOR+SHIF(j))/MU(1)
         p_eps    = HRO*j/ROC*ABC/(RTOR+SHIF(j))
         p_grphi  = ER(j)
         p_gr2phi = grad(ER,j)
C Warning: The NCLASS version 1.2 returns NaN resistivity
C          if p_eb = 0.
C         p_eb     = max(.00001,-ULON(j)*BTOR/(GP2*RTOR))
         p_eb     = -ULON(j)*BTOR/(GP2*RTOR)
         if (abs(p_eb) .lt. 1.d-3)	p_eb = sign(0.001,p_eb)
         YTH2     = RHO(j)*G22(j)*(MU(j)/RTOR)**2
         p_b2     = (1.+YTH2)*G33(j)*(BTOR*IPOL(j))**2

C Approximate expressions are used:
         p_bm2   = (1.0+1.5*p_eps**2)/BTOR**2
         p_fhat	 = -1./(MU(j)*p_eps)
	 include 'fml/tpf'			! large R/a, low pressure
         p_ft	 = TPF
         p_grbm2 = 1.0/BTOR**2
         p_ngrth = -MU(j)/(RTOR+SHIF(j))
!  p_fm(3)-poloidal moments of geometric factor for PS viscosity [-]
         do i=1,3
            p_fm(i)=0.0
         enddo
         do i=1,3
            p_fm(i)=i*((1.0-SQRT(1.0-p_eps**2))/p_eps)**(2.0*i)
     #           *(1.0+i*SQRT(1.0-p_eps**2))/((1.0-p_eps**2)**1.5
     #           )*p_ngrth**2
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
         if (j .lt. NA)	then
            YH = HRO
         else
            YH = HROA
         endif

         m_i = 1                ! Reserved for electrons 
         m_z = 1
         amu_i(m_i)  = 5.4463e-4
         temp_i(m_i) = TE(j)
         grt_i(m_i)  = (TE(j+1)-TE(j))/YH
         grti	     = (TI(j+1)-TI(j))/YH
         den_iz(m_i,m_z) = 1.e19*NE(j)
         grp_iz(m_i,m_z) = 1.e19*(TE(j+1)*NE(j+1)-TE(j)*NE(j))/YH

         if (NHYDR(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = 1.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            den_iz(m_i,m_z)=1.e19*NHYDR(j)
            grp_iz(m_i,m_z)=1.e19*(TI(j+1)*NHYDR(j+1)-TI(j)*NHYDR(j))/YH
         endif
         if (NDEUT(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = 2.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            den_iz(m_i,m_z) = 1.e19*NDEUT(j)
            grp_iz(m_i,m_z) = 1.e19*(TI(j+1)*NDEUT(j+1)-TI(j)*NDEUT(j))
     &           /YH
         endif
         if (NTRIT(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = 3.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            den_iz(m_i,m_z) = 1.e19*NTRIT(j)
            grp_iz(m_i,m_z) = 1.e19*(TI(j+1)*NTRIT(j+1)-TI(j)*NTRIT(j))
     &           /YH
         endif

         if (m_i .eq. 1 .and. mx_mz .lt. 2 .and. ZMJ .ge. 2 )	then
            write(*,*)">>> NEOCL: ion composition error:"
            write(*,*)" mx_mz < 2 and no hydrogen. Call ignored."
            return
         endif
         if (NHE3(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = 3.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            m_z = 2
            den_iz(m_i,m_z) = 1.e19*NHE3(j)
            grp_iz(m_i,m_z) = 1.e19*(TI(j+1)*NHE3(j+1)-TI(j)*NHE3(j))/YH
         endif
         if (NALF(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = 4.
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            m_z = 2
            den_iz(m_i,m_z) = 1.e19*NALF(j)
            grp_iz(m_i,m_z) = 1.e19*(TI(j+1)*NALF(j+1)-TI(j)*NALF(j))/YH
         endif

C It is assumed here that the main ion species is one of H,D,T,He3,He4
         if (m_i .eq. 1)	then ! No specification for main ions
            m_i = 2 
            amu_i(m_i) = AMAIN(j)
            temp_i(m_i) = TI(j)
            grt_i(m_i) = grti
            den_iz(m_i,m_z) = 1.e19*NI(j)
            grp_iz(m_i,m_z) = 1.e19*(TI(j+1)*NI(j+1)-TI(j)*NI(j))/YH
         endif

         if (m_i .eq. mx_mi)	goto	5

C consider up to three impurity species
         jj = ZIM1(j)+0.5
         if (jj .gt. mx_mz)	goto	3
         if (NIZ1(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = AIM1
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj .gt. m_z)	m_z = jj
            den_iz(m_i,jj) = 1.e19*NIZ1(j)
            grp_iz(m_i,jj) = 1.e19*(TI(j+1)*NIZ1(j+1)-TI(j)*NIZ1(j))/YH
         endif
         if (m_i .eq. mx_mi)	goto	5
 3       jj = ZIM2(j)+0.5
         if (jj .gt. mx_mz)	goto	4
         if (NIZ2(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = AIM2
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj .gt. m_z)	m_z = jj
            den_iz(m_i,jj) = 1.e19*NIZ2(j)
            grp_iz(m_i,jj) = 1.e19*(TI(j+1)*NIZ2(j+1)-TI(j)*NIZ2(j))/YH
         endif
         if (m_i .eq. mx_mi)	goto	5
 4       jj = ZIM3(j)+0.5
         if (jj .gt. mx_mz)	goto	5
         if (NIZ3(j) .gt. y_den)	then
            m_i = m_i+1
            amu_i(m_i) = AIM3
            grt_i(m_i) = grti
            temp_i(m_i) = TI(j)
            if (jj .gt. m_z)	m_z = jj
            den_iz(m_i,jj) = 1.e19*NIZ3(j)
            grp_iz(m_i,jj) = 1.e19*(TI(j+1)*NIZ3(j+1)-TI(j)*NIZ3(j))/YH
         endif
 5       continue

         CALL NCLASS(k_order,k_potato,m_i,m_z,c_den,c_potb,c_potl,p_b2,
     #        p_bm2,p_eb,p_fhat,p_fm,p_ft,p_grbm2,p_grphi,p_gr2phi,
     #        p_ngrth,amu_i,grt_i,temp_i,den_iz,fex_iz,grp_iz,
C output:
     #        m_s,jm_s,jz_s,p_bsjb,p_etap,p_exjb,calm_i,caln_ii,capm_ii,
     #        capn_ii,bsjbp_s,bsjbt_s,dn_s,gfl_s,qfl_s,sqz_s,upar_s,
     #        utheta_s,vn_s,veb_s,qeb_s,xi_s,ymu_s,chip_ss,chit_ss,
     #        dp_ss,dt_ss,iflag)

         IF (iflag.gt.0 .and. k_out.gt.0) THEN
            nout=6
C            OPEN(unit=nout,file='out_nclass_pt.dat',status='unknown')
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
            close(nout)
            GOTO 1000
         ENDIF
C Note! All quantities are returned on the auxiliary grid
C                (wrong for BS & CC)
C            (j) = -1.e-6*p_exjb/BTOR	! Driven current density
         YBS(j) = -1.e-6*p_bsjb/BTOR	! Bootstrap current density
         YCC(j) = 1.e-6/p_etap		! Conductivity

!  jm_s(s)-isotope number of s (-)
!  jz_s(s)-charge state of s (-)
!  upar_s(3,m,s)-parallel flow of s from force m (T*m/s)
!                m=1, p', T', Phi'
!                m=2, <E.B>
!                m=3, fex_iz
!  utheta_s(3,m,s)-poloidal flow of s from force m (m/s/T)
!                  m=1, p', T'
!                  m=2, <E.B>
!                  m=3, fex_iz
	nout = 6
!  Flow velocities on outside midplane
C        label='     *** Flow Velocities on Outside Midplane ***'
C        CALL WRITE_LINE(nout,label,2,1)
C        label='     Species       v-tor       v-pol       v-par'//
C     #        '      v-perp'
C        CALL WRITE_LINE(nout,label,0,0)
C        label='           -         m/s         m/s         m/s'//
C     #        '         m/s'
C        CALL WRITE_LINE(nout,label,0,0)                        
        ybtor=BTOR/(1.0+p_eps)
        ybpol=ybtor/p_fhat
        ybtot=SQRT(ybtor**2+ybpol**2)*ybtor/ABS(ybtor)
      z_coulomb=1.6022e-19
      z_j7kv=1.6022e-16
        DO i=1,m_s
          idum(1)=i
          im=jm_s(i)
          iz=jz_s(i)
          iza=IABS(iz)
          yppr=p_fhat*grp_iz(im,iza)*z_j7kv
     #        /(z_coulomb*iz*den_iz(im,iza))+p_fhat*p_grphi
          yuthai=utheta_s(1,1,i)+utheta_s(1,2,i)+utheta_s(1,3,i)
!         Toroidal
          ddum(1)=yuthai*ybtor-yppr/ybtor
!         Poloidal
          ddum(2)=yuthai*ybpol
!         Parallel
          ddum(3)=yuthai*ybtot-yppr/ybtot
!         Perpendicular
          ddum(4)=yppr*ybpol/ybtot/ybtor
	if (i .eq. 2)	then	! Ion component only
C----------------------------------------------------------------------|
        if(amu_i(i).gt.4.5.or.amu_i(i).lt.0.5)write(*,*)"Wrong species"
	if (iz .ne. 1) write(*,*)"Wrong species"
        YHE(j) = ddum(1) ! Toroidal ion velocity on outside midplane
        YXI(j) = ddum(2) ! Poloidal ion velocity on outside midplane
	endif
        ENDDO

         YHE(j) = chit_ss(1,1)*VRS(j)/G11(j) ! El. heat conductivity
         YXI(j) = chit_ss(2,2)*VRS(j)/G11(j) ! Ion heat conductivity
      enddo

      YBS(NA1) = YBS(NA)*(ROC/HRO-(NA-1))+YBS(NA-1)*(NA-ROC/HRO)
      YCC(NA1) = YCC(NA)*(ROC/HRO-(NA-1))+YCC(NA-1)*(NA-ROC/HRO)
      YHE(NA1) = YHE(NA)*(ROC/HRO-(NA-1))+YHE(NA-1)*(NA-ROC/HRO)
      YXI(NA1) = YXI(NA)*(ROC/HRO-(NA-1))+YXI(NA-1)*(NA-ROC/HRO)
 1000 end
C======================================================================|
