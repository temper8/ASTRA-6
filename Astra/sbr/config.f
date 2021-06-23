      subroutine config()
      implicit double precision (a-h,o-z)
      double precision con,tem,abtor,rm,x0,z0,rh1,rh,rha
     &,delta,ell,gamma,cdl,cly,cgm,cmy,coeffs,zero,drhodr
     &,amy,shift0,ell0,gamma0,con0,tem0
      include 'for/parameter.inc'
      include 'for/const.inc'
      include 'for/status.inc'
      external polin,polin1
      dimension con(NRD),tem(NRD)
      dimension rh(NRD),delta(NRD),ell(NRD),gamma(NRD)
      dimension rha(NRD),drhodr(NRD),amy(NRD)
      dimension cdl(10),cly(10),cgm(10),cmy(10),coeffs(10)

      parameter(zero=0.d0,ipsy=5)
cc*********************************************************************
cc   ipsy = number of polinomial decomposition coefficients
cc           used for interpolation of Zakharov's moments.
cc*********************************************************************
!!!!!   Arrays for testing:
c      parameter(np=11, npp=101)
c      double precision arho,adelta,aell,agamma,ddelta,dell,dgamma,fmy,dfmy
c     &,tt,xx,zz,xxx,zzz,bpzakh,tetfiz,fiztet,fdf
c      dimension arho(np),adelta(np),aell(np),agamma(np)
c      dimension ddelta(np),dell(np),dgamma(np),fmy(np)
cc*********************************************************************

cc*********************************************************************
cc    Co-ordinates used in ray-tracing:
cc         (x-x0)/rm=r*cos(teta)-delta+gamma*sin^2(teta)
cc         (z-z0)/rm=ell*r*sin(teta)
cc    Definitions:
cc    (x0,z0) - magnetic axis position, centimeters
cc    rm      - minor radius in mid-plane, cenrimeters
cc    r(rho_ASTRA),delta(r),gamma(r),ell(r) - dimensionless functions
cc    rho_ASTRA=sqrt(Phi_tor/GP/BTOR)
cc    Interval for r:  0.<= r <=1.
cc*********************************************************************

      UPDWN=0.
      ipsy1=ipsy-1
      inpt=NA1          ! ASTRA radial grid number
      do i=1,inpt
       rh(i)=dble(AMETR(i)/ABC)
       rha(i)=dble((i-0.5)*HRO/ABC) ! /ABC instead of /ROC is not a mistake!
       delta(i)=dble(SHIF(i))
       ell(i)=dble(ELON(i))
       gamma(i)=dble(TRIA(i))
       con(i)=dble(NE(i))
       tem(i)=dble(TE(i))
      end do

      call approx(rh,delta,inpt,polin,ipsy,coeffs)
      shift0=coeffs(1)  !SHIF extrapolation to the magnetic axis
      call approx(rh,ell,inpt,polin,ipsy,coeffs)
      ell0=coeffs(1)    !ELON extrapolation to the magnetic axis
      call approx(rh,gamma,inpt,polin,ipsy,coeffs)
      gamma0=coeffs(1)  !TRIA extrapolation to the magnetic axis
      call approx(rh,con,inpt,polin,ipsy,coeffs)
      con0=coeffs(1)    !NE extrapolation to the magnetic axis
      call approx(rh,tem,inpt,polin,ipsy,coeffs)
      tem0=coeffs(1)    !TE extrapolation to the magnetic axis

      do i=2,inpt
       delta(i)=dble((shift0-SHIF(i))/ABC)  !LHCD definition of Shafr. shift
      end do

      rh1=rh(1)          !saving the first ASTRA radial grid element
      rh(1)=zero         !shifting the first element to zero
      rha(1)=zero        !shifting the first element to zero
      delta(1)=zero      !putting delta(rh=0.)=0.
      ell(1)=ell0        !putting ell(rh=0.)=ell0
      gamma(1)=gamma0    !putting gamma(rh=0.)=gamma0
      con(1)=con0        !putting con(rh=0.)=con0
      tem(1)=tem0        !putting tem(rh=0.)=tem0

      abtor=1.d4*dble(BTOR*RTOR/(RTOR+shift0)) !B_tor_(magnetic axis), Gauss
      rm=dble(100.*ABC)                      !minor radius in mid-plane, cm
      x0=dble(100.*(RTOR+shift0))     !x-coordinate of the magnetic axis, cm
      z0=dble(100.*UPDWN)             !z-coordinate of the magnetic axis, cm

cccc   shift as a function of "minor radius":
        call approx(rh,delta,inpt,polin1,ipsy1,coeffs)
          cdl(1)=zero
          do k=2,ipsy
            cdl(k)=coeffs(k-1)
          end do

cccc   triangularity as a function of "minor radius":
        call approx(rh,gamma,inpt,polin1,ipsy1,coeffs)
          cgm(1)=zero
          do k=2,ipsy
            cgm(k)=coeffs(k-1)
          end do

cccc   ellipticity as a function of "minor radius":
        call approx(rh,ell,inpt,polin,ipsy,cly)

cccc   "poloidal magnetic field":
       call diff(rh,rha,inpt,drhodr)
       do i=2,inpt
        amy(i)=1.d4*dble(BTOR*MU(i))*rha(i)*drhodr(i)
       end do
       amy(1)=zero
!! amy=(btor/q)*rho*(drho/dr) is a function of "minor radius" r=rh(i).
!! Poloidal magnetic field: B_pol=amy(r)*sqrt(g22/g), where g is
!! determinant of 3D metric tensor and g22 is the (22) element of
!! the tensor, normalized on ABC^4 and ABC^2, correspondingly.
!!
!!  Polinomial approximation of the amy(r):
       inpt2=inpt-3
       call approx(rh,amy,inpt2,polin1,ipsy1,coeffs)
       cmy(1)=zero
        do k=2,ipsy
         cmy(k)=coeffs(k-1)
        end do

!!TEST!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c       open(97,file='dat/uhr/drhodr', status='unknown')                  c
c       do i=1,inpt                                                       c
c        write(97,*) rh(i),drhodr(i),rha(i)                               c
c       end do                                                            c
c       close(97)                                                         c
c                                                                         c
c      open(97, file='dat/uhr/coeffs', status='unknown')                  c
c       do i=1,ipsy                                                       c
c        write(97,666) cdl(i),cly(i),cgm(i),cmy(i)                        c
c       end do                                                            c
c      close(97)                                                          c
c                                                                         c
c      open(97, file='dat/uhr/moments', status='unknown')                 c
c      do i=1,np                                                          c
c        arho(i)=(i-1.)/(np-1.)                                           c
c        adelta(i)=fdf(arho(i),cdl,ipsy,ddelta(i))                        c
c        aell(i)=fdf(arho(i),cly,ipsy,dell(i))                            c
c        agamma(i)=fdf(arho(i),cgm,ipsy,dgamma(i))                        c
c        fmy(i)=fdf(arho(i),cmy,ipsy,dfmy)                                c
c        write(97,666) arho(i),adelta(i),aell(i),agamma(i),fmy(i)         c
c      end do                                                             c
c      close(97)                                                          c
c                                                                         c
c      open(97, file='dat/uhr/surfzakh', status='unknown')                c
c      open(98, file='dat/uhr/bpzakh', status='unknown')                  c
c      open(99, file='dat/uhr/bpol', status='unknown')                    c
c      do ii=2,np                                                         c
c        if(ii.gt.2)  write(97,*)                                         c
c        if(ii.gt.2)  write(98,*)                                         c
c        do j=1,npp                                                       c
c         tt=2.d0*dble(gp*(j-1)/(npp-1.))                                 c
c         xx=arho(ii)*dcos(tt)-adelta(ii)-agamma(ii)*dsin(tt)**2          c
c         zz=aell(ii)*arho(ii)*dsin(tt)                                   c
c         xxx=(x0+xx*rm)/1.d2                                             c
c         zzz=(z0+zz*rm)/1.d2                                             c
c         write(97,*) xxx,zzz                                             c
c          call bpol(arho(ii),tt,bpzakh,x0,rm,                            c
c     *     adelta(ii),ddelta(ii),aell(ii),dell(ii),                      c
c     *             agamma(ii),dgamma(ii),fmy(ii))                        c
c          tetfiz=fiztet(xx,zz)/gp                                        c
c         write(98,*) tetfiz,bpzakh                                       c
c          if(j.eq.1) then                                                c
c            write(99,*) arho(ii),bpzakh                                  c
c          end if                                                         c
c        end do                                                           c
c      end do                                                             c
c      close(97)                                                          c
c      close(98)                                                          c
c      close(99)                                                          c
cc************************************************************************

cc************************************************************************
c   Saving data necessary for ray-tracing
c   This file can be used for running ray-tracing without ASTRA
c   for a given plasma slice
cc************************************************************************
       open(97, file='dat/uhr/input', status='unknown')
        write(97,*) inpt       !number of ASTRA grid points
        write(97,*) ipsy       !decomposition coeffs. number
        write(97,*) rm         !minor radius in mid-plane, cm
        write(97,*) x0         !magnetic axis x-coordinate, cm
        write(97,*) z0         !magnetic axis z-coordinate, cm
        write(97,*) rh1        !smallest ASTRA  rh=AMETR(1)/ABC
        write(97,*) abtor      !B_tor at magnetic axis, Gauss
       do i=1,ipsy
        write(97,*) cdl(i)     !power decomp. coeffs. for Shafr. shift
       end do
       do i=1,ipsy
        write(97,*) cly(i)     !power decomp. coeffs. for ellipticity
       end do
       do i=1,ipsy
        write(97,*) cgm(i)     !power decomp. coeffs. for triangularity
       end do
       do i=1,ipsy
        write(97,*) cmy(i)     !power decomp. coeffs. for "B_poloidal"
       end do
         do i=1,inpt
          write(97,*) rh(i)    !rh grid points, rh(1)=0.
         end do
         do i=1,inpt
          write(97,*) con(i)   !plasma density in the rh(i) grid points
         end do
         do i=1,inpt
          write(97,*) tem(i)   !el. temperature in the rh(i) grid points
         end do
      close(97)
666   format(6(e12.5,1x))
      end
