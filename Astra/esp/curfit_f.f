        subroutine prefit(rk,zk,ncpfc,NECON,WECON,rax,zax,alp_b,psi_bnd)

         include 'double.inc'
         include 'parcur.inc'
          INCLUDE 'prm.inc'
         include 'comevl.inc'
       common/com_flag/kastr


       real*8  rk(*),zk(*)
       DIMENSION   WECON(*)
       integer  necon(*),ncpfc

c   -----------------------------
        alph=alp_b
        r_ax=rax
        z_ax=zax

c   -----------------------------
         if(kastr.eq.0) then 

        write(fname,'(a,a)') path(1:kname),'bonfit.dat'
        open(1,file=fname,form='formatted')
          !open(1,file='bonfit.dat')

           read(1,*) wwl,wdk,wsig

           read(1,*) (d_wght(ik),ik=1,NEQUI)

           read(1,*) Lfit
             if(Lfit.ne.0) then
            do l=1,Lfit
           read(1,*) rfit(l),zfit(l)
            enddo
             endif
           read(1,*) Lpre
             if(Lpre.ne.0) then
            do l=1,Lpre
           read(1,*) rpre(l),zpre(l)
            enddo
             endif
           read(1,*) rx_p,zx_p
           read(1,*) psi_bnd
          close(1)
!temporary
!            do l=1,Lfit
!            rfit(l)=rfit(l)+0.02d0
!            enddo
!            do l=1,Lpre
!            rpre(l)=rpre(l)+0.02d0
!            enddo
!            rx_p=rx_p+0.02d0
!temporary

         elseif(kastr.eq.1) then 

          wwl= 1.d0
		  wdk= 1.d0
		  wsig=1.d-2

            do ik=1,NEQUI
           d_wght(ik)=1.d0
            enddo

        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
        open(1,file=fname,form='formatted')
          !open(1,file='tab_bnd.dat')

           read(1,*) Lfit

            do l=1,Lfit
           read(1,*) rfit(l),zfit(l)
            enddo

           Lpre=0

           rx_p=rax !dummy
		 zx_p=zax !dummy

         endif

              do l=1,Lfit
              w_wght(l)=wwl*2.d0
              enddo
              s_wght=wsig

           call grimat(rk,zk,ncpfc,NECON, WECON)

!        do iq=1,NEQUI
!          curref(iq)=PFCEQW(iq)
!        enddo

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine curfit(rk,zk,nk,NECON,WECON,psi_bnd)

         include 'double.inc'
         include 'parcur.inc'
          INCLUDE 'prm.inc'
         include 'comevl.inc'


c   -----------------------------
       real*8  rk(*),zk(*)
       DIMENSION   WECON(*)
       integer  necon(*),nk

           call psi_cpn

           call cursol(NEQUI,PFCEQW,psi_bnd)

            write(*,*) 'coil currents' 

            do ik=1,nequi
            write(*,*) PFCEQW(ik),ik 
            enddo

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine curfit_L(rk,zk,nk,NECON,WECON,psi_bnd)

         include 'double.inc'
         include 'parcur.inc'
          INCLUDE 'prm.inc'
         include 'comevl.inc'


c   -----------------------------
       real*8  rk(*),zk(*)
       DIMENSION   WECON(*)
       integer  necon(*),nk

           call psi_cpn

           call cursol_(NEQUI,PFCEQW,psi_bnd)

            write(*,*) 'coil currents' 

            do ik=1,nequi
            write(*,*) PFCEQW(ik),ik 
            enddo

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine curfit_(rk,zk,nk,NECON,WECON,psi_bnd)

         include 'double.inc'
         include 'dim.inc'
         include 'parcur.inc'
          INCLUDE 'prm.inc'
         include 'comevl.inc'

c   -----------------------------

       real*8  rk(*),zk(*)
       DIMENSION   WECON(*)
       integer  necon(*),nk

           call bonpsi

           call cursol_(NEQUI,PFCEQW,psi_bnd)

            write(*,*) 'coil currents' 

            do ik=1,nequi
            write(*,*) PFCEQW(ik),ik 
            enddo

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine grimat(rk,zk,ncpfc,NECON, WECON )

         include 'double.inc'
         include 'parcur.inc'
          INCLUDE 'prm.inc'
         include 'comevl.inc'

         real*8 rk(*),zk(*), WECON(*)
         integer nk
         integer NECON(*)

         pi = 3.14159265359d0 

         !ncpfc=nloc(npfc)

         do 1010 iq=1,NEQUI

            do 1876 l=1,Lfit
               Gindk(iq,l)=0.d0
 1876       continue

            do l=1,Lpre
               Gindp(iq,l)=0.d0
            enddo

               Gx_r(iq)=0.d0
               Gx_z(iq)=0.d0


             if(Lfit.ne.0) then
            do 110 l=1,Lfit

               do 112 ik=1,ncpfc
                 if( necon(ik) .eq. iq )  then
                  zgindk=greeni(rfit(l),zfit(l),rk(ik),zk(ik))/pi
                  Gindk(iq,l)=Gindk(iq,l)+zgindk*wecon(ik)
                 endif
 112           continue

 110        continue
             endif

             if(Lpre.ne.0) then
            do l=1,Lpre

               do ik=1,ncpfc
                 if( necon(ik) .eq. iq )  then
                  zgindp=greeni(rpre(l),zpre(l),rk(ik),zk(ik))/pi
                  Gindp(iq,l)=Gindp(iq,l)+zgindp*wecon(ik)
                 endif
               enddo

            enddo
             endif

               do ik=1,ncpfc
                 if( necon(ik) .eq. iq )  then

               call grGREN(rk(ik),zk(ik),Rx_p,Zx_p, dGdr,dGdz)
               Gx_r(iq)=Gx_r(iq)+dGdr*wecon(ik)/pi
               Gx_z(iq)=Gx_z(iq)+dGdz*wecon(ik)/pi

                 endif
               enddo

 1010    continue

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine precal(rk,zk,nk,NECON, WECON )

         include 'double.inc'
         include 'parcur.inc'
          INCLUDE 'prm.inc'
         include 'comevl.inc'
          parameter(nshp=10)
         include 'param_.inc'
         include 'comblc.inc'
          real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)

         real*8 rk(*),zk(*), WECON(*)
         integer nk
         integer NECON(*)

         !ncpfc=nloc(npfc)
         ncpfc=nk

         do 1010 iq=1,NEQUI

               Ginda(iq)=0.d0
               Gindx(iq)=0.d0

               do ik=1,ncpfc
                 if( necon(ik) .eq. iq )  then
                  zginda=greeni(r_ax,z_ax,rk(ik),zk(ik))/pi
                  Ginda(iq)=Ginda(iq)+zginda*wecon(ik)
                  zgindx=greeni(rx_p,zx_p,rk(ik),zk(ik))/pi
                  Gindx(iq)=Gindx(iq)+zgindx*wecon(ik)
                 endif
               enddo

 1010    continue

         r0=r_ax
         z0=z_ax

         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1

c---definition of the nearest knote

         sdmin=rmax

        do 320 k=0,1
         rr=r(ic+k)
        do 325 l=0,1
         zz=z(jc+l)
         dlx=dsqrt( (rr-r_ax)**2+(zz-z_ax)**2 )
       if(dlx.lt.sdmin) then
        sdmin=dlx
         ik=ic+k
         jk=jc+l
       endif
 325    continue
 320    continue

          nsh=1
          xs(nsh)=r(ik)
          ys(nsh)=z(jk)
         fun(nsh)=ui(ik,jk)

        do 400 k=-1,1

          i= ik+k

        do 410 l=-1,1

          j= jk+l

          if(k.eq.0  .AND. l.eq.0 ) go to 410
          nsh=nsh+1
          xs(nsh)=r(i)
          ys(nsh)=z(j)
         fun(nsh)=ui(i,j)

 410    continue
 400    continue

         rra=xs(1)
         zza=ys(1)

          call deriv5(xs,ys,fun,nsh,5,dp)

       psip_a=fun(1)+ dp(1)*(r_ax-rra) + dp(2)*(z_ax-zza)
     +          + 0.5d0*dp(3)*(r_ax-rra)*(r_ax-rra)
     +          +       dp(4)*(r_ax-rra)*(z_ax-zza)
     +          + 0.5d0*dp(5)*(z_ax-zza)*(z_ax-zza)

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine cursol(NEQUI,PFCEQW,psi_bnd)

         include 'double.inc'
         include 'parcur.inc'

        DIMENSION A(ncf_p,ncf_p),X(ncf_p),Y(ncf_p),IP(ncf_p)
        DIMENSION psictr(ncf_p)
        DIMENSION PFCEQW(*)

!      a(j,k)*I(k)+a(j,NEQUI+l)*Lam(l)+a(j,NEQUI+Lpre+1)*Psi_b=y(j)

         pi=3.14159265359d0
         amu0=0.4d0*pi
         !amu0=1.0d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(I_j)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          a(j,k)*I(k) ,j-number of equation

         do k=1,NEQUI 
         do j=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+Gindk(k,l)*Gindk(j,l)*w_wght(l)
          enddo
          a(j,k)=asum
        enddo
        enddo

         do j=1,NEQUI 
          a(j,j)=a(j,j)+d_wght(j)*s_wght
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation

         do j=1,NEQUI 
         do l=1,Lpre 
          a(j,NEQUI+l)=Gindp(j,l)
        enddo
        enddo

!          a(j,NEQUI+Lpre+1)*Lam_r(l) ,j-number of equation

         do j=1,NEQUI 
          a(j,NEQUI+Lpre+1)=Gx_r(j)
        enddo

!          a(j,NEQUI+Lpre+2)*Lam_z(l) ,j-number of equation

         do j=1,NEQUI 
          a(j,NEQUI+Lpre+2)=Gx_z(j)
        enddo




!          a(j,NEQUI+Lpre+3)*Lamx ,j-number of equation

         do j=1,NEQUI 
          a(j,NEQUI+Lpre+3)=(1.d0-alph)*Ginda(j)+alph*Gindx(j)
        enddo




!          a(j,NEQUI+Lpre+4)*Psi_b  ,j-number of equation

         do j=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*Gindk(j,l)
          enddo
          a(j,NEQUI+Lpre+4)=-asum
        enddo

!        right hand side

         do j=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*Gindk(j,l)*psifit(l)
          enddo
          y(j)=-asum+d_wght(j)*s_wght*curref(j)*amu0
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(Lam_l)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          a(j,k)*I(k) ,j-number of equation

         do l=1,Lpre 
          j=NEQUI+l
         do k=1,NEQUI 
          a(j,k)=Gindp(k,l)
        enddo
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation
           
        do l=1,Lpre 
          j=NEQUI+l
        do ll=1,Lpre 
          k=NEQUI+ll
          a(j,k)=0.d0
        enddo
        enddo

!          a(j,NEQUI+Lpre+1)*Lam_r ,j-number of equation
!          a(j,NEQUI+Lpre+2)*Lam_z ,j-number of equation
!          a(j,NEQUI+Lpre+3)*Lamx ,j-number of equation

        do l=1,Lpre 
          j=NEQUI+l
          k=NEQUI+Lpre+1
          a(j,k)=0.d0
          k=NEQUI+Lpre+2
          a(j,k)=0.d0
          k=NEQUI+Lpre+3
          a(j,k)=0.d0
        enddo

!          a(j,NEQUI+Lpre+4)*Psi_b  ,j-number of equation
 
        do l=1,Lpre 
          j=NEQUI+l
          k=NEQUI+Lpre+4 
          a(j,k)=-1.d0
        enddo
!        right hand side
          do l=1,Lpre 
           j=NEQUI+l
           y(j)=-psipre(l)
          enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(Lam_r)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          a(j,k)*I(k) ,j-number of equation

          j=NEQUI+Lpre+1

         do k=1,NEQUI 
          a(j,k)=Gx_r(k)
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation
           
        do ll=1,Lpre 
          k=NEQUI+ll
          a(j,k)=0.d0
        enddo

!          a(j,NEQUI+Lpre+1)*Lam_r ,j-number of equation
!          a(j,NEQUI+Lpre+2)*Lam_z ,j-number of equation
!          a(j,NEQUI+Lpre+3)*Lamx ,j-number of equation

          k=NEQUI+Lpre+1
          a(j,k)=0.d0
          k=NEQUI+Lpre+2
          a(j,k)=0.d0
          k=NEQUI+Lpre+3
          a(j,k)=0.d0

!          a(j,NEQUI+Lpre+4)*Psi_b  ,j-number of equation
 
          k=NEQUI+Lpre+4 
          a(j,k)=0.d0
!        right hand side
           y(j)=-dpsx_r

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(Lam_z)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          a(j,k)*I(k) ,j-number of equation

          j=NEQUI+Lpre+2

         do k=1,NEQUI 
          a(j,k)=Gx_z(k)
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation
           
        do ll=1,Lpre 
          k=NEQUI+ll
          a(j,k)=0.d0
        enddo

!          a(j,NEQUI+Lpre+1)*Lam_r ,j-number of equation
!          a(j,NEQUI+Lpre+2)*Lam_z ,j-number of equation
!          a(j,NEQUI+Lpre+3)*Lamx ,j-number of equation

          k=NEQUI+Lpre+1
          a(j,k)=0.d0
          k=NEQUI+Lpre+2
          a(j,k)=0.d0
          k=NEQUI+Lpre+3
          a(j,k)=0.d0

!          a(j,NEQUI+Lpre+3)*Psi_b  ,j-number of equation
 
          k=NEQUI+Lpre+4 
          a(j,k)=0.d0
!        right hand side
           y(j)=-dpsx_z






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(Lamx)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          a(j,k)*I(k) ,j-number of equation

          j=NEQUI+Lpre+3

         do k=1,NEQUI 
          a(j,k)=(1.d0-alph)*Ginda(k)+alph*Gindx(k)
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation
           
        do ll=1,Lpre 
          k=NEQUI+ll
          a(j,k)=0.d0
        enddo

!          a(j,NEQUI+Lpre+1)*Lam_r ,j-number of equation
!          a(j,NEQUI+Lpre+2)*Lam_z ,j-number of equation
!          a(j,NEQUI+Lpre+3)*Lamx ,j-number of equation

          k=NEQUI+Lpre+1
          a(j,k)=0.d0
          k=NEQUI+Lpre+2
          a(j,k)=0.d0
          k=NEQUI+Lpre+3
          a(j,k)=0.d0

!          a(j,NEQUI+Lpre+4)*Psi_b  ,j-number of equation
 
          k=NEQUI+Lpre+4 
          a(j,k)=-1.d0

!        right hand side
           y(j)=-(1.d0-alph)*psip_a-alph*psip_x







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(-psi_b)=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!          a(j,k)*I(k) ,j-number of equation

          j=NEQUI+Lpre+4

         do k=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*Gindk(k,l)
          enddo
          a(j,k)=asum
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation

         do l=1,Lpre 
          k=NEQUI+l
          a(j,k)=1.d0
        enddo

!          a(j,NEQUI+Lpre+1)*Lam_r ,j-number of equation
          a(j,NEQUI+Lpre+1)=0.d0
!          a(j,NEQUI+Lpre+2)*Lam_z ,j-number of equation
          a(j,NEQUI+Lpre+2)=0.d0
!          a(j,NEQUI+Lpre+3)*Lamx ,j-number of equation
          a(j,NEQUI+Lpre+3)=1.d0

!          a(j,NEQUI+Lpre+4)*Psi_b  ,j-number of equation

           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)
          enddo
          a(j,NEQUI+Lpre+4)=-asum
!!!
!        right hand side
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*psifit(l)
          enddo
          y(j)=-asum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call GE(NEQUI+Lpre+4,ncf_p,A,Y,X,IP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do iq=1,NEQUI
         PFCEQW(iq)=x(iq)
        enddo

!!!!!!checking!!!!!
        do l=1,Lfit
          psictr(l)=psifit(l)
        do iq=1,NEQUI
          psictr(l)=psictr(l)+Gindk(iq,l)*PFCEQW(iq)
        enddo
        enddo

          gr_xp=dpsx_r
          gz_xp=dpsx_z
          ps_xp=psip_x
          ps_ma=psip_a
        do iq=1,NEQUI
          gr_xp=gr_xp+Gx_r(iq)*PFCEQW(iq)
          gz_xp=gz_xp+Gx_z(iq)*PFCEQW(iq)

          ps_xp=ps_xp+Gindx(iq)*PFCEQW(iq)
          ps_ma=ps_ma+Ginda(iq)*PFCEQW(iq)
        enddo

        do l=1,Lpre
          psictr(Lfit+l)=psipre(l)
        do iq=1,NEQUI
          psictr(Lfit+l)=psictr(Lfit+l)+Gindp(iq,l)*PFCEQW(iq)
        enddo
        enddo

        do iq=1,NEQUI
         PFCEQW(iq)=x(iq)/amu0
        enddo

        psi_bnd=x(NEQUI+Lpre+4)

        alp_out=(ps_ma-psi_bnd)/(ps_ma-ps_xp)

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine psi_cpn

         include 'double.inc'
          parameter(nshp=10)
         include 'parcur.inc'
         include 'param.inc'
         include 'comblc.inc'
          real*8 xs(nshp),ys(nshp),fun(nshp),dp(5)

        do l=1,Lfit

         r0=rfit(l)
         z0=zfit(l)

         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1

         psifit(l)=blin(ic,jc,r0,z0)

        enddo 

        do l=1,Lpre

         r0=rpre(l)
         z0=zpre(l)

         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1

         psipre(l)=blin(ic,jc,r0,z0)

        enddo 

         r0=rx_p
         z0=zx_p

         ic=(r0-rmin)/dr(1)+1
         jc=(z0-zmin)/dz(1)+1


c---definition of the nearest knote

         sdmin=rmax

        do 320 k=0,1
         rr=r(ic+k)
        do 325 l=0,1
         zz=z(jc+l)
         dlx=dsqrt( (rr-rx_p)**2+(zz-zx_p)**2 )
       if(dlx.lt.sdmin) then
        sdmin=dlx
         ik=ic+k
         jk=jc+l
       endif
 325    continue
 320    continue

          nsh=1
          xs(nsh)=r(ik)
          ys(nsh)=z(jk)
         fun(nsh)=ui(ik,jk)

        do 400 k=-1,1

          i= ik+k

        do 410 l=-1,1

          j= jk+l

          if(k.eq.0  .AND. l.eq.0 ) go to 410
          nsh=nsh+1
          xs(nsh)=r(i)
          ys(nsh)=z(j)
         fun(nsh)=ui(i,j)

 410    continue
 400    continue

         rrx=xs(1)
         zzx=ys(1)

          call deriv5(xs,ys,fun,nsh,5,dp)

         dpsx_r = dp(1) + dp(3)*(rx_p-rrx) + dp(4)*(zx_p-zzx)
         dpsx_z = dp(2) + dp(5)*(zx_p-zzx) + dp(4)*(rx_p-rrx)

       psip_x=fun(1)+ dp(1)*(rx_p-rrx) + dp(2)*(zx_p-zzx)
     +          + 0.5d0*dp(3)*(rx_p-rrx)*(rx_p-rrx)
     +          +       dp(4)*(rx_p-rrx)*(zx_p-zzx)
     +          + 0.5d0*dp(5)*(zx_p-zzx)*(zx_p-zzx)

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine cursol_(NEQUI,PFCEQW,psi_bnd)

         include 'double.inc'
         include 'parcur.inc'

        DIMENSION A(ncf_p,ncf_p),X(ncf_p),Y(ncf_p),IP(ncf_p)
        DIMENSION psictr(ncf_p)
        DIMENSION PFCEQW(*)


!          a(j,k)*I(k)+a(j,NEQUI+l)*Lam(l)+a(j,NEQUI+Lpre+1)*Psi_b=y(j)

         pi=3.14159265359d0
         amu0=0.4d0*pi

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(I_j)=0
!
!          a(j,k)*I(k) ,j-number of equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do k=1,NEQUI 
         do j=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+Gindk(k,l)*Gindk(j,l)*w_wght(l)
          enddo
          a(j,k)=asum
        enddo
        enddo

         do j=1,NEQUI 
          a(j,j)=a(j,j)+d_wght(j)*s_wght
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation

         do j=1,NEQUI 
         do l=1,Lpre 
          a(j,NEQUI+l)=Gindp(j,l)
        enddo
        enddo

!          a(j,NEQUI+Lpre+1)*Psi_b  ,j-number of equation

         do j=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*Gindk(j,l)
          enddo
          a(j,NEQUI+Lpre+1)=-asum
        enddo

         do j=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*Gindk(j,l)*psifit(l)
          enddo
          y(j)=-asum+d_wght(j)*s_wght*curref(j)*amu0
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!        equation (.)*d(Lam_l)=0
!
!          a(j,k)*I(k) ,j-number of equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do l=1,Lpre 
          j=NEQUI+l
         do k=1,NEQUI 
          a(j,k)=Gindp(k,l)
        enddo
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation
           
        do l=1,Lpre 
          j=NEQUI+l
        do ll=1,Lpre 
          k=NEQUI+ll
          a(j,k)=0.d0
        enddo
        enddo

!          a(j,NEQUI+Lpre+1)*Psi_b  ,j-number of equation
 
        do l=1,Lpre 
          j=NEQUI+l
          k=NEQUI+Lpre+1 
          a(j,k)=-1.d0
        enddo

          do l=1,Lpre 
           j=NEQUI+l
           y(j)=-psipre(l)
          enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        equation (.)*d(-psi_b)=0
!
!          a(j,k)*I(k) ,j-number of equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         do k=1,NEQUI 
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*Gindk(k,l)
          enddo
          a(NEQUI+Lpre+1,k)=asum
        enddo

!          a(j,NEQUI+l)*Lam(l) ,j-number of equation

         do l=1,Lpre 
          k=NEQUI+l
          a(NEQUI+Lpre+1,k)=1.d0
        enddo

!          a(j,NEQUI+Lpre+1)*Psi_b  ,j-number of equation

           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)
          enddo
          a(NEQUI+Lpre+1,NEQUI+Lpre+1)=-asum
!!!
           asum=0.d0
          do l=1,Lfit 
           asum=asum+w_wght(l)*psifit(l)
          enddo
          y(NEQUI+Lpre+1)=-asum

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        call GE(NEQUI+Lpre+1,ncf_p,A,Y,X,IP)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do iq=1,NEQUI
         PFCEQW(iq)=x(iq)
        enddo


        do l=1,Lfit
          psictr(l)=psifit(l)
        do iq=1,NEQUI
          psictr(l)=psictr(l)+Gindk(iq,l)*PFCEQW(iq)
        enddo
        enddo

        do l=1,Lpre
          psictr(Lfit+l)=psipre(l)
        do iq=1,NEQUI
          psictr(Lfit+l)=psictr(Lfit+l)+Gindp(iq,l)*PFCEQW(iq)
        enddo
        enddo

        do iq=1,NEQUI
         PFCEQW(iq)=x(iq)/amu0
        enddo

        psi_bnd=x(NEQUI+Lpre+1)

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine bonpsi

         include 'double.inc'
         include 'dim.inc'
         include 'parcur.inc'
         include 'compol.inc'

         common/comaaa/ a12(nrp,ntp),a23(nrp,ntp),a34(nrp,ntp),
     +                  a14(nrp,ntp),a13(nrp,ntp),a24(nrp,ntp)

         dimension binadg(ntp,ntp),dgdn(ntp)

          i=iplas        !!!!!!!!!!!!!!!!!!!!!!1

         do 10 j=2,nt1

       a1=a13(i-1,j-1)
       a2=a34(i-1,j-1)+a12(i-1,j)
       a3=a24(i-1,j)

       g1=psi(i-1,j-1)
       g2=psi(i-1,j)
       g3=psi(i-1,j+1)

        dltk=(dlt(i,j-1)+dlt(i,j))*0.5d0
        dgdnl=a1*g1+a2*g2+a3*g3
        dgdn(j)=dgdnl/dltk

 10      continue

        dgdn(1)=dgdn(nt1)
        dgdn(nt)=dgdn(2)


        do l=1,Lfit

        rr=rfit(l)
        zz=zfit(l)

        do jb=1,nt1

        r0=r(iplas,jb)
        z0=z(iplas,jb)

        r1=r(iplas,jb+1)
        z1=z(iplas,jb+1)

        call bint(rr,zz,R0,Z0,r1,z1,Fint,1)

        binadg(jb,l)=fint

        enddo 
        enddo 

       do l=1,Lfit

        psb=0.d0

       do jb=2,nt1

        psb=psb+binadg(jb,l)*(dgdn(jb)+dgdn(jb+1))*0.5d0

       enddo

        psifit(l)=-psb

       enddo

        do l=1,Lpre

        rr=rpre(l)
        zz=zpre(l)

        do jb=1,nt1

        r0=r(iplas,jb)
        z0=z(iplas,jb)

        r1=r(iplas,jb+1)
        z1=z(iplas,jb+1)

        call bint(rr,zz,R0,Z0,r1,z1,Fint,1)

        binadg(jb,l)=fint

        enddo 
        enddo 

       do l=1,Lpre

        psb=0.d0

       do jb=2,nt1

        psb=psb+binadg(jb,l)*(dgdn(jb)+dgdn(jb+1))*0.5d0

       enddo

        psipre(l)=-psb

       enddo

        return
        end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

