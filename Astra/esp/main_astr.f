          include 'double.inc'
          include 'dim.inc'
          parameter(ni_p=nrp,nj_p=ntp)
          parameter(np=1000,nb1=np+1) 
          parameter(nbtabp=1000,nbtabp2=nbtabp*2)
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)

        real*8 rbtab(nbtabp),zbtab(nbtabp),nbtab
        real*8 pstab(nursp),pptab(nursp),fptab(nursp),nutab
        real*8 Cc(np),CUbs(np),Te(np),Cd(np)
        real*8 Cu(np)

        real*8 rzbnd(nbtabp2),eqpf(np),eqff(np),fp(np),ipl,ybetpl,yli3,
     *       rtor,btor,rho(np),roc,
     *       g11(np),g22(np),g33(np),                          
     *       vr(np),vrs(np),slat(np),gradro(np),rocnew, 
     *       mu(np),ipol(np),bmaxt(np),bmint(np),bdb02(np),bdb0(np),
     *	   b0db2(np),droda(np)                          

        real*8 rout(ni_p,nj_p),zout(ni_p,nj_p)
              character*40 prename

           
! !!!!!!!!!!!!+++++test++++++++++!!!!!!!!!!!!!!!!!!!!!!
!
!          call blic_d(r0,z0, r1,r2,r3,r4, z1,z2,z3,z4,
!     *                             u1,u2,u3,u4, u0, dudr,dudz)
!
!
! !!!!!!!!!!!!+++++test++++++++++!!!!!!!!!!!!!!!!!!!!!!

 !!!!!!!!!!!!+++++temporary++++++++++!!!!!!!!!!!!!!!!!!!!!!

        ! open(1,file='spidat.dat.fail')
        ! open(1,file='spidat.dat-0')
         open(1,file='spidat.dat')
	     read(1,*) 
	     read(1,*) neql,nteta,nbnd
	     read(1,*) 
	     read(1,*) (rzbnd(i),i=1,nbnd*2)
	     read(1,*) 
	     read(1,*) na1
	     read(1,*) 
	     read(1,*) (eqpf(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (eqff(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (fp(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (rho(i),i=1,na1)
	     read(1,*) 
	     read(1,*) ipl,rtor,btor,roc,nstep
	     read(1,*) 
 	     read(1,*) (fp(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (g11(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (g22(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (g33(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (vr(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (vrs(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (slat(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (gradro(i),i=1,na1)
	     read(1,*) 
 	     read(1,*) (mu(i),i=1,na1)
	     read(1,*) 
	     read(1,*) (ipol(j),j=1,na1)

	read(1,*) 	 
	read(1,*) (cc(j),j=1,na1)
	read(1,*) 	 
	read(1,*) (Te(j),j=1,na1)
	read(1,*) 	 
	read(1,*) (cubs(j),j=1,na1)
	read(1,*) 	 
	read(1,*) (cd(j),j=1,na1)
        close(1)

           hro=rho(2)-rho(1)

!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!
        na=na1-1
        nutab=na1

        eqff(na1)=eqff(na)
                  !neql=98
                  !nteta=110
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!

	    kpr=1   !print in spider
          call  kpr_calc(kpr)
	    prename=''
          kname=4
      call  put_name(prename,kname)

          k_grid= 0   ! k_grid= 0   rect. grid
                      ! k_grid= 1   adap. grid

             k_auto= 1   ! k_auto= 1->   full initialization, 
                         ! k_auto= 0-> preinitialization is assumed to be done
             nstep=0     !nstep=0->initialization,#0->time stepping
             key_dmf=0   !=1->diff.mag.field, =0->without
             k_fixfree=1   !=0->only fixed boundary spider 

             dt=1.d-2
             time=0.d0

        call aspid_flag(1)


        call astra2spider(neql,nteta,nbnd,rzbnd,key_dmf,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,nstep,yreler,mu,
     *                    cc,Te,cubc,cd )


        call spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree)  

        call cur_avg
!                     stop
!!!!!!!!!!!!!!!!!!!!!!gap test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!         call get_psix(r_xp,z_xp,psi_xp)
!         rg=5.575d0             !8.2806d0
!         zg=-3.891d0            !0.4665d0
!         rt=rg
!         zt=zg
!         ut=psi_xp
!        
!         call gap(ut,rt,zt,rg,zg)
!
!        D_gap=dsqrt((rg-rt)**2+(zg-zt)**2)
!
!!!!!!!!!!!!!!!!!!!!!!gap test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        call spider2astra(rout,zout,rtor,btor,rho,roc,na1,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *				    yreler,yli3,ni_p,nj_p,platok,cu,fp)


        do nstep=1,300

              time=time+dt
              call put_tim(dt,time)

        call astra2spider(neql,nteta,nbnd,rzbnd,key_dmf,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,nstep,yreler,mu,
     *                    cc,Te,cubs,cd )


        call spider(nstep,time,dt,key_dmf,k_grid,k_auto,k_fixfree)  
         
         !write(*,*) 'before NEWGRD_s'
         !write(*,*) 'roc hro nb1 na1 btor hroa '
         !write(*,*) roc, hro, nb1, na1, btor, hroa 

	  call	NEWGRD_s( ROC, HRO, NB1, NA1, RHO,btor, HROA  )
         
         !write(*,*) 'after NEWGRD_s'
         !write(*,*) 'roc hro nb1 na1 btor hroa '
         !write(*,*) roc, hro, nb1, na1, btor, hroa 
         !write(*,*) 'rho'
         !write(*,*) (rho(i),i=1,na1)
 
        call spider2astra(rout,zout,rtor,btor,rho,roc,na1,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *				    yreler,yli3,ni_p,nj_p,platok,cu,fp)




        enddo


        stop
        end

