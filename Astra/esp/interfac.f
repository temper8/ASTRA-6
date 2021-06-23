        subroutine astra2spider(neql,nteta,nbnd,rzbnd,key_dmf,
     *                    na1,eqpf,eqff,fp,ipl,
     *                    rtor,btor,rho,roc,nstep,yreler,mu,
     *                    cc,Te,cubs,cd,key_ini,eqdfn )
                               
          include 'double.inc'
	    include 'dimpl1.inc'
          parameter(nbtabp=1000)
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
          common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab
          common/com_pas/ pstab(nursp),pptab(nursp),fptab(nursp),nutab
	    common /creler/ relerr

          real*8 rout(nrp,ntp),zout(nrp,ntp)

          real   t_start, t_finish

        real*8 rzbnd(*),eqpf(*),eqff(*),fp(*),ipl,ybetpl,yli3,
     *         rtor,btor,rho(*),roc,yreler,mu(*),
     *         cc(*),Te(*),cubs(*),cd(*) 
                              
          character*40 eqdfn

          pi=3.14159265359d0
          amu0=0.4d0*pi

		relerr	=yreler
 				!write(*,*) 'b_eqb, rz', (rzbnd(j),j=1,nbnd)
 				!write(*,*) 'nstep', nstep
             if(nbnd.gt.nbtabp) then
              write(*,*) 'spider:nbnd exeeds maximum value',nbtabp		   	
              write(*,*) 'nbnd=',nbnd		   				
              write(*,*) 'program is interrupted'
              stop					   				
             endif				
				
             if(na1.gt.nursp) then
              write(*,*) 'spider:na1 exeeds maximum value',nursp		   			
              write(*,*) 'na1=',na1		   				
              write(*,*) 'program is interrupted'
              stop					   				
             endif
		   
             if( neql.gt.nrp .or. nteta.gt.ntp) then
              write(*,*) 'spider: neql or nteta exeeds maximum value'		
              write(*,*) 'neql=',neql,'nrp=',nrp		   				
              write(*,*) 'nteta=',nteta,'ntp=',ntp		   				
              write(*,*) 'program is interrupted'
              stop					   				
             endif
				          
             nurs=-399
             i_betp=0
             i_bsh=1
             igdf=2
             i_eqdsk=0
             epsro=1.d-6
             maxit=250
            ! nstep=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
             n_tht=nteta
             n_psi=neql
             platok=ipl
             tokf=ipl

           call put_Ipl(tokf)

        if(nstep.eq.0 .AnD. key_ini.eq.0) then
               call tab_efit(tokf,psax,eqdfn,rax,zax,b0,r0)
             keyctr=0 
        write(fname,'(a,a)') path(1:kname),'durs_d.dat'
        open(1,file=fname,form='formatted')
          !open(1,file='durs_d.dat')
              write(1,*) n_tht,n_psi,igdf,nurs,keyctr,i_eqdsk,i_betp
              write(1,*) epsro,betplx,tokf,psax,b0,r0,rax,zax
              write(1,*) alf0,alf1,alf2,bet0,bet1,bet2
          close(1)
         return
        endif

!!!!!!!!!!!!+++++temporary++++++++++!!!!!!!!!!!!!!!!!!!!!!
!
!         open(1,file='tab_bnd.dat')
!	     read(1,*) nbtab
!	     do ib=1,nbtab
!	        read(1,*) rbtab(ib),zbtab(ib)
!	     enddo
!         close(1)
!
!         nbnd=nbtab
!
!	     do ib=1,nbtab
!	        rzbnd(ib)=rbtab(ib)
!              rzbnd(ib+nbnd)=zbtab(ib)
!	     enddo
!
!
!      open(1,file='tabppf.dat')
!
!              read(1,*) nutab
!              do i=1,nutab
!                 read(1,*) pstab(i), pptab(i), fptab(i)
!              enddo
!
!	        do i=1,nutab
!                 pstab(i)=pstab(i)/pstab(nutab)
!              enddo
!
!      close(1)
!
!!!!!!!!!!!!+++++++++++++++++++++++++++++++!!!!!!!!!!!!!!!!!!!!!!

         nbtab=nbnd

	     do ib=1,nbtab

	        rbtab(ib)=rzbnd(ib)
              zbtab(ib)=rzbnd(ib+nbnd)

	     enddo
!__CH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        write(fname,'(a,a)') path(1:kname),'tab_bnd.dat'
        open(1,file=fname,form='formatted')
         !open(1,file='tab_bnd.dat')
	     write(1,*) nbtab
	     do ib=1,nbtab
	        write(1,*) rbtab(ib),zbtab(ib)
	     enddo
         close(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!__CH

             rbmax=rbtab(1)
             rbmin=rbtab(1)
             zbmax=zbtab(1)
             zbmin=zbtab(1)

	     do ib=1,nbtab

            if(rbtab(ib).ge.rbmax) rbmax=rbtab(ib)
            if(rbtab(ib).le.rbmin) rbmin=rbtab(ib)
            if(zbtab(ib).ge.zbmax) zbmax=zbtab(ib)
            if(zbtab(ib).le.zbmin) zbmin=zbtab(ib)

	     enddo

          rc0=0.5d0*(rbmax+rbmin)           
          zc0=0.5d0*(zbmax+zbmin)           

          rax=rc0
          zax=zc0
!__CH
          b0=btor
          r0=rtor
!__CH

!! p' and ff' profiles

              do i=1,na1
                eqpf(i)=eqpf(i)*amu0  
                eqff(i)=eqff(i)*amu0
              enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ps1=fp(1)           
             ps2=fp(2)           
             ps3=fp(3)           

             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)

             ffx1=eqff(1)
             ffx2=eqff(2)
             ffx3=eqff(3)

       call EXTRP2(rh0,ps0, rh1,rh2,rh3, ps1,ps2,ps3)
       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
       call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

           psb=fp(na1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	if(mu(1).gt.0.) then
		dps0=.5d0*(mu(1)*rh1*rh1)
		ps0=ps1-dps0
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                pptab(1)=ppx0
                fptab(1)=ffx0
                pstab(1)=0.d0


        nutab=na1+1

              do i=2,nutab
                pptab(i)=eqpf(i-1)
                fptab(i)=eqff(i-1)
              enddo
               pstab(2)=(ps1-ps0)/(psb-ps0)
               pstab(3)=(ps2-ps0)/(psb-ps0)
               pstab(4)=(ps3-ps0)/(psb-ps0)

              do i=2,nutab
        if(i.gt.4)        pstab(i)=(fp(i-1)-ps0)/(psb-ps0)
	   if(pstab(i).le.pstab(i-1))   then
          write(*,*) 'psi is nonmonotonic!!! ', pstab(i),i
         endif
              enddo

!normalization        

              do i=1,nutab
                pptab(i)=pptab(i)/rtor
                fptab(i)=fptab(i)*rtor
              enddo
!__CH

        if(nstep.eq.0) then

        write(fname,'(a,a)') path(1:kname),'tabppf.dat'
        open(1,file=fname,form='formatted')
          !open(1,file='tabppf.dat')

              write(1,*) nutab

           do i=1,nutab
              write(1,*) pstab(i),pptab(i),fptab(i)
           enddo

          close(1)                                                 
        endif

!__CH

         !write(*,*)'nstep,key_dmf',nstep,key_dmf

          if(nstep.ne.0) then
        
            call retab_L_ast(pstab,pptab,fptab,nutab) !
        
           if(key_dmf.ne.0) then
            call profro_p(na1,neql,rtor,eqpf,rho,roc)
            call profro_sbc(na1,neql,rtor,btor,cc,Te,cubs,cd,rho,roc)
           else
            call profro(na1,neql,rtor,eqpf,eqff,rho,roc)
           endif
 
          endif
!__CH
              !do i=1,nutab
              !  pptab(i)=eqpf(i)
              !  fptab(i)=eqff(i)
              !  pstab(i)=fp(i)
              !enddo

          !open(17,file='out.pr')



 1149     format(a40)

        if(nstep.eq.0) then
             keyctr=0 
        write(fname,'(a,a)') path(1:kname),'durs_d.dat'
        open(1,file=fname,form='formatted')
          !open(1,file='durs_d.dat')
              write(1,*) n_tht,n_psi,igdf,nurs,keyctr,i_eqdsk,i_betp
              write(1,*) epsro,betplx,tokf,psax,b0,r0,rax,zax
              write(1,*) alf0,alf1,alf2,bet0,bet1,bet2
          close(1)
        endif

         return
         end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subroutine spider2astra(rout,zout,rtor,btor,rho,roc,na1,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *				    yreler,yli3,ni_p,nj_p,platok,cu,fp,pres,W_Dj)
                               
          include 'double.inc'
	    include 'dimpl1.inc'
          parameter(nbtabp=1000)
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
          common /com_bas/ rbtab(nbtabp),zbtab(nbtabp),nbtab
          common/com_pas/ pstab(nursp),pptab(nursp),fptab(nursp),nutab
      	common /creler/ relerr
	integer na1,ni_p,nj_p
          real*8 rout(ni_p,nj_p),zout(ni_p,nj_p)

          real   t_start, t_finish

        real*8 
     *       g11(*),g22(*),g33(*),                          
     *       vr(*),vrs(*),slat(*),gradro(*),rocnew,ybetpl,betpol,platok,
     *       mu(*),ipol(*),bmaxt(*),bmint(*),bdb02(*),bdb0(*),
     *	     b0db2(*),droda(*),rtor,btor,rho(*),roc,cu(*),fp(*),pres(*),
     *	     W_Dj(*)
                              
          character*40 eqdfn


         call get_rz(rout,zout,ni_p,nj_p)

cw	write(*,*) 'after eqb'

             ybetpl=betpol

           call comet(na1,platok,
     *                    rtor,btor,rho,roc,nstep,
     *                    g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                    mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *                    ybetpl,yli3,cu,fp,pres,W_Dj)



       return
       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!last change: 27.05.03
 
           subroutine comet(na1,platok,
     *                      rtor,btor,rho,roc,nstep,
     *                      g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                      mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *					  ybetpl,yli3,cu,fp,pres,W_Dj)

         include 'double.inc'
         include 'dim.inc'
         parameter(nrpl=nrp+1) 
          parameter(nrpl4=nrpl+4,nrpl6=nrpl4*6)
         include 'compol.inc'
       common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh
        common /com_jb/ BJ_av(nrp),curfi_av(nrp)
        common/com_heat_Dj/ WDj(nrp)

          real*8 rhos(nrp),vols(nrpl),g11s(nrpl),g22s(nrpl),g33s(nrpl)
          real*8 gradrs(nrpl),amus(nrpl),fnors(nrpl)
          real*8 dvdro(nrpl),sa(nrp),rhocs(nrpl),rhocn(nrpl),rhosn(nrpl)
          real*8 b_maxt(nrpl),b_mint(nrpl),b_db02(nrpl),b_db0(nrpl),
     *		     b_0db2(nrpl),d_roda(nrpl)
          real*8 psi_1(nrp)
          real*8 prs_1(nrp)
     	
          real*8 RRK(nrpl4),CCK(nrpl4),WRK(nrpl6)
          real*8 CWK(4)
          real*8 rhowr(1000)         
        real*8 rtor,btor,rho(*),roc,ybetpl,yli3,
     *       g11(*),g22(*),g33(*),                          
     *       vr(*),vrs(*),slat(*),gradro(*),rocnew,droda(*), 
     *       mu(*),ipol(*),bmaxt(*),bmint(*),bdb02(*),b0db2(*),bdb0(*),
     *       cu(*),fp(*),pres(*),W_Dj(*)
          
          !common /fp_sav/ fp_0(1000)
          !common /fp_dot/ dfpdt(1000),nna1
                          
           sqrt(xx)=dsqrt(xx)
       
         !do i=1,na1  
          !fp_0(i)=fp(i)
         !enddo
         !nna1=na1
       
       
         na=na1-1
         platok=tokp       
         do i=1,iplas
          !rhos(i)=sqrt(flucf(i)/(pi*btor))
          rhos(i)=sqrt(flx_fi(i)/(pi*btor))
         enddo
         
         call get_psibon(psi_bn1)
         
         do i=1,iplas
          psi_1(i)=(psi(i,2)+psi_bn1)*2.d0*pi
         enddo
         
         do i=1,iplas
          psn=psia(i)
          prs_1(i)=funppp(psn)*(psim-psip)*cnor/amu0
         enddo
CP
cp
         rhocs(1)=0.d0
          rhos(1)=0.d0

         do i=2,iplas
          rhocs(i)=0.5d0*(rhos(i-1)+rhos(i))
         enddo

          rocnew=rhos(iplas)
          rhocs(iplas+1)=rhos(iplas)

c=============================andrey=02.03.2003vvvvvvvvvvvvv
          ARC0=rhocs(iplas+1)
          ARC1=rhocs(iplas)
          ARC2=rhocs(iplas-1)
          ARC3=rhocs(iplas-2)
          SPW1=(ARC0-ARC2)*(ARC0-ARC3)
     &        /(ARC1-ARC2)/(ARC1-ARC3)
          SPW2=(ARC0-ARC1)*(ARC0-ARC3)
     &        /(ARC2-ARC1)/(ARC2-ARC3)
          SPW3=(ARC0-ARC1)*(ARC0-ARC2)
     &        /(ARC3-ARC1)/(ARC3-ARC2)
          q(iplas)=SPW1*q(iplas-1)+SPW2*q(iplas-2)+SPW3*q(iplas-3)

c=============================andrey=02.03.2003^^^^^^^^^^^^^^
 


         do i=1,iplas+1
          rhocn(i)=rhocs(i)/rhos(iplas)
         enddo
          rhocn(1)=0.d0

         do i=1,iplas
          rhosn(i)=rhos(i)/rhos(iplas)
	    droda(i)=rhosn(i)
         enddo
          rhosn(1)=0.d0
cp

          rhocn(iplas+1)=1.d0
          rhosn(iplas)=1.d0

         do i=1,na
          rhowr(i)=rho(i)/roc
         enddo
          rhowr(na1)=1.d0-1.d-9

          sa(1)=0.d0
cw	write(*,*) (rhocn(j),j=1,iplas)
cw	write(*,*) (rhosn(j),j=1,iplas)

         do i=1,iplas-1

           dsqi=0.d0
           dvoli=0.d0
           av_r2=0.d0
           av_gr1=0.d0
           av_gr2=0.d0
           av_grr2=0.d0

          u1=rhos(i)
          u2=rhos(i+1)
          u3=rhos(i+1)
          u4=rhos(i)

          do j=2,nt1

           dvoli=dvoli+vol(i,j)
           dsqi=dsqi+sr(i+1,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

          z1=z(i,j)
          z2=z(i+1,j)
          z3=z(i+1,j+1)
          z4=z(i,j+1)

           sss=s(i,j)
           
           dadr= 0.5d0*((u1+u2)*(z2-z1)+(u2+u3)*(z3-z2)+
     *                   (u3+u4)*(z4-z3)+(u4+u1)*(z1-z4))/sss

           dadz=-0.5d0*((u1+u2)*(r2-r1)+(u2+u3)*(r3-r2)+
     *                   (u3+u4)*(r4-r3)+(u4+u1)*(r1-r4))/sss

           gr2=dadr**2+dadz**2

           av_gr1=av_gr1+vol(i,j)*sqrt(gr2)
           av_gr2=av_gr2+vol(i,j)*gr2
           av_grr2=av_grr2+vol(i,j)*gr2/r0**2
           av_r2=av_r2+vol(i,j)/r0**2

          enddo
          sa(i+1)=dsqi*2.d0*pi
c          vols(i+1)=dvoli
          gradrs(i+1)=av_gr1/dvoli
          amus(i+1)=(2.d0*pi)/q(i)
          fnors(i+1)=f(i)/(rtor*btor)

          dvdro(i+1)=2.d0*pi*dvoli/(rhos(i+1)-rhos(i))
cc          g11s(i+1)=av_gr2/dvoli
cc          g22s(i+1)=av_grr2/(dvoli*4.d0*pi*pi)*rtor
          g11s(i+1)=av_gr2*dvdro(i+1)/dvoli
        g22s(i+1)=av_grr2*dvdro(i+1)/(dvoli*4.d0*pi*pi)*rtor/fnors(i+1)
          g33s(i+1)=av_r2*rtor**2/dvoli
         enddo

cccc      magnetic field

CW	write(*,*) 'before b'

           bp2av=0.d0
		   bp2=0.d0
cc          do i=2,iplas-1 
          do i=1,iplas-1 

           dvoli=0.d0
           bmx=0.d0
           bmn=1.d9*btor
           bavr=0.d0
           bavr2=0.d0
           bavrd2=0.d0

           do j=2,nt1

           dvoli=dvoli+vol(i,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

            bm_pol=psim*(psia(i+1)-psia(i))/st(i,j)
            bp_pol=psim*(psia(i+1)-psia(i))/st(i,j+1)

        bpol2v=bm_pol**2*( vol1(i,j)/sin1(i,j)+vol2(i,j)/sin2(i,j))+ 
     +         bp_pol**2*( vol3(i,j)/sin3(i,j)+vol4(i,j)/sin4(i,j)) 

            bpol2=bpol2v /vol(i,j)
            b_fi2=(f(i)/r0)**2

            b_tot2=bpol2+b_fi2
            b_tot=dsqrt(b_tot2)

            bmx=dmax1(b_tot,bmx)
            bmn=dmin1(b_tot,bmn)

           bavr2=bavr2+b_tot2*vol(i,j)
           bavr =bavr +b_tot *vol(i,j)
	   bavrd2=bavrd2+vol(i,j)/(b_tot2+1.d-8)
           bp2av=bp2av+bpol2v
		bp2=bp2+bpol2v
          enddo            
           b_maxt(i+1)=bmx
           b_mint(i+1)=bmn
           b_db02(i+1)=bavr2/(dvoli*btor**2)
           b_0db2(i+1)=(bavrd2*btor**2)/(dvoli)
           b_db0(i+1)=bavr/(dvoli*btor)


        enddo 
           b_maxt(1)=.5d0*(b_maxt(2)+b_mint(2))
           b_mint(1)=b_maxt(1)
           b_db02(1)=b_db02(2)
           b_0db2(1)=b_0db2(2)
           b_db0(1)=b_db0(2)


c            yli3  =4.d0*pi*bp2av/(rtor*tokp*tokp)
            yli3  =  4.d0*pi*bp2/(rtor*(.4d0*pi*tokp)**2)     
C*NEW
 
c	goto 9991 
	ip0=1
	ip1=2
	ip2=3
	ip3=4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x3=rhocs(ip3)
          x2=rhocs(ip2)
          x1=rhocs(ip1)
          x0=rhocs(ip0)
!!!!!!!!!!!!g11(1)
call	b_extrp(x0, X1,X2,X3, ip1,ip2,ip3,ip0,g11s,jerr)
!!!!!!!!!!!!g22(1)
call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g33s,jerr)
          
!!!!!!!!!!!!fnors(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,fnors,jerr)
 
!!!!!!!!!!!!gradrs(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,gradrs,jerr)

!!!!!!!!!!!!amus(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,amus,jerr)


!          x0=rhocs(1)
!          x1=rhocs(2)
!          x2=rhocs(3)
!          x3=rhocs(4)

!
 9991	continue

	ip0=iplas+1
	ip1=iplas-2
	ip2=iplas-1
	ip3=iplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x1=rhocs(ip3)
          x2=rhocs(ip2)
          x3=rhocs(ip1)
          x0=rhocs(ip0)


cw	write(*,*) 'passed 1'
!!!!!!!!!!!gradrs(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,gradrs,jerr)
 
!!!!!!!!!!!amus(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,amus,jerr)

        
!!!!!!!!!!!dvdro(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,dvdro,jerr)

 
!!!!!!!!!!!!g11(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g11s,jerr)

!!!!!!!!!!!!g22(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(iplas+1)
	call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g33s,jerr)

!!!!!!!!!!!!b_maxt(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_maxt,jerr)

!!!!!!!!!!!!b_mint(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_mint,jerr)

!!!!!!!!!!!!b_db02(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db02,jerr)


!!!!!!!!!!!!b_0db2(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_0db2,jerr)

 
!!!!!!!!!!!!b_db0(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db0,jerr)


!!!!!!!!!!!!!!!!!!!!!splining to astra grids

	rhocn(1)=0.d0
	rhocn(iplas+1)=1.d0
	dvdro(1)=0.d0
	sa(1)=0.d0
	g22s(1)=0.d0
	g11s(1)=0.d0
        fnors(iplas+1)=fvac/(rtor*btor)

cp        n3spl=iplas
        n3spl=iplas+1

cw	write(*,*) (b_mint(i), i=1,iplas)
!!!!!!!!!!!!!!!!!!!!! bmint !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_mint,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmint(i)=CWk(1)
           enddo
               bmint(na1)=b_mint(iplas)

cw	write(*,*) 'bmint'

!!!!!!!!!!!!!!!!!!!!! bmaxt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_maxt,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmaxt(i)=CWk(1)
           enddo
               bmaxt(na1)=b_maxt(iplas)

cw	write(*,*) 'bmaxt'

!!!!!!!!!!!!!!!!!!!!! bdb02 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db02,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb02(i)=CWk(1)
           enddo
               bdb02(na1)=b_db02(iplas)

cw	write(*,*) 'bdb02'

!!!!!!!!!!!!!!!!!!!!! b0db2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_0db2,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               b0db2(i)=CWk(1)
           enddo
               b0db2(na1)=b_0db2(iplas)

cw	write(*,*) 'b0db2'

!!!!!!!!!!!!!!!!!!!!! bdb0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db0,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
cp               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb0(i)=CWk(1)
          enddo
              bdb0(na1)=b_db0(iplas)

cw	write(*,*) 'bdb0'

!!!!!!!!!!!!!!!!!!!!!! cu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,BJ_av,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               cu(i)=CWk(1)/btor/amu0
           enddo
               cu(na1)=BJ_av(iplas)/btor/amu0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! pres !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,prs_1,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               pres(i)=CWk(1)
           enddo
               pres(na1)=prs_1(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! fp !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,psi_1,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               fp(i)=-CWk(1)
          enddo
               fp(na1)=-psi_1(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!! W_Dj !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,WDj,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               W_Dj(i)=-CWk(1)
          enddo
               W_Dj(na1)=WDj(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          !do i=1,na1
          ! dfpdt(i)=(fp_0(i)-fp(i))/dtim
          !enddo

!!!!!!!!!!!!!!!!!!!!!! slat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,sa,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               slat(i)=CWk(1)
           enddo
               slat(na1)=sa(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         n3spl=iplas+1

!!!!!!!!!!!!!!!!!!!!!!! ipol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,fnors,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               ipol(i)=CWk(1)
           enddo
               ipol(na1)=fnors(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!! mu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,amus,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               mu(i)=CWk(1)
           enddo
               mu(na1)=amus(iplas+1)



!!!!!!!!!!!!!!!! dvdro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhocn,dvdro,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vr(i)=CWk(1)
           enddo

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vrs(i)=CWk(1)
           enddo
               vrs(na1)=dvdro(iplas)

          
!!!!!!!!!!!!!!!!!! g11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cw	write(*,*) 'g11s',(g11s(j),j=1,iplas)

        CALL E01BAF(n3spl,rhocn,g11s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g11(i)=CWk(1)
           enddo
		g11(na1)=g11s(iplas+1)



!!!!!!!!!!!!!!!!!!!! g22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g22s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
              g22(i)=CWk(1)

          enddo

		g22(na1)=g22s(iplas+1)

!!!!!!!!!!!!!!!!!!!!!! g33 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g33s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g33(i)=CWk(1)
           enddo


!!!!!!!!!!!!!!!!!!!!!!! gradro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        CALL E01BAF(n3spl,rhocn,gradrs,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               gradro(i)=CWk(1)
           enddo
               gradro(na1)=gradrs(iplas+1)

CW		write(*,*) 'end'


           return
           end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           subroutine comet_(na1,platok,
     *                      rtor,btor,rho,roc,nstep,
     *                      g11,g22,g33,vr,vrs,slat,gradro,rocnew,
     *                      mu,ipol,bmaxt,bmint,bdb02,b0db2,bdb0,droda,
     *					  ybetpl,yli3)

         include 'double.inc'
         include 'dim.inc'
         parameter(nrpl=nrp+1) 
          parameter(nrpl4=nrpl+4,nrpl6=nrpl4*6)
         include 'compol.inc'
       common /combsh/ rm0,zm0,rc0,zc0,asp0,el_up,el_lw,tr_up,tr_lw,nbsh

          real*8 rhos(nrp),vols(nrpl),g11s(nrpl),g22s(nrpl),g33s(nrpl)
          real*8 gradrs(nrpl),amus(nrpl),fnors(nrpl)
          real*8 dvdro(nrpl),sa(nrp),rhocs(nrpl),rhocn(nrpl),rhosn(nrpl)
          real*8 b_maxt(nrpl),b_mint(nrpl),b_db02(nrpl),b_db0(nrpl),
     *		   b_0db2(nrpl),d_roda(nrpl)	
          real*8 RRK(nrpl4),CCK(nrpl4),WRK(nrpl6)
          real*8 CWK(4)
          real*8 rhowr(1000)         
        real*8 rtor,btor,rho(*),roc,ybetpl,yli3,
     *       g11(*),g22(*),g33(*),                          
     *       vr(*),vrs(*),slat(*),gradro(*),rocnew,droda(*), 
     *       mu(*),ipol(*),bmaxt(*),bmint(*),bdb02(*),b0db2(*),bdb0(*)
                          
           sqrt(xx)=dsqrt(xx)

         na=na1-1
       
         platok=tokp       
         do i=1,iplas
          !rhos(i)=sqrt(flucf(i)/(pi*btor))
          rhos(i)=sqrt(flx_fi(i)/(pi*btor))
         enddo
CP
cp
         rhocs(1)=0.d0
          rhos(1)=0.d0

         do i=2,iplas
          rhocs(i)=0.5d0*(rhos(i-1)+rhos(i))
         enddo

          rocnew=rhos(iplas)
          rhocs(iplas+1)=rhos(iplas)

c=============================andrey=02.03.2003vvvvvvvvvvvvv
          ARC0=rhocs(iplas+1)
          ARC1=rhocs(iplas)
          ARC2=rhocs(iplas-1)
          ARC3=rhocs(iplas-2)
          SPW1=(ARC0-ARC2)*(ARC0-ARC3)
     &        /(ARC1-ARC2)/(ARC1-ARC3)
          SPW2=(ARC0-ARC1)*(ARC0-ARC3)
     &        /(ARC2-ARC1)/(ARC2-ARC3)
          SPW3=(ARC0-ARC1)*(ARC0-ARC2)
     &        /(ARC3-ARC1)/(ARC3-ARC2)
          q(iplas)=SPW1*q(iplas-1)+SPW2*q(iplas-2)+SPW3*q(iplas-3)

c=============================andrey=02.03.2003^^^^^^^^^^^^^^
 


         do i=1,iplas+1
          rhocn(i)=rhocs(i)/rhos(iplas)
         enddo
          rhocn(1)=0.d0

         do i=1,iplas
          rhosn(i)=rhos(i)/rhos(iplas)
	    droda(i)=rhosn(i)
         enddo
          rhosn(1)=0.d0
cp

          rhocn(iplas+1)=1.d0
          rhosn(iplas)=1.d0

         do i=1,na
          rhowr(i)=rho(i)/roc
         enddo
          rhowr(na1)=1.d0-1.d-9

          sa(1)=0.d0
cw	write(*,*) (rhocn(j),j=1,iplas)
cw	write(*,*) (rhosn(j),j=1,iplas)

         do i=1,iplas-1

           dsqi=0.d0
           dvoli=0.d0
           av_r2=0.d0
           av_gr1=0.d0
           av_gr2=0.d0
           av_grr2=0.d0

          u1=rhos(i)
          u2=rhos(i+1)
          u3=rhos(i+1)
          u4=rhos(i)

          do j=2,nt1

           dvoli=dvoli+vol(i,j)
           dsqi=dsqi+sr(i+1,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

          z1=z(i,j)
          z2=z(i+1,j)
          z3=z(i+1,j+1)
          z4=z(i,j+1)

           sss=s(i,j)
           
           dadr= 0.5d0*((u1+u2)*(z2-z1)+(u2+u3)*(z3-z2)+
     *                   (u3+u4)*(z4-z3)+(u4+u1)*(z1-z4))/sss

           dadz=-0.5d0*((u1+u2)*(r2-r1)+(u2+u3)*(r3-r2)+
     *                   (u3+u4)*(r4-r3)+(u4+u1)*(r1-r4))/sss

           gr2=dadr**2+dadz**2

           av_gr1=av_gr1+vol(i,j)*sqrt(gr2)
           av_gr2=av_gr2+vol(i,j)*gr2
           av_grr2=av_grr2+vol(i,j)*gr2/r0**2
           av_r2=av_r2+vol(i,j)/r0**2

          enddo
          sa(i+1)=dsqi*2.d0*pi
c          vols(i+1)=dvoli
          gradrs(i+1)=av_gr1/dvoli
          amus(i+1)=(2.d0*pi)/q(i)
          fnors(i+1)=f(i)/(rtor*btor)

          dvdro(i+1)=2.d0*pi*dvoli/(rhos(i+1)-rhos(i))
cc          g11s(i+1)=av_gr2/dvoli
cc          g22s(i+1)=av_grr2/(dvoli*4.d0*pi*pi)*rtor
          g11s(i+1)=av_gr2*dvdro(i+1)/dvoli
        g22s(i+1)=av_grr2*dvdro(i+1)/(dvoli*4.d0*pi*pi)*rtor/fnors(i+1)
          g33s(i+1)=av_r2*rtor**2/dvoli
         enddo

cccc      magnetic field

CW	write(*,*) 'before b'

           bp2av=0.d0
		bp2=0.d0
cc          do i=2,iplas-1 
          do i=1,iplas-1 

           dvoli=0.d0
           bmx=0.d0
           bmn=1.d9*btor
           bavr=0.d0
           bavr2=0.d0
           bavrd2=0.d0

           do j=2,nt1

           dvoli=dvoli+vol(i,j)

          r1=r(i,j)
          r2=r(i+1,j)
          r3=r(i+1,j+1)
          r4=r(i,j+1)
          r0=(r1+r2+r3+r4)*0.25d0

            bm_pol=psim*(psia(i+1)-psia(i))/st(i,j)
            bp_pol=psim*(psia(i+1)-psia(i))/st(i,j+1)

        bpol2v=bm_pol**2*( vol1(i,j)/sin1(i,j)+vol2(i,j)/sin2(i,j))+ 
     +         bp_pol**2*( vol3(i,j)/sin3(i,j)+vol4(i,j)/sin4(i,j)) 

            bpol2=bpol2v /vol(i,j)
            b_fi2=(f(i)/r0)**2

            b_tot2=bpol2+b_fi2
            b_tot=dsqrt(b_tot2)

            bmx=dmax1(b_tot,bmx)
            bmn=dmin1(b_tot,bmn)

           bavr2=bavr2+b_tot2*vol(i,j)
           bavr =bavr +b_tot *vol(i,j)
	   bavrd2=bavrd2+vol(i,j)/(b_tot2+1.d-8)
           bp2av=bp2av+bpol2v
		bp2=bp2+bpol2v
          enddo            
           b_maxt(i+1)=bmx
           b_mint(i+1)=bmn
           b_db02(i+1)=bavr2/(dvoli*btor**2)
           b_0db2(i+1)=(bavrd2*btor**2)/(dvoli)
           b_db0(i+1)=bavr/(dvoli*btor)


        enddo 
           b_maxt(1)=.5d0*(b_maxt(2)+b_mint(2))
           b_mint(1)=b_maxt(1)
           b_db02(1)=b_db02(2)
           b_0db2(1)=b_0db2(2)
           b_db0(1)=b_db0(2)


c            yli3  =4.d0*pi*bp2av/(rtor*tokp*tokp)
            yli3  =  4.d0*pi*bp2/(rtor*(.4d0*pi*tokp)**2)     
C*NEW
 
c	goto 9991 
	ip0=1
	ip1=2
	ip2=3
	ip3=4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x3=rhocs(ip3)
          x2=rhocs(ip2)
          x1=rhocs(ip1)
          x0=rhocs(ip0)
!!!!!!!!!!!!g11(1)
call	b_extrp(x0, X1,X2,X3, ip1,ip2,ip3,ip0,g11s,jerr)
!!!!!!!!!!!!g22(1)
call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,g33s,jerr)
          
!!!!!!!!!!!!fnors(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,fnors,jerr)
 
!!!!!!!!!!!!gradrs(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,gradrs,jerr)

!!!!!!!!!!!!amus(1)
	call	b_extrp(x0,X1,X2,X3, ip1,ip2,ip3,ip0,amus,jerr)


!          x0=rhocs(1)
!          x1=rhocs(2)
!          x2=rhocs(3)
!          x3=rhocs(4)

!
 9991	continue

	ip0=iplas+1
	ip1=iplas-2
	ip2=iplas-1
	ip3=iplas
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          x1=rhocs(ip3)
          x2=rhocs(ip2)
          x3=rhocs(ip1)
          x0=rhocs(ip0)


cw	write(*,*) 'passed 1'
!!!!!!!!!!!gradrs(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,gradrs,jerr)
 
!!!!!!!!!!!amus(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,amus,jerr)

        
!!!!!!!!!!!dvdro(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,dvdro,jerr)

 
!!!!!!!!!!!!g11(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g11s,jerr)

!!!!!!!!!!!!g22(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,g22s,jerr)
          
!!!!!!!!!!!!g33(iplas+1)
	call	b_extrp(x0,X1,X2,X3, ip3,ip2,ip1,ip0,g33s,jerr)

!!!!!!!!!!!!b_maxt(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_maxt,jerr)

!!!!!!!!!!!!b_mint(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_mint,jerr)

!!!!!!!!!!!!b_db02(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db02,jerr)


!!!!!!!!!!!!b_0db2(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_0db2,jerr)

 
!!!!!!!!!!!!b_db0(iplas+1)
	call	b_extrp(x0, X1,X2,X3, ip3,ip2,ip1,ip0,b_db0,jerr)


!!!!!!!!!!!!!!!!!!!!!splining to astra grids

	rhocn(1)=0.d0
	rhocn(iplas+1)=1.d0
	dvdro(1)=0.d0
	sa(1)=0.d0
	g22s(1)=0.d0
	g11s(1)=0.d0
        fnors(iplas+1)=fvac/(rtor*btor)

cp        n3spl=iplas
        n3spl=iplas+1

cw	write(*,*) (b_mint(i), i=1,iplas)
!!!!!!!!!!!!!!!!!!!!! bmint !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_mint,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmint(i)=CWk(1)
           enddo
               bmint(na1)=b_mint(iplas)

cw	write(*,*) 'bmint'

!!!!!!!!!!!!!!!!!!!!! bmaxt !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_maxt,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bmaxt(i)=CWk(1)
           enddo
               bmaxt(na1)=b_maxt(iplas)

cw	write(*,*) 'bmaxt'

!!!!!!!!!!!!!!!!!!!!! bdb02 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db02,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb02(i)=CWk(1)
           enddo
               bdb02(na1)=b_db02(iplas)

cw	write(*,*) 'bdb02'

!!!!!!!!!!!!!!!!!!!!! b0db2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_0db2,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               b0db2(i)=CWk(1)
           enddo
               b0db2(na1)=b_0db2(iplas)

cw	write(*,*) 'b0db2'

!!!!!!!!!!!!!!!!!!!!! bdb0 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,b_db0,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
cp               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               bdb0(i)=CWk(1)
          enddo
              bdb0(na1)=b_db0(iplas)

cw	write(*,*) 'bdb0'


!!!!!!!!!!!!!!!!!!!!!! slat !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        n3spl=iplas

        CALL E01BAF(n3spl,rhosn,sa,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               slat(i)=CWk(1)
           enddo
               slat(na1)=sa(iplas)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         n3spl=iplas+1

!!!!!!!!!!!!!!!!!!!!!!! ipol !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,fnors,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               ipol(i)=CWk(1)
           enddo
               ipol(na1)=fnors(iplas+1)
!!!!!!!!!!!!!!!!!!!!!!! mu !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,amus,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               mu(i)=CWk(1)
           enddo
               mu(na1)=amus(iplas+1)



!!!!!!!!!!!!!!!! dvdro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CALL E01BAF(n3spl,rhocn,dvdro,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vr(i)=CWk(1)
           enddo

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               vrs(i)=CWk(1)
           enddo
               vrs(na1)=dvdro(iplas)

          
!!!!!!!!!!!!!!!!!! g11 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cw	write(*,*) 'g11s',(g11s(j),j=1,iplas)

        CALL E01BAF(n3spl,rhocn,g11s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g11(i)=CWk(1)
           enddo
		g11(na1)=g11s(iplas+1)



!!!!!!!!!!!!!!!!!!!! g22 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g22s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
		zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
              g22(i)=CWk(1)

          enddo

		g22(na1)=g22s(iplas+1)

!!!!!!!!!!!!!!!!!!!!!! g33 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        CALL E01BAF(n3spl,rhocn,g33s,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na1
               zspl=rhowr(i)
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               g33(i)=CWk(1)
           enddo


!!!!!!!!!!!!!!!!!!!!!!! gradro !!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        CALL E01BAF(n3spl,rhocn,gradrs,RRK,CCK,
     *                          n3spl+4,WRK,6*n3spl+16,IFAIL)

          do i=1,na
               zspl=(rhowr(i)+rhowr(i+1))*0.5d0
              CALL E02BCF(n3spl+4,RRK,CCK,zspl,0,CWk,IFAIL)
               gradro(i)=CWk(1)
           enddo
               gradro(na1)=gradrs(iplas+1)

CW		write(*,*) 'end'


           return
           end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro(na1,neql,rtor,eqpf,eqff,rho,roc)
c	29 JAN 2003
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

       real*8 rosn(nrp),ppr(nursp),ffpr(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 eqpf(*),eqff(*),rho(*),rtor,roc

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)

             ffx1=eqff(1)
             ffx2=eqff(2)
             ffx3=eqff(3)

       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
       call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

           ppr(1)=ppx0
           ffpr(1)=ffx0
           rhow(1)=0.d0

          do i=2,na1+1
           ppr(i)=eqpf(i-1)
           ffpr(i)=eqff(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           dpdpsi(1)=ppr(1)
           dfdpsi(1)=ffpr(1)

           dpdpsi(iplas)=ppr(na1+1)
           dfdpsi(iplas)=ffpr(na1+1)


           CALL E01BAF(nspl,rhow,ppr,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dpdpsi(i)=cwk(1)

           enddo


           CALL E01BAF(nspl,rhow,ffpr,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dfdpsi(i)=cwk(1)

           enddo

              do i=1,iplas
                dpdpsi(i)=dpdpsi(i)/rtor
                dfdpsi(i)=dfdpsi(i)*rtor
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_p(na1,neql,rtor,eqpf,rho,roc)
c	march 2007
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'

       real*8 rosn(nrp),rhot(nursp),ppr(nursp),ffpr(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 eqpf(*),rho(*),rtor,roc

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ppx1=eqpf(1)
             ppx2=eqpf(2)
             ppx3=eqpf(3)

             !ffx1=eqff(1)
             !ffx2=eqff(2)
             !ffx3=eqff(3)

       call EXTRP2(rh0,ppx0, rh1,rh2,rh3, ppx1,ppx2,ppx3)
       !call EXTRP2(rh0,ffx0, rh1,rh2,rh3, ffx1,ffx2,ffx3)

           ppr(1)=ppx0
           !ffpr(1)=ffx0
           rhow(1)=0.d0

          do i=2,na1+1
           ppr(i)=eqpf(i-1)
           !ffpr(i)=eqff(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           dpdpsi(1)=ppr(1)
           !dfdpsi(1)=ffpr(1)

           dpdpsi(iplas)=ppr(na1+1)
           !dfdpsi(iplas)=ffpr(na1+1)


           CALL E01BAF(nspl,rhow,ppr,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              dpdpsi(i)=cwk(1)

           enddo


!           CALL E01BAF(nspl,rhow,ffpr,RRK,CCK,
!     *                        nspl+4,WRK,6*nspl+16,IFAIL)
!
!           do i=2,iplas-1
!
!              zrho=rosn(i)
!              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
!              dfdpsi(i)=cwk(1)
!
!           enddo

              do i=1,iplas
                dpdpsi(i)=dpdpsi(i)/rtor
                !dfdpsi(i)=dfdpsi(i)*rtor
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       subroutine profro_sbc(na1,neql,rtor,btor,cc,Te,cubs,cd,rho,roc)
c	03 JAN 2008
           include 'double.inc'
          parameter(nursp=1000,nursp4=nursp+4,nursp6=nursp4*6)
           include 'dim.inc'
           include 'compol.inc'
	    common /com_sigcd/ C_sig(nrp),T_el(nrp),C_bts(nrp),C_driv(nrp)

          real*8 rosn(nrp),funw(nursp),rhow(nursp)
          real*8 rrk(nursp4),cck(nursp4),wrk(nursp6)
          real*8 cwk(4)

        real*8 rho(*),rtor,btor,roc,cc(*),Te(*),cubs(*),cd(*)

        iplas=neql
        nspl=na1+1
         do i=1,iplas
          rosn(i)=(i-1.d0)/(iplas-1.d0)
         enddo

!! exstrapolation on magn. axis
       
             rh0=0.d0
		              
             rh1=rho(1)           
             rh2=rho(2)           
             rh3=rho(3)           

             ccx1=cc(1)
             ccx2=cc(2)
             ccx3=cc(3)

             cdx1=cd(1)
             cdx2=cd(2)
             cdx3=cd(3)

             cbx1=cubs(1)
             cbx2=cubs(2)
             cbx3=cubs(3)

             cex1=te(1)
             cex2=te(2)
             cex3=te(3)

       call EXTRP2(rh0,ccx0, rh1,rh2,rh3, ccx1,ccx2,ccx3)
       call EXTRP2(rh0,cdx0, rh1,rh2,rh3, cdx1,cdx2,cdx3)
       call EXTRP2(rh0,cbx0, rh1,rh2,rh3, cbx1,cbx2,cbx3)
       call EXTRP2(rh0,cex0, rh1,rh2,rh3, cex1,cex2,cex3)

!!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=ccx0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cc(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_sig(1)=funw(1)
           C_sig(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_sig(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cc!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!!!!!!!!!!!!!!!cd!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cdx0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cd(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_driv(1)=funw(1)
           C_driv(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_driv(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cd!!!!!!!!!!!!!!!!!!!!!!!
  
 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cbx0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=cubs(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           C_bts(1)=funw(1)
           C_bts(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              C_bts(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!cubs!!!!!!!!!!!!!!!!!!!!!!!
 
 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!

           funw(1)=cex0
           rhow(1)=0.d0

          do i=2,na1+1
           funw(i)=te(i-1)
           rhow(i)=rho(i-1)/roc
          enddo		

           T_el(1)=funw(1)
           T_el(iplas)=funw(na1+1)

           CALL E01BAF(nspl,rhow,funw,RRK,CCK,
     *                        nspl+4,WRK,6*nspl+16,IFAIL)

           do i=2,iplas-1

              zrho=rosn(i)
              CALL E02BCF(nspl+4,RRK,CCK,zrho,0,CWk,IFAIL)
              T_el(i)=cwk(1)

           enddo

 !!!!!!!!!!!!!!!!!!!te!!!!!!!!!!!!!!!!!!!!!!!
  
             do i=1,iplas
                C_bts(i)=C_bts(i)*btor
                C_driv(i)=C_driv(i)*btor
                T_el(i)=T_el(i)*1.d3
              enddo

            return
            end
		 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         subroutine get_rz(rout,zout,ni_p,nj_p)

           include 'double.inc'
           include 'dim.inc'
           include 'compol.inc'

          real*8 rout(ni_p,nj_p),zout(ni_p,nj_p)

        do i=1,iplas
        do j=1,nt
         rout(i,j)=r(i,j)
         zout(i,j)=z(i,j)
        enddo
        enddo

      return
      end


	subroutine	NEWGRD_ss( ROC, HRO, NB1, NA1, RHO,btor, HROA  )
C----------------------------------------------------------------------|
C Define size of the edge cell to be:   .6 <= HROA/HRO < 1.8
C	with a hysteresis of 0.2*HRO, so that
C		HROA<=0.6*HRO  jumps to  HROA<=1.6*HRO
C		HROA>=1.8*HRO  jumps to  HROA>=0.8*HRO
C	write(*,*)NA1,TIME,RHO(NA1)," -> ",ROC
C----------------------------------------------------------------------|
         include 'double.inc'
         include 'dim.inc'
         include 'compol.inc'
         dimension RHO(*)

           na=na1-1
          roc=sqrt(flx_fi(iplas)/(pi*btor))

C	write(*,*)NA1,ROC,HRO*(NA+0.1),HRO*(NA1+0.3),HRO*(NB1-0.49)
C	write(*,'(1P,6E13.5)')(RHO(j),j=NA-2,NA1+1)
C	write(*,'(1P,6E13.5)')(AMETR(j),j=NA-2,NA1+1)
	if (ROC .gt. HRO*(NB1-0.49))	then
C rho_edge gets out of the grid
	   write(*,*)'>>> WARNING: the allocated grid is too small'
	   write(*,*)'Time =',TIME,'   Rho_b =',ROC,
     >		'   Rho_max = ',HRO*(NB1-0.5)
	   write(*,*)'Try to increase AWALL'
	   write(*,*)'If this does not help inspect equilibrium input'
	   NA = NB1-1
	elseif (ROC .le. HRO*(NA+0.1))	then	! 0.1 <=> 0.6-0.5
!	    if (HROA < 0.6*HRO) then   reduce NA
C	   write(*,*)" <- ",ROC,HRO*NA,RHO(NA1),HRO*NA1,NA1
	   do	j=NA-1,1,-1
		NA = j
C	   	write(*,*)j,NA*HRO,ROC,(NA+1)*HRO,ROC/HRO-(NA-0.5)
		if (ROC .gt. HRO*(j+0.1))	goto	10
	   enddo
	elseif (ROC .gt. HRO*(NA1+0.3))	then	! 0.3 <=> 1.8-1.5
!	    if (HROA > 1.8*HRO) then   increase NA
C	   write(*,*)" -> ",ROC,HRO*NA,RHO(NA1),HRO*NA1,NA1
	   do	j=NA1,NB1
		NA = j
C	   	write(*,*)j,NA*HRO,ROC,(NA+1)*HRO
		if (ROC .le. HRO*(j+1.3))	goto	10
	   enddo
	endif
 10	continue

C HRO*(NA+0.6) < ROC <= HRO*(NA+1.8)
C	write(*,*)AB,ABC,TIME,TSTART
C	write(*,*)NA+1,NA+1-NA1,HRO*(NA+0.6),ROC,HRO*(NA+1.8)

	NA1 = NA+1
	do	J=1,NB1
	   RHO(J)=(J-0.5)*HRO
	enddo
	HROA = ROC-RHO(NA)
	RHO(NA1) = ROC
	
      return
	end
C======================================================================|

	subroutine	get_roc(ROC,btor)

c	implicit none
         include 'double.inc'
         include 'dim.inc'
         include 'compol.inc'

	double precision ROC,btor
          roc=sqrt(flx_fi(iplas)/(pi*btor))
 
	return
	end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine pla_volt(dt)
        include 'double.inc'
        include 'dim.inc'
        include 'compol.inc'
          common /com_volt/ upls(nrp)
          common /com_psisave/ psi_0(nrp)
 
          call get_psibon(psi_bn1)
          
        do i=1,iplas
         !upls(i)=(psia(i)*psim-psi_0(i))/dt
         upls(i)=(psi(i,2)+psi_bn1-psi_0(i))/dt
         enddo
          
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine savepsi
        include 'double.inc'
        include 'dim.inc'
        include 'compol.inc'
          common /com_psisave/ psi_0(nrp)
          
         call get_psibon(psi_bn1)
          
         do i=1,iplas
          !psi_0(i)=psim*psia(i)
          psi_0(i)=psi(i,2)+psi_bn1
         enddo
          
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         subroutine get_psibon(psi_bnd)
        include 'double.inc'
        include 'dim.inc'
      common/selcon/ psi_d(nrp),fi_d(nrp),f_d(nrp),ri_d(nrp),
     *               ps_pnt(nrp),del_psb,psi_bn1
          
          psi_bnd=psi_bn1
          
        return
        end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	