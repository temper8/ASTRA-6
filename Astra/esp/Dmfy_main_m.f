C---------------------------------------------------------------
C   MAIN  PROGRAM  OF  THE EVOLUTION CODE  "PET"
C---------------------------------------------------------------

       IMPLICIT REAL*8( A-H, O-Z )
        dimension  contvals_mat(2500),voltpf(500),d_pf_mat(500)
        dimension  d_cam_mat(500)
        dimension  Rcp(10),Zcp(10)
         character*40 prename
c  kpr=1 for debugging, kpr=0 no printing

	    kpr=1
	    prename=''
          kname=1
      call aspid_flag(0)
      call  kpr_calc(kpr)
      call  put_name(prename,kname)
  
       key_dmf=-2 

       nstop = 300

       nstep=0
       time=0.d0
       dt=1.0d-2
       tau_con=2.0d-2
       k_con=(tau_con+1.d-8)/dt
       !k_con=100000
       !k_con=1
              call put_tim(dt,time)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       KLUCH = 0
       k_fixfree = 1

!----------------only fix boundary adaptive grid case
!
!
         if(k_fixfree.eq.0) then    
       k_auto= 1
          call B_STEPON( KLUCH, k_auto, nstep, dt, time,
     *                   rax,zax ,key_dmf)
          call cur_avg
          call  wrb
          call field_c
        write(*,*) 'nstep,time',nstep,time
        write(*,*) '**'
       KLUCH = 1
       k_auto= 0
 
       do nstep = 1,nstop
               time=time+dt
              call put_tim(dt,time)
      
          call B_STEPON( KLUCH, k_auto, nstep, dt, time,
     *                   rax,zax ,key_dmf)
          call cur_avg
          call field_c
          call  wrb
        write(*,*) '**'
        write(*,*) 'nstep,time',nstep,time
        write(*,*) '**'

       enddo
       stop
        endif
!        
!        
!---------------------only fix boundary adaptive grid case
       k_auto= 1
c----------------
!!!       k_grid= 0   ! rect. grid
          k_grid= 1   ! adap. grid

          !call tab_build
!
!!!!! basic free bound rectan, equilibrium ( KLUCH=0 )
!
        key_ini=0

       if(key_ini .eq. 0) then
          call  sstepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
        else
          call  cf_init( k_auto, nstep, dt, time,
     *                     voltpf, d_pf_mat,d_tcam_mat )
        endif
          !call  wrfb
          call  wrrec

       !stop 

          !call get_bouL_bra(contvals_mat,g1r,g1z,g2r,g2z,g5r,g5z)

       k_auto= 1

        if(k_grid.eq.1) then
!
!!!!! basic free bound, adaptive equilibrium ( KLUCH=0 )
!

          call  f_stepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
          call cur_avg
          call f_wrd
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         pause 'pause:initial equilibrium '
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       KLUCH = 1

       do nstep = 1,nstop
         

         nstep1=nstep-1
      if(nstep1/k_con*k_con .eq. nstep1) then
         do i=1,19
	 voltpf(i)=0.d0
         enddo
	 call cntrlr_iter(k_grid,voltpf)
      endif

      if(nstep1 .eq. 240) then
	 pause 'pause '
      endif

              time=time+dt
              call put_tim(dt,time)
	print *,' Next step: stepon',nstep

        if(k_grid.eq.0) then

          call  sstepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,key_dmf)
	print *,'  stepon done'
         call wrd
	print *,'  wrd done'
        elseif(k_grid.eq.1) then

          call  f_stepon( KLUCH, k_auto,nstep,dt,time,
     *                  voltpf, d_pf_mat ,d_cam_mat,rax,zax,key_dmf)
         call cur_avg
         call f_wrd
        endif

	print *,' After stepon'


          !call get_bouL_bra(contvals_mat,g1r,g1z,g2r,g2z,g5r,g5z)

        if(nstep/5*5.eq.nstep) then 
          call wrd_tim
          !pause 'pause'
        endif
        write(*,*) '**'
        write(*,*) 'nstep,time',nstep,time
        write(*,*) '**'
       enddo

       stop
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



