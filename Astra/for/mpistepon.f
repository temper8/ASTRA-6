      program   main
      implicit none
      logical   i_am_parent
      integer max_grid, ll_my_id, idum

! Init parallel mode use mpi ! 1st call in prog for common mpi prog
      print *, 'START'
      call init_ll("mpi")			! sbr/my_mpi.f90
!      call chdir("/afs/ipp-garching.mpg.de/home/g/grp/a6/")
      call get_i_am_parent(i_am_parent)
      if (.not.i_am_parent) then
        print *, 'child starts '
        call glf_trans_ll(1, 1)
      endif
      call astra
      call get_ll_my_id(ll_my_id)
      call ll_end(ll_my_id)
      end program    main
      subroutine astra
      include 'for/stepon.f'
