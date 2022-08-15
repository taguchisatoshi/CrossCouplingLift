!******************************************************************************************
    PROGRAM main
!******************************************************************************************
    USE constants; USE coordinates; USE molecular_velocity; USE parameters
    USE counter; USE VDF; USE global; USE specifier; USE discontinuity
    USE infinity; USE trajectory; USE IO
    IMPLICIT NONE
    INTEGER :: sz
    CHARACTER :: name*99
    REAL(8) :: Ia, Ib, Ic, force_lift
    real(8) :: kn_tmp, kn0
    integer :: InputStatus
    character :: file_u*128, file_s*128, char_tmp*999
    logical :: exists
    namelist /nlist/ kn0, n_unit

!******************************************************************************************
      !PRINT*, 'ID for U and S'
      !READ(*,*) ID_u
      !ID_u = 278; ID_s = 13 ! kn=0.5
      !ID_u = 273; ID_s = 9 ! kn=1
      !ID_u = 327; ID_s = 14  ! kn=2
      !ID_u = 280; ID_s = 16  ! kn=5
      !PRINT*, 'ID for S'
      !READ(*,*) ID_s
      !ID_s = 13
      !ID_s = 9
      !ID_s = 14


    open(unit=10, file='input.txt')
      read(10, nml=nlist)
    close(10)

      !Print*, 'kn ?'
      !read(*,*) kn0, n_unit
      print*, kn0
      
      !print*, 'n_unit ?'
      !read(*,*) n_unit
      print*, n_unit

      WRITE(ID_u_char,*) n_unit
      file_u = 'H:/work2/SPBGK/20220604_SPES20U1M/ID_list_summary_n'//TRIM(ADJUSTL(ID_u_char))//'.txt'
      inquire (file = file_u, exist = exists)
      if (.not. exists) then
        write (*,'(2a/)') ' >> cannot find file ', file_u
        stop
      end if
    
      !print*, file_u
      open(unit=10, file=file_u, status='old')
        read(10,*) char_tmp
        read(10,*) char_tmp
        do
          read(10,*,IOSTAT=InputStatus) ID_u, kn_tmp
          if (kn_tmp == kn0) then
            WRITE(ID_u_char,*) ID_u
            dir_u = 'H:/work2/SPBGK/20220604_SPES20U1M/ID'//TRIM(ADJUSTL(ID_u_char))
            print*, dir_u
            exit
          elseif (InputStatus < 0) then
            print*, 'No more data'
            stop
          endif
        enddo
      close(10)

      WRITE(ID_s_char,*) n_unit
      file_s = 'H:/work2/SPBGK/SPBGK20S1M/ID_list_summary_n'//TRIM(ADJUSTL(ID_s_char))//'.txt'
      inquire (file = file_s, exist = exists)
      if (.not. exists) then
        write (*,'(2a/)') ' >> cannot find file ', file_s
        stop
      end if

      !print*, trim(file_s)
      open(unit=10, file=file_s, status='old')
        read(10,*) char_tmp
        read(10,*) char_tmp
        do
          read(10,*,IOSTAT=InputStatus) ID_s, kn_tmp
          if (kn_tmp == kn0) then
            WRITE(ID_s_char,*) ID_s
            dir_s = 'H:/work2/SPBGK/SPBGK20S1M/ID'//TRIM(ADJUSTL(ID_s_char))
            print*, dir_s
            exit
          elseif (InputStatus < 0) then
            print*, 'No more data'
            stop
          endif
        enddo
      close(10)

      pi = 4.d0*ATAN(1.d0)

    !WRITE(ID_u_char,*) ID_u
    !dir_u = 'H:/work2/SPBGK/20220604_SPES20U1M/ID'//TRIM(ADJUSTL(ID_u_char))
    !WRITE(ID_s_char,*) ID_s
    !dir_s = 'H:/work2/SPBGK/SPBGK20S1M/ID'//TRIM(ADJUSTL(ID_s_char))
    
    name = TRIM(dir_u)//'/log_BIN.dat'
    name = trim(name)
    OPEN(UNIT=10, STATUS='OLD', FILE=name, FORM='UNFORMATTED')

      INQUIRE (FILE = name, SIZE = sz) ! get the file size
      READ(10) kn_u, number_unit_u, L_r_u, imax1_u, M_ref_u, M_z_u, M_z_th_u, &
              order_ipl_u, order_fd_u, step_out_u, IDfm_u, ID_dsc_u, ID_inf_u, conc_u, &
              n_crit_conv_u, step_add_conv_u, ID_dev_u, ID_dev_level_u, &
              coeff_under_inf_u, d1g_u, d2g_u, d3g_u, &
              ID_es_u, nu_u, ID_dev_noneq_u, accommo_u

    WRITE(*,*) 'Data read from ... ', trim(name)

    name = TRIM(dir_s)//'/log_s_BIN.dat'
    name = trim(name)
    OPEN(UNIT=10, STATUS='OLD', FILE=name, FORM='UNFORMATTED')

      INQUIRE (FILE = name, SIZE = sz) ! get the file size
      if (sz == 138) then
        READ(10) kn_s, number_unit_s, L_r_s, imax1_s, M_ref_s, M_z_s, M_z_th_s, &
                order_ipl_s, order_fd_s, step_out_s, IDfm_s, ID_dsc_s, ID_dsc_analyt_s, &
                ID_inf_s, conc_s, n_crit_conv_s, step_add_conv_s, ID_dev_s, ID_dev_level_s, &
                coeff_under_inf_s, d1g_s, d2g_s, d3g_s, &
                ID_es_s, nu_s, ID_dev_noneq_s, accommo_s, grid_conc_z_th_nw
      endif
      if (sz < 138) then
        READ(10) kn_s, number_unit_s, L_r_s, imax1_s, M_ref_s, M_z_s, M_z_th_s, &
                order_ipl_s, order_fd_s, step_out_s, IDfm_s, ID_dsc_s, ID_dsc_analyt_s, &
                ID_inf_s, conc_s, n_crit_conv_s, step_add_conv_s, ID_dev_s, ID_dev_level_s, &
                coeff_under_inf_s, d1g_s, d2g_s, d3g_s, &
                ID_es_s, nu_s, ID_dev_noneq_s, accommo_s
        grid_conc_z_th_nw = 0.d0
      endif
    CLOSE(10)

    WRITE(*,*) 'Data read from ... ', trim(name)

    CALL summary
    
! ---- Grid 
    ! region 1
    ALLOCATE( r1(0:imax1), dr1(0:imax1) )
! ---- Velocity space
    ALLOCATE( z(0:M_z),  dz(0:M_z) )
    ALLOCATE( z_th_ne(0:M_z_th), dz_th_ne(0:M_z_th) )
    ALLOCATE( z_th_nw(0:M_z_th), dz_th_nw(0:M_z_th) )

! ---- Allocation
! ---- U
    ALLOCATE( G1_ne_u(0:M_z,0:M_z_th,0:imax1) )
    ALLOCATE( G1_nw_u(0:M_z,0:M_z_th,0:imax1) )
    ALLOCATE( H1_ne_u(0:M_z,0:M_z_th,0:imax1) )
    ALLOCATE( H1_nw_u(0:M_z,0:M_z_th,0:imax1) )
    ALLOCATE( G1_dsc_u(0:M_z,0:1,0:imax1), H1_dsc_u(0:M_z,0:1,0:imax1) )
    ALLOCATE( dG1_dsc_u(0:M_z,0:imax1), dH1_dsc_u(0:M_z,0:imax1) )

    ALLOCATE( rho1_u(0:imax1,0:1), u1_r_u(0:imax1,0:1), u1_th_u(0:imax1,0:1), &
               t1_u(0:imax1,0:1), p1_u(0:imax1,0:1) )

    ALLOCATE( st1_rr_u(0:imax1,0:1), st1_rth_u(0:imax1,0:1), st1_thth_u(0:imax1,0:1), &
              st1_phph_u(0:imax1,0:1), q1_r_u(0:imax1,0:1), q1_th_u(0:imax1,0:1) )
    
    ALLOCATE( sigma1_u(0:1) )
       
    ALLOCATE( c1(0:imax1,0:1), c2(0:imax1,0:1), c3(0:imax1,0:1) )

    n_trj_max = imax1 + M_z_th-1
    ALLOCATE(r_trj(0:n_trj_max), z_th_trj(0:n_trj_max), s_trj(0:n_trj_max))
    ALLOCATE(r_map(0:imax1), z_th_map(1:M_z_th-1))
    ALLOCATE( G1_dsc_trj_u(0:M_z,0:1,0:n_trj_max), H1_dsc_trj_u(0:M_z,0:1,0:n_trj_max) )
    ALLOCATE( dG1_dsc_trj_u(0:M_z,0:n_trj_max), dH1_dsc_trj_u(0:M_z,0:n_trj_max) )

! ---- S
    ALLOCATE( G1_ne_s(0:M_z,0:M_z_th,0:imax1) )
    ALLOCATE( G1_nw_s(0:M_z,0:M_z_th,0:imax1) )
    ALLOCATE( G1_dsc_s(0:M_z,0:1,0:imax1) )
    ALLOCATE( dG1_dsc_s(0:M_z,0:imax1) )

    ALLOCATE( u1_ph_s(0:imax1,0:1) )

    ALLOCATE( st1_rph_s(0:imax1,0:1), q1_ph_s(0:imax1,0:1) )

    ALLOCATE( c4(0:imax1,0:1) )

    ALLOCATE( G1_dsc_trj_s(0:M_z,0:1,0:n_trj_max) )
    ALLOCATE( dG1_dsc_trj_s(0:M_z,0:n_trj_max) )
! ----
    CALL grid_cyl_r1
    CALL uniform
    CALL ordinates

! ---- Data
    CALL read_macro_U
    CALL read_const_inf_U
    CALL read_VDF_U
      
    CALL read_macro_S
    CALL read_const_inf_S
    CALL read_VDF_S
! ---- Trajectory
     CALL get_pts_trj
!*****************************************************************************************
! ---- Main iteration
!*****************************************************************************************
10 CONTINUE

    !CALL integral_a0(Ia)
    !print*, Ia
    !CALL integral_a01(Ia)
    !print*, Ia
    CALL integral_a(Ia)
    print*, Ia
    CALL integral_b(Ib)
    print*, Ib
    !CALL integral_b0(Ib)
    !print*, Ib
    CALL term_c(Ic)
    
    force_lift = Ia + Ib + Ic

    write(*,*)
    print*, 'hl: ', real(force_lift), 'Ia+Ib: ', real(Ia + Ib), 'Ic: ', real(Ic)
    write(*,*)

    name = TRIM(dir_u)//'/force_lift.dat'
    name = trim(name)
    
    WRITE(*,*) 'Data written to ... ', TRIM(name)
    OPEN(UNIT=10, STATUS='unknown', FILE=name)
      WRITE(10,*) force_lift
    CLOSE(10)
    
    name = TRIM(dir_u)//'/force_lift_BIN.dat'
    name = TRIM(name)
    
    WRITE(*,*) 'Data written to ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, FORM='UNFORMATTED')
      WRITE(10) force_lift, Ia, Ib, Ic
    CLOSE(10)
    
    name = TRIM(dir_u)//'/force_lift.txt'
    name = TRIM(name)
    
    WRITE(*,*) 'Data written to ... ', TRIM(name)
    OPEN(UNIT=10, STATUS='unknown', FILE=name)
      WRITE(10,*) real(kn), real(force_lift), real(Ia+Ib), real(Ic)
    CLOSE(10)
!******************************************************************************************
    END PROGRAM main
!******************************************************************************************
