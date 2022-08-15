!*****************************************************************************************
    MODULE global
!*****************************************************************************************
    IMPLICIT NONE
! ---- Region 1
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: rho1_u, u1_r_u, u1_th_u, t1_u, p1_u
    REAL(8), ALLOCATABLE :: sigma1_u(:)
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: st1_rr_u, st1_rth_u, st1_thth_u, st1_phph_u, &
                                                q1_r_u, q1_th_u

    REAL(8), DIMENSION(:,:), ALLOCATABLE :: u1_ph_s
    
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: st1_rph_s, q1_ph_s

!*****************************************************************************************
    END MODULE global
!*****************************************************************************************



!*****************************************************************************************
    MODULE coordinates
!*****************************************************************************************
    IMPLICIT NONE
! ---- Space variables

    INTEGER :: number_unit ! number of subdivisions in the unit length = 4
    INTEGER :: L_r ! length of the region 1 in the r direction in the unit 1
    INTEGER :: imax1
    REAL(8), ALLOCATABLE :: r1(:), dr1(:)
    REAL(8), ALLOCATABLE :: x1(:,:), y1(:,:)

!    REAL(8) :: dx, dy ! grid spacing
    REAL(8) :: r_d
    REAL(8) :: conc, d1g, d2g, d3g, d1g0
    REAL(8) :: dr1_max, dr1_min

    INTEGER :: number_unit_u, number_unit_s
    INTEGER :: L_r_u, L_r_s
    INTEGER :: imax1_u, imax1_s
    REAL(8) :: conc_u, d1g_u, d2g_u, d3g_u
    REAL(8) :: conc_s, d1g_s, d2g_s, d3g_s
!*****************************************************************************************
    END MODULE coordinates
!*****************************************************************************************



!*****************************************************************************************
    MODULE molecular_velocity
!*****************************************************************************************
    IMPLICIT NONE
! ---- Molecular velocity
    INTEGER :: M_z, M_z_th, M_ref
    REAL(8), DIMENSION(:), ALLOCATABLE :: z, dz
    REAL(8), DIMENSION(:), ALLOCATABLE :: z_th_ne, dz_th_ne
    REAL(8), DIMENSION(:), ALLOCATABLE :: z_th_nw, dz_th_nw
    REAL(8) :: dz_th
    real(8) :: grid_conc_z_th_nw

    INTEGER :: M_z_u, M_z_th_u, M_ref_u, M_z_s, M_z_th_s, M_ref_s
!*****************************************************************************************
    END MODULE molecular_velocity
!*****************************************************************************************



!*****************************************************************************************
    MODULE VDF
!*****************************************************************************************
    IMPLICIT NONE

    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_ne_u, G1_nw_u, H1_ne_u, H1_nw_u
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_ne_s, G1_nw_s
!*****************************************************************************************
    END MODULE VDF
!*****************************************************************************************


!*****************************************************************************************
    MODULE constants
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8) :: pi
    REAL(8) :: kn ! -- Knudsen number
    REAL(8) :: nu ! -- Controle parameter of the Prandtl number in ES model
    REAL(8) :: prandtl
    REAL(8) :: accommo

    REAL(8) :: kn_u, kn_s
    REAL(8) :: nu_u, nu_s
    REAL(8) :: accommo_u, accommo_s
!*****************************************************************************************
  END MODULE constants
!*****************************************************************************************



  
!*****************************************************************************************
    MODULE parameters
!*****************************************************************************************
    IMPLICIT NONE
! ---- Parameter
    INTEGER :: order_ipl ! -- order of interpolation
    INTEGER :: order_fd ! -- order of finite difference approximation
    INTEGER :: n_crit_conv   ! convergence criteria 10^n_crit (n_crit < 0)
    INTEGER :: step_add_conv ! steps to add the variation

    INTEGER :: order_ipl_u, order_ipl_s
    INTEGER :: order_fd_u, order_fd_s
    INTEGER :: n_crit_conv_u, n_crit_conv_s
    INTEGER :: step_add_conv_u, step_add_conv_s
!*****************************************************************************************
    END MODULE parameters
!*****************************************************************************************


!*****************************************************************************************
    MODULE counter
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER :: step ! step counter
    INTEGER :: step_out ! step to output data

    INTEGER :: step_u, step_s
    INTEGER :: step_out_u, step_out_s
!*****************************************************************************************
    END MODULE counter
!*****************************************************************************************


!*****************************************************************************************
    MODULE specifier
!*****************************************************************************************
    IMPLICIT NONE
    CHARACTER :: ID_dsc
    INTEGER :: ID_fmf = 0
    CHARACTER :: IDfm
    CHARACTER :: ID_inf
    CHARACTER :: ID_dev
    INTEGER :: ID_dev_level
    CHARACTER:: ID_aitken
    character :: ID_dsc_analyt
    integer :: ID_es
    integer :: ID_VDF_file
    integer :: ID_const_inf_fix
    CHARACTER :: ID_dev_noneq

    CHARACTER :: ID_dsc_u, ID_dsc_s
    CHARACTER :: IDfm_u, IDfm_s
    CHARACTER :: ID_inf_u, ID_inf_s
    CHARACTER :: ID_dev_u, ID_dev_s
    INTEGER :: ID_dev_level_u, ID_dev_level_s
    CHARACTER:: ID_aitken_u, ID_aitken_s
    character :: ID_dsc_analyt_u, ID_dsc_analyt_s
    integer :: ID_es_u, ID_es_s
    integer :: ID_VDF_file_u, ID_VDF_file_s
    integer :: ID_const_inf_fix_u, ID_const_inf_fix_s
    CHARACTER :: ID_dev_noneq_u, ID_dev_noneq_s
!*****************************************************************************************
    END MODULE specifier
!*****************************************************************************************



!*****************************************************************************************
    MODULE dsc_limit
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8) :: inf_limit = 1d-12
!*****************************************************************************************
    END MODULE dsc_limit
!*****************************************************************************************


!*****************************************************************************************
    MODULE trajectory
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), ALLOCATABLE :: r_trj(:), z_th_trj(:), s_trj(:)
    INTEGER :: n_trj, n_trj_max
    INTEGER, ALLOCATABLE :: r_map(:), z_th_map(:)
!*****************************************************************************************
    END MODULE trajectory
!*****************************************************************************************



!*****************************************************************************************
    MODULE infinity
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(0:1) :: c1_inf, c2_inf, c3_inf, c4_inf
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: c1, c2, c3, c4
    REAL(8) :: coeff_under_inf
    type struc_const_inf_ls
      real(8) :: c1, c2, c3, c4
      real(8) :: pw1, pw2, pw3
      real(8) :: dev
      real(8) :: rmin, rmax
      character :: label*6
    end type struc_const_inf_ls
    type(struc_const_inf_ls), allocatable :: const_inf_ls(:)
    real(8), allocatable :: x(:), y_c1(:), y_c2(:), y_c3(:), y_c4(:)

    REAL(8) :: coeff_under_inf_u, coeff_under_inf_s
!*****************************************************************************************
    END MODULE infinity
!*****************************************************************************************




!!*****************************************************************************************
!    MODULE MPI_global
!!*****************************************************************************************
!    IMPLICIT NONE
!    INTEGER :: ierr, nprocs, myrank
!    INTEGER, ALLOCATABLE :: k_bgn_MPI(:), k_end_MPI(:)
!    REAL(8), ALLOCATABLE :: QSUM1(:)
!    INTEGER :: ISUMR8
!    INCLUDE 'mpif.h'
!!*****************************************************************************************
!    END MODULE MPI_global
!!*****************************************************************************************
!
!
!!*****************************************************************************************
!    MODULE MPI_input
!!*****************************************************************************************
!    IMPLICIT NONE
!    CHARACTER :: filename_MPI
!!*****************************************************************************************
!    END MODULE MPI_input
!!*****************************************************************************************



!*****************************************************************************************
    MODULE discontinuity
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8) :: tor_dsc = 1d-8
    
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_dsc_u, H1_dsc_u
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: dG1_dsc_u, dH1_dsc_u
    ! G1_dsc(k,l,i) : l = 0, 1 (0 for the limit from below, 1 for the limit from above)
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_dsc_trj_u, H1_dsc_trj_u
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: dG1_dsc_trj_u, dH1_dsc_trj_u

    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_dsc_s
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: dG1_dsc_s
    ! G1_dsc(k,l,i) : l = 0, 1 (0 for the limit from below, 1 for the limit from above)
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_dsc_trj_s
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: dG1_dsc_trj_s
!*****************************************************************************************
    END MODULE discontinuity
!*****************************************************************************************



!*****************************************************************************************
    MODULE misc
!*****************************************************************************************
    IMPLICIT NONE
    
    CONTAINS

    FUNCTION search(xp,istart,iend,f) RESULT(i_low)
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: xp
    INTEGER, INTENT(IN) :: istart, iend
    REAL(8), INTENT(IN) :: f(istart:iend)
    INTEGER :: i_mid, i_low, i_high
    REAL(8) :: epsilon = 0.d0

    IF (ABS(xp-f(istart)) <= epsilon ) THEN
      i_low = istart
      RETURN
    ELSEIF (ABS(xp-f(iend)) <= epsilon) THEN
      i_low = iend-1
      RETURN
    ELSEIF (xp < f(istart)-epsilon .OR. xp > f(iend)+epsilon) THEN
      WRITE(*,*) 'Error -- (search)'
      WRITE(*,*) xp
      STOP
    ENDIF
          
    i_low = istart
    i_high = iend
    
    DO
      i_mid = (i_low + i_high)/2

      IF (xp >= f(i_mid)) THEN
        i_low = i_mid
      ELSEIF (xp < f(i_mid)) THEN
        i_high = i_mid
      ENDIF
      
      IF (i_high - i_low == 1) THEN
        EXIT
      ENDIF
      
    ENDDO

    END FUNCTION search

    END MODULE misc
!*****************************************************************************************




  
!*****************************************************************************************
    MODULE IO
!*****************************************************************************************
    IMPLICIT NONE
    CHARACTER :: dir_u*128, ID_u_char*128, dir_s*128, ID_s_char*128
    INTEGER :: ID_u, ID_s
    integer :: n_unit
!*****************************************************************************************
  END MODULE IO
!*****************************************************************************************
