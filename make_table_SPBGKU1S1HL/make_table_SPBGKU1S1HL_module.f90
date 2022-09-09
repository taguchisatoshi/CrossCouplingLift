!*****************************************************************************************
    MODULE global
!*****************************************************************************************
    IMPLICIT NONE
! ---- Region 1
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: rho1, u1_r, u1_th, t1, p1
    REAL(8), DIMENSION(:), ALLOCATABLE :: rho1_ne, u1_r_ne, u1_th_ne, t1_ne, p1_ne
    REAL(8), DIMENSION(:), ALLOCATABLE :: rho1_nw, u1_r_nw, u1_th_nw, t1_nw, p1_nw
    REAL(8), ALLOCATABLE :: sigma1(:)
    
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: st1_rr, st1_rth, st1_thth, st1_phph, &
                                                q1_r, q1_th
    REAL(8), DIMENSION(:), ALLOCATABLE :: st1_rr_ne, st1_rth_ne, st1_thth_ne, &
                                            st1_phph_ne, q1_r_ne, q1_th_ne
    REAL(8), DIMENSION(:), ALLOCATABLE :: st1_rr_nw, st1_rth_nw, st1_thth_nw, &
                                            st1_phph_nw, q1_r_nw, q1_th_nw 

! ---- Force on the sphere
    REAL(8) :: force_lift
    REAL(8), DIMENSION(:), ALLOCATABLE :: h_d
! ---- Variables for the convergence check
    REAL(8), DIMENSION(:), ALLOCATABLE :: delrho1, delu1_r, delu1_th, delt1, delp1, &
                                            delst1_rr, delst1_rth, delst1_thth, delst1_phph, &
                                            delq1_r, delq1_th
    REAL(8) :: maxdelrho1, maxdelu1_r, maxdelu1_th, maxdelt1, maxdelp1, &
              maxdelst1_rr, maxdelst1_rth, maxdelst1_thth, maxdelst1_phph, &
              maxdelq1_r, maxdelq1_th

    REAL(8) :: maxdelmacro1
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
!*****************************************************************************************
    END MODULE molecular_velocity
!*****************************************************************************************



!*****************************************************************************************
    MODULE VDF
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_ne, G1_nw, H1_ne, H1_nw
!    REAL(8), DIMENSION(:,:), ALLOCATABLE :: G1_ne_b, G1_nw_b, H1_ne_b, H1_nw_b
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
    real(8) :: accommo ! -- Accommodation coefficient in the Maxwell's BC
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
!*****************************************************************************************
    END MODULE parameters
!*****************************************************************************************


!*****************************************************************************************
    MODULE counter
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER :: step ! step counter
    INTEGER :: step_out ! step to output data
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
    REAL(8), DIMENSION(0:1) :: c1_inf, c2_inf, c3_inf
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: c1, c2, c3
    REAL(8) :: coeff_under_inf
    type struc_const_inf_ls
      real(8) :: c1, c2, c3
      real(8) :: pw1, pw2, pw3
      real(8) :: dev
      real(8) :: rmin, rmax
      character :: label*6
    end type struc_const_inf_ls
    type(struc_const_inf_ls), allocatable :: const_inf_ls(:)
    real(8), allocatable :: x(:), y_c1(:), y_c2(:), y_c3(:)
!*****************************************************************************************
    END MODULE infinity
!*****************************************************************************************




!*****************************************************************************************
!    MODULE MPI_global
!*****************************************************************************************
!    IMPLICIT NONE
!    INTEGER :: ierr, nprocs, myrank
!    INTEGER, ALLOCATABLE :: k_bgn_MPI(:), k_end_MPI(:)
!    REAL(8), ALLOCATABLE :: QSUM1(:)
!    INTEGER :: ISUMR8
!    INCLUDE 'mpif.h'
!*****************************************************************************************
!    END MODULE MPI_global
!*****************************************************************************************
!
!
!*****************************************************************************************
!    MODULE MPI_input
!*****************************************************************************************
!    IMPLICIT NONE
!    CHARACTER :: filename_MPI
!*****************************************************************************************
!    END MODULE MPI_input
!*****************************************************************************************



!*****************************************************************************************
    MODULE discontinuity
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_dsc, H1_dsc
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: dG1_dsc, dH1_dsc
    ! G1_dsc(k,l,i) : l = 0, 1 (0 for the limit from below, 1 for the limit from above)
    REAL(8), DIMENSION(:,:,:), ALLOCATABLE :: G1_dsc_trj, H1_dsc_trj
    REAL(8), DIMENSION(:,:), ALLOCATABLE :: dG1_dsc_trj, dH1_dsc_trj
    REAL(8) :: tor_dsc = 1d-8
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
    MODULE fmf
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), DIMENSION(:), ALLOCATABLE :: rho_fmf, u_r_fmf, u_t_fmf, t_fmf, p_fmf
    REAL(8), DIMENSION(:), ALLOCATABLE :: st_rr_fmf, st_rt_fmf, st_tt_fmf, &
                                            st_pp_fmf, q_r_fmf, q_t_fmf

    CONTAINS

    FUNCTION integral_fmf(m,n,r) RESULT(S)
    ! integral of g(t) = t^m * SQRT(1-t^2)^(n/2) * SQRT(1-r^2*t^2) from t=0 to t=1/r
    ! m>=1: INTEGER, n: INTEGER
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: m, n
    REAL(8), INTENT(IN) :: r
    REAL(8) :: S
    INTEGER, parameter :: imax = 1600; INTEGER :: i
    REAL(8) :: x(0:imax), dx, f(0:imax), x_min, x_max
    
!    IF (r == 1.d0) THEN
!      x_min = 0.d0; x_max = 1.d0/r
!      dx = (x_max-x_min)/imax
!      DO i = 0, imax
!        x(i) = MIN(x_min + dx*i, x_max)
!        f(i) = x(i)**m * SQRT(MAX(1.d0-x(i)**2,0.d0))**m * SQRT(MAX(1.d0-r**2*x(i)**2,0.d0))
!      ENDDO
!
!    ELSE
      x_min = SQRT(1.d0-1.d0/r**2); x_max = 1.d0
      dx = (x_max-x_min)/imax
      DO i = 0, imax
        x(i) = MIN(x_min + dx*i, x_max)
        f(i) = SQRT(MAX(1.d0-x(i)**2,0.d0))**(m-1) * x(i)**(n+1) * SQRT(MAX(1.d0-r**2+r**2*x(i)**2,0.d0))
      ENDDO
    
!    ENDIF
    
    S = 0.d0
    DO i = imax-1, 1, -2
      S = S + f(i-1) + 4.d0*f(i) + f(i+1)
    ENDDO
    S = S/3.d0 * dx
    
    END function integral_fmf

!*****************************************************************************************
    END MODULE fmf
!*****************************************************************************************






!*****************************************************************************************
    module type_list
!*****************************************************************************************
    implicit none
    type ID_unit
      character :: ID_name*4
      ! log_BIN
      real(8) :: kn, conc, dr1_min, dr1_max, coeff_under_inf, d1g, d2g, d3g, accommo
      integer :: L_r, number_unit, imax1, M_ref, M_z, M_z_th, ID_dev_level, order_fd, &
                 order_ipl, n_crit_conv, step_add_conv, step_aitken, step_out, &
                 ID_es
      character :: ID_dsc, ID_inf, ID_dev, ID_aitken, IDfm, ID_dsc_analyt, ID_dev_noneq
      ! force_lift_BIN
      real(8) :: force_lift
      ! constants
      real(8) :: c1, c2, c3
      real(8) :: c1_ls, c2_ls, c3_ls
      real(8) :: c1_pw, c2_pw, c3_pw, c_dev
      integer :: n_posi
      real(8) :: nu

      type(ID_unit), pointer :: next
    end type ID_unit
    
    integer :: n_ID

    integer :: n_unit
    character :: ID_u_char*99, accommo_char*99
    
    real(8) :: accommo_IO
!*****************************************************************************************
    END MODULE type_list
!*****************************************************************************************



