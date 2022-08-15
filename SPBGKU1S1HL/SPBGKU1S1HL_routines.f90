!*****************************************************************************************
    SUBROUTINE ordinates
!*****************************************************************************************
    USE molecular_velocity
    IMPLICIT NONE
    INTEGER :: k
    REAL(8) :: a0, a, b
!*****************************************************************************************
   
    a0 = 0.002d0   
    a = a0 * 60.d0
    b = 5.d0 - a

    z = 0.d0;  dz = 0.d0

    DO k = 0, M_z
      z(k) = a * DBLE(k)/DBLE(M_ref) + b * (DBLE(k)/DBLE(M_ref))**3
      dz(k) = a/DBLE(M_ref) + 3.d0/DBLE(M_ref) * b * (DBLE(k)/DBLE(M_ref))**2
    ENDDO
!*****************************************************************************************
    END SUBROUTINE ordinates
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE uniform
!*****************************************************************************************
    USE molecular_velocity; USE constants
    IMPLICIT NONE
    INTEGER :: j
!*****************************************************************************************
    z_th_ne = 0.d0; z_th_nw = 0.d0
    dz_th_ne = 0.d0; dz_th_nw = 0.d0
    
    dz_th = pi/2.d0/DBLE(M_z_th)
    
    DO j = 0, M_z_th
      z_th_ne(j) = dz_th * DBLE(j)
      dz_th_ne(j) = dz_th
    ENDDO
    
    DO j = 0, M_z_th
      z_th_nw(j) = pi/2.d0 + z_th_ne(j)
      dz_th_nw(j) = dz_th
    ENDDO
    z_th_nw(M_z_th) = MIN(z_th_nw(M_z_th),pi)

!*****************************************************************************************
    END SUBROUTINE uniform
!*****************************************************************************************



!*****************************************************************************************
      SUBROUTINE grid_cyl_r1
!*****************************************************************************************
      USE coordinates; USE constants; USE MPI_global
      IMPLICIT NONE
      INTEGER :: i
      REAL(8) :: r_min, r_max, radius, dr_tmp
!*****************************************************************************************
      r_min = 1.d0
      r_max = DBLE(L_r)
      radius = 1.d0
      dr_tmp = (r_max - radius)/DBLE(imax1)
!      d1g = d1g0/DBLE(number_unit)/conc
!      d3g = d1g / conc
!      d2g = (d3g*(r_max - r_min) - d1g * (1.d0-EXP(-d3g*imax1))) &
!          / (d3g*imax1 - 1.d0 + EXP(-d3g*imax1))

      DO i = 0, imax1
        IF (conc == 0.d0) THEN
          r1(i) = radius + dr_tmp * DBLE(i)
          dr1(i) = dr_tmp
        ELSEIF (conc > 0.d0) THEN
          !CALL non_uniform(imax1,r_min,r_max,conc,r1(:),dr1(:))
!          CALL non_uniform_b(imax1,r_min,r_max,d1g,d2g,d3g,r1(:),dr1(:))
          !CALL non_uniform_b2(imax1,d1g,d2g,d3g,r1(:),dr1(:))
          CALL non_uniform_a(imax1,d1g,d2g,d3g,r1(:),dr1(:))
          !CALL non_uniform_c(imax1,r_min,r_max,1.d0/conc,r1(:),dr1(:))
        ENDIF
      ENDDO

      r1 = r_min + (r_max-r_min) * r1
      dr1 = (r_max-r_min) * dr1

!      DO j = 0, jmax1
!        DO i = 0, imax1
!          x1(i,j) = r1(i) * COS(th1(j))
!          y1(i,j) = r1(i) * SIN(th1(j))
!        ENDDO
!      ENDDO


!*****************************************************************************************
      END SUBROUTINE grid_cyl_r1
!*****************************************************************************************






!*****************************************************************************************
    SUBROUTINE non_uniform(imax,xmin,xmax,conc,x,dx)
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: xmin, xmax, conc
    INTEGER, INTENT(IN) :: imax
    REAL(8), INTENT(OUT) :: x(0:imax), dx(0:imax)
    INTEGER :: i
!*****************************************************************************************
    DO i = 0, imax
      IF (i == 0) THEN ! to avoid numerical error
        x(0) = xmin
      ELSEIF (i == imax) THEN ! to avoid numerical error
        x(imax) = xmax
      ELSE
        x(i) = xmin + (xmax-xmin)/2.d0 &
             * (TANH(conc*(DBLE(i)/DBLE(imax)-0.5d0))/TANH(conc/2.d0) + 1.d0)
      ENDIF
    ENDDO

    DO i = 0, imax
      dx(i) = (xmax-xmin)/2.d0 &
            / COSH(conc*(DBLE(i)/DBLE(imax)-0.5d0))**2 &
            / TANH(conc/2.d0) &
            * conc/DBLE(imax)
    ENDDO
!*****************************************************************************************
    END SUBROUTINE non_uniform
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE non_uniform_b(imax,xmin,xmax,d1,d2,d3,x,dx)
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: xmin, xmax, d1, d2, d3
    INTEGER, INTENT(IN) :: imax
    REAL(8), INTENT(OUT) :: x(0:imax), dx(0:imax)
    INTEGER :: i
!*****************************************************************************************
    DO i = 0, imax
      IF (i == 0) THEN ! to avoid numerical error
        x(0) = xmin
      ELSEIF (i == imax) THEN ! to avoid numerical error
        x(imax) = xmax
      ELSE
        x(i) = xmin + d2 * DBLE(i) &
             - (d2-d1)/d3 * (1.d0 - EXP(-d3*DBLE(i)))
      ENDIF
    ENDDO

    DO i = 0, imax
      dx(i) = d2 - (d2-d1) * EXP(-d3*DBLE(i))
    ENDDO
!*****************************************************************************************
    END SUBROUTINE non_uniform_b
!*****************************************************************************************





!*****************************************************************************************
    SUBROUTINE non_uniform_b2(imax,d1,d2,d3,x,dx)
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: d1, d2, d3
    INTEGER, INTENT(IN) :: imax
    REAL(8), INTENT(OUT) :: x(0:imax), dx(0:imax)
    REAL(8) :: x0
    INTEGER :: i
!*****************************************************************************************
    DO i = 0, imax
      x(i) = d2 * DBLE(i) &
           - (d2-d1)/d3 * (1.d0 - EXP(-d3*DBLE(i)))
    ENDDO

    DO i = 0, imax
      dx(i) = d2 - (d2-d1) * EXP(-d3*DBLE(i))
    ENDDO
    
    x0 = x(imax)
    x = x/x0
    dx = dx/x0
!*****************************************************************************************
    END SUBROUTINE non_uniform_b2
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE non_uniform_a(imax,d1,d2,d3,x,dx)
!*****************************************************************************************
    USE coordinates
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: d1, d2, d3
    INTEGER, INTENT(IN) :: imax
    REAL(8), INTENT(OUT) :: x(0:imax), dx(0:imax)
    REAL(8) :: x0
    INTEGER :: i
!*****************************************************************************************
    DO i = 0, imax
      x(i) = 1.d0/d3 * LOG(1.d0 - d1/d2 + d1/d2 * EXP(d2*d3*DBLE(i)/DBLE(number_unit)))
    ENDDO

    DO i = 0, imax
      dx(i) = (d2 - (d2-d1) * EXP(-d3*x(i)))/DBLE(number_unit)
    ENDDO
    
    x0 = x(imax)
    x = x/x0
    dx = dx/x0
!*****************************************************************************************
    END SUBROUTINE non_uniform_a
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE non_uniform_c(imax,xmin,xmax,a,x,dx)
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: xmin, xmax, a
    INTEGER, INTENT(IN) :: imax
    REAL(8), INTENT(OUT) :: x(0:imax), dx(0:imax)
    INTEGER :: i
    REAL(8) :: d
!*****************************************************************************************
    d = DBLE(imax) - (1.d0 - EXP(-a*DBLE(imax)))/a
    d = 1.d0/d
    
    DO i = 0, imax
      x(i) = DBLE(i) - (1.d0 - EXP(-a*DBLE(i)))/a
      dx(i) = 1.d0 - EXP(-a*DBLE(i))
    ENDDO
    
    x = d * x
    dx = d * dx
    
    x = xmin + (xmax - xmin) * x
    dx = (xmax - xmin) * dx
!*****************************************************************************************
    END SUBROUTINE non_uniform_c
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE dist_equilibrium(z,z_th,rho,u_r,u_th,t,g,h)
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z, z_th, rho, u_r, u_th, t
    REAL(8), INTENT(OUT) :: g, h
!*****************************************************************************************
    g = rho + 2.d0 * u_r * z * COS(z_th) + (z**2 - 1.5d0) * t
    h = 2.d0 * u_th * z * SIN(z_th)
!*****************************************************************************************
    END SUBROUTINE dist_equilibrium
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE dist_gain_es(z,z_th,rho,u_r,u_th,t,p,st_rr,st_rth,st_thth,g,h)
!*****************************************************************************************
    USE constants
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: z, z_th, rho, u_r, u_th, t, p, st_rr, st_thth, st_rth
    REAL(8), INTENT(OUT) :: g, h
!*****************************************************************************************
    g = rho + 2.d0 * u_r * z * COS(z_th) + (z**2 - 1.5d0) * t &
      + nu * z**2 * (st_rr * COS(z_th)**2 + st_thth * SIN(z_th)**2 - p)
    h = 2.d0 * u_th * z * SIN(z_th) &
      + 2.d0 * nu * st_rth * z**2 * COS(z_th) * SIN(z_th)
!*****************************************************************************************
    END SUBROUTINE dist_gain_es
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE simpson_z2(M_z0,M_z1,M_z_th,weight,f,S)
!*****************************************************************************************
    USE MPI_global
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M_z0, M_z1, M_z_th
    REAL(8), INTENT(IN) :: weight(M_z0:M_z1), f(M_z0:M_z1,0:M_z_th)
    REAL(8), INTENT(OUT) :: S(0:M_z_th)
    INTEGER :: i, j
!*****************************************************************************************
    IF (MOD(M_z1-M_z0,2) /= 0) THEN
      PRINT*, 'Error -- (simpson_z)'
      PRINT*, 'Invalid number of subdivisions'
      CALL MPI_ABORT(MPI_COMM_WORLD,9,ierr)
      CALL MPI_FINALIZE(ierr)
      STOP
    ENDIF

    S = 0.d0

    DO j = 0, M_z_th
      DO i = M_z1-1, M_z0+1, -2
        S(j) = S(j) + weight(i-1)*f(i-1,j) + 4.d0*weight(i)*f(i,j) + weight(i+1)*f(i+1,j)
      ENDDO
    ENDDO
    S = S/3.d0
!*****************************************************************************************
    END SUBROUTINE simpson_z2
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE simpson_z_th(M_z_th,weight,f,S)
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M_z_th
    REAL(8), INTENT(IN) :: weight(0:M_z_th), f(0:M_z_th)
    REAL(8), INTENT(OUT) :: S
    INTEGER :: i
!*****************************************************************************************
    S = 0.d0

    DO i = 1, M_z_th-1, 2
      S = S + weight(i-1)*f(i-1) + 4.d0*weight(i)*f(i) + weight(i+1)*f(i+1)
    ENDDO
    S = S/3.d0
!*****************************************************************************************
    END SUBROUTINE simpson_z_th
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE simpson(M_bgn,M_end,weight,f,S)
!*****************************************************************************************
    USE MPI_global
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M_bgn, M_end
    REAL(8), INTENT(IN) :: weight(M_bgn:M_end), f(M_bgn:M_end)
    REAL(8), INTENT(OUT) :: S
    INTEGER :: i
!*****************************************************************************************
    IF (MOD(M_end-M_bgn,2) /= 0) THEN
      PRINT*, 'Error -- (simpson)'
      PRINT*, 'Invalid number of subdivisions'
      CALL MPI_ABORT(MPI_COMM_WORLD,9,ierr)
      CALL MPI_FINALIZE(ierr)
      STOP
    ENDIF

    S = 0.d0

    DO i = M_end-1, M_bgn+1, -2
      S = S + weight(i-1)*f(i-1) + 4.d0*weight(i)*f(i) + weight(i+1)*f(i+1)
    ENDDO
    S = S/3.d0
!*****************************************************************************************
    END SUBROUTINE simpson
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE interpolation_1d_macro(order,x,rho,u,v,t,xp, &
                                       res_rho,res_u,res_v,res_t)
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: order
    REAL(8), INTENT(IN) :: x(0:order)
    REAL(8), INTENT(IN) :: rho(0:order), u(0:order), v(0:order), t(0:order)
    REAL(8), INTENT(IN) :: xp
    REAL(8), INTENT(OUT) :: res_rho, res_u, res_v, res_t
!*****************************************************************************************
    CALL interpolation_1d(order,x,rho,xp,res_rho)
    CALL interpolation_1d(order,x,u,xp,res_u)
    CALL interpolation_1d(order,x,v,xp,res_v)
    CALL interpolation_1d(order,x,t,xp,res_t)
!*****************************************************************************************
    END SUBROUTINE interpolation_1d_macro
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE interpolation_1d_macro_full(order,x,rho,u,v,stxx,stxy,styy,stzz,qx,qy,xp, &
                                           res_rho,res_u,res_v,res_stxx,res_stxy, &
                                           res_styy,res_stzz,res_qx,res_qy)
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: order
    REAL(8), INTENT(IN) :: x(0:order)
    REAL(8), INTENT(IN) :: rho(0:order), u(0:order), v(0:order), stxx(0:order), &
                           stxy(0:order), styy(0:order), stzz(0:order), &
                           qx(0:order), qy(0:order)
    REAL(8), INTENT(IN) :: xp
    REAL(8), INTENT(OUT) :: res_rho, res_u, res_v, res_stxx, res_stxy, res_styy, &
                            res_stzz, res_qx, res_qy
!*****************************************************************************************
    CALL interpolation_1d(order,x,rho,xp,res_rho)
    CALL interpolation_1d(order,x,u,xp,res_u)
    CALL interpolation_1d(order,x,v,xp,res_v)
    CALL interpolation_1d(order,x,stxx,xp,res_stxx)
    CALL interpolation_1d(order,x,stxy,xp,res_stxy)
    CALL interpolation_1d(order,x,styy,xp,res_styy)
    CALL interpolation_1d(order,x,stzz,xp,res_stzz)
    CALL interpolation_1d(order,x,qx,xp,res_qx)
    CALL interpolation_1d(order,x,qy,xp,res_qy)
!*****************************************************************************************
    END SUBROUTINE interpolation_1d_macro_full
!*****************************************************************************************
