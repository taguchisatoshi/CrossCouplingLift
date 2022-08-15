!*****************************************************************************************
    SUBROUTINE integral_f_bgk(i,G_s,H_u,ID_dir,fa,fb,fc)
!*****************************************************************************************
    USE specifier; USE molecular_velocity; USE discontinuity; USE coordinates
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    REAL(8), INTENT(IN) :: G_s(0:M_z,0:M_z_th), H_u(0:M_z,0:M_z_th)
    CHARACTER, INTENT(IN) :: ID_dir*2
    REAL(8), INTENT(OUT) :: fa, fb, fc
    REAL(8) :: z_th_dsc, g_dsc_loc(0:M_z,0:1), h_dsc_loc(0:M_z,0:1)
!*****************************************************************************************
    IF (ID_dir == 'ne') THEN
      IF (i /= 0) THEN
      z_th_dsc = ASIN(1.d0/r1(i))
      g_dsc_loc(0:M_z,0:1) = G1_dsc_s(0:M_z,0:1,i)
      h_dsc_loc(0:M_z,0:1) = H1_dsc_u(0:M_z,0:1,i)
        IF (MAXVAL( ABS(g_dsc_loc(0:M_z,0)*EXP(-z(0:M_z)**2 & 
                  -g_dsc_loc(0:M_z,1)*EXP(-z(0:M_z)**2)))) > tor_dsc) THEN
          CALL integral_f_bgk_loc_dsc(z_th_ne(0:M_z_th),dz_th_ne(0:M_z_th), &
                fa,fb,fc, &
                G_s(0:M_z,0:M_z_th),H_u(0:M_z,0:M_z_th), &
                z_th_dsc, &
                g_dsc_loc(0:M_z,0:1),h_dsc_loc(0:M_z,0:1))
        ELSE
          CALL integral_f_bgk_loc_nocor(z_th_ne(0:M_z_th),dz_th_ne(0:M_z_th), &
                fa,fb,fc,G_s(0:M_z,0:M_z_th),H_u(0:M_z,0:M_z_th))
        ENDIF

      ELSE
        CALL integral_f_bgk_loc_nocor(z_th_ne(0:M_z_th),dz_th_ne(0:M_z_th), &
              fa,fb,fc,G_s(0:M_z,0:M_z_th),H_u(0:M_z,0:M_z_th))
    ENDIF

    ELSEIF (ID_dir == 'nw') THEN
      CALL integral_f_bgk_loc_nocor(z_th_nw(0:M_z_th),dz_th_nw(0:M_z_th), &
            fa,fb,fc,G_s(0:M_z,0:M_z_th),H_u(0:M_z,0:M_z_th))
    ENDIF
!*****************************************************************************************
  END SUBROUTINE integral_f_bgk
!*****************************************************************************************

  
  
  
!*****************************************************************************************
    SUBROUTINE integral_f_bgk_loc_nocor(z_th,dz_th_part,fa,fb,fc,G_s,H_u)
!*****************************************************************************************
    USE molecular_velocity; USE constants
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: fa, fb, fc
    REAL(8), INTENT(IN) :: z_th(0:M_z_th), dz_th_part(0:M_z_th)
    REAL(8), INTENT(IN) :: G_s(0:M_z,0:M_z_th), H_u(0:M_z,0:M_z_th)
    REAL(8) :: w_z(0:M_z), w_z_th(0:M_z_th), B=2.d0
!*****************************************************************************************
! --- fa
    w_z = z**3 * dz * EXP(-z**2)/SQRT(pi)/2.d0
    w_z_th = dz_th_part * SIN(z_th)**2

    CALL get_macro_simpson(0,M_z,M_z_th,G_s(0:M_z,:),w_z,w_z_th,fa)
    
! --- fb
    w_z = z**3 * dz * EXP(-z**2)/SQRT(pi)/2.d0
    w_z_th = dz_th_part * SIN(z_th)**2

    CALL get_macro_simpson(0,M_z,M_z_th,H_u(:,:),w_z,w_z_th,fb)

! --- fc
    w_z = z**4 * dz * EXP(-z**2)/SQRT(pi)/2.d0 * B
    w_z_th = dz_th_part * COS(z_th) * SIN(z_th)**2

    CALL get_macro_simpson(0,M_z,M_z_th,H_u(:,:),w_z,w_z_th,fc)

!*****************************************************************************************
  END SUBROUTINE integral_f_bgk_loc_nocor
!*****************************************************************************************


  
  
  

!*****************************************************************************************
    SUBROUTINE integral_f_bgk_loc_dsc(z_th,dz_th_part,fa,fb,fc, &
                                     G_s,H_u, &
                                     z_th_dsc,G_dsc_s,H_dsc_u)
!*****************************************************************************************
    USE molecular_velocity; USE constants
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: fa, fb, fc
    REAL(8), INTENT(IN) :: z_th(0:M_z_th), dz_th_part(0:M_z_th)
    REAL(8), INTENT(IN) :: G_s(0:M_z,0:M_z_th), H_u(0:M_z,0:M_z_th)
    REAL(8), INTENT(IN) :: z_th_dsc, G_dsc_s(0:M_z,0:1), H_dsc_u(0:M_z,0:1)
    REAL(8) :: w_z(0:M_z), w_z_th(0:M_z_th), B=2.d0
    INTEGER :: ID_side, l
    REAL(8) :: f_th_dsc(0:1), f_th(0:M_z_th)
    REAL(8) :: w_z_th_dsc
!*****************************************************************************************
! --- fa
    w_z = z**3 * dz * EXP(-z**2)/SQRT(pi)/2.d0
    w_z_th = dz_th_part * SIN(z_th)**2

    DO l = 0, M_z_th
      CALL simpson(0,M_z,w_z(0:M_z),G_s(0:M_z,l),f_th(l))
    ENDDO

    DO ID_side = 0, 1
      CALL simpson(0,M_z,w_z(0:M_z),G_dsc_s(0:M_z,ID_side),f_th_dsc(ID_side))
    ENDDO
    
    w_z_th_dsc = dz_th * SIN(z_th_dsc)**2
    CALL get_macro_dsc(M_z_th,z_th,w_z_th,f_th,dz_th,z_th_dsc,w_z_th_dsc,f_th_dsc,fa)

! --- fb
    w_z = z**3 * dz * EXP(-z**2)/SQRT(pi)/2.d0
    w_z_th = dz_th_part * SIN(z_th)**2

    DO l = 0, M_z_th
      CALL simpson(0,M_z,w_z(0:M_z),H_u(0:M_z,l),f_th(l))
    ENDDO

    DO ID_side = 0, 1
      CALL simpson(0,M_z,w_z(0:M_z),H_dsc_u(0:M_z,ID_side),f_th_dsc(ID_side))
    ENDDO
    
    w_z_th_dsc = dz_th * SIN(z_th_dsc)**2
    CALL get_macro_dsc(M_z_th,z_th,w_z_th,f_th,dz_th,z_th_dsc,w_z_th_dsc,f_th_dsc,fb)
    
! --- fc
    w_z = z**4 * dz * EXP(-z**2)/SQRT(pi)/2.d0 * B
    w_z_th = dz_th_part * COS(z_th) * SIN(z_th)**2

    DO l = 0, M_z_th
      CALL simpson(0,M_z,w_z(0:M_z),H_u(0:M_z,l),f_th(l))
    ENDDO

    DO ID_side = 0, 1
      CALL simpson(0,M_z,w_z(0:M_z),H_dsc_u(0:M_z,ID_side),f_th_dsc(ID_side))
    ENDDO
    
    w_z_th_dsc = dz_th * COS(z_th_dsc) * SIN(z_th_dsc)**2
    CALL get_macro_dsc(M_z_th,z_th,w_z_th,f_th,dz_th,z_th_dsc,w_z_th_dsc,f_th_dsc,fc)

!*****************************************************************************************
  END SUBROUTINE integral_f_bgk_loc_dsc
!*****************************************************************************************
  
  
  
  
  
!*****************************************************************************************
    SUBROUTINE integral_f_bgk_d(i,G_ne_s,G_nw_s,H_ne_u,H_nw_u,fd)
!*****************************************************************************************
    USE specifier; USE molecular_velocity; USE discontinuity; USE coordinates
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: i
    REAL(8), INTENT(IN) :: G_ne_s(0:M_z,0:M_z_th), H_ne_u(0:M_z,0:M_z_th)
    REAL(8), INTENT(IN) :: G_nw_s(0:M_z,0:M_z_th), H_nw_u(0:M_z,0:M_z_th)
    REAL(8), INTENT(OUT) :: fd
    REAL(8) :: z_th_dsc, g_dsc_loc(0:M_z,0:1), h_dsc_loc(0:M_z,0:1)
!*****************************************************************************************
    IF (i /= 0) THEN
    z_th_dsc = ASIN(1.d0/r1(i))
    g_dsc_loc(0:M_z,0:1) = G1_dsc_s(0:M_z,0:1,i)
    h_dsc_loc(0:M_z,0:1) = H1_dsc_u(0:M_z,0:1,i)
      IF (MAXVAL( ABS(g_dsc_loc(0:M_z,0)*EXP(-z(0:M_z)**2 & 
                -g_dsc_loc(0:M_z,1)*EXP(-z(0:M_z)**2)))) > tor_dsc) THEN
        CALL integral_f_bgk_loc_dsc_d(z_th_ne(0:M_z_th),dz_th_ne(0:M_z_th), &
                                      fd, &
                                      G_ne_s(0:M_z,0:M_z_th),G_nw_s(0:M_z,0:M_z_th), &
                                      H_ne_u(0:M_z,0:M_z_th),H_nw_u(0:M_z,0:M_z_th), &
                                      z_th_dsc, &
                                      g_dsc_loc(0:M_z,0:1),h_dsc_loc(0:M_z,0:1))
      ELSE
        CALL integral_f_bgk_loc_nocor_d(z_th_ne(0:M_z_th),dz_th_ne(0:M_z_th), &
                    fd,G_ne_s(0:M_z,0:M_z_th),G_nw_s(0:M_z,0:M_z_th), &
                    H_ne_u(0:M_z,0:M_z_th),H_nw_u(0:M_z,0:M_z_th))
      ENDIF

    ELSE
      CALL integral_f_bgk_loc_nocor_d(z_th_ne(0:M_z_th),dz_th_ne(0:M_z_th), &
                  fd,G_ne_s(0:M_z,0:M_z_th),G_nw_s(0:M_z,0:M_z_th), &
                  H_ne_u(0:M_z,0:M_z_th),H_nw_u(0:M_z,0:M_z_th))
    ENDIF

!*****************************************************************************************
  END SUBROUTINE integral_f_bgk_d
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE integral_f_bgk_loc_nocor_d(z_th,dz_th_part,fd,G_ne_s,G_nw_s,H_ne_u,H_nw_u)
!*****************************************************************************************
    USE molecular_velocity; USE constants
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: fd
    REAL(8), INTENT(IN) :: z_th(0:M_z_th), dz_th_part(0:M_z_th)
    REAL(8), INTENT(IN) :: G_ne_s(0:M_z,0:M_z_th), H_ne_u(0:M_z,0:M_z_th)
    REAL(8), INTENT(IN) :: G_nw_s(0:M_z,0:M_z_th), H_nw_u(0:M_z,0:M_z_th)
    REAL(8) :: w_z(0:M_z), w_z_th(0:M_z_th), f(0:M_z,0:M_z_th)
    INTEGER :: l
!*****************************************************************************************
! --- fd
    w_z = z**2 * dz * EXP(-z**2)/SQRT(pi)/2.d0
    w_z_th = dz_th_part * SIN(z_th)
    
    DO l = 0, M_z_th
      f(:,l) = H_ne_u(:,l) * G_nw_s(:,M_z_th-l) + G_ne_s(:,l) * H_nw_u(:,M_z_th-l)
    ENDDO

    CALL get_macro_simpson(0,M_z,M_z_th,f(:,:),w_z,w_z_th,fd)

!*****************************************************************************************
  END SUBROUTINE integral_f_bgk_loc_nocor_d
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE integral_f_bgk_loc_dsc_d(z_th,dz_th_part,fd, &
                                 G_ne_s,G_nw_s,H_ne_u,H_nw_u, &
                                 z_th_dsc,G_dsc_s,H_dsc_u)
!*****************************************************************************************
    USE molecular_velocity; USE constants; USE parameters
    IMPLICIT NONE
    REAL(8), INTENT(INOUT) :: fd
    REAL(8), INTENT(IN) :: z_th(0:M_z_th), dz_th_part(0:M_z_th)
    REAL(8), INTENT(IN) :: G_ne_s(0:M_z,0:M_z_th), H_ne_u(0:M_z,0:M_z_th)
    REAL(8), INTENT(IN) :: G_nw_s(0:M_z,0:M_z_th), H_nw_u(0:M_z,0:M_z_th)
    REAL(8), INTENT(IN) :: z_th_dsc, G_dsc_s(0:M_z,0:1), H_dsc_u(0:M_z,0:1)
    REAL(8) :: w_z(0:M_z), w_z_th(0:M_z_th), f(0:M_z,0:M_z_th), f_dsc(0:M_z,0:1)
    INTEGER :: ID_side, l, k
    REAL(8) :: f_th_dsc(0:1), f_th(0:M_z_th), G_s_ipl, H_u_ipl
    REAL(8) :: w_z_th_dsc, z_th_tmp
!*****************************************************************************************
! --- fa
    w_z = z**2 * dz * EXP(-z**2)/SQRT(pi)/2.d0
    w_z_th = dz_th_part * SIN(z_th)

    DO l = 0, M_z_th
      f(:,l) = H_ne_u(:,l) * G_nw_s(:,M_z_th-l) + G_ne_s(:,l) * H_nw_u(:,M_z_th-l)
    ENDDO

    DO l = 0, M_z_th
      CALL simpson(0,M_z,w_z(0:M_z),f(0:M_z,l),f_th(l))
    ENDDO

    z_th_tmp = pi-z_th_dsc
    DO k = 0, M_z
      CALL interpolation_VDF(z_th_tmp,order_ipl,G_nw_s(k,:),H_nw_u(k,:),G_s_ipl,H_u_ipl)
      f_dsc(k,0:1) = H_dsc_u(k,0:1) * G_s_ipl + G_dsc_s(k,0:1) * H_u_ipl
    ENDDO
    DO ID_side = 0, 1
      CALL simpson(0,M_z,w_z(0:M_z),f_dsc(0:M_z,ID_side),f_th_dsc(ID_side))
    ENDDO
    
    w_z_th_dsc = dz_th * SIN(z_th_dsc)
    CALL get_macro_dsc(M_z_th,z_th,w_z_th,f_th,dz_th,z_th_dsc,w_z_th_dsc,f_th_dsc,fd)

!*****************************************************************************************
    END SUBROUTINE integral_f_bgk_loc_dsc_d
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE interpolation_VDF(z_th_tmp,order,G_s,H_u,G_s_ipl,H_u_ipl)
!*****************************************************************************************
    USE molecular_velocity; USE misc
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: order
    REAL(8), INTENT(IN) :: z_th_tmp, G_s(0:M_z_th), H_u(0:M_z_th)
    REAL(8), INTENT(OUT) :: G_s_ipl, H_u_ipl
    INTEGER :: l_tmp
    REAL(8) :: z_th_ipl(0:order), G_ipl(0:order), H_ipl(0:order)
!*****************************************************************************************
    l_tmp = search(z_th_tmp,0,M_z_th,z_th_nw(0:M_z_th))
    l_tmp = MIN(MAX(l_tmp-order/2,0), M_z_th-order)
    z_th_ipl(0:order) = z_th_nw(l_tmp:l_tmp+order)
    G_ipl(0:order) = G_s(l_tmp:l_tmp+order)
    H_ipl(0:order) = H_u(l_tmp:l_tmp+order)
    CALL interpolation_1d(order,z_th_ipl(:),G_ipl(:),z_th_tmp,G_s_ipl)
    CALL interpolation_1d(order,z_th_ipl(:),H_ipl(:),z_th_tmp,H_u_ipl)
!*****************************************************************************************
  END SUBROUTINE interpolation_VDF
!*****************************************************************************************



  
  !*****************************************************************************************
    SUBROUTINE integral_a(Ia)
!*****************************************************************************************
    USE global; USE coordinates; USE constants; USE infinity
    IMPLICIT NONE
    REAL(8) :: f(0:imax1), S, c1_tmp, c4_tmp, gamma1
    REAL(8), INTENT(OUT) :: Ia
!*****************************************************************************************
    gamma1 = 1.d0
    c1_tmp = c1_inf(1)
    c4_tmp = c4_inf(1)
    f(:) = -c4_tmp * (1.d0 + c1_tmp/2.d0/r1(:)) * (st1_rr_u(:,1)-p1_u(:,1)) &
         + r1(:)**2 * (2.d0 * kn * gamma1 * c1_tmp/r1(:)**2 + st1_rr_u(:,1)-p1_u(:,1)) &
         * ( &
            - (1.d0 + c1_tmp/2.d0/r1(:)) * u1_ph_s(:,1) &
            + c4_tmp/r1(:)**2 * u1_th_u(:,1) + u1_th_u(:,1) * u1_ph_s(:,1))

    !CALL trapezoid(imax1,r1(:),f(:),S)
    CALL simpson(0,imax1,dr1(:),f(:),S)

    !do i = 0, imax1, 2
    !  print*, i, f(i)
    !enddo

    Ia = -0.5d0 * S + kn * c1_tmp * c4_tmp * (1.d0 + c1_tmp/4.d0)
    Ia = -8.d0/3.d0 * pi/kn * Ia
!*****************************************************************************************
  END SUBROUTINE integral_a
!*****************************************************************************************


!*****************************************************************************************
    SUBROUTINE integral_a0(Ia)
!*****************************************************************************************
    USE global; USE coordinates; USE constants; USE infinity
    IMPLICIT NONE
    REAL(8) :: f(0:imax1), S, c1_tmp, c4_tmp, gamma1
    REAL(8), INTENT(OUT) :: Ia
    INTEGER :: i
    REAL(8) :: u_th_tmp, u_ph_tmp, p_tmp, st_rr_tmp
!*****************************************************************************************
    gamma1 = 1.d0
    c1_tmp = c1_inf(1)
    c4_tmp = c4_inf(1)

    do i = 0, imax1
      u_th_tmp = u1_th_u(i,1)
      u_ph_tmp = u1_ph_s(i,1)
      p_tmp = p1_u(i,1)
      st_rr_tmp = st1_rr_u(i,1)

      u_th_tmp = u_th_tmp - c1_inf(1)/r1(i)/2.d0 - 1.d0
      u_ph_tmp = u_ph_tmp + c4_inf(1)/r1(i)**2
      p_tmp = p_tmp + kn * gamma1 * c1_inf(1)/r1(i)**2
      st_rr_tmp = st_rr_tmp + 3.d0 * kn * gamma1 * c1_inf(1)/r1(i)**2

      f(i) = r1(i)**2 * u_th_tmp * u_ph_tmp * (st_rr_tmp - p_tmp)
    enddo

    !CALL trapezoid(imax1,r1(:),f(:),S)
    CALL simpson(0,imax1,dr1(:),f(:),S)


    !do i = 0, imax1, 2
    !  print*, i, f(i)
    !enddo

    Ia = -0.5d0 * S
    Ia = -8.d0/3.d0 * pi/kn * Ia
    
!*****************************************************************************************
  END SUBROUTINE integral_a0
!*****************************************************************************************

  
  !*****************************************************************************************
    SUBROUTINE integral_a01(Ia)
!*****************************************************************************************
    USE global; USE coordinates; USE constants; USE infinity
    IMPLICIT NONE
    REAL(8) :: f(0:imax1), S, c1_tmp, c4_tmp, gamma1
    REAL(8), INTENT(OUT) :: Ia
    INTEGER :: i
    REAL(8) :: u_th_tmp, u_ph_tmp, p_tmp, st_rr_tmp
!*****************************************************************************************
    gamma1 = 1.d0
    c1_tmp = c1_inf(1)
    c4_tmp = c4_inf(1)

    do i = 0, imax1
      u_th_tmp = u1_th_u(i,1)
      u_ph_tmp = u1_ph_s(i,1)
      p_tmp = p1_u(i,1)
      st_rr_tmp = st1_rr_u(i,1)

      u_th_tmp = u_th_tmp - c1_inf(1)/r1(i)/2.d0 - 1.d0
      u_ph_tmp = u_ph_tmp + c4_inf(1)/r1(i)**2
      p_tmp = p_tmp + kn * gamma1 * c1_inf(1)/r1(i)**2
      st_rr_tmp = st_rr_tmp + 3.d0 * kn * gamma1 * c1_inf(1)/r1(i)**2

      f(i) = r1(i)**2 * u_th_tmp * u_ph_tmp * (st_rr_tmp - p_tmp)
    enddo

    !CALL trapezoid(imax1,r1(:),f(:),S)
    CALL simpson(0,imax1,dr1(:),f(:),S)


    !do i = 0, imax1, 2
    !  print*, i, f(i)
    !enddo

    Ia = -0.5d0 * S
    Ia = Ia + kn * c1_tmp * c4_tmp/r1(imax1)  ! correction
    Ia = -8.d0/3.d0 * pi/kn * Ia
    
!*****************************************************************************************
  END SUBROUTINE integral_a01
!*****************************************************************************************



  
!*****************************************************************************************
    SUBROUTINE integral_test
!*****************************************************************************************
    USE global; USE coordinates; USE constants; USE infinity
    IMPLICIT NONE
    REAL(8) :: f(0:imax1), S
    INTEGER :: i
!*****************************************************************************************
    
    do i = 0, imax1

      f(i) = 1/r1(i)**2
    enddo

    CALL simpson(0,imax1,dr1(:),f(:),S)
    print*, s
    CALL trapezoid(imax1,r1(:),f(:),S)
    
    !do i = 0, imax1, 2
    !  print*, i, f(i)
    !enddo

    print*, S, 1.d0-1.d0/r1(imax1)
    
!*****************************************************************************************
  END SUBROUTINE integral_test
!*****************************************************************************************



!*****************************************************************************************
  SUBROUTINE trapezoid(imax,x,f,S)
!*****************************************************************************************
  IMPLICIT NONE
  INTEGER :: imax
  REAL(8), INTENT(IN) :: x(0:imax), f(0:imax)
  REAL(8), INTENT(OUT) :: S
  INTEGER :: i
!*****************************************************************************************
  S = 0.d0

  DO i = 0, imax-1
    S = S + 0.5d0 * (f(i) + f(i+1)) * (x(i+1)-x(i)) 
  ENDDO

!*****************************************************************************************
  END SUBROUTINE trapezoid
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE integral_b(Ib)
!*****************************************************************************************
    USE global; USE coordinates; USE constants; USE infinity; USE VDF
    IMPLICIT NONE
    REAL(8) :: f(0:imax1), ff(0:imax1), S, c1_tmp, c4_tmp, c3_tmp, gamma1
    REAL(8), INTENT(OUT) :: Ib
    REAL(8) :: fa, fb, fc, fd
    REAL(8) :: fa_ne, fb_ne, fc_ne, fa_nw, fb_nw, fc_nw
    CHARACTER :: ID_dir*2
    INTEGER :: i
!*****************************************************************************************
    gamma1 = 1.d0
    c1_tmp = c1_inf(1)
    c4_tmp = c4_inf(1)
    c3_tmp = c3_inf(1)

    DO i = 0, imax1
      ID_dir = 'nw'
      CALL integral_f_bgk(i,G1_nw_s(:,:,i),H1_nw_u(:,:,i),ID_dir,fa_nw,fb_nw,fc_nw)
      ID_dir = 'ne'
      CALL integral_f_bgk(i,G1_ne_s(:,:,i),H1_ne_u(:,:,i),ID_dir,fa_ne,fb_ne,fc_ne)
      fa = fa_ne + fa_nw
      fb = fb_ne + fb_nw
      fc = fc_ne + fc_nw

      CALL integral_f_bgk_d(i,G1_ne_s(:,:,i),G1_nw_s(:,:,i),H1_ne_u(:,:,i),H1_nw_u(:,:,i),fd)
      
      f(i) = 2.d0 * (1.d0 + c1_tmp/2.d0/r1(i)) * fa &
           - 2.d0 * c4_tmp/r1(i)**2 * fb &
           + 3.d0 * kn /r1(i)**3 * c4_tmp * fc &
           - fd
    ENDDO

    ff(:) = r1(:)**2 * (c4_tmp/r1(:)**2 + u1_ph_s(:,1)) &
                     * ( -(1.d0 + c1_tmp/r1(:) + u1_r_u(:,1)) * st1_rth_u(:,1) &
                         +(c3_tmp/r1(:)**2 + t1_u(:,1)) * Q1_th_u(:,1) ) &
          + r1(:)**2 * ((kn * gamma1 * c1_tmp - c3_tmp)/r1(:)**2 + rho1_u(:,1)) &
                     * ( -(1.d0 + c1_tmp/2.d0/r1(:)) * u1_ph_s(:,1) &
                         + u1_th_u(:,1) * (c4_tmp/r1(:)**2 + u1_ph_s(:,1)) + f(:))

    !CALL trapezoid(imax1,r1(:),ff(:),S)
    CALL simpson(0,imax1,dr1(:),ff(:),S)

    Ib = -8.d0/3.d0 * pi/kn * S
    
!*****************************************************************************************
  END SUBROUTINE integral_b
!*****************************************************************************************



  
  !*****************************************************************************************
    SUBROUTINE integral_b0(Ib)
!*****************************************************************************************
    USE global; USE coordinates; USE constants; USE infinity; USE VDF
    USE molecular_velocity; USE discontinuity
    IMPLICIT NONE
    REAL(8) :: f(0:imax1), ff(0:imax1), S, c1_tmp, c4_tmp, c3_tmp, gamma1
    REAL(8), INTENT(OUT) :: Ib
    REAL(8) :: fd
    REAL(8) :: B, z_th_dsc
    INTEGER :: i, l, k
    REAL(8) :: rho_tmp, u_r_tmp, u_th_tmp, u_ph_tmp, t_tmp, p_tmp, st_rr_tmp, st_rth_tmp, Q_th_tmp
!*****************************************************************************************
    gamma1 = 1.d0
    c1_tmp = c1_inf(1)
    c4_tmp = c4_inf(1)
    c3_tmp = c3_inf(1)
    B = 2.d0

    DO i = 0, imax1
      DO l = 0, M_z_th
        DO k = 0, M_z
          H1_ne_u(k,l,i) = H1_ne_u(k,l,i) - 2.d0 * (1.d0 + c1_tmp/2.d0/r1(i)) * z(k) * SIN(z_th_ne(l))
          H1_nw_u(k,l,i) = H1_nw_u(k,l,i) - 2.d0 * (1.d0 + c1_tmp/2.d0/r1(i)) * z(k) * SIN(z_th_nw(l))
          G1_ne_s(k,l,i) = G1_ne_s(k,l,i) + c4_tmp * z(k) * SIN(z_th_ne(l)) &
                         * (2.d0/r1(i)**2 + 3.d0 * kn/r1(i)**3 * z(k) * COS(z_th_ne(l)) * B)
          G1_nw_s(k,l,i) = G1_nw_s(k,l,i) + c4_tmp * z(k) * SIN(z_th_nw(l)) &
                         * (2.d0/r1(i)**2 + 3.d0 * kn/r1(i)**3 * z(k) * COS(z_th_nw(l)) * B)
        ENDDO
      ENDDO
    ENDDO
      
    DO i = 0, imax1
      z_th_dsc = ASIN(1.d0/r1(i))
      DO k = 0, M_z
        H1_dsc_u(k,0:1,i) = H1_dsc_u(k,0:1,i) - 2.d0 * (1.d0 + c1_tmp/2.d0/r1(i)) * z(k) * SIN(z_th_dsc)
        G1_dsc_s(k,0:1,i) = G1_dsc_s(k,0:1,i) + c4_tmp * z(k) * SIN(z_th_dsc) &
                          * (2.d0/r1(i)**2 + 3.d0 * kn/r1(i)**3 * z(k) * COS(z_th_dsc) * B)
      ENDDO
    ENDDO

    DO i = 0, imax1

      CALL integral_f_bgk_d(i,G1_ne_s(:,:,i),G1_nw_s(:,:,i),H1_ne_u(:,:,i),H1_nw_u(:,:,i),fd)
      
      f(i) = -fd

    ENDDO

    DO i = 0, imax1
      rho_tmp = rho1_u(i,1)
      u_r_tmp = u1_r_u(i,1)
      u_th_tmp = u1_th_u(i,1)
      u_ph_tmp = u1_ph_s(i,1)
      t_tmp = t1_u(i,1)
      p_tmp = p1_u(i,1)
      st_rr_tmp = st1_rr_u(i,1)
      st_rth_tmp = st1_rth_u(i,1)
      Q_th_tmp = Q1_th_u(i,1)

      rho_tmp = rho_tmp + (kn * gamma1 * c1_inf(1) - c3_inf(1))/r1(i)**2
      u_r_tmp = u_r_tmp + c1_inf(1)/r1(i) + 1.d0
      u_th_tmp = u_th_tmp - c1_inf(1)/r1(i)/2.d0 - 1.d0
      u_ph_tmp = u_ph_tmp + c4_inf(1)/r1(i)**2
      t_tmp = t_tmp + c3_inf(1)/r1(i)**2
      p_tmp = p_tmp + kn * gamma1 * c1_inf(1)/r1(i)**2
      st_rr_tmp = st_rr_tmp + 3.d0 * kn * gamma1 * c1_inf(1)/r1(i)**2

      ff(i) = r1(i)**2 * u_ph_tmp * (-u_r_tmp * st_rth_tmp + t_tmp * Q_th_tmp) &
            + r1(i)**2 * rho_tmp * ( u_th_tmp * u_ph_tmp + f(i) )
    ENDDO

    !CALL trapezoid(imax1,r1(:),ff(:),S)
    CALL simpson(0,imax1,dr1(:),ff(:),S)

    Ib = -8.d0/3.d0 * pi/kn * S
    
!*****************************************************************************************
  END SUBROUTINE integral_b0
!*****************************************************************************************


  !*****************************************************************************************
    SUBROUTINE term_c(Ic)
!*****************************************************************************************
    USE global; USE coordinates; USE constants; USE infinity
    IMPLICIT NONE
    REAL(8) :: c1_tmp, c3_tmp, gamma1
    REAL(8), INTENT(OUT) :: Ic
    REAL(8) :: sigma_tmp, st_rth_tmp
!*****************************************************************************************
    gamma1 = 1.d0/(1.d0-nu)
    c1_tmp = c1_inf(1)
    c3_tmp = c3_inf(1)
    sigma_tmp = sigma1_u(1)
    st_rth_tmp = st1_rth_u(0,1)

    sigma_tmp = sigma_tmp - SQRT(pi) * (1.d0 + c1_tmp) &
              - c3_tmp/2.d0 &
              + 2.d0 * gamma1 * c1_tmp * kn

    Ic = 4.d0/3.d0 * pi * sigma_tmp * st_rth_tmp
!*****************************************************************************************
  END SUBROUTINE term_c
!*****************************************************************************************



