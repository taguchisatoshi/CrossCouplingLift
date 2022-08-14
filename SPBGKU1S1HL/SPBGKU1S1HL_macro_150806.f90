!*****************************************************************************************
    SUBROUTINE get_macro_simpson(M_z0,M_z1,M_th,f,w,w_th,res)
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M_z0, M_z1, M_th
    REAL(8), INTENT(IN) :: f(M_z0:M_z1,0:M_th), w(M_z0:M_z1), w_th(0:M_th)
    REAL(8), INTENT(OUT) :: res
    REAL(8) :: f_imd(0:M_th)
!*****************************************************************************************
    CALL simpson_z2(M_z0,M_z1,M_th,w,f(:,:),f_imd)
    CALL simpson_z_th(M_th,w_th,f_imd,res)
!*****************************************************************************************
    END SUBROUTINE get_macro_simpson
!*****************************************************************************************






!*****************************************************************************************
    SUBROUTINE get_macro_dsc(M_z_th,z_th,w_z_th,f_th,dz_th_base, &
                              z_th_dsc,w_z_th_dsc,f_th_dsc,macro)
!*****************************************************************************************
    USE misc
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: M_z_th
    REAL(8), INTENT(IN) :: z_th(0:M_z_th), w_z_th(0:M_z_th), f_th(0:M_z_th)
    REAL(8), INTENT(IN) :: dz_th_base
    REAL(8), INTENT(IN) :: z_th_dsc, f_th_dsc(0:1), w_z_th_dsc
    REAL(8), INTENT(OUT) :: macro
    INTEGER :: l, m
    REAL(8) :: S, dS, dz_th
    REAL(8) :: z_th_set(0:3), w_z_th_set(0:3), f_th_set(0:1,0:3)
    INTEGER :: l_tmp, l_tmp_simp
!*****************************************************************************************
    S = 0.d0

    DO l = 1, M_z_th-1, 2
      dS = w_z_th(l-1)*f_th(l-1) + 4.d0*w_z_th(l)*f_th(l) + w_z_th(l+1)*f_th(l+1)
      S = S + dS/3.d0
    ENDDO

! ---- Correction
    l_tmp = search(ABS(z_th_dsc),0,M_z_th,ABS(z_th(0:M_z_th)))
    IF (MOD(l_tmp,2) == 0) THEN
      l_tmp_simp = l_tmp
    ELSE
      l_tmp_simp = l_tmp-1
    ENDIF

    l = l_tmp_simp
    dS = w_z_th(l)*f_th(l) + 4.d0*w_z_th(l+1)*f_th(l+1) + w_z_th(l+2)*f_th(l+2)
    S = S - dS/3.d0

    IF (l_tmp == l_tmp_simp) THEN
      z_th_set(0:3) = (/ z_th(l), z_th_dsc, z_th(l+1:l+2) /)
      w_z_th_set(0:3) = (/ w_z_th(l), w_z_th_dsc, w_z_th(l+1:l+2) /)
      f_th_set(0:1,0) = f_th(l)
      f_th_set(0:1,1) = f_th_dsc(0:1)
      f_th_set(0:1,2) = f_th(l+1)
      f_th_set(0:1,3) = f_th(l+2)

    ELSEIF (l_tmp == l_tmp_simp+1) THEN
      z_th_set(0:3) = (/ z_th(l:l+1), z_th_dsc, z_th(l+2) /)
      w_z_th_set(0:3) = (/ w_z_th(l:l+1), w_z_th_dsc, w_z_th(l+2) /)
      f_th_set(0:1,0) = f_th(l)
      f_th_set(0:1,1) = f_th(l+1)
      f_th_set(0:1,2) = f_th_dsc(0:1)
      f_th_set(0:1,3) = f_th(l+2)

    ENDIF
    
    DO m = 1, 3
      dz_th = ABS(z_th_set(m) - z_th_set(m-1))/dz_th_base
      dS = 0.5d0*(w_z_th_set(m) * f_th_set(0,m) + w_z_th_set(m-1) * f_th_set(1,m-1)) * dz_th
      S = S + dS
    ENDDO

    macro = S

!*****************************************************************************************
    END SUBROUTINE get_macro_dsc
!*****************************************************************************************





!*****************************************************************************************
    SUBROUTINE interpolation(xp,x,f,order,fp)
!*****************************************************************************************
    USE misc
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: xp, x(0:), f(0:SIZE(x)-1)
    INTEGER, INTENT(IN) :: order
    REAL(8), INTENT(OUT) :: fp
    INTEGER :: i_tmp, n
    REAL(8) :: x_ipl(0:order), f_ipl(0:order)
!*****************************************************************************************
    n = SIZE(x)-1
    i_tmp = search(xp,0,n,x(0:n))
    i_tmp = MIN(MAX(i_tmp-order/2,0), n-order)
    x_ipl(0:order) = x(i_tmp:i_tmp+order)
    f_ipl(0:order) = f(i_tmp:i_tmp+order)
    CALL interpolation_1d(order,x_ipl(:),f_ipl(:),xp,fp)
!*****************************************************************************************
  END SUBROUTINE interpolation
!*****************************************************************************************

