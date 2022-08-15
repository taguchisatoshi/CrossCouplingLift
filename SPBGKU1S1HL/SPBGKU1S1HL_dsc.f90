!*****************************************************************************************
    SUBROUTINE get_pts_trj
!*****************************************************************************************
    USE trajectory; USE molecular_velocity; USE coordinates
    IMPLICIT NONE
    INTEGER :: i, m, l
    LOGICAL :: exist
    REAL(8) :: epsilon
    REAL(8) :: r_add(M_z_th-1), z_th_add(M_z_th-1)
!*****************************************************************************************
    epsilon = 1d-12
    
    r_trj(0:imax1) = r1(0:imax1)
    z_th_trj(0:imax1) = ASIN(1.d0/r1(0:imax1))

    r_add(1:M_z_th-1) = 1.d0/SIN(z_th_ne(M_z_th-1:1:-1))
    z_th_add(1:M_z_th-1) = z_th_ne(M_z_th-1:1:-1)
    
    n_trj = imax1
    DO i = 1, M_z_th-1

      exist = .FALSE.

      DO m = 0, n_trj
        IF (ABS(r_trj(m)-r_add(i)) < epsilon) THEN
          exist = .TRUE.
          EXIT
        ENDIF
      ENDDO
      
      IF ((exist .eqv. .FALSE.) .AND. (r_add(i) < MAXVAL(r1))) THEN
        n_trj = n_trj + 1
        r_trj(n_trj) = r_add(i)
        z_th_trj(n_trj) = z_th_add(i)
      ENDIF
      
    ENDDO

    s_trj(0) = 0.d0
    DO m = 1, n_trj
      s_trj(m) = r_trj(m)*COS(z_th_trj(m))
    ENDDO

    CALL sort_is(s_trj(0:n_trj),r_trj(0:n_trj),z_th_trj(0:n_trj),n_trj)

    DO i = 0, imax1
      DO m = 0, n_trj
        IF (r_trj(m) == r1(i)) THEN
          r_map(i) = m
          EXIT
        ENDIF
      ENDDO
    ENDDO
    
    DO l = 1, M_z_th-1
      DO m = 0, n_trj
        IF (ABS(z_th_trj(m) - z_th_ne(l)) < epsilon) THEN
          z_th_map(l) = m
          EXIT
        ENDIF
      ENDDO
    ENDDO
!*****************************************************************************************
    END SUBROUTINE get_pts_trj
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE sort_is(func,func2,func3,nn)
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER :: i, j
    INTEGER, INTENT(IN) :: nn
    REAL(8), INTENT(INOUT) :: func(0:nn), func2(0:nn), func3(0:nn)
    REAL(8) :: tmp
!*****************************************************************************************
    DO i = nn, 1, -1
      DO j = 0, i-1
        IF (func(j) > func(j+1)) THEN
          tmp = func(j)
          func(j) = func(j+1)
          func(j+1) = tmp
          
          tmp = func2(j)
          func2(j) = func2(j+1)
          func2(j+1) = tmp

          tmp = func3(j)
          func3(j) = func3(j+1)
          func3(j+1) = tmp
        ENDIF
      ENDDO
    ENDDO
!*****************************************************************************************
    END SUBROUTINE sort_is
!*****************************************************************************************
