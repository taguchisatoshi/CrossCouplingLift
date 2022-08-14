!*****************************************************************************************
    SUBROUTINE summary
!*****************************************************************************************
    USE global; USE constants; USE coordinates; USE molecular_velocity; USE parameters
    USE counter; USE specifier; USE infinity; USE accel
    IMPLICIT NONE
!*****************************************************************************************
    PRINT*, REPEAT('=', 72)
    PRINT*, '     SUMMARY'
    PRINT*, REPEAT('=', 72)

    IF (kn_u == kn_s) THEN; kn = kn_u; ELSE; STOP; ENDIF
    IF (number_unit_u == number_unit_s) THEN; number_unit = number_unit_u; ELSE; STOP; ENDIF
    IF (L_r_u == L_r_s) THEN; L_r = L_r_u; ELSE; STOP; ENDIF
    IF (imax1_u == imax1_s) THEN; imax1 = imax1_u; ELSE; STOP; ENDIF
    IF (M_ref_u == M_ref_s) THEN; M_ref = M_ref_u; ELSE; STOP; ENDIF
    IF (M_z_u == M_z_s) THEN; M_z = M_z_u; ELSE; STOP; ENDIF
    IF (M_z_th_u == M_z_th_s) THEN; M_z_th = M_z_th_u; ELSE; STOP; ENDIF
    IF (order_ipl_u == order_ipl_s) THEN; order_ipl = order_ipl_u; ELSE; STOP; ENDIF
    IF (order_fd_u == order_fd_s) THEN; order_fd = order_fd_u; ELSE; STOP; ENDIF
    IF (IDfm_u == IDfm_s) THEN; IDfm = IDfm_u; ELSE; STOP; ENDIF
    IF (ID_dsc_u == ID_dsc_s) THEN; ID_dsc = ID_dsc_u; ELSE; STOP; ENDIF
    IF (ID_inf_u == ID_inf_s) THEN; ID_inf = ID_inf_u; ELSE; STOP; ENDIF
    IF (conc_u == conc_s) THEN; conc = conc_u; ELSE; STOP; ENDIF
    IF (ID_dev_u == ID_dev_s) THEN; ID_dev = ID_dev_u; ELSE; STOP; ENDIF
    IF (d1g_u == d1g_s) THEN; d1g = d1g_u; ELSE; STOP; ENDIF
    IF (d2g_u == d2g_s) THEN; d2g = d2g_u; ELSE; STOP; ENDIF
    IF (d3g_u == d3g_s) THEN; d3g = d3g_u; ELSE; STOP; ENDIF
    IF (ID_es_u == ID_es_s) THEN; ID_es = ID_es_u; ELSE; STOP; ENDIF
    IF (nu_u == nu_s) THEN; nu = nu_u; ELSE; STOP; ENDIF
    IF (ID_dev_noneq_u == ID_dev_noneq_s) THEN; ID_dev_noneq = ID_dev_noneq_u; ELSE; STOP; ENDIF

    ID_dsc_analyt = ID_dsc_analyt_s

    PRINT*, REPEAT('*', 36)
    PRINT*, '     * Parameters'
    PRINT*, REPEAT('*', 36)
    
    IF (IDfm /= 'y') WRITE(*,*) '* Knudsen number Kn:', kn

    PRINT*, REPEAT('*', 36)
    PRINT*, '     * Physical space'
    PRINT*, REPEAT('*', 36)
 
    PRINT*, '* Number of subdivisions in the unit length (unit=1): ', number_unit
    PRINT*, '* Length of the computational domain (unit=1): ', L_r
    PRINT*, '* Number of subdivisions in the radial direction: ', imax1
       
    PRINT*, REPEAT('*', 36)
    PRINT*, '     * Velocity space'
    PRINT*, REPEAT('*', 36)
    PRINT*, '* Number of subdivsions in zeta between [0, 5]: ', M_ref
    PRINT*, '* Number of actual subdivsions in zeta: ', M_z
    PRINT*, '* Number of subdivsions in theta_z for [0, pi/2]: ', M_z_th

    PRINT*, REPEAT('*', 36)
    PRINT*, '     * Other parameters'
    PRINT*, REPEAT('*', 36)
    PRINT*, '* Order of interpolation: ', order_ipl
    PRINT*, '* Order of the finite difference scheme: ', order_fd
    WRITE(*,*) '* Discontinuity: ', ID_dsc
    WRITE(*,*) '* Connecting to the asymptotic solution for large r: ', ID_inf
    WRITE(*,*) '* Solve the deviation: ', ID_dev
    WRITE(*,*) '* d1, d2, d3 : ', SNGL(d1g), SNGL(d2g), SNGL(d3g)
    WRITE(*,*) '* ES MODEL: ', ID_es
    IF (ID_es == 1) WRITE(*,*) '* nu: ', nu
    WRITE(*,*)
!*****************************************************************************************
    END SUBROUTINE summary
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE read_macro_U
!*****************************************************************************************
    USE global; USE counter; USE IO
    IMPLICIT NONE
    CHARACTER :: name*99
!*****************************************************************************************
    name = TRIM(dir_u)//'/macro_BIN.dat'
    name = TRIM(name)

    OPEN(UNIT=10, FILE = name, STATUS = 'OLD', &
      FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) rho1_u, u1_r_u, u1_th_u, t1_u, p1_u, &
               st1_rr_u, st1_rth_u, st1_thth_u, st1_phph_u, q1_r_u, q1_th_u, &
               sigma1_u, &
               step_u
    CLOSE(10)
    WRITE(*,*) 'Data read from ... ', TRIM(name)
!*****************************************************************************************
    END SUBROUTINE read_macro_U
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE read_const_inf_U
!*****************************************************************************************
    USE infinity; USE IO
    IMPLICIT NONE
    CHARACTER :: name*99
!*****************************************************************************************
    name = trim(dir_u)//'/const_inf_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, STATUS='OLD', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) c1_inf, c2_inf, c3_inf
    CLOSE(10)
    
    WRITE(*,'(X,A,3ES14.5)') '(c1,c2,c3)=', c1_inf(1), c2_inf(1), c3_inf(1)

    name = trim(dir_u)//'/const_inf_dist_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, STATUS='OLD', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) c1, c2, c3
    CLOSE(10)

!*****************************************************************************************
  END SUBROUTINE read_const_inf_U
!*****************************************************************************************


  
  !*****************************************************************************************
    SUBROUTINE read_VDF_U
!*****************************************************************************************
    USE global; USE VDF; USE counter; USE discontinuity; USE IO
    IMPLICIT NONE
    CHARACTER :: name*99
!*****************************************************************************************
    name = trim(dir_u)//'/VDF_u_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) G1_ne_u, G1_nw_u, H1_ne_u, H1_nw_u
    CLOSE(10)

    name = trim(dir_u)//'/VDF_dsc_u_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) G1_dsc_u, H1_dsc_u, dG1_dsc_u, dH1_dsc_u
    CLOSE(10)

    name = trim(dir_u)//'/VDF_dsc_trj_u_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) G1_dsc_trj_u, H1_dsc_trj_u, dG1_dsc_trj_u, dH1_dsc_trj_u
    CLOSE(10)

!*****************************************************************************************
    END SUBROUTINE read_VDF_U
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE read_macro_S
!*****************************************************************************************
    USE global; USE counter; USE IO
    IMPLICIT NONE
    CHARACTER :: name*99
!*****************************************************************************************
    name = TRIM(dir_s)//'/macro_s_BIN.dat'
    name = TRIM(name)

    OPEN(UNIT=10, FILE=name, STATUS = 'OLD', &
         FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) u1_ph_s, st1_rph_s, q1_ph_s, &
               step_s
    CLOSE(10)
    WRITE(*,*) 'Data read from ... ', trim(name)
!*****************************************************************************************
    END SUBROUTINE read_macro_S
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE read_const_inf_S
!*****************************************************************************************
    USE infinity; USE IO
    IMPLICIT NONE
    CHARACTER :: name*99
!*****************************************************************************************
    name = trim(dir_s)//'/const_inf_s_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, STATUS='OLD', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) c4_inf
    CLOSE(10)
    
    WRITE(*,'(X,A,1ES14.5)') 'c4 =', c4_inf(1)

    name = trim(dir_s)//'/const_inf_dist_s_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, STATUS='OLD', FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) c4
    CLOSE(10)

!*****************************************************************************************
    END SUBROUTINE read_const_inf_S
!*****************************************************************************************





!*****************************************************************************************
    SUBROUTINE read_VDF_S
!*****************************************************************************************
    USE global; USE VDF; USE counter; USE discontinuity; USE IO
    IMPLICIT NONE
    CHARACTER :: name*99
!*****************************************************************************************
    name = trim(dir_s)//'/VDF_s_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) G1_ne_s, G1_nw_s
    CLOSE(10)

    name = trim(dir_s)//'/VDF_dsc_s_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) G1_dsc_s, dG1_dsc_s
    CLOSE(10)

    name = trim(dir_s)//'/VDF_dsc_trj_s_BIN.dat'
    name = TRIM(name)

    WRITE(*,*) 'Data read from ... ', TRIM(name)
    OPEN(UNIT=10, FILE=name, FORM='UNFORMATTED', CONVERT='LITTLE_ENDIAN')
      READ(10) G1_dsc_trj_s, dG1_dsc_trj_s
    CLOSE(10)

!*****************************************************************************************
    END SUBROUTINE read_VDF_S
!*****************************************************************************************
