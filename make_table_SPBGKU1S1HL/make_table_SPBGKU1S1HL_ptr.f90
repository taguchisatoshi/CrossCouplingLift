!*****************************************************************************************
    SUBROUTINE get_list_ID(ID)
!*****************************************************************************************
    use type_list
    implicit none
    !character :: ID_u_tmp*5, ID_s_tmp*5
    character :: ID_u_tmp*5
    type(ID_unit), pointer :: ID
    type(ID_unit), pointer :: tmpptr
    integer :: InputStatus, AllocateStatus
    real(8) :: kn_tmp, conc_tmp, coeff_under_inf_tmp, d1g_tmp, d2g_tmp, d3g_tmp, nu_tmp, &
               accommo_tmp
    integer :: number_unit_tmp, L_r_tmp, imax1_tmp, M_ref_tmp, M_z_tmp, M_z_th_tmp, &
               order_ipl_tmp, order_fd_tmp, step_out_tmp, &
               n_crit_conv_tmp, step_add_conv_tmp, ID_dev_level_tmp, &
               step_aitken_tmp, ID_es_tmp
    character :: IDfm_tmp, ID_dsc_tmp, ID_inf_tmp, ID_dev_tmp, ID_aitken_tmp, &
                 ID_dsc_analyt_tmp, ID_dev_noneq_tmp
    real(8) :: force_lift_tmp, tmp
    !real(8) :: c1_tmp(0:1), c2_tmp(0:1), c3_tmp(0:1)
    character :: name*10, name1*99, file*99, dir*99, char_accommo*99
    !real(8) :: c1_pw_tmp, c2_pw_tmp, c3_pw_tmp, c_dev_tmp
    !integer :: n_posi_tmp
    INTEGER :: sz
!*****************************************************************************************
!    ALLOCATE(ID)
!    PRINT*, ASSOCIATED(ID)
    
    print*, 'accommo ?'
    read(*,*) accommo_IO

    if (accommo_IO == 0.8d0) then
      char_accommo = '0.8'
    elseif (accommo_IO == 0.5d0) then
      char_accommo = '0.5'
    else
      print*, 'Invalid value of the accommodation coefficient'
      stop
    endif
    
    print*, 'n_unit ?'
    read(*,*) n_unit
    write(ID_u_char,*) n_unit

    dir = 'H:/work2/SPBGK/20220604_SPES20U1M/ID'
    file = 'H:/work2/SPBGK/20220604_SPES20U1M/ID_list_accommo'//TRIM(ADJUSTL(char_accommo))//'_n'//TRIM(ADJUSTL(ID_u_char))//'.txt'
    !file = 'ID_list.txt'
    
    n_ID = 0
    OPEN(UNIT=10, STATUS ='OLD', FILE = TRIM(file))
      DO
        !READ(10,*,IOSTAT = InputStatus) kn_tmp, ID_u_tmp, ID_s_tmp
        READ(10,*,IOSTAT = InputStatus) ID_u_tmp
        IF (InputStatus < 0) EXIT
        ALLOCATE(tmpptr, STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
        tmpptr%ID_name = ID_u_tmp
        tmpptr%next => ID
        ID => tmpptr
        n_ID = n_ID+1
      ENDDO
    CLOSE(10)

    tmpptr => ID
    DO WHILE (ASSOCIATED(tmpptr))
      name = TRIM(tmpptr%ID_name)
      name1 = trim(dir)//TRIM(name)//'/log_BIN.dat'
      OPEN(UNIT=10, STATUS='OLD', FILE=name1, FORM='UNFORMATTED')
        INQUIRE (FILE = name1, SIZE = sz) ! get the file size

        READ(10) kn_tmp, number_unit_tmp, L_r_tmp, imax1_tmp, M_ref_tmp, M_z_tmp, M_z_th_tmp, &
                  order_ipl_tmp, order_fd_tmp, step_out_tmp, IDfm_tmp, ID_dsc_tmp, ID_inf_tmp, conc_tmp, &
                  n_crit_conv_tmp, step_add_conv_tmp, ID_dev_tmp, ID_dev_level_tmp, &
                  coeff_under_inf_tmp, d1g_tmp, d2g_tmp, d3g_tmp, &
                  ID_es_tmp, nu_tmp, ID_dev_noneq_tmp, accommo_tmp
        ID_dsc_analyt_tmp = 'y'

      CLOSE(10)
      tmpptr%kn = kn_tmp
      tmpptr%number_unit = number_unit_tmp
      tmpptr%L_r = L_r_tmp
      tmpptr%imax1 = imax1_tmp
      tmpptr%M_ref = M_ref_tmp
      tmpptr%M_z = M_z_tmp
      tmpptr%M_z_th = M_z_th_tmp
      tmpptr%order_ipl = order_ipl_tmp
      tmpptr%order_fd = order_fd_tmp
      tmpptr%step_out = step_out_tmp
      tmpptr%IDfm = IDfm_tmp
      tmpptr%ID_dsc = ID_dsc_tmp
      tmpptr%ID_dsc_analyt = ID_dsc_analyt_tmp
      tmpptr%ID_inf = ID_inf_tmp
      tmpptr%conc = conc_tmp
      tmpptr%n_crit_conv = n_crit_conv_tmp
      tmpptr%step_add_conv = step_add_conv_tmp
      tmpptr%ID_dev = ID_dev_tmp
      tmpptr%ID_dev_level = ID_dev_level_tmp
      tmpptr%coeff_under_inf = coeff_under_inf_tmp
      tmpptr%ID_aitken = ID_aitken_tmp
      tmpptr%step_aitken = step_aitken_tmp
      tmpptr%d1g = d1g_tmp
      tmpptr%d2g = d2g_tmp
      tmpptr%d3g = d3g_tmp
      tmpptr%ID_es = ID_es_tmp
      tmpptr%nu = nu_tmp
      tmpptr%ID_dev_noneq = ID_dev_noneq_tmp
      tmpptr%accommo = accommo_tmp

      name1 = trim(dir)//TRIM(name)//'/force_lift_BIN.dat'
      OPEN(UNIT=10, STATUS='OLD', FILE=name1, FORM='UNFORMATTED')
        READ(10) force_lift_tmp, tmp, tmp, tmp
      CLOSE(10)
      tmpptr%force_lift = force_lift_tmp

      tmpptr => tmpptr%next
    ENDDO

    write(accommo_char,'(1f3.1)') accommo_tmp
    print*, accommo_char

!*****************************************************************************************
    END SUBROUTINE get_list_ID
!*****************************************************************************************


!*****************************************************************************************
    RECURSIVE SUBROUTINE merge_sort_ID(x)
!*****************************************************************************************
    USE type_list
    TYPE(ID_unit), POINTER :: a, b, x, y, z
    
    INTERFACE
      SUBROUTINE merge_list_ID(x,y,z)
        USE type_list
        TYPE(ID_unit), POINTER :: x, y, z
      END SUBROUTINE merge_list_ID
    END INTERFACE
!*****************************************************************************************
    ! no need for sorting
    IF ((ASSOCIATED(x) .eqv. .FALSE.) .OR. (ASSOCIATED(x%next) .eqv. .FALSE.)) THEN
      z => x
      RETURN
    ENDIF
    
    ! divide the list into two parts
    a => x
    b => x%next
    IF (ASSOCIATED(b)) b => b%next
       
    DO WHILE(ASSOCIATED(b))
      a => a%next
      b => b%next
      IF (ASSOCIATED(b)) b => b%next
    ENDDO

    y => a%next      ! make the head of the second list
    NULLIFY(a%next) ! make the end of the first list

    CALL merge_sort_ID(x)
    CALL merge_sort_ID(y)
    CALL merge_list_ID(x,y,z)
    x => z
!    print*, ''

!*****************************************************************************************
    END SUBROUTINE merge_sort_ID
!*****************************************************************************************


!*****************************************************************************************
    SUBROUTINE merge_list_ID(x,y,z1)
!*****************************************************************************************
    USE type_list
    IMPLICIT NONE
    TYPE(ID_unit), POINTER :: x, y, z1, curr, z
!*****************************************************************************************
    !ALLOCATE(z,z1)
    ALLOCATE(z)
    curr => z

    DO WHILE (ASSOCIATED(x) .AND. ASSOCIATED(y))
      IF (x%kn <= y%kn) THEN
        curr%next => x
        curr => x
        x => x%next
      ELSE
        curr%next => y
        curr => y 
        y => y%next
      ENDIF
    ENDDO

    IF (ASSOCIATED(x) .eqv. .FALSE.) THEN
      curr%next => y
    ELSE
      curr%next => x
    ENDIF
    
    z1 => z%next ! the head of z is dummy
!*****************************************************************************************
    END SUBROUTINE merge_list_ID
!*****************************************************************************************


!*****************************************************************************************
    SUBROUTINE printlist_ID(x)
!*****************************************************************************************
    USE type_list
    IMPLICIT NONE
    TYPE(ID_unit), POINTER :: x, currptr
!*****************************************************************************************
    currptr => x
    DO WHILE (ASSOCIATED(currptr))
      PRINT*, currptr%ID_name, real(currptr%kn)
      currptr => currptr%next
    ENDDO
!*****************************************************************************************
    END SUBROUTINE printlist_ID
!*****************************************************************************************


!*****************************************************************************************
    SUBROUTINE write_list_ID(x)
!*****************************************************************************************
    USE type_list
    IMPLICIT NONE
    TYPE(ID_unit), POINTER :: x, currptr
    CHARACTER :: name*128
!*****************************************************************************************
    currptr => x

    !name = 'ID_list_summary.txt'
    
    name = 'ID_list_summary_accommo'//trim(adjustl(accommo_char))//'_n'//TRIM(ADJUSTL(ID_u_char))//'.txt'

    !    WRITE(*,*) 'Data written to ...', TRIM(name)
    OPEN(UNIT=10, STATUS='unknown', FILE=name)
      WRITE(10,'(1X,A2)',advance='no') 'ID'
      write(10,'(1X,A6)',advance='no') 'kn'
      write(10,'(3X,A)',advance='no') 'L_r'
      write(10,'(2X,A)',advance='no') 'n_unit'
      write(10,'(3X,A)',advance='no') 'imax1'
      write(10,'(3X,A)',advance='no') 'M_z'
      write(10,'(3X,A)',advance='no') 'M_z_th'
      write(10,'(2X,A)',advance='no') 'd1g'
      write(10,'(8X,A)',advance='no') 'd2g'
      write(10,'(8X,A)',advance='no') 'd3g'
      write(10,'(8X,A)',advance='no') 'IP'
      write(10,'(3X,A)',advance='no') 'FD'
      !write(10,'(3X,A)',advance='no') 'conc'
      write(10,'(1X,A)',advance='no') 'IDfm'
      write(10,'(1X,A)',advance='no') 'ID_dsc'
      write(10,'(1X,A)',advance='no') 'ID_dsc_analyt'
      write(10,'(1X,A)',advance='no') 'ID_inf'
      write(10,'(1X,A)',advance='no') 'n_crit_conv'
      write(10,'(1X,A)',advance='no') 'step_add_conv'
      write(10,'(1X,A)',advance='no') 'ID_dev'
      write(10,'(1X,A)',advance='no') 'ID_dev_level'
      write(10,'(1X,A)',advance='no') 'ID_dev_noneq'
      write(10,'(1X,A)',advance='no') 'coeff_under_inf'
      write(10,'(3X,A)',advance='no') 'accommo'
      write(10,'(1X,A)',advance='no') 'ID_es'
      write(10,'(3X,A)',advance='yes') 'nu'

      WRITE(10,'(A)') REPEAT('-',212)
      DO WHILE (ASSOCIATED(currptr))
        write(10,'(1X,A4)',advance='no') currptr%ID_name
        write(10,'(F6.3)',advance='no') currptr%kn
        write(10,'(I5)',advance='no') currptr%L_r
        write(10,'(2X,I5)',advance='no') currptr%number_unit
        write(10,'(3X,I5)',advance='no') currptr%imax1
        write(10,'(1X,I6)',advance='no') currptr%M_z
        write(10,'(2X,I5)',advance='no') currptr%M_z_th
        write(10,'(2X,ES11.3)',advance='no') currptr%d1g
        write(10,'(ES11.3)',advance='no') currptr%d2g
        write(10,'(ES11.3)',advance='no') currptr%d3g
        write(10,'(I4)',advance='no') currptr%order_ipl
        write(10,'(I4)',advance='no') currptr%order_fd
        !write(10,'(F7.1)',advance='no') currptr%conc
        write(10,'(1X,A4)',advance='no') currptr%IDfm
        write(10,'(1X,A4)',advance='no') currptr%ID_dsc
        write(10,'(8X,A4)',advance='no') currptr%ID_dsc_analyt
        write(10,'(2X,A8)',advance='no') currptr%ID_inf
        write(10,'(7X,I5)',advance='no') currptr%n_crit_conv
        write(10,'(5X,I5)',advance='no') currptr%step_add_conv
        write(10,'(10X,A)',advance='no') currptr%ID_dev
        write(10,'(5X,I3)',advance='no') currptr%ID_dev_level
        write(10,'(10X,A4)',advance='no') currptr%ID_dev_noneq
        write(10,'(8X,F11.8)',advance='no') currptr%coeff_under_inf
        write(10,'(6X,F6.3)',advance='no') currptr%accommo
        write(10,'(1X,I5)',advance='no') currptr%ID_es
        write(10,'(3X,F6.3)',advance='yes') currptr%nu
        currptr => currptr%next
      ENDDO
    CLOSE(10)

    PRINT*, 'File created : ', TRIM(name)
!*****************************************************************************************
    END SUBROUTINE write_list_ID
!*****************************************************************************************




!*****************************************************************************************
    subroutine write_force_lift_ID(x)
!*****************************************************************************************
    use type_list; use constants
    implicit none
    type(ID_unit), pointer :: x, currptr
    character :: name*128
    integer :: i, imax, id_es
    real(8) :: kn_max, kn_min, dkn_tmp, kn_tmp, gamma1, k0, a1, a2, a3, b1, Ad, hl_tmp
    real(8) :: gamma1_hs
    real(8) :: kn_th_tmp, accommo_tmp
!*****************************************************************************************
    currptr => x

    if (currptr%ID_es == 0) then
    ! BGK
      gamma1 = 1.d0
      k0 = -1.01619d0
      a1 = 0.76632d0
      a2 = 0.5d0
      a3 = -0.26632d0
      b1 = 0.11684d0
    elseif (currptr%ID_es == 1) then
    ! ES
      nu = currptr%nu
      gamma1 = 1.d0/(1.d0-nu)
      k0 = -1.01619d0/(1.d0-nu)
      a1 = 0.76632d0/(1.d0-nu)
      a2 = 0.5d0/(1.d0-nu)
      a3 = -0.26632d0/(1.d0-nu)
      b1 = 0.11684d0/(1.d0-nu)
    endif
    Ad = 3.d0*k0**2 - 4.d0*a1 + a2 + a3 + 2.d0*b1

    
    ! HSS
!    k0 = -1.2540d0
    gamma1_hs = 1.270042427d0
    accommo_tmp = currptr%accommo

    if (currptr%ID_es == 0) then
      !name = 'hl_vs_kn_bgk.dat'
      name = 'hl_vs_kn_bgk_accommo'//trim(adjustl(accommo_char))//'_n'//TRIM(ADJUSTL(ID_u_char))//'.dat'
      id_es = 0
    elseif (currptr%ID_es == 1) then
      !name = 'hl_vs_kn_es.dat'
      name = 'hl_vs_kn_es_acccommo'//trim(adjustl(accommo_char))//'_n'//TRIM(ADJUSTL(ID_u_char))//'.dat'
      id_es = 1
    endif
    open(unit=10, status='unknown', file=name)
      write(10,'(1X,A)') 'VARIABLES = "kn", "hl"'

      if (currptr%ID_es == 0) then
        write(10,*) 'ZONE T = "num (BGK) : Accommo ='//trim(adjustl(accommo_char))//'"'
      elseif (currptr%ID_es == 1) then
        write(10,*) 'ZONE T = "num (ES) : Accommo ='//trim(adjustl(accommo_char))//'"'
      endif
      !write(10,*) 'ZONE T = "kn vs hl"'
      do while (associated(currptr))
        write(10,'(2ES16.8)') currptr%kn, currptr%force_lift
        currptr => currptr%next
      enddo
      write(10,*)
      write(10,*)

      imax = 50
      kn_min = 0.001d0; kn_max = 0.1d0
      dkn_tmp = (kn_max-kn_min)/dble(imax)
      if (id_es == 0) then; write(10,*) 'ZONE T = "asy (BGK 0)"'
      elseif (id_es == 0) then; write(10,*) 'ZONE T = "asy (ES 0)"'
      endif
      hl_tmp = 2.d0 * pi
      write(10,'(2ES16.8)') kn_min, hl_tmp
      write(10,'(2ES16.8)') kn_max, hl_tmp
      write(10,*)
      write(10,*)
      
      !imax = 60
      !kn_min = 0.001d0; kn_max = 0.15d0
      !dkn_tmp = (kn_max-kn_min)/dble(imax)
      !if (id_es == 0) then; write(10,*) 'ZONE T = "asy (BGK 1)"'
      !elseif (id_es == 0) then; write(10,*) 'ZONE T = "asy (ES 2)"'
      !endif
      !do i = 0, imax
      !  kn_tmp = kn_min + dkn_tmp * i
      !  hl_tmp = 2.d0 * pi * (1.d0 + 3.d0 *  k0 * kn_tmp)
      !  write(10,'(2ES16.8)') kn_tmp, hl_tmp
      !enddo
      !write(10,*)
      !write(10,*)

      write(10,*) 'ZONE T = "asy fmf : Accommo ='//trim(adjustl(accommo_char))//'"'
      hl_tmp = -4.d0/3.d0 * pi * accommo_tmp
      write(10,'(2ES16.8)') 10.d0, hl_tmp
      write(10,'(2ES16.8)') 100.d0, hl_tmp

      write(10,*)
      write(10,*)

      write(10,*) 'ZONE T = "hl=0"'
      hl_tmp = 0.d0
      write(10,'(2ES16.8)') kn_min, hl_tmp
      write(10,'(2ES16.8)') 100.d0, hl_tmp

      write(10,*)
      write(10,*)

      write(10,*) 'ZONE T = "hl=0 (2)"'
      hl_tmp = 0.d0
      kn_th_tmp = 0.710d0
      write(10,'(2ES16.8)') kn_min, hl_tmp
      write(10,'(2ES16.8)') kn_th_tmp, hl_tmp

      write(10,*)
      write(10,*)

      write(10,*) 'ZONE T = "k_th"'
      write(10,'(2ES16.8)') kn_th_tmp, 0.d0
      write(10,'(2ES16.8)') kn_th_tmp, -6.d0

      write(10,*)
      write(10,*)
      !currptr => x
      !write(10,*) 'ZONE T = "ES (num:converted to HS)"'
      !do while (associated(currptr))
      !  write(10,'(2ES16.8)') currptr%kn/gamma1_hs/(1.d0-nu), currptr%force_lift
      !  currptr => currptr%next
      !enddo

    close(10)

    print*, 'File created : ', trim(name)

! ---- ! For data base
    
    currptr => x

    if (currptr%ID_es == 0) then
      !name = 'hl_vs_kn_db_bgk.dat'
      name = 'hl_vs_kn_db_bgk_accommo'//trim(adjustl(accommo_char))//'_n'//TRIM(ADJUSTL(ID_u_char))//'.dat'
    elseif (currptr%ID_es == 1) then
      !name = 'hl_vs_kn_db_es.dat'
      name = 'hl_vs_kn_db_es_accommo'//trim(adjustl(accommo_char))//'_n'//TRIM(ADJUSTL(ID_u_char))//'.dat'
    endif

    open(unit=10, status='unknown', file=name)
      do while (associated(currptr))
        write(10,'(2e24.16)') currptr%kn, currptr%force_lift
        currptr => currptr%next
      enddo
    close(10)

    print*, 'File created : ', trim(name)


!*****************************************************************************************
  end subroutine write_force_lift_ID
!*****************************************************************************************


