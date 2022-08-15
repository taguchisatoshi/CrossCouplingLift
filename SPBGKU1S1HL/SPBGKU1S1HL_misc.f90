!*****************************************************************************************
      SUBROUTINE para_range(n1,n2,nprocs,irank,ibgn,iend)
!*****************************************************************************************
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n1, n2, nprocs, irank
      INTEGER, INTENT(OUT) :: ibgn, iend
      INTEGER :: iwork1, iwork2
!*****************************************************************************************
      iwork1 = (n2-n1+1)/nprocs
      iwork2 = MOD(n2-n1+1,nprocs)
      ibgn = irank*iwork1 + n1 + MIN(irank,iwork2)
      iend = ibgn + iwork1 - 1
      IF (iwork2 > irank) iend = iend + 1
!*****************************************************************************************
      END SUBROUTINE para_range
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE MYSUM(RIN,RINOUT,LEN,ITYPE)
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8) :: RIN(*), RINOUT(*)
    INTEGER :: LEN, ITYPE
    INTEGER :: i
!*****************************************************************************************
    DO i = 1, LEN
      RINOUT(i) = RIN(i) + RINOUT(i)
    ENDDO
!*****************************************************************************************
    END SUBROUTINE MYSUM
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE display_date_time
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER :: date_time(8)
    CHARACTER :: date*8
    CHARACTER :: char_time*10
!*****************************************************************************************
    CALL DATE_AND_TIME (DATE=date,TIME=char_time,values=date_time)
    WRITE(*,11) ' Time: ', date_time(5), ':', date_time(6), ':', date_time(7), &
                 &' (',date_time(2), '/', date_time(3), ')'

11  FORMAT(a7,i2,a1,i2,a1,i2,a2,i2,a1,i2,a1)
!*****************************************************************************************
    END SUBROUTINE display_date_time
!*****************************************************************************************






!*****************************************************************************************
    SUBROUTINE interpolation_1d(order,x,f,xp,res)
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: order
    REAL(8), INTENT(IN) :: x(0:order), f(0:order), xp
    REAL(8), INTENT(OUT) :: res
    INTEGER :: i, j
    REAL(8) :: slope(order), add
!*****************************************************************************************
    res = 0.d0
    CALL Newton_difference(x,f,order,slope)
    DO i = order, 1, -1  ! "i" is the order of polynomial
      add = slope(i)
      DO j = 0, i-1
        add = add * ( xp-x(j) )
      ENDDO
        res = res + add
    ENDDO

    res = res + f(0)
!*****************************************************************************************
    END SUBROUTINE interpolation_1d
!*****************************************************************************************




!*****************************************************************************************
    SUBROUTINE Newton_difference(x,f,n,slope)
!*****************************************************************************************
    IMPLICIT NONE
    INTEGER :: i,j
    INTEGER, INTENT(IN) :: n  !    n-th polynomial is to be used.
    REAL(8) :: M(0:n,0:n)
    REAL(8), INTENT(IN) :: x(0:n), f(0:n)
    REAL(8), INTENT(OUT) :: slope(n)
!*****************************************************************************************    
    M(:,0) = f(:)
    DO j = 1, n
      DO i = 0, n-j
        M(i,j) = ( M(i+1,j-1)-M(i,j-1) )/( x(i+j)-x(i) )
      ENDDO
      slope(j) = M(0,j)
    ENDDO
!*****************************************************************************************
    END SUBROUTINE Newton_difference
!*****************************************************************************************



!*****************************************************************************************
    SUBROUTINE aitken_pred(a,pred)
!*****************************************************************************************
    IMPLICIT NONE
    REAL(8), INTENT(IN) :: a(0:2)
    REAL(8), INTENT(OUT) :: pred
    REAL(8) :: tmp
!*****************************************************************************************
    tmp = a(2) - 2.d0*a(1) + a(0)
    IF (tmp /= 0.d0) THEN
      pred = a(2) - (a(2)-a(1))**2/tmp
    ELSE
      pred = a(2)
    ENDIF
    
!*****************************************************************************************
    END SUBROUTINE aitken_pred
!*****************************************************************************************
