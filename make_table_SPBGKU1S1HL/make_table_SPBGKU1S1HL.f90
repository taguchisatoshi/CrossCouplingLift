!******************************************************************************************
    PROGRAM make_table
!******************************************************************************************
    USE type_list; USE constants
    implicit none
    type(ID_unit), pointer :: IDL

    interface
      subroutine get_list_ID(ID)
        use type_list
        type(ID_unit), pointer :: ID
      end subroutine get_list_ID
    end interface

    INTERFACE
      RECURSIVE SUBROUTINE merge_sort_ID(x)
        USE type_list
        TYPE(ID_unit), POINTER :: x
      END SUBROUTINE merge_sort_ID
    END INTERFACE

    INTERFACE
      SUBROUTINE printlist_ID(x)
        USE type_list
        TYPE(ID_unit), POINTER :: x
      END SUBROUTINE printlist_ID
    END INTERFACE

    INTERFACE
      SUBROUTINE write_list_ID(x)
        USE type_list
        TYPE(ID_unit), POINTER :: x
      END SUBROUTINE write_list_ID
    END INTERFACE

    INTERFACE
      SUBROUTINE write_force_lift_ID(x)
        USE type_list
        TYPE(ID_unit), POINTER :: x
      END SUBROUTINE write_force_lift_ID
    END INTERFACE

    !INTERFACE
    !  SUBROUTINE write_const_ID(x)
    !    USE type_list
    !    TYPE(ID_unit), POINTER :: x
    !  END SUBROUTINE write_const_ID
    !END INTERFACE

    !INTERFACE
    !  SUBROUTINE visualization_force_vs_v(x)
    !    USE type_list
    !    TYPE(ID_unit), POINTER :: x
    !  END SUBROUTINE visualization_force_vs_v
    !END INTERFACE

    !INTERFACE
    !  SUBROUTINE write_table_ID(x)
    !    USE type_list
    !    TYPE(ID_unit), POINTER :: x
    !  END SUBROUTINE write_table_ID
    !END INTERFACE
!******************************************************************************************
    pi = 4.d0*ATAN(1.d0)

    CALL get_list_ID(IDL)
    CALL merge_sort_ID(IDL)
    CALL printlist_ID(IDL)
    CALL write_list_ID(IDL)
    CALL write_force_lift_ID(IDL)


!******************************************************************************************
    END PROGRAM make_table
!******************************************************************************************


