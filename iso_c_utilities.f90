!
! Copyright 2014 Aleksandar Donev
! Public Domain
! From http://cims.nyu.edu/~donev/Fortran/DLL/DLL.Forum.txt
!
MODULE ISO_C_UTILITIES
  USE ISO_C_BINDING
  CHARACTER(C_CHAR), DIMENSION(1), SAVE, TARGET, PRIVATE :: dummy_string="?"
CONTAINS
  FUNCTION C_F_STRING(CPTR) RESULT(FPTR)
    TYPE(C_PTR), INTENT(IN) :: CPTR
    CHARACTER(KIND=C_CHAR), DIMENSION(:), POINTER :: FPTR

    INTERFACE
       FUNCTION strlen(string) RESULT(len) BIND(C,NAME="strlen")
         USE ISO_C_BINDING
         TYPE(C_PTR), VALUE :: string
       END FUNCTION strlen
    END INTERFACE

    IF(C_ASSOCIATED(CPTR)) THEN
       CALL C_F_POINTER(FPTR=FPTR, CPTR=CPTR, SHAPE=[strlen(CPTR)])
    ELSE
       FPTR=>dummy_string
    END IF
  END FUNCTION C_F_STRING
END MODULE ISO_C_UTILITIES
