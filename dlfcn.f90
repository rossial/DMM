!
! Copyright 2014 Aleksandar Donev
! Public Domain
! From http://cims.nyu.edu/~donev/Fortran/DLL/DLL.Forum.txt
!
MODULE DLFCN
  USE ISO_C_BINDING
  USE ISO_C_UTILITIES
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: DLOpen, DLSym, DLClose, DLError

  INTEGER(C_INT), PARAMETER, PUBLIC :: RTLD_LAZY=1, RTLD_NOW=2, RTLD_GLOBAL=8, RTLD_LOCAL=4

  INTERFACE
     FUNCTION DLOpen(file,mode) RESULT(handle) BIND(C,NAME="dlopen")
       USE ISO_C_BINDING
       CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: file
       INTEGER(C_INT), VALUE :: mode
       TYPE(C_PTR) :: handle
     END FUNCTION DLOpen
     FUNCTION DLSym(handle,name) RESULT(funptr) BIND(C,NAME="dlsym")
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE :: handle
       CHARACTER(C_CHAR), DIMENSION(*), INTENT(IN) :: name
       TYPE(C_FUNPTR) :: funptr
     END FUNCTION DLSym
     FUNCTION DLClose(handle) RESULT(status) BIND(C,NAME="dlclose")
       USE ISO_C_BINDING
       TYPE(C_PTR), VALUE :: handle
       INTEGER(C_INT) :: status
     END FUNCTION DLClose
     FUNCTION DLError() RESULT(error) BIND(C,NAME="dlerror")
       USE ISO_C_BINDING
       TYPE(C_PTR) :: error
     END FUNCTION DLError
  END INTERFACE
END MODULE DLFCN
