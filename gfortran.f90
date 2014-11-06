!
! Copyright (C) 2014 Houtan Bastani
!
! This code is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This code is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
!

MODULE GFORTRAN
  IMPLICIT NONE
CONTAINS
  FUNCTION LOGICAL2INTEGER(V)
#if defined(OCTAVE_MEX_FILE)
    LOGICAL*4, DIMENSION(:), INTENT(IN) :: V
#else
    LOGICAL, DIMENSION(:), INTENT(IN) :: V
#endif
    INTEGER :: I
    INTEGER, DIMENSION(SIZE(V)) :: LOGICAL2INTEGER
    DO I=1,SIZE(V)
       IF (V(I)) THEN
          LOGICAL2INTEGER(I) = 1
       ELSE
          LOGICAL2INTEGER(I) = 0
       END IF
    END DO
  END FUNCTION LOGICAL2INTEGER
END MODULE GFORTRAN
