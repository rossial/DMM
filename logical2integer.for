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
      SUBROUTINE LOGICAL2INTEGER(INV, dim, OUTV)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: dim
      LOGICAL, DIMENSION(dim), INTENT(IN) :: INV
      INTEGER, DIMENSION(dim), INTENT(OUT) :: OUTV
      INTEGER :: I
      DO I=1,dim
         IF (INV(I)) THEN
            OUTV(I) = 1
         ELSE
            OUTV(I) = 0
         END IF
      END DO
      END SUBROUTINE LOGICAL2INTEGER
