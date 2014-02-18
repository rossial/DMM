C ----------------------------------------------------------------------
C TNORMI generates a truncated normal random number 
C through the inversion method        
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C         
C Copyright (C) 2010-2014 European Commission 
C
C This file is part of Program DMM
C
C DMM is free software developed at the Joint Research Centre of the 
C European Commission: you can redistribute it and/or modify it under 
C the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C DMM is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with DMM.  If not, see <http://www.gnu.org/licenses/>.   
C
C In addition, as a special exception, the copyright holders give
C permission to link the code of portions of this program with the
C NAG Fortran library under certain conditions as described in each
C individual source file, and distribute linked combinations including
C the two.
C
C You must obey the GNU General Public License in all respects for all
C of the code used other than NAG Fortran library. If you modify file(s)
C with this exception, you may extend this exception to your
C version of the file(s), but you are not obligated to do so. If
C you do not wish to do so, delete this exception statement from
C your version. If you delete this exception statement from all
C source files in the program, then also delete it here.      
C ----------------------------------------------------------------------
	DOUBLE PRECISION FUNCTION TNORMI(PHIL,PHIP)
      INTEGER IFAIL
	DOUBLE PRECISION PHIL,PHIP,T,G01FAF

C	T=G05CAF(T) ! Sampling from U(0,1)
      T = ranf()  ! Sampling from U(0,1)
      T=PHIL+T*(PHIP-PHIL) ! Rescaling U(PHIL,PHIP)
      TNORMI=G01FAF('L',T,IFAIL) ! INVERSE of N(0,1)

      RETURN
      END