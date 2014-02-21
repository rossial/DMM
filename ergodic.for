C ---------------------------------------------------------------------
C ERGODIC solves the linear system: PE*(I-P') = 0 
C Developed by A.Rossi, C.Planas and G.Fiorentini           
C PE=(PE1,PE2); PE1=PE(1),...,PE(n-1) and PE2=1-sum(PE(1:n-1))
C       A   = I-P' = [a  b
C     n x n           c  d]  
C     where a = A(1:n-1,1:n-1); c = A(n,1:n-1) (1 x n-1)  
C     => PE1 = -c*(a-1*c)**(-1), 1 = (1,1,...1)' (n-1 x 1)
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
C ---------------------------------------------------------------------
	SUBROUTINE ERGODIC(n,P,PE)
C INPUT
	INTEGER n
	DOUBLE PRECISION P(n,n)
C OUTPUT
	DOUBLE PRECISION PE(n)
C LOCALS
      INTEGER I,IFAIL,LWORK,IPIV(n-1)
	DOUBLE PRECISION A(n,n),AC(n-1,n-1)
	DOUBLE PRECISION, ALLOCATABLE:: WORK(:)

	LWORK = 64*(n-1)
	ALLOCATE(WORK(LWORK))

C A = I - P'	
	DO 10 I = 1,n
	A(I,:) = -P(:,I)
10	A(I,I) = A(I,I) + 1.D0

C AC = a - 1*c
	DO 20 I =1,n-1
20	AC(I,1:n-1) = A(I,1:n-1) - A(n,1:n-1)

C inv of AC	
C	CALL F07ADF(n-1,n-1,AC,n-1,IPIV,IFAIL)
C	CALL F07AJF(n-1,AC,n-1,IPIV,WORK,LWORK,IFAIL)
      CALL DGETRF(n-1,n-1,AC,n-1,IPIV,IFAIL)
      CALL DGETRI(n-1,AC,n-1,IPIV,WORK,LWORK,IFAIL)

C PE = -c*(A-1*c)**(-1)
	DO 30 I=1,n-1
30	PE(I) = -SUM(A(n,1:n-1)*AC(:,I))
	PE(n) = 1.D0-SUM(PE(1:n-1))

	DEALLOCATE(WORK)

	RETURN
	END