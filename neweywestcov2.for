C ------------------------------------------------------------------------
C NEWEYWESTCOV2 implements Newey and West 1987, A simple positive semi-definite 
C hetheroscedasticity and autocorrelation consistent covariance matrix,
C Econometrica, 55, 703-08.
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C      
C OUTPUT:
C OMEGA = OMEGA0 + SUM(is=1,nq) (1-is/(nq+1))*(Omegas+Omegas')  
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
C ------------------------------------------------------------------------	
	SUBROUTINE NEWEYWESTCOV2(N,nvar,NQ,MAT,MEAN,OMEGA)
! INPUT
	INTEGER N,nvar,NQ
	DOUBLE PRECISION MAT(N,nvar),MEAN(nvar)
! OUTPUT
	DOUBLE PRECISION OMEGA(nvar,nvar)
! LOCALS
	INTEGER I,J,K,is
      DOUBLE PRECISION OMEGAS(nvar,nvar)
	DOUBLE PRECISION, ALLOCATABLE::MAT0(:,:)
	DOUBLE PRECISION ZERO,ONE
	DATA ZERO/0.0D0/, ONE/1.0D0/

      ALLOCATE (MAT0(N,nvar))
	DO 10 I = 1,nvar 
10	MAT0(:,I) = MAT(:,I) - MEAN(I) ! Remove the mean

	OMEGA(:,:) = ZERO  ! lag-0 covariance matrix
	DO 30 I =1,nvar
	DO 30 J =1,I
	OMEGA(I,J) = sum(MAT0(:,I)*MAT0(:,J))
30    OMEGA(J,I) = OMEGA(I,J)
	OMEGA(:,:) = OMEGA(:,:)/dfloat(N)

	DO 100 is = 1, NQ
	 OMEGAS(:,:) = ZERO  ! lag-s covariance matrix
	 DO 50 I =1,nvar
	 DO 50 J =1,nvar
 	 DO 50 K =1,N-is
50	 OMEGAS(I,J) = OMEGAS(I,J)+MAT0(K+is,I)*MAT0(K,J)/dfloat(N-is)

	 DO 60 I =1,nvar
	 DO 60 J =1,I
	 OMEGA(I,J) = OMEGA(I,J) + 
     +              (ONE-is/dfloat(NQ+1))*(OMEGAS(J,I)+OMEGAS(I,J))
60     OMEGA(J,I) = OMEGA(I,J)

100   CONTINUE
	DEALLOCATE (MAT0)
	RETURN
	END