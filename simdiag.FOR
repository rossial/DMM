C -------------------------------------------------------------------------
C SIMDIAG returns the Inverse of   A + k*M, A, M NxN, A psd, M pd
C Write          M  = PM*DM*PM'
C                MA = (PM*DM^-.5)'*A*(PM*DM^-.5)
C and            MA = PMA * DMA * PMA'
C then           Q = PM*DM^-.5*PMA verifies:
C                Q'*M*Q = I and Q'*A*Q = DMA
C implying:      A + k*M = inv(Q)'*DMA*inv(Q) + k inv(Q)'*inv(Q)
C                        = inv(Q)'*(DMA + k I) inv(Q)
C so        inv(A + k*M) = Q * inv(DMA + k I) * Q'
C           
C CARE!!!! inv (k*A + (1-k)*M) = (1/k) * Q *inv(DMA + (1-k)/k)) Q'
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
C -------------------------------------------------------------------------
	SUBROUTINE SIMDIAG(A,M,N,DMA,Q)
! INPUT
	INTEGER N
	DOUBLE PRECISION A(N,N),M(N,N)
! OUTPUT
	DOUBLE PRECISION DMA(N),Q(N,N)
! LOCALS
      INTEGER I,J,K,IFAIL
	DOUBLE PRECISION DM(N),PM(N,N),PMA(N,N),COM(N,N),WORK(3*N)
	DOUBLE PRECISION ZERO
	DATA ZERO/0.0D0/
! EXTERMAL SUBROUTINES
      EXTERNAL DSYEV

C  [COM DM] = eig(M) 
	IFAIL = -1
	PM(:,:) = M(:,:)
C	CALL F02FAF('V','U',N,PM,N,DM,WORK,3*N,IFAIL)  ! COM = P
      CALL DSYEV('V','U',N,PM,N,DM,WORK,3*N,IFAIL)  ! COM = P

C  PM = PM * DM^-.5
      DO 10 I = 1,N
	 DO 10 J = 1,N
10	  PM(I,J) = PM(I,J)/dsqrt(DM(J))

C  COM = (PM*DM^-.5)'*A
	COM(:,:) = ZERO
	DO 20 I = 1,N
	 DO 20 J = 1,N
        DO 20 K = 1,N
20	   COM(I,J) = COM(I,J) + PM(K,I)*A(K,J)
C  Q = (PM*DM^-.5)'*A*(PM*DM^-.5)
	PMA(:,:) = ZERO
	DO 30 I = 1,N
	 DO 30 J = 1,I
        DO 31 K = 1,N
31	   PMA(I,J) = PMA(I,J) + COM(I,K)*PM(K,J)
30      PMA(J,I)=PMA(I,J)
C [PMA,DMA] = eig(A2)    
	IFAIL = -1
C	CALL F02FAF('V','U',N,PMA,N,DMA,WORK,3*N,IFAIL)  ! COM = P
      CALL DSYEV('V','U',N,PMA,N,DMA,WORK,3*N,IFAIL)  ! COM = P
      
C Q = PM*DM^-.5*PMA      
      Q(:,:)=ZERO
      DO 40 I = 1,N
	 DO 40 J = 1,N
        DO 40 K = 1,N
40	   Q(I,J) = Q(I,J) + PM(I,K)*PMA(K,J)
	RETURN
	END