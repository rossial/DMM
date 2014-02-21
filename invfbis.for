C ----------------------------------------------------------------------
C INVFBIS computes the inverse of the np x np matrix (A + k*B) when k->inf.
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C
C CARE!!  A and B must have the following form 
C   A = [a11  a12
C        a12' a22], where a22 is a (np-nq)x(np-nq) matrix of full rank
C
C   B psd with rank(B) = nq <= np   
C
C
C  OUTPUT:  Am, Bm: inv(A+kB) = Am + (1/k)Bm - (1/k^2)Bm*A*Bm + O(1/k^3) 
C           FFF: Bm*A*Bm
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
C -----------------------------------------------------------------------                                 
	SUBROUTINE INVFBIS(A,B,np,nq,Am,Bm,FFF)  
C INPUT
	INTEGER NP,NQ
	DOUBLE PRECISION A(np,np),B(np,np)
C OUTPUT	
	DOUBLE PRECISION Am(np,np),Bm(np,np),FFF(np,np)
C LOCALS	
	INTEGER IFAIL,i,j
	DOUBLE PRECISION ZERO
	
	DOUBLE PRECISION Q(np,np)
	DOUBLE PRECISION DM(np),PM(np,np),PMA(np,np),COM(np,np),
	1 WORK(3*np),W(np)
	
	EXTERNAL DSYEV
	DATA ZERO/0.0D0/ 		 
	
C Inverse of     A + k*M, A, M NxN, A psd, M pd
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
C ---------------------------------------------
C  [COM DM] = eig(M) 
	IFAIL = -1
	PM(:,:) = A(:,:)
C	CALL F02FAF('V','U',np,PM,np,DM,WORK,3*np,IFAIL)  ! COM = P
      CALL DSYEV('V','U',np,PM,np,DM,WORK,3*np,IFAIL)  ! COM = P

C  PM = PM * DM^-.5
      DO 10 J = 1,np
10	PM(:,J) = PM(:,J)/dsqrt(DM(J))

C  COM = (PM*DM^-.5)'*FI
	COM(:,:) = ZERO
	DO 20 I = 1,np
	DO 20 J = 1,np
20	COM(I,J) = SUM(PM(1:np,I)*B(1:np,J))

C  Q = (PM*DM^-.5)'*FI*(PM*DM^-.5)
	PMA(:,:) = ZERO
	DO 30 I = 1,np
	PMA(I,I) = SUM(COM(I,1:np)*PM(1:np,I))
	DO 30 J = 1,I-1
       PMA(I,J) = SUM(COM(I,1:np)*PM(1:np,J))
30     PMA(J,I)=PMA(I,J)

C [PMA,W] = eig(A2)    
	IFAIL = -1
C	CALL F02FAF('V','U',np,PMA,np,W,WORK,3*np,IFAIL)
      CALL DSYEV('V','U',np,PMA,np,W,WORK,3*np,IFAIL)
C Q = PM*DM^-.5*PMA      
      Q(:,:)=ZERO
      DO 40 I = 1,np
	DO 40 J = 1,np
40	Q(I,J) = SUM(PM(I,1:np)*PMA(1:np,J))

C  Am = Q1*Q1' Q1 p x p-q
	DO 110 I=1,np
      Am(I,I) = SUM(Q(I,1:np-nq)*Q(I,1:np-nq))    
	DO 110 J=1,I-1
      Am(I,J) = SUM(Q(I,1:np-nq)*Q(J,1:np-nq))
110   Am(J,I) = Am(I,J)      

C  Bm = Q2*Q2'/W(np-nq+1:np) Q2 p x q
	DO 120 I=1,np
	DO 120 J=1,np
120   Bm(I,J) = SUM(Q(I,np-nq+1:np)*Q(J,np-nq+1:np)/W(np-nq+1:np))

C FFF  is Q2 * Q2 / W^2
	DO 130 I=1,np
	DO 130 J=1,np
130   FFF(I,J) =SUM(Q(I,np-nq+1:np)*Q(J,np-nq+1:np)/W(np-nq+1:np)**2.D0)

	RETURN
	END