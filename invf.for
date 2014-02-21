C ----------------------------------------------------------------------
C INVF computes the inverse of the np x np matrix (A + k*B) when k->inf.
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C
C CARE!!  A and B must have the following structure:
C  A = [a11  a12
C	  a12' a22], where a22 is (np-nq) x (np-nq) of full rank
C
C  B = [b1 0
C	  0 0],      where b1 is  nq x nq of full rank
C
C  OUTPUT:  Am, Bm: inv(A+kB) = Am+(1/k)Bm-(1/k^2)Bm*A*Bm+O(1/k^3) 
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
	SUBROUTINE INVF(A,B,np,nq,Am,Bm,FFF)  
C INPUT
	INTEGER np,nq
	DOUBLE PRECISION A(np,np),B(np,np)
C OUTPUT
	DOUBLE PRECISION Am(np,np),Bm(np,np),FFF(np,np)
C LOCALS	
	INTEGER IFAIL,IPIV(np),I,J
	DOUBLE PRECISION Q1(np,nq),Q2(np,np-nq),V(np-nq,np-nq),
     1 VV(nq,nq),M(nq,np-nq),Q1A(nq,np-nq),Q2Ainv(np-nq,np-nq),
     3 AQ2inv(np-nq,np-nq),AQ2(np-nq,np-nq),BmA(np,np)
	DOUBLE PRECISION W(np),WORK(64*np)
	DOUBLE PRECISION ZERO
	DATA ZERO/0.0D0/ 
      EXTERNAL DSYEV,DGETRF,DGETRI
      
	Q1(:,:) = ZERO
	Q2(:,:) = ZERO
	V(:,:)  = A(nq+1:np,nq+1:np) ! V = A22 
      
C [V,W] = eig(A(nq+1:np,nq+1:np))  
	IFAIL = 0 
c	CALL F02FAF('V','U',np-nq,V,np-nq,W(1:np-nq),WORK,64*np,IFAIL)	
      CALL DSYEV('V','U',np-nq,V,np-nq,W(1:np-nq),WORK,64*np,IFAIL)	

C  Q2(nq+1:np,:) = V*W^-.5  		
      DO 10 J = 1,np-nq
10    Q2(nq+1:np,J) = V(:,J)/dsqrt(W(J))

c  [VV,W]=eig(B(1:nq,1:nq));
	VV(:,:) = B(1:nq,1:nq)
	IFAIL = 0
c	CALL F02FAF('V','U',nq,VV,nq,W(1:nq),WORK,64*np,IFAIL)	
      CALL DSYEV('V','U',nq,VV,nq,W(1:nq),WORK,64*np,IFAIL)	

C  Q1(1:nq,1:nq) = VV*W^-.5;
	DO 20 I = 1,nq
20	Q1(1:nq,I) = VV(1:nq,I)/dsqrt(W(I)) 

c M =-Q1(1:nq,:)'*A(1:nq,nq+1:np)*Q2(nq+1:np,:)*inv(A(nq+1:np,nq+1:np)*Q2(nq+1:np,:))
	DO 30 I = 1, np-nq
	DO 30 J = 1, np-nq	
30	AQ2(I,J) = SUM(A(I+nq,1+nq:np)*Q2(1+nq:np,J))

	IFAIL = 0
C	CALL F07ADF(np-nq,np-nq,AQ2,np-nq,IPIV(1:np-nq),IFAIL)
C	CALL F07AJF(np-nq,AQ2,np-nq,IPIV(1:np-nq),WORK,64*np,IFAIL)
	CALL DGETRF(np-nq,np-nq,AQ2,np-nq,IPIV(1:np-nq),IFAIL)
	CALL DGETRI(np-nq,AQ2,np-nq,IPIV(1:np-nq),WORK,64*np,IFAIL)
      
      AQ2inv(:,:) = AQ2(:,:)  

	DO 40 I=1,nq
	DO 40 J=1,np-nq	
40    Q1A(I,J) = SUM(Q1(:,I)*A(1:nq,J+nq))

C  Q2(nq+1:np,:)*inv(A(nq+1:np,nq+1:np)*Q2(nq+1:np,:))
	DO 50 I=1,np-nq
	DO 50 J=1,np-nq	
50    Q2Ainv(I,J) = SUM(Q2(I+nq,1:np-nq)*AQ2inv(1:np-nq,J))

	DO 60 I=1,nq
	DO 60 J=1,np-nq	
60    M(I,J) = -SUM(Q1A(I,1:np-nq)*Q2Ainv(1:np-nq,J))

c  Q1(nq+1:np,:) = M'
	DO 70 I = 1,nq
	DO 70 J = 1,np-nq  
70    Q1(nq+J,I) = M(I,J)		

C  Am = Q2*Q2'
      Am(:,:) = ZERO
	DO 80 I=nq+1,np
      Am(I,I) = SUM(Q2(I,1:np-nq)*Q2(I,1:np-nq))
	DO 80 J=nq+1,I-1
      Am(I,J) = SUM(Q2(I,1:np-nq)*Q2(J,1:np-nq))
80    Am(J,I) = Am(I,J)
      
C  Bm = Q1*Q1'
	DO 90 I=1,np
	DO 90 J=1,I
      Bm(I,J) = SUM(Q1(I,1:nq)*Q1(J,1:nq))
90    Bm(J,I) = Bm(I,J)

C  FFF = Bm*A*Bm
	DO 100 I=1,np
	DO 100 J=1,np	
100   BmA(I,J) = SUM(Bm(I,1:np)*A(1:np,J))

	DO 110 I=1,np
	DO 110 J=1,I	
      FFF(I,J) = SUM(BmA(I,1:np)*Bm(1:np,J))
110   FFF(J,I) = FFF(I,J)

	RETURN
	END