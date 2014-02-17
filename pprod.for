C -------------------------------------------------------------
C PPROD computes PALL = P1 x P2 x ...x Pnv where 
C Pk(i,j)= Pr[Sk(t+1)=i|Sk(t)=j], k = 1,2,...,min(6,nv)
C P(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv  
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
C -------------------------------------------------------------
	SUBROUTINE PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,P)
C INPUT
	INTEGER nv,nstot
	INTEGER INFOS(9,6)
	DOUBLE PRECISION P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6))
C OUTPUT
	DOUBLE PRECISION P(nstot,nstot)
C LOCALS
      INTEGER n1,n2,n3,n4,n5,n6,I,J
	DOUBLE PRECISION PC(nstot,nstot)
	
	n1 = INFOS(8,1)
	IF (nv.EQ.1) THEN
	 P(1:n1,1:n1) = P1(1:n1,1:n1)
	ELSEIF (nv.GT.1) THEN
	 n2 = INFOS(8,2)
	 DO 20 I =1,n1
	 DO 20 J =1,n1
20	 P(n2*(I-1)+1:n2*I,n2*(J-1)+1:n2*J)=P1(I,J)*P2(1:n2,1:n2)
	ELSEIF (nv.GT.2) THEN
	 PC(:,:) = P(:,:) 
	 n3 = INFOS(8,3)
	 DO 30 I =1,n1*n2
	 DO 30 J =1,n1*n2
30	 P(n3*(I-1)+1:n3*I,n3*(J-1)+1:n3*J)=PC(I,J)*P3(1:n3,1:n3)
      ELSEIF (nv.GT.3) THEN
 	 PC(:,:) = P(:,:) 
	 n4 = INFOS(8,4)
	 DO 40 I =1,n1*n2*n3
	 DO 40 J =1,n1*n2*n3
40	 P(n4*(I-1)+1:n4*I,n4*(J-1)+1:n4*J)=PC(I,J)*P4(1:n4,1:n4)
	ELSEIF (nv.GT.4) THEN
 	 PC(:,:) = P(:,:) 
	 n5 = INFOS(8,5)
	 DO 50 I =1,n1*n2*n3*n4
	 DO 50 J =1,n1*n2*n3*n4
50     P(n5*(I-1)+1:n5*I,n5*(J-1)+1:n5*J)=PC(I,J)*P5(1:n5,1:n5)
	ELSEIF (nv.GT.5) THEN
 	 PC(:,:) = P(:,:) 
	 n6 = INFOS(8,6)
	 DO 60 I =1,n1*n2*n3*n4*n5
	 DO 60 J =1,n1*n2*n3*n4*n5
60     P(n6*(I-1)+1:n6*I,n6*(J-1)+1:n6*J)=PC(I,J)*P6(1:n6,1:n6)
	ENDIF

	RETURN
	END