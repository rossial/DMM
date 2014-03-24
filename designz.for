C --------------------------------------------------------------------
C DESIGNZ sets transition probs for latent variables using INFOS
C Developed by A.Rossi, C.Planas and G.Fiorentini
C
C by cols: S1,S2,...,Snv; with nv <=6
C by row: the 1st contains the # of matrices affected by Si
C         the 2nd-3rd etc point to c (1),H (2),G (3),a (4),F (5),R (6)
C         the 8-th row contains the # of states
C         the 9-th row spec the dynamics for Sj
C
C nstot: total # of states i.e. ns1 x ns2 x ...x nsv
C
C OUTPUT: P1,P2,...,P6 where
C         Pk(i,j)   = Pr[Sk(t+1)=i|Sk(t)=j], k = 1,2,...,min(6,nv)
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
	SUBROUTINE DESIGNZ(nv,np,psi,INFOS,P1,P2,P3,P4,P5,P6)
C INPUT
	INTEGER nv,np,INFOS(9,6)
	DOUBLE PRECISION psi(np)
C OUTPUT
	DOUBLE PRECISION P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6))

C LOCALS
	INTEGER I,J,K

C Transition probability matrix
	K = 0
	IF (nv.EQ.1) THEN

	 IF (INFOS(9,1).EQ.1) THEN      ! S~IID
	  DO I = 1,INFOS(8,1)-1
	   P1(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,1)-1
	 ELSEIF (INFOS(9,1).EQ.2) THEN  ! S~Markov
	  DO J = 1,INFOS(8,1)
	   DO I = 1,INFOS(8,1)-1
	    P1(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,1)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,1)
	  P1(INFOS(8,1),J) = 1.D0-SUM(P1(1:INFOS(8,1)-1,J))
	 ENDDO

	ELSEIF (nv.EQ.2) THEN

	 IF (INFOS(9,1).EQ.1) THEN
	  DO I = 1,INFOS(8,1)-1
	   P1(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,1)-1
	 ELSEIF (INFOS(9,1).EQ.2) THEN
	  DO J = 1,INFOS(8,1)
	   DO I = 1,INFOS(8,1)-1
	    P1(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,1)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,1)
	  P1(INFOS(8,1),J) = 1.D0-SUM(P1(1:INFOS(8,1)-1,J))
	 ENDDO
	 IF (INFOS(9,2).EQ.1) THEN
	  DO I = 1,INFOS(8,2)-1
	   P2(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,2)-1
	 ELSEIF (INFOS(9,2).EQ.2) THEN
	  DO J = 1,INFOS(8,2)
	   DO I = 1,INFOS(8,2)-1
	    P2(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,2)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,2)
	  P2(INFOS(8,2),J) = 1.D0-SUM(P2(1:INFOS(8,2)-1,J))
	 ENDDO

	ELSEIF (nv.EQ.3) THEN

	 IF (INFOS(9,1).EQ.1) THEN
	  DO I = 1,INFOS(8,1)-1
	   P1(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,1)-1
	 ELSEIF (INFOS(9,1).EQ.2) THEN
	  DO J = 1,INFOS(8,1)
	   DO I = 1,INFOS(8,1)-1
	    P1(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,1)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,1)
	  P1(INFOS(8,1),J) = 1.D0-SUM(P1(1:INFOS(8,1)-1,J))
	 ENDDO
	 IF (INFOS(9,2).EQ.1) THEN
	  DO I = 1,INFOS(8,2)-1
	   P2(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,2)-1
	 ELSEIF (INFOS(9,2).EQ.2) THEN
	  DO J = 1,INFOS(8,2)
	   DO I = 1,INFOS(8,2)-1
	    P2(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,2)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,2)
	  P2(INFOS(8,2),J) = 1.D0-SUM(P2(1:INFOS(8,2)-1,J))
	 ENDDO
	 IF (INFOS(9,3).EQ.1) THEN
	  DO I = 1,INFOS(8,3)-1
	   P3(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,3)-1
	 ELSEIF (INFOS(9,3).EQ.2) THEN
	  DO J = 1,INFOS(8,3)
	   DO I = 1,INFOS(8,3)-1
	    P3(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,3)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,3)
	  P3(INFOS(8,3),J) = 1.D0-SUM(P3(1:INFOS(8,3)-1,J))
	 ENDDO

	ELSEIF (nv.EQ.4) THEN

	 IF (INFOS(9,1).EQ.1) THEN
	  DO I = 1,INFOS(8,1)-1
	   P1(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,1)-1
	 ELSEIF (INFOS(9,1).EQ.2) THEN
	  DO J = 1,INFOS(8,1)
	   DO I = 1,INFOS(8,1)-1
	    P1(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,1)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,1)
	  P1(INFOS(8,1),J) = 1.D0-SUM(P1(1:INFOS(8,1)-1,J))
	 ENDDO
	 IF (INFOS(9,2).EQ.1) THEN
	  DO I = 1,INFOS(8,2)-1
	   P2(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,2)-1
	 ELSEIF (INFOS(9,2).EQ.2) THEN
	  DO J = 1,INFOS(8,2)
	   DO I = 1,INFOS(8,2)-1
	    P2(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,2)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,2)
	  P2(INFOS(8,2),J) = 1.D0-SUM(P2(1:INFOS(8,2)-1,J))
	 ENDDO
	 IF (INFOS(9,3).EQ.1) THEN
	  DO I = 1,INFOS(8,3)-1
	   P3(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,3)-1
	 ELSEIF (INFOS(9,3).EQ.2) THEN
	  DO J = 1,INFOS(8,3)
	   DO I = 1,INFOS(8,3)-1
	    P3(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,3)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,3)
	  P3(INFOS(8,3),J) = 1.D0-SUM(P3(1:INFOS(8,3)-1,J))
	 ENDDO
	 IF (INFOS(9,4).EQ.1) THEN
	  DO I = 1,INFOS(8,4)-1
	   P4(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,4)-1
	 ELSEIF (INFOS(9,4).EQ.2) THEN
	  DO J = 1,INFOS(8,4)
	   DO I = 1,INFOS(8,4)-1
	    P4(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,4)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,4)
	  P4(INFOS(8,4),J) = 1.D0-SUM(P4(1:INFOS(8,4)-1,J))
	 ENDDO

	ELSEIF (nv.EQ.5) THEN

	 IF (INFOS(9,1).EQ.1) THEN
	  DO I = 1,INFOS(8,1)-1
	   P1(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,1)-1
	 ELSEIF (INFOS(9,1).EQ.2) THEN
	  DO J = 1,INFOS(8,1)
	   DO I = 1,INFOS(8,1)-1
	    P1(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,1)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,1)
	  P1(INFOS(8,1),J) = 1.D0-SUM(P1(1:INFOS(8,1)-1,J))
	 ENDDO
	 IF (INFOS(9,2).EQ.1) THEN
	  DO I = 1,INFOS(8,2)-1
	   P2(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,2)-1
	 ELSEIF (INFOS(9,2).EQ.2) THEN
	  DO J = 1,INFOS(8,2)
	   DO I = 1,INFOS(8,2)-1
	    P2(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,2)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,2)
	  P2(INFOS(8,2),J) = 1.D0-SUM(P2(1:INFOS(8,2)-1,J))
	 ENDDO
	 IF (INFOS(9,3).EQ.1) THEN
	  DO I = 1,INFOS(8,3)-1
	   P3(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,3)-1
	 ELSEIF (INFOS(9,3).EQ.2) THEN
	  DO J = 1,INFOS(8,3)
	   DO I = 1,INFOS(8,3)-1
	    P3(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,3)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,3)
	  P3(INFOS(8,3),J) = 1.D0-SUM(P3(1:INFOS(8,3)-1,J))
	 ENDDO
	 IF (INFOS(9,4).EQ.1) THEN
	  DO I = 1,INFOS(8,4)-1
	   P4(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,4)-1
	 ELSEIF (INFOS(9,4).EQ.2) THEN
	  DO J = 1,INFOS(8,4)
	   DO I = 1,INFOS(8,4)-1
	    P4(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,4)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,4)
	  P4(INFOS(8,4),J) = 1.D0-SUM(P4(1:INFOS(8,4)-1,J))
	 ENDDO
	 IF (INFOS(9,5).EQ.1) THEN
	  DO I = 1,INFOS(8,5)-1
	   P5(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,5)-1
	 ELSEIF (INFOS(9,5).EQ.2) THEN
	  DO J = 1,INFOS(8,5)
	   DO I = 1,INFOS(8,5)-1
	    P5(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,5)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,5)
	  P5(INFOS(8,5),J) = 1.D0-SUM(P5(1:INFOS(8,5)-1,J))
	 ENDDO

	ELSEIF (nv.EQ.6) THEN

	 IF (INFOS(9,1).EQ.1) THEN
	  DO I = 1,INFOS(8,1)-1
	   P1(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,1)-1
	 ELSEIF (INFOS(9,1).EQ.2) THEN
	  DO J = 1,INFOS(8,1)
	   DO I = 1,INFOS(8,1)-1
	    P1(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,1)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,1)
	  P1(INFOS(8,1),J) = 1.D0-SUM(P1(1:INFOS(8,1)-1,J))
	 ENDDO
	 IF (INFOS(9,2).EQ.1) THEN
	  DO I = 1,INFOS(8,2)-1
	   P2(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,2)-1
	 ELSEIF (INFOS(9,2).EQ.2) THEN
	  DO J = 1,INFOS(8,2)
	   DO I = 1,INFOS(8,2)-1
	    P2(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,2)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,2)
	  P2(INFOS(8,2),J) = 1.D0-SUM(P2(1:INFOS(8,2)-1,J))
	 ENDDO
	 IF (INFOS(9,3).EQ.1) THEN
	  DO I = 1,INFOS(8,3)-1
	   P3(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,3)-1
	 ELSEIF (INFOS(9,3).EQ.2) THEN
	  DO J = 1,INFOS(8,3)
	   DO I = 1,INFOS(8,3)-1
	    P3(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,3)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,3)
	  P3(INFOS(8,3),J) = 1.D0-SUM(P3(1:INFOS(8,3)-1,J))
	 ENDDO
	 IF (INFOS(9,4).EQ.1) THEN
	  DO I = 1,INFOS(8,4)-1
	   P4(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,4)-1
	 ELSEIF (INFOS(9,4).EQ.2) THEN
	  DO J = 1,INFOS(8,4)
	   DO I = 1,INFOS(8,4)-1
	    P4(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,4)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,4)
	  P4(INFOS(8,4),J) = 1.D0-SUM(P4(1:INFOS(8,4)-1,J))
	 ENDDO
	 IF (INFOS(9,5).EQ.1) THEN
	  DO I = 1,INFOS(8,5)-1
	   P5(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,5)-1
	 ELSEIF (INFOS(9,5).EQ.2) THEN
	  DO J = 1,INFOS(8,5)
	   DO I = 1,INFOS(8,5)-1
	    P5(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,5)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,5)
	  P5(INFOS(8,5),J) = 1.D0-SUM(P5(1:INFOS(8,5)-1,J))
	 ENDDO
	 IF (INFOS(9,6).EQ.1) THEN
	  DO I = 1,INFOS(8,6)-1
	   P6(I,:) = psi(K+I)
	  ENDDO
	  K = K + INFOS(8,6)-1
	 ELSEIF (INFOS(9,6).EQ.2) THEN
	  DO J = 1,INFOS(8,6)
	   DO I = 1,INFOS(8,6)-1
	    P6(I,J) = psi(K+I)
	   ENDDO
	   K = K + INFOS(8,6)-1
	  ENDDO
 	 ENDIF
	 DO J = 1,INFOS(8,6)
	  P6(INFOS(8,6),J) = 1.D0-SUM(P6(1:INFOS(8,6)-1,J))
	 ENDDO
	ENDIF

	RETURN
	END
