C ---------------------------------------------------------------
C The Cholesky decomposition of the MxM real symmetric positive 
C semi-definite matrix A=LU 
C U is the transpose of L, is performed and stored in the lower
C triangle of the array L.
C A is retained so that the solution obtained can be subsequently
C improved. The procedure will fail if A, modified by the rounding 
C errors, is not positive semi-definite.
C
C A input matrix
C L lower triangular matrix of the Cholesky decomposition
C M actual dimension of A
C NULL nullity of A, i.e. # of "zeros" diagonal elements of L
C TINY used as tolerance
C IFAIL 1 if A not positive semi-definite, 0 otherwise
C Modified from Wilkinson & Reinsch (1971) p.21 and Healy AS6
C
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
C ---------------------------------------------------------------
      SUBROUTINE SCHOLLU(A,L,M,NULL,TINY,IFAIL)
C INPUT	
	INTEGER M,NULL
	DOUBLE PRECISION A(M,M),TINY,X
C OUTPUT
	DOUBLE PRECISION L(M,M)
	INTEGER IFAIL
C LOCALS 
	INTEGER I,J,K
      
      NULL=0
      DO 10 I=1,M
       X=A(I,I)
       DO 20 K=I-1,1,-1
   20   X=X-L(I,K)*L(I,K)
       IF (X.GE.TINY) THEN
        L(I,I)=DSQRT(X)
        DO 21 J=I+1,M
         X=A(I,J)
         DO 22 K=I-1,1,-1
   22    X=X-L(J,K)*L(I,K)
   21   L(J,I)=X/L(I,I)
       ELSE IF (X.GE.0.D0) THEN
        NULL=NULL+1
        L(I,I)=0.D0
        DO 23 J=I+1,M
   23   L(J,I)=0.D0
       ELSE
        IFAIL=1
        RETURN
       ENDIF
   10 CONTINUE

      IFAIL=0
      RETURN
      END