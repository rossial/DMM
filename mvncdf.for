C ----------------------------------------------------------------------
C MVNCDF Multivariate normal cumulative distribution function
C computes the multivariate Normal cumulative distribution 
C function with mean vector MU, variance matrix SIGMA inside 
C LB and UB. Algorithm due to Alan Genz (1992): 
C "Numerical Computation of Multivariate Normal Probabilities", 
C Journal of Computational and Graphical Statistics, pp. 141-149.
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
C ----------------------------------------------------------------------
	SUBROUTINE MVNCDF(LB,UB,MU,SIGMA,K,EPS,NMAX,INTSUM,ERROR,N)
	
! INPUT
	INTEGER K,NMAX 
	DOUBLE PRECISION LB(K),UB(K),MU(K),SIGMA(K,K),EPS	

! OUTPUT
	INTEGER N 
	DOUBLE PRECISION INTSUM,ERROR	

! LOCALS
	INTEGER IFAIL,I,J
	DOUBLE PRECISION LBC(K),UBC(K),C(K,K),F(K),E(K),D(K),Y(K)
      DOUBLE PRECISION ALPHA,VARSUM,W,DELTA,SUMCY      

! EXTERNAL SUBROUTINES
      EXTERNAL DPOTRF
! EXTERNAL FUNCTIONS      
      DOUBLE PRECISION CUMNORM,GENUNF,PPND16
      

	LBC(:) = LB(:) - MU(:)
	UBC(:) = UB(:) - MU(:)
	ALPHA =  2.326347874040841 ! 99-TH PERCENTILE FOR A N(0,1)
	C(1:K,1:K) = SIGMA(1:K,1:K)
	IFAIL = -1
      CALL DPOTRF('L',K,C,K,IFAIL) ! SIGMA = C*C'
      	
C INITIALIZATIONS
	F(:) = 0.D0
	D(:) = 0.D0
	E(:) = 0.D0
	Y(:) = 0.D0
	IFAIL = 0
      E(1) = CUMNORM(UBC(1)/C(1,1)) ! CDF
	IFAIL = 0
      D(1) = CUMNORM(LBC(1)/C(1,1)) ! CDF     
	F(1) = E(1) - D(1)
    	INTSUM = 0.D0
	N = 0
      VARSUM = 0.D0	
	ERROR = EPS+1.D0
	DO WHILE ((ERROR.GT.EPS).AND.(N.LT.NMAX))
	 DO 100 I = 2,K	
         W = GENUNF(0.D0,1.D0)
	   IF((D(I-1)+W*F(I-1)).LE.0.D0) THEN
	    Y(I-1) = -10.D10
	   ELSEIF ((D(I-1)+W*F(I-1)).GE.1.D0) THEN
	    Y(I-1) = 10.D10
	   ELSE
          Y(I-1) = PPND16(D(I-1)+W*F(I-1),IFAIL) 	
	   ENDIF
	   SUMCY = 0.D0
	   DO 50 J = 1, I-1
50	   SUMCY = SUMCY + C(I,J)*Y(J)
         E(I) =  CUMNORM((UBC(I) - SUMCY) / C(I,I))
	   IFAIL = 0
         D(I) =  CUMNORM((LBC(I) - SUMCY) / C(I,I))
100      F(I) = (E(I) - D(I))*F(I-1)
       N = N + 1
	 DELTA = (F(K) - INTSUM)/DFLOAT(N)
       INTSUM = INTSUM + DELTA
       VARSUM = (N-2)*VARSUM/DFLOAT(N) + DELTA**2
       ERROR = ALPHA * DSQRT(VARSUM)
      END DO  
	
	RETURN
	END