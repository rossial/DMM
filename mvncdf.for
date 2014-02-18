C ----------------------------------------------------------------------
C   MVNCDF Multivariate normal cumulative distribution function
C   P = MVNCDF(LB,UB,MU,SIGMA) computes the multivariate normal cdf 
C   F(LB,UB) with mean vector MU and variance matrix SIGMA 
C
C   P = MVNCDF(LB, UB, MU,SIGMA,ERRMAX,CI,NMAX) uses additional control 
C   parameters. The difference between P and the true value of the
C   cdf is less than ERRMAX CI percent of the time. NMAX is the 
C   maximum number of iterations that the algorithm makes. By 
C   default, ERRMAX is 0.01, CI is 99, and NMAX is 300.
C
C   [P,ERR,N] = MVNCDF(...) also returns the estimated error and the
C   number of iterations made.
C
C   Algorithm from Alan Genz (1992) Numerical Computation of 
C   Multivariate Normal Probabilities, Journal of Computational and 
C   Graphical Statistics, pp. 141-149.
C   Copyright 2005 Alex Strashny (alex@strashny.org)
C   version 1, April 29, 2005
C   Recoded in Fortran by A.Rossi, C.Planas and G.Fiorentini 	
C      
C  In addition, as a special exception, the copyright holders give
C  permission to link the code of portions of this program with the
C  NAG Fortran library under certain conditions as described in each
C  individual source file, and distribute linked combinations including
C  the two.
C
C  You must obey the GNU General Public License in all respects for all
C  of the code used other than NAG Fortran library. If you modify file(s)
C  with this exception, you may extend this exception to your
C  version of the file(s), but you are not obligated to do so. If
C  you do not wish to do so, delete this exception statement from
C  your version. If you delete this exception statement from all
C  source files in the program, then also delete it here.      
C ----------------------------------------------------------------------
	SUBROUTINE mvncdf(LB,UB,MU,SIGMA,K,errMax,Nmax,P,err,N)
	
! INPUT
	INTEGER K,Nmax 
	DOUBLE PRECISION LB(K),UB(K),MU(K),SIGMA(K,K),errMax	

! OUTPUT
	INTEGER N 
	DOUBLE PRECISION P,err	

! LOCALS
	INTEGER IFAIL,I,J
	DOUBLE PRECISION LBC(K),UBC(K),alph,AP((K*(K+1))/2),
	1 C(K,K),varSum,F(K),E(K),D(K),Y(K),W,DEL,Q
	DOUBLE PRECISION S15ABF,genunf,G01FAF

	LBC(:) = LB(:) - MU(:)
	UBC(:) = UB(:) - MU(:)
	alph = 2.32634787D0  ! deviate at 99% of N(0,1)
	
      DO 10 I = 1,K
	DO 10 J = 1,I
10     AP(I+(2*K-J)*(J-1)/2) = SIGMA(I,J)
      
	IFAIL = 0
	CALL F07GDF('L',K,AP,IFAIL) ! SIGMA = C*C'
	C(:,:) = 0.D0
      DO 20 I=1,K
	DO 20 J=1,I
20    C(I,J) = AP(I+(2*K-J)*(J-1)/2)
	
C d is always zero
	F(:) = 0.D0
	D(:) = 0.D0
	E(:) = 0.D0
	Y(:) = 0.D0
	IFAIL = 0
	E(1) = S15ABF(UBC(1)/C(1,1), IFAIL) ! CDF
	IFAIL = 0
	D(1) = S15ABF(LBC(1)/C(1,1), IFAIL) ! CDF
	F(1) = E(1) - D(1)
	
	err = 2.D0*errMax
	P = 0.D0
	N = 0
	varSum = 0.D0 
	DO WHILE ((err.GT.errMax).AND.(N.LT.Nmax))
	 DO 100 I =2,K	
C	   W = G05CAF(W)
         W = genunf(0.D0,1.D0)
	   IF((D(I-1)+W*F(I-1)).LE.0.D0) THEN
	    Y(I-1) = -10.D10
	   ELSEIF ((D(I-1)+W*F(I-1)).GE.1.D0) THEN
	    Y(I-1) = 10.D10
	   ELSE
	    Y(I-1) = G01FAF('L', D(I-1)+W*F(I-1), IFAIL) 	
	   ENDIF
	   Q = 0.D0
	   DO 50 J = 1, I-1
50	   Q = Q + C(I,J)*Y(J)
         E(I) =  S15ABF((UBC(I) - Q) / C(I,I),IFAIL)
		IFAIL = 0
         D(I) =  S15ABF((LBC(I) - Q) / C(I,I),IFAIL)
100      F(I) = (E(I) - D(I))*F(I-1)
       N = N + 1
	 del = (F(K) - P) / dfloat(N)
       P = P + del
       varSum = (N-2)*varSum /dfloat(N) + del**2
       err = alph * dsqrt(varSum)
      END DO  
	
	RETURN
	END