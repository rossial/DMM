C ***********************************************************************
C     ALGORITHM AS7, APPLIED STATISTICS, VOL.17, 1968.
C
C     ARGUMENTS:-
C     A()     = INPUT, THE SYMMETRIC MATRIX TO BE INVERTED, STORED IN
C               LOWER TRIANGULAR FORM
C     N       = INPUT, ORDER OF THE MATRIX
C     C()     = OUTPUT, THE INVERSE OF A (A GENERALIZED INVERSE IF C IS
C               SINGULAR), ALSO STORED IN LOWER TRIANGULAR.
C               C AND A MAY OCCUPY THE SAME LOCATIONS.
C     W()     = WORKSPACE, DIMENSION AT LEAST N.
C     NULLTY  = OUTPUT, THE RANK DEFICIENCY OF A.
C     IFAULT  = OUTPUT, ERROR INDICATOR
C                     = 1 IF N < 1
C                     = 2 IF A IS NOT +VE SEMI-DEFINITE
C                     = 0 OTHERWISE
C     RMAX    = OUTPUT, APPROXIMATE BOUND ON THE ACCURACY OF THE DIAGONAL
C               ELEMENTS OF C.  E.G. IF RMAX = 1.E-04 THEN THE DIAGONAL
C               ELEMENTS OF C WILL BE ACCURATE TO ABOUT 4 DEC. DIGITS.
C
C     LATEST REVISION - 18 October 1985
C
C*************************************************************************
      SUBROUTINE SYMINV(A,N,C,W,NULLTY,IFAULT,RMAX)
      implicit doubleprecision (a-h, o-z)
	INTEGER N,NROW,IFAULT
      DIMENSION A(*),C(*),W(N)
      DATA ZERO/0.D0/, ONE/1.D0/
C
      NROW=N
      IFAULT=1
      IF(NROW.LE.0) GO TO 100
      IFAULT=0
C
C     CHOLESKY FACTORIZATION OF A, RESULT IN C
C
      CALL CHOLA(A,NROW,C,NULLTY,IFAULT,RMAX,W)
      IF(IFAULT.NE.0) GO TO 100
C
C     INVERT C & FORM THE PRODUCT (CINV)'*CINV, WHERE CINV IS THE INVERSE
C     OF C, ROW BY ROW STARTING WITH THE LAST ROW.
C     IROW = THE ROW NUMBER, NDIAG = LOCATION OF LAST ELEMENT IN THE ROW.
C
      NN=NROW*(NROW+1)/2
      IROW=NROW
      NDIAG=NN
 10   IF(C(NDIAG).EQ.ZERO) GO TO 60
      L=NDIAG
      DO 20 I=IROW,NROW
        W(I)=C(L)
        L=L+I
 20   CONTINUE
      ICOL=NROW
      JCOL=NN
      MDIAG=NN
 30   L=JCOL
      X=ZERO
      IF(ICOL.EQ.IROW) X=ONE/W(IROW)
      K=NROW
 40   IF(K.EQ.IROW) GO TO 50
      X=X-W(K)*C(L)
      K=K-1
      L=L-1
      IF(L.GT.MDIAG) L=L-K+1
      GO TO 40
 50   C(L)=X/W(IROW)
      IF(ICOL.EQ.IROW) GO TO 80
      MDIAG=MDIAG-ICOL
      ICOL=ICOL-1
      JCOL=JCOL-1
      GO TO 30
c
c     Special case, zero diagonal element.
c
 60   L=NDIAG
      DO 70 J=IROW,NROW
        C(L)=ZERO
        L=L+J
 70   CONTINUE
c
c      End of row.
c
 80   NDIAG=NDIAG-IROW
      IROW=IROW-1
      IF(IROW.NE.0) GO TO 10
 100  RETURN
      END

       SUBROUTINE CHOLA(A, N, U, NULLTY, IFAULT, RMAX, R)
C
C     ALGORITHM AS6, APPLIED STATISTICS, VOL.17, 1968, WITH
C     MODIFICATIONS BY A.J.MILLER
C
C     ARGUMENTS:-
C     A()     = INPUT, A +VE DEFINITE MATRIX STORED IN LOWER-TRIANGULAR
C               FORM.
C     N       = INPUT, THE ORDER OF A
C     U()     = OUTPUT, A LOWER TRIANGULAR MATRIX SUCH THAT U*U' = A.
C               A & U MAY OCCUPY THE SAME LOCATIONS.
C     NULLTY  = OUTPUT, THE RANK DEFICIENCY OF A.
C     IFAULT  = OUTPUT, ERROR INDICATOR
C                     = 1 IF N < 1
C                     = 2 IF A IS NOT +VE SEMI-DEFINITE
C                     = 0 OTHERWISE
C     RMAX    = OUTPUT, AN ESTIMATE OF THE RELATIVE ACCURACY OF THE
C               DIAGONAL ELEMENTS OF U.
C     R()     = OUTPUT, ARRAY CONTAINING BOUNDS ON THE RELATIVE ACCURACY
C               OF EACH DIAGONAL ELEMENT OF U.
C
C     LATEST REVISION - 18 October 1985
C
C*************************************************************************
C
      implicit double precision (a-h, o-z)
      DIMENSION A(*),U(*),R(N)
C
C     ETA SHOULD BE SET EQUAL TO THE SMALLEST +VE VALUE SUCH THAT
C     1.0 + ETA IS CALCULATED AS BEING GREATER THAN 1.0 IN THE ACCURACY
C     BEING USED.
C
      DATA ETA/1.D-16/, ZERO/0.D0/, FIVE/5.D0/
C
      IFAULT=1
      IF(N.LE.0) GO TO 100
      IFAULT=2
      NULLTY=0
      RMAX=ETA
      R(1)=ETA
      J=1
      K=0
C
C     FACTORIZE COLUMN BY COLUMN, ICOL = COLUMN NO.
C
      DO 80 ICOL=1,N
        L=0
C
C     IROW = ROW NUMBER WITHIN COLUMN ICOL
C
        DO 40 IROW=1,ICOL
          K=K+1
          W=A(K)
          IF(IROW.EQ.ICOL) RSQ=(W*ETA)**2
          M=J
          DO 10 I=1,IROW
            L=L+1
            IF(I.EQ.IROW) GO TO 20
            W=W-U(L)*U(M)
            IF(IROW.EQ.ICOL) RSQ=RSQ+(U(L)**2*R(I))**2
            M=M+1
 10       CONTINUE
 20       IF(IROW.EQ.ICOL) GO TO 50
          IF(U(L).EQ.ZERO) GO TO 30
          U(K)=W/U(L)
          GO TO 40
 30       U(K)=ZERO
          IF(ABS(W).GT.ABS(RMAX*A(K))) GO TO 100
 40     CONTINUE
C
C     END OF ROW, ESTIMATE RELATIVE ACCURACY OF DIAGONAL ELEMENT.
C
 50     RSQ=SQRT(RSQ)
        IF(ABS(W).LE.FIVE*RSQ) GO TO 60
        IF(W.LT.ZERO) GO TO 100
        U(K)=SQRT(W)
        R(I)=RSQ/W
        IF(R(I).GT.RMAX) RMAX=R(I)
        GO TO 70
 60     U(K)=ZERO
        NULLTY=NULLTY+1
 70     J=J+ICOL
 80   CONTINUE
      IFAULT=0
 100  RETURN
      END
