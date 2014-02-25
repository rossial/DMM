      DOUBLE PRECISION FUNCTION PPCHI2(P,V,IFAULT)
C
C ALGORITHM AS 91 APPL. STATIST, (1975) VOL.24, No.3 
C      
C TO EVALUATE THE PERCENTAGE POINTS OF THE CHI-SQUARED 
C PROABILITY DISTRIBUTION FUNCTION. 
C P MUST LIE IN THE RANGE 0.000002 TO 0.999998, V MUST BE POSITIVE, 
C
      INTEGER IFAULT,IF1
      DOUBLE PRECISION P,V
      DOUBLE PRECISION PPND16,GAMAIN,gammln
      DOUBLE PRECISION A,AA,G,C,CH,P1,P2,XX,E,T
      DATA E,AA/.5D-6, .6931471805D0/ 
C 
C AFTER DEFINING ACCURACY AND LN(2), TEST ARGUMENTS AND INITIALIZE
C
      PPCHI2 = -1.D0 
      IF (P.LT..000002D0) THEN 
       PPCHI2 = 0.D0
       IFAIL = 0
       RETURN
      ENDIF
      IFAULT = 1 
      IF (P.GT..999998D0) RETURN 
      IFAULT = 2 
      IF (V.LE.0.D0) RETURN 
      IFAULT = 0 
      G  = gammln(V/2.D0)
      XX = 0.5D0 * V 
      C  = XX - 1.D0
C 
C STARTING APPROXIMIATION FOR SMALL CHI-SQUARED 
C 
      IF (V.GE.-1.24D0*DLOG(P)) GOTO 1
       CH = (P * XX * DEXP(G + XX * AA)) ** (1.D0 / XX) 
       IF (CH - E) 6, 4, 4 
C 
C STARTING APPROXIMATION FOR V LESS THAN OR EQUAL TO 0.32 
C      
1      IF (V.GT.0.32D0) GOTO 3 
        CH = 0.4D0 
        A = DLOG(1.D0 - P) 
2       Q = CH 
        P1 = 1.D0 + CH * (4.67D0 + CH) 
        P2 = CH * (6.73D0 + CH * (6.66D0 + CH)) 
        T = -0.5D0 + (4.67D0 + 2.D0 * CH) / P1 - 
     *      (6.73D0 + CH * (13.32D0 + 3.D0 * CH)) / P2 
        CH = CH - (1.D0 - DEXP(A + G + 0.5D0 * CH + C * AA) * P2 /P1)/T
        IF (DABS(Q / CH - 1.D0) - 0.01D0) 4, 4, 2 
C 
C CALL TO ALGORITHM AS241
C 
3       X = PPND16(P, IF1) 
C 
C STARTING APPROXIMATION USING WILSON AND HIILFERTY ESTIMATE 
C        
        P1 = 0.222222D0 / V 
        CH = V * (X * DSQRT(P1) + 1.D0 - P1)**3 
C 
C STARTING APPROXIMATION FOR P TENDING TO i 
C
        IF (CH.GT.2.2D0 * V + 6.D0) CH = -2.D0 * (DLOG(1.D0 - P) - C 
     *    * DLOG(.5D0 * CH) + G)
C        
C CALL TO ALGORITHM AS 32 AND CALCULATION OF SEVEN TERM TAYLOR SERIES 
C 
4       Q = CH 
        P1 = 0.5D0 * CH 
        P2 = P - GAMAIN(P1, XX, G, IF1) 
        IF (IF1.EQ.0.D0) GOTO 5 
        IFAULT = 3 
        RETURN        
5       T = P2 * DEXP(XX * AA + G + P1 - C * DLOG(CH)) 
        B = T / CH 
        A = 0.5D0 * T - B * C 
        S1 = (210.D0+A*(140.D0+A*(105.D0+A*(84.D0+A*(70.D0+60.D0*A)))))
     *     / 420.D0 
        S2 = (420.D0+A*(735.D0+A*(966.D0+A*(1141.D0+1278.D0*A)))) 
     *     / 2520.D0 
        S3 = (210.D0 + A * (462.D0 + A * (707.D0 + 932.D0 * A))) 
     *     / 2520.D0 
        S4=(252.D0+A*(672.D0+1182.D0*A)+C*(294.D0+A*(889.D0+1740.D0*A)))
     *     / 5040.D0 
        S5 = (84.D0 + 264.D0 * A + C * (175.D0 + 606.D0 * A)) 
     *     / 2520.D0 
        S6 = (120.D0 + C * (346.D0 + 127.D0 * C)) 
     *     / 5040.D0 
        CH = CH + T*(1.D0+.5D0*T*S1-B*C  
     *     * (S1-B*(S2-B*(S3-B*(S4-B*(S5-B*S6)))))) 
        IF (DABS(Q / CH - 1.D0).GT.E) GOTO 4         
6       PPCHI2 = CH 
       RETURN 
      END

        
        
      DOUBLE PRECISION FUNCTION GAMAIN(X,P,G,IFAULT) 
      DOUBLE PRECISION X,P,G
      INTEGER IFAULT
C 
C ALGORITHM AS 32 J.R.STATIST.SOC. C. (1970) VOL.19 NO.3 
C 
C COMPUTES INCOMPLETE GAMMA RATIO FOR POSITIVE VALUES OF 
C ARGUMENTS X AND P. G MUST BE SUPPLIED AND SHOULD BE EOUAL TO 
C LN(GAMMA(P)). C IFAULT = 1 IF P.LE.O ELSE 2 IF X.LT.O ELSE 0. 
C USES SERIES EXPANSION IF P.GT.X OR X.LE.1, OTHERWISE A 
C CONTINUED FRACTION APPROXIMATION. 
C DIMENSION PN(6) 
      DOUBLE PRECISION gammln
      DOUBLE PRECISION ACU,OFLO,GIN,FACTOR,A,B,TERM,AN,RN,DIF
      DOUBLE PRECISION PN(6)
      INTEGER I
      
C 
C DEFINE ACCURACY AND INITIALIZE 
C 
      ACU    = 1.D-8 
      OFLO   = 1.0D30 
      GIN    = 0.D0 
      IFAULT = 0 
C 
C TEST FOR ADMISSIBILITY OF ARGUMENTS 
C 
      IF(P.LE.0.D0) IFAULT=1 
      IF(X.LT.0.D0) IFAULT=2 
      IF(IFAULT.GT.0.OR.X.EQ.0.D0) GO TO 50 
      FACTOR = DEXP(P*DLOG(X)-X-G) 
      IF(X.GT.1.D0.AND.X.GE.P) GO TO 30 
C 
C CALCULATION BY SERIES EXPANSION 
C 
      GIN  = 1.D0 
      TERM = 1.D0 
      RN   = P 
20    RN   = RN+1.D0 
      TERM = TERM*X/RN 
      GIN  = GIN+TERM 
      IF(TERM.GT.ACU) GO TO 20 
      GIN=GIN*FACTOR/P 
      GO TO 50        
C 
C CALCULATION BY CONTINUED FRACTION 
C 
30    A     = 1.D0-P 
      B     = A+X+1.D0 
      TERM  = 0.D0 
      PN(1) = 1.D0 
      PN(2) = X
      PN(3) = X+1.D0 
      PN(4) = X*B
      GIN   = PN(3)/PN(4)
 32   A     = A+1.D0 
      B     = B+2.D0 
      TERM  = TERM+1.D0 
      AN    = A*TERM 
      DO 33 I=1,2 
33    PN(I+4)=B*PN(I+2)-AN*PN(I) 
      IF(PN(6).EQ.0.D0) GO TO 35 
      RN  = PN(5)/PN(6) 
      DIF = DABS(GIN-RN) 
      IF(DIF.GT.ACU) GO TO 34 
      IF(DIF.LE.ACU*RN) GO TO 42 
34    GIN = RN 
35    DO 36 I=1,4 
36    PN(I) = PN(I+2) 
      IF(DABS(PN(5)).LT.OFLO) GO TO 32 
      DO 41 I=1,4 
41    PN(I) = PN(I)/OFLO 
      GO TO 32 
42    GIN = 1.D0-FACTOR*GIN 
      
50    GAMAIN = GIN 
      RETURN 
      END      