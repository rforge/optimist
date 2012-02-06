      SUBROUTINE MTC2(N,W,C,Z,X,JDN,JDL,JFO,JCK,XX,WR,PR,M,L)
C
C THIS SUBROUTINE SOLVES THE UNBOUNDED CHANGE-MAKING PROBLEM
C
C MINIMIZE  Z = X(1) + ... + X(N)
C
C SUBJECT TO:   W(1)*X(1) + ... + W(N)*X(N) = C ,
C               X(J) .GE. 0 AND INTEGER  FOR J=1,...,N.
C
C THE PROGRAM IS INCLUDED IN THE VOLUME
C   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS
C   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990
C AND IMPLEMENTS THE ENUMERATIVE ALGORITHM DESCRIBED IN
C SECTION  5.6 .
C
C THE INPUT PROBLEM MUST SATISFY THE CONDITIONS
C
C   1) 2 .LE. N .LE. JDN - 1 ;
C   2) W(J), C  POSITIVE INTEGERS;
C   3) MAX (W(J)) .LT. C .
C
C MTC2 CALLS  5  PROCEDURES: CHMTC2, COREC, MAXT, MTC1 AND SORTI.
C
C THE PROGRAM IS COMPLETELY SELF-CONTAINED AND COMMUNICATION TO IT IS
C ACHIEVED SOLELY THROUGH THE PARAMETER LIST OF MTC2.
C NO MACHINE-DEPENDENT CONSTANT IS USED.
C THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN
C AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE
C SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY).
C THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P.
C 9000/840.
C
C MTC2 NEEDS
C   5  ARRAYS ( W ,  X ,  XX ,  WR  AND  PR ) OF LENGTH AT LEAST
C               JDN ;
C   2  ARRAYS ( M  AND  L ) OF LENGTH AT LEAST  JDL .
C
C MEANING OF THE INPUT PARAMETERS:
C N     = NUMBER OF ITEM TYPES;
C W(J)  = WEIGHT OF EACH ITEM OF TYPE  J  (J=1,...,N);
C C     = CAPACITY;
C JDN   = DIMENSION OF ARRAYS  W ,  X ,  XX ,  WR  AND  PR ;
C JDL   = DIMENSION OF ARRAYS  M  AND  L ( SUGGESTED VALUE
C                   JDL = MAX (W(J)) - 1 ;
C         IF THE CORE MEMORY IS NOT ENOUGH,  JDL  SHOULD BE SET
C         TO THE LARGEST POSSIBLE VALUE);
C JFO   = 1 IF OPTIMAL SOLUTION IS REQUIRED,
C       = 0 IF APPROXIMATE SOLUTION IS REQUIRED (AT MOST 100000
C         BACKTRACKINGS ARE PERFORMED);
C JCK   = 1 IF CHECK ON THE INPUT DATA IS DESIRED,
C       = 0 OTHERWISE.
C
C MEANING OF THE OUTPUT PARAMETERS:
C Z    = VALUE OF THE SOLUTION FOUND IF  Z .GT. 0 ,
C      = NO FEASIBLE SOLUTION EXISTS IF Z .EQ. 0 ,
C      = ERROR IN THE INPUT DATA (WHEN JCK=1) IF Z .LT. 0 : CONDI-
C        TION  - Z  IS VIOLATED;
C X(J) = NUMBER OF ITEMS OF TYPE  J  IN THE SOLUTION FOUND.
C
C ARRAYS XX, M, L, WR AND PR ARE DUMMY.
C
C ALL THE PARAMETERS ARE INTEGER. ON RETURN OF MTC2 ALL THE INPUT
C PARAMETERS ARE UNCHANGED.
C
      INTEGER W(JDN),X(JDN),C,Z
      INTEGER XX(JDN),WR(JDN),PR(JDN)
      INTEGER M(JDL),L(JDL)
      INTEGER PRJ,CWS,S1,S2
      Z = C + 1
      IF ( JCK .EQ. 1 ) CALL CHMTC2(N,W,C,Z,JDN)
      IF ( Z .LT. 0 ) RETURN
      MAXBCK = - 1
      IF ( JFO .EQ. 0 ) MAXBCK = 100000
C LOWER BOUND COMPUTATION.
      CALL MAXT(N,W,I1,I2,I3,JDN)
      S1 = C/W(I1)
      S2 = (C - S1*W(I1))/W(I2)
      IP = S1 + S2
      CWS = C - S1*W(I1) - S2*W(I2)
      IF ( CWS .NE. 0 ) GO TO 20
      Z = IP
      DO 10 J=1,N
        X(J) = 0
   10 CONTINUE
      X(I1) = S1
      X(I2) = S2
      RETURN
   20 LB = IP + (CWS + W(I3) - 1)/W(I3)
      L1 = IP - 1 + (CWS + W(I1) + W(I2) - 1)/W(I2)
      IF ( L1 .LT. LB ) LB = L1
      IF ( N .LE. 500 ) GO TO 90
C DETERMINE AND SOLVE THE CORE PROBLEM.
      CALL COREC(N,W,I1,I2,I3,JDN,NC,PR)
      CALL SORTI(NC,W,PR,JDN)
      DO 30 J=1,NC
        PRJ = PR(J)
        WR(J) = W(PRJ)
   30 CONTINUE
      CALL MTC1(NC,WR,C,LB,Z,XX,JDN,JDL,MAXBCK,X,M,L)
      IF ( Z .GT. 0 ) GO TO 40
      Z = C + 1
      GO TO 90
   40 DO 50 J=1,N
        X(J) = 0
   50 CONTINUE
      DO 60 J=1,NC
        PRJ = PR(J)
        X(PRJ) = XX(J)
   60 CONTINUE
      IF ( Z .EQ. LB ) GO TO 80
      DO 70 J=1,N
        XX(J) = X(J)
   70 CONTINUE
      GO TO 90
C THE CORE PROBLEM SOLUTION IS OPTIMAL.
   80 RETURN
C SOLVE THE COMPLETE PROBLEM.
   90 DO 100 J=1,N
        PR(J) = J
  100 CONTINUE
      CALL SORTI(N,W,PR,JDN)
      DO 110 J=1,N
        PRJ = PR(J)
        WR(J) = W(PRJ)
  110 CONTINUE
      IZH = Z
      CALL MTC1(N,WR,C,LB,Z,XX,JDN,JDL,MAXBCK,X,M,L)
      IF ( Z .EQ. 0 ) RETURN
C STORE IN X THE FINAL SOLUTION.
      IF ( Z .EQ. IZH ) GO TO 130
      DO 120 J=1,N
        PRJ = PR(J)
        X(PRJ) = XX(J)
  120 CONTINUE
      RETURN
  130 DO 140 J=1,N
        X(J) = XX(J)
  140 CONTINUE
      RETURN
      END
      SUBROUTINE CHMTC2(N,W,C,Z,JDN)
C
C CHECK THE INPUT DATA.
C
      INTEGER W(JDN),C,Z
      IF ( N .GE. 2 .AND. N .LE. JDN - 1 ) GO TO 10
      Z = - 1
      RETURN
   10 IF ( C .GT. 0 ) GO TO 30
   20 Z = - 2
      RETURN
   30 DO 40 J=1,N
        IF ( W(J) .LE. 0 ) GO TO 20
        IF ( W(J) .GE. C ) GO TO 50
   40 CONTINUE
      RETURN
   50 Z = - 3
      RETURN
      END
      SUBROUTINE COREC(N,W,I1,I2,I3,JDN,NC,PR)
C
C DETERMINE THE CORE PROBLEM.
C
      INTEGER W(JDN),PR(JDN)
      INTEGER IV(3),IVS(3)
      NC = N/20
      IF ( NC .LT. 500 ) NC = 500
      IV(1) = I1
      PR(1) = 1
      IV(2) = I2
      PR(2) = 2
      IV(3) = I3
      PR(3) = 3
      CALL SORTI(3,IV,PR,3)
      IVS(1) = - IV(PR(3))
      IVS(2) = - IV(PR(2))
      IVS(3) = - IV(PR(1))
      JTEST = 1
      JVTEST = - IVS(1)
      LLI = 1
      DO 30 J=1,NC
   10   IF ( LLI .LT. JVTEST ) GO TO 20
          IF ( LLI .EQ. JVTEST ) IVS(JTEST) = - IVS(JTEST)
          JTEST = JTEST + 1
          JVTEST = N + 1
          IF ( JTEST .LE. 3 ) JVTEST = - IVS(JTEST)
        GO TO 10
   20   PR(J) = LLI
        IF ( J .LT. NC ) LLI = LLI + (N - LLI)/(NC - J)
   30 CONTINUE
      NT = NC + 1
      IF ( IVS(1) .GT. 0 ) GO TO 40
      NT = NT - 1
      IF ( IVS(2) .EQ. PR(NT) .OR. IVS(3) .EQ. PR(NT) ) NT = NT - 1
      IF ( IVS(2) .EQ. PR(NT) .OR. IVS(3) .EQ. PR(NT) ) NT = NT - 1
      PR(NT) = - IVS(1)
   40 IF ( IVS(2) .GT. 0 ) GO TO 50
      NT = NT - 1
      IF ( IVS(1) .EQ. PR(NT) .OR. IVS(3) .EQ. PR(NT) ) NT = NT - 1
      IF ( IVS(1) .EQ. PR(NT) .OR. IVS(3) .EQ. PR(NT) ) NT = NT - 1
      PR(NT) = - IVS(2)
   50 IF ( IVS(3) .GT. 0 ) RETURN
      NT = NT - 1
      IF ( IVS(1) .EQ. PR(NT) .OR. IVS(2) .EQ. PR(NT) ) NT = NT - 1
      IF ( IVS(1) .EQ. PR(NT) .OR. IVS(2) .EQ. PR(NT) ) NT = NT - 1
      PR(NT) = - IVS(3)
      RETURN
      END
      SUBROUTINE MAXT(N,W,I1,I2,I3,JDN)
C
C DETERMINE THE THREE ITEMS OF MAXIMUM WEIGHT.
C
      INTEGER W(JDN)
      MAX1 = -1
      MAX2 = -1
      MAX3 = -1
      I1 = 0
      I2 = 0
      DO 30 I=1,N
        IF ( W(I) .LE. MAX3 ) GO TO 30
        IF ( W(I) .GT. MAX1 ) GO TO 20
        IF ( W(I) .GT. MAX2 ) GO TO 10
        MAX3 = W(I)
        I3 = I
        GO TO 30
   10   MAX3 = MAX2
        MAX2 = W(I)
        I3 = I2
        I2 = I
        GO TO 30
   20   MAX3 = MAX2
        MAX2 = MAX1
        MAX1 = W(I)
        I3 = I2
        I2 = I1
        I1 = I
   30 CONTINUE
      RETURN
      END
      SUBROUTINE MTC1(N,W,C,LB,Z,X,JDN,JDL,MAXBCK,XX,M,L)
C
C THIS SUBROUTINE SOLVES A CHANGE-MAKING PROBLEM THROUGH  THE
C BRANCH-AND-BOUND ALGORITHM PRESENTED IN
C  S. MARTELLO, P. TOTH, "OPTIMAL AND CANONICAL SOLUTIONS OF THE
C  CHANGE-MAKING PROBLEM", EUROPEAN JOURNAL OF OPERATIONAL RESEARCH,
C  1980.
C
      INTEGER W(JDN),X(JDN),C,Z
      INTEGER XX(JDN),CWF,S,PROFIT
      INTEGER M(JDL),L(JDL)
      KBCK = 0
      CWF = C
      W(N+1) = 1
      JDOM = JDL
      IF ( JDOM .GE. W(1) ) JDOM = W(1) - 1
      IF ( JDOM .LT. W(N) ) GO TO 30
      N2 = N + 2
      DO 20 JJ=2,N
        J = N2 - JJ
        K1 = W(J)
        K2 = W(J-1) - 1
        IF ( K2 .GT. JDOM ) K2 = JDOM
        DO 10 K=K1,K2
          M(K) = J
          L(K) = 0
   10   CONTINUE
      IF ( K2 .EQ. JDOM ) GO TO 30
   20 CONTINUE
   30 K1 = W(N) - 1
      IF ( K1 .GT. JDOM ) K1 = JDOM
      IF ( K1 .EQ. 0 ) GO TO 50
      DO 40 K=1,K1
        L(K) = Z
   40 CONTINUE
   50 XX(1) = C/W(1)
      DO 60 J=2,N
        XX(J) = 0
   60 CONTINUE
      PROFIT = XX(1)
      C = C - XX(1)*W(1)
      II = 2
C
C STEP (2.A).
C
   70 IF ( C .LE. JDOM ) GO TO 90
      IF ( C .LT. W(N) ) GO TO 150
      IIOLD = II
   80 IF ( C .GE. W(II) ) GO TO 100
      II = II + 1
      GO TO 80
   90 IF ( L(C) .GE. Z - PROFIT ) GO TO 150
      IIOLD = II
      II = M(C)
C
C STEP 2.
C
  100 S = C/W(II)
      ICT = C - S*W(II)
      IF ( Z .LE. PROFIT + S + (ICT + W(II+1) - 1)/W(II+1) ) GO TO 110
      IF ( ICT .EQ. 0 ) GO TO 130
      IF ( II .LT. N ) GO TO 120
  110 II = IIOLD
      GO TO 150
C
C STEP 3.
C
  120 C = ICT
      PROFIT = PROFIT + S
      XX(II) = S
      II = II + 1
      GO TO 70
C
C STEP 4.
C
  130 Z = PROFIT + S
      DO 140 J=1,N
        X(J) = XX(J)
  140 CONTINUE
      X(II) = S
      IF ( Z .NE. LB ) GO TO 150
      C = CWF
      RETURN
C
C STEP 5.
C
  150 KBCK = KBCK + 1
      IF ( KBCK .EQ. MAXBCK ) GO TO 170
      IB = II - 1
      DO 160 J=1,IB
        IIMJ = II - J
        IF ( XX(IIMJ) .GT. 0 ) GO TO 180
  160 CONTINUE
  170 C = CWF
      IF ( Z .GT. C ) Z = 0
      RETURN
  180 KK = II - J
      IF ( C .GE. W(KK) ) GO TO 190
      IF ( C .GT. JDOM ) GO TO 190
      IF ( Z - PROFIT .GT. L(C) ) L(C) = Z - PROFIT
  190 C = C + W(KK)
      PROFIT = PROFIT - 1
      XX(KK) = XX(KK) - 1
      IF ( Z .GT. PROFIT + (C + W(KK+1) - 1)/W(KK+1) ) GO TO 200
      C = C + XX(KK)*W(KK)
      PROFIT = PROFIT - XX(KK)
      XX(KK) = 0
      II = KK + 1
      GO TO 150
  200 II = KK + 1
      IIOLD = II
      IF ( C .GT. JDOM ) GO TO 210
      IF ( L(C) .GE. Z - PROFIT ) GO TO 150
  210 IF ( C - W(KK) .GE. W(N) ) GO TO 100
      IH = KK
C
C STEP 6.
C
  220 IH = IH + 1
      IF ( Z .LE. PROFIT + (C + W(IH) - 1)/W(IH) ) GO TO 150
      IF ( IH .GT. N ) GO TO 150
      IF ( C - W(IH) .LT. W(N) ) GO TO 220
      II = IH
      IIOLD = II
      GO TO 100
      END
      SUBROUTINE SORTI(N,A,V,JDA)
C
C SORT THE INTEGER ARRAY A BY DECREASING VALUES (DERIVED FROM
C SUBROUTINE SORTZV OF THE C.E.R.N. LIBRARY).
C
C JDA           = LENGTH OF ARRAY A;
C N             = NUMBER OF ELEMENTS OF A TO BE SORTED;
C V(I) (INPUT)  = POINTER TO THE I-TH ELEMENT TO BE SORTED;
C V(I) (OUTPUT) = POINTER TO THE I-TH ELEMENT OF THE SORTED ARRAY.
C
C ON RETURN, ARRAY A IS UNCHANGED.
C
      INTEGER V(N),IU(20),IL(20)
      INTEGER A(JDA),T
      II = 1
      JJ = N
      IF ( N .LE. 1 ) RETURN
      M = 1
      I = II
      J = JJ
   10 IF ( I .GE. J ) GO TO 80
   20 K = I
      IJ = (J + I)/2
      IV = V(IJ)
      T = A(IV)
      KI = V(I)
      IF ( A(KI) .GE. T ) GO TO 30
      V(IJ) = KI
      V(I) = IV
      IV = V(IJ)
      T = A(IV)
   30 L = J
      KI = V(J)
      IF ( A(KI) .LE. T ) GO TO 50
      V(IJ) = KI
      V(J) = IV
      IV = V(IJ)
      T = A(IV)
      KI = V(I)
      IF ( A(KI) .GE. T ) GO TO 50
      V(IJ) = KI
      V(I) = IV
      IV = V(IJ)
      T = A(IV)
      GO TO 50
   40 V(L) = V(K)
      V(K) = IVT
   50 L = L - 1
      KI = V(L)
      IF ( A(KI) .LT. T ) GO TO 50
      IVT = KI
   60 K = K + 1
      KI = V(K)
      IF ( A(KI) .GT. T ) GO TO 60
      IF ( K .LE. L ) GO TO 40
      IF ( L - I .LE. J - K ) GO TO 70
      IL(M) = I
      IU(M) = L
      I = K
      M = M + 1
      GO TO 90
   70 IL(M) = K
      IU(M) = J
      J = L
      M = M + 1
      GO TO 90
   80 M = M - 1
      IF ( M .EQ. 0 ) RETURN
      I = IL(M)
      J = IU(M)
   90 IF ( J - I .GE. II ) GO TO 20
      IF ( I .EQ. II ) GO TO 10
      I = I - 1
  100 I = I + 1
      IF ( I .EQ. J ) GO TO 80
      IV = V(I+1)
      T = A(IV)
      KI = V(I)
      IF ( A(KI) .GE. T ) GO TO 100
      K = I
  110 V(K+1) = V(K)
      K = K - 1
      KI = V(K)
      IF ( T .GT. A(KI) ) GO TO 110
      V(K+1) = IV
      GO TO 100
      END
