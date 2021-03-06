      SUBROUTINE MTCB(N,W,B,C,Z,X,JDN,JDL,JFO,JCK,XX,WR,BR,PR,M,L)
C THIS SUBROUTINE SOLVES THE BOUNDED CHANGE-MAKING PROBLEM
C
C MINIMIZE  Z = X(1) + ... + X(N)
C
C SUBJECT TO:   W(1)*X(1) + ... + W(N)*X(N) = C ,
C               0 .LE. X(J) .LE. 0 FOR J=1,...,N,
C               X(J) INTEGER       FOR J=1,...,N.
C
C THE PROGRAM IS INCLUDED IN THE VOLUME
C   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS
C   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990
C AND IMPLEMENTS THE BRANCH-AND-BOUND ALGORITHM DESCRIBED
C IN SECTION  5.8 .
C
C THE INPUT PROBLEM MUST SATISFY THE CONDITIONS
C
C   1) 2 .LE. N .LE. JDN - 1 ;
C   2) W(J), B(J), C  POSITIVE INTEGERS;
C   3) MAX (W(J)) .LT. C ;
C   4) B(J)*W(J) .LE. C FOR J=1,...,N;
C   5) B(1)*W(1) + ...+ B(N)*W(N) .GT. C .
C
C MTCB CALLS  3  PROCEDURES: CHMTCB, CMPB AND SORTI.
C
C THE PROGRAM IS COMPLETELY SELF-CONTAINED AND COMMUNICATION TO IT IS
C ACHIEVED SOLELY THROUGH THE PARAMETER LIST OF MTCB.
C NO MACHINE-DEPENDENT CONSTANT IS USED.
C THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN
C AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE
C SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY).
C THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P.
C 9000/840.
C
C MTCB NEEDS
C   7  ARRAYS ( W ,  B ,  X ,  XX ,  WR ,  BR  AND  PR ) OF LENGTH
C               AT LEAST  JDN ;
C   2  ARRAYS ( M  AND  L ) OF LENGTH AT LEAST  JDL .
C
C MEANING OF THE INPUT PARAMETERS:
C N     = NUMBER OF ITEM TYPES;
C W(J)  = WEIGHT OF EACH ITEM OF TYPE  J  (J=1,...,N);
C B(J)  = NUMBER OF AVAILABLE ITEMS OF TYPE J  (J=1,...,N);
C C     = CAPACITY;
C JDN   = DIMENSION OF ARRAYS  W ,  B ,  X ,  XX ,  WR ,  BR  AND  PR ;
C JDL   = DIMENSION OF ARRAYS  M  AND  L (SUGGESTED VALUE
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
C ARRAYS XX, M, L, WR, BR AND PR ARE DUMMY.
C
C ALL THE PARAMETERS ARE INTEGER. ON RETURN OF MTCB ALL THE INPUT
C PARAMETERS ARE UNCHANGED.
C
      INTEGER W(JDN),B(JDN),X(JDN),C,Z
      INTEGER XX(JDN),WR(JDN),BR(JDN),PR(JDN)
      INTEGER M(JDL),L(JDL)
      Z = C + 1
      IF ( JCK .EQ. 1 ) CALL CHMTCB(N,W,B,C,Z,JDN)
      IF ( Z .LT. 0 ) RETURN
      MAXBCK = - 1
      IF ( JFO .EQ. 0 ) MAXBCK = 100000
C SORTING.
      DO 10 J=1,N
        PR(J) = J
   10 CONTINUE
      CALL SORTI(N,W,PR,JDN)
      DO 20 J=1,N
        JPR = PR(J)
        WR(J) = W(JPR)
        BR(J) = B(JPR)
   20 CONTINUE
C SOLUTION.
      CALL CMPB(N,WR,BR,C,Z,XX,JDN,JDL,MAXBCK,X,M,L)
C STORE IN X THE FINAL SOLUTION.
      DO 30 J=1,N
        JPR  = PR(J)
        X(JPR) = XX(J)
   30 CONTINUE
      RETURN
      END
      SUBROUTINE CHMTCB(N,W,B,C,Z,JDN)
C
C CHECK THE INPUT DATA.
C
      INTEGER W(JDN),B(JDN),C,Z
      IF ( N .GE. 2 .AND. N .LE. JDN - 1 ) GO TO 10
      Z = - 1
      RETURN
   10 IF ( C .GT. 0 ) GO TO 30
   20 Z = - 2
      RETURN
   30 JSUM = 0
      DO 40 J=1,N
        IF ( W(J) .LE. 0 ) GO TO 20
        IF ( B(J) .LE. 0 ) GO TO 20
        IF ( W(J) .GE. C ) GO TO 50
        IF ( B(J)*W(J) .GT. C ) GO TO 60
        JSUM = JSUM + B(J)*W(J)
   40 CONTINUE
      IF ( JSUM .LE. C ) GO TO 70
      RETURN
   50 Z = - 3
      RETURN
   60 Z = - 4
      RETURN
   70 Z = - 5
      RETURN
      END
      SUBROUTINE CMPB(N,W,B,C,Z,X,JDN,JDL,MAXBCK,XX,M,L)
C
C THIS SUBROUTINE SOLVES A BOUNDED CHANGE-MAKING PROBLEM THROUGH THE
C BRANCH-AND-BOUND ALGORITHM PRESENTED IN
C  S. MARTELLO, P. TOTH, "SOLUTION OF THE BOUNDED AND UNBOUNDED CHANGE-
C  MAKING PROBLEM", TIMS/ORSA JOINT NATIONAL MEETING, SAN FRANCISCO,
C  1977.
C
      INTEGER W(JDN),B(JDN),X(JDN),C,Z
      INTEGER XX(JDN),CWF,PROFIT,CWS
      INTEGER M(JDL),L(JDL)
C
C STEP 1.
C
      KBCK = 0
      CWF = C
      W(N+1) = 1
      B(N+1) = C + 1
C LOWER BOUND COMPUTATION.
      CWS = C
      JSB = 0
      DO 10 J=1,N
        IF ( B(J)*W(J) .GT. CWS )GO TO 20
        CWS = CWS - B(J)*W(J)
        JSB = JSB + B(J)
   10 CONTINUE
   20 JS = J
      JZP = JSB + CWS/W(JS)
      JCP = CWS - (CWS/W(JS))*W(JS)
      LB = JZP + ( JCP + W(JS+1) - 1)/W(JS+1)
      LB1 = JZP - 1 + ( JCP + W(JS-1) + W(JS) - 1 )/W(JS)
      IF ( LB1 .LT. LB ) LB = LB1
      IF ( JCP .GT. 0 ) GO TO 50
      Z = JZP
      DO 30 J=1,JS
        X(J) = B(J)
   30 CONTINUE
      DO 40 J =JS,N
        X(J) = 0
   40 CONTINUE
      X(JS) = CWS/W(JS)
      RETURN
   50 JDOM = JDL
      IF ( JDOM .GE. W(1) ) JDOM = W(1) - 1
      IF ( JDOM .LT. W(N) ) GO TO 80
      N2 = N + 2
      DO 70 JJ=2,N
        J = N2 - JJ
        K1 = W(J)
        K2 = W(J-1) - 1
        IF ( K2 .GT. JDOM ) K2 = JDOM
        DO 60 K=K1,K2
          M(K) = J
          L(K) = 0
   60   CONTINUE
      IF ( K2 .EQ. JDOM ) GO TO 80
   70 CONTINUE
   80 K1 = W(N) - 1
      IF ( K1 .GT. JDOM ) K1 = JDOM
      IF ( K1 .EQ. 0 ) GO TO 100
      DO 90 K=1,K1
        L(K) = Z
   90 CONTINUE
  100 XX(1) = C/W(1)
      IF ( B(1) .LT. XX(1) ) XX(1) = B(1)
      DO 110 J=2,N
        XX(J) = 0
  110 CONTINUE
      PROFIT = XX(1)
      C = C - XX(1)*W(1)
      II = 2
      GO TO 150
C
C STEP (2.A).
C
  120 IF ( C .LE. JDOM ) GO TO 140
      IF ( C .LT. W(N) ) GO TO 230
      IIOLD = II
  130 IF ( C .GE. W(II) ) GO TO 150
      II = II + 1
      GO TO 130
  140 IF ( L(C) .GE. Z - PROFIT ) GO TO 230
      IIOLD = II
      II = M(C)
C
C STEP 2.
C
  150 JYP = 0
      JCT = C
      DO 160 I =II,N
        JY = JCT/W(I)
        IF ( JY .GT. B(I) ) JY = B(I)
        JYP = JYP + JY
        JCT = JCT - JY*W(I)
        JW = W(I+1)
        IF ( JY .LT. B(I) ) JW = W(I)
        IF ( Z .LE. PROFIT + JYP + (JCT + JW - 1)/JW ) GO TO 230
        IF ( JCT .EQ. 0 ) GO TO 200
        IF ( JY .LT. B(I) ) GO TO 170
  160 CONTINUE
      GO TO 230
  170 PROFIT = PROFIT + ( JYP - JY )
      I1 = I - 1
      IF ( I1 .LT.II ) GO TO 190
      DO 180 K=II,I1
        XX(K) = B(K)
  180 CONTINUE
  190 II = I
C
C STEP 3.
C
      C = JCT
      PROFIT = PROFIT + JY
      XX(II) = JY
      II = II + 1
      GO TO 120
C
C STEP 4.
C
  200 Z = PROFIT + JYP
      DO 210 J=1,N
        X(J) = XX(J)
  210 CONTINUE
      DO 220 J=II,I
        X(J) = B(J)
  220 CONTINUE
      X(I) = JY
      IF ( Z .NE. LB ) GO TO 230
      C = CWF
      RETURN
C
C STEP 5.
C
  230 KBCK = KBCK + 1
      IF ( KBCK .EQ. MAXBCK ) GO TO 250
      IB = II - 1
      DO 240 J=1,IB
        IIJ = II - J
        IF ( XX(IIJ) .GT. 0 ) GO TO 260
  240 CONTINUE
  250 C = CWF
      IF ( Z .GT. C ) Z = 0
      RETURN
  260 KK = II - J
      IF ( C .GE. W(KK) ) GO TO 270
      IF ( C .GT. JDOM ) GO TO 270
      IF ( Z - PROFIT .GT. L(C) ) L(C) = Z - PROFIT
  270 C = C + W(KK)
      PROFIT = PROFIT - 1
      XX(KK) = XX(KK) - 1
      IF ( Z .GT. PROFIT + (C + W(KK+1) - 1)/W(KK+1) ) GO TO 280
      C = C + XX(KK)*W(KK)
      PROFIT = PROFIT - XX(KK)
      XX(KK) = 0
      II = KK + 1
      GO TO 230
  280 II = KK + 1
      IIOLD = II
      IF ( C .GT. JDOM ) GO TO 290
      IF ( L(C) .GE. Z - PROFIT ) GO TO 230
  290 IF ( C - W(KK) .GE. W(N) ) GO TO 150
      IH = KK
C
C STEP 6.
C
  300 IH = IH + 1
      IF ( Z .LE. PROFIT + (C + W(IH) - 1)/W(IH) ) GO TO 230
      IF ( IH .GT. N ) GO TO 230
      IF ( C - W(IH) .LT. W(N) ) GO TO 300
      II = IH
      IIOLD = II
      GO TO 150
      END
