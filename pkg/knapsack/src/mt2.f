      SUBROUTINE MT2(N,P,W,C,Z,X,JDIM,JFO,JFS,JCK,JUB,
     1               IA1,IA2,IA3,IA4,RA)
C
C THIS SUBROUTINE SOLVES THE 0-1 SINGLE KNAPSACK PROBLEM
C
C MAXIMIZE  Z = P(1)*X(1) + ... + P(N)*X(N)
C
C SUBJECT TO:   W(1)*X(1) + ... + W(N)*X(N) .LE. C ,
C               X(J) = 0 OR 1  FOR J=1,...,N.
C
C THE PROGRAM IS INCLUDED IN THE VOLUME
C   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS
C   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990
C AND IMPLEMENTS THE ENUMERATIVE ALGORITHM DESCRIBED IN
C SECTION  2.9.3 .
C
C THE INPUT PROBLEM MUST SATISFY THE CONDITIONS
C
C   1) 2 .LE. N .LE. JDIM - 3 ;
C   2) P(J), W(J), C  POSITIVE INTEGERS;
C   3) MAX (W(J)) .LE. C ;
C   4) W(1) + ... + W(N) .GT. C ;
C
C AND, IF  JFS = 1 ,
C
C   5) P(J)/W(J) .GE. P(J+1)/W(J+1) FOR J=1,...,N-1.
C
C MT2 CALLS  9  PROCEDURES: CHMT2, CORE, CORES, FMED, KP01M, NEWB,
C                           REDNS, REDS AND SORTR.
C
C THE PROGRAM IS COMPLETELY SELF-CONTAINED AND COMMUNICATION TO IT IS
C ACHIEVED SOLELY THROUGH THE PARAMETER LIST OF MT2.
C NO MACHINE-DEPENDENT CONSTANT IS USED.
C THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN
C AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE
C SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY).
C THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P.
C 9000/840.
C
C MT2 NEEDS  8  ARRAYS ( P ,  W ,  X ,  IA1 ,  IA2 ,  IA3 ,  IA4  AND
C                        RA ) OF LENGTH AT LEAST  N + 3 .
C
C MEANING OF THE INPUT PARAMETERS:
C N    = NUMBER OF ITEMS;
C P(J) = PROFIT OF ITEM  J  (J=1,...,N);
C W(J) = WEIGHT OF ITEM  J  (J=1,...,N);
C C    = CAPACITY OF THE KNAPSACK;
C JDIM = DIMENSION OF THE 8 ARRAYS;
C JFO  = 1 IF OPTIMAL SOLUTION IS REQUIRED,
C      = 0 IF APPROXIMATE SOLUTION IS REQUIRED;
C JFS  = 1 IF THE ITEMS ARE ALREADY SORTED ACCORDING
C          TO DECREASING PROFIT PER UNIT WEIGHT,
C      = 0 OTHERWISE;
C JCK  = 1 IF CHECK ON THE INPUT DATA IS DESIRED,
C      = 0 OTHERWISE.
C
C MEANING OF THE OUTPUT PARAMETERS:
C Z    = VALUE OF THE SOLUTION FOUND IF  Z .GT. 0 ,
C      = ERROR IN THE INPUT DATA (WHEN JCK=1) IF Z .LT. 0 : CONDI-
C        TION  - Z  IS VIOLATED;
C X(J) = 1 IF ITEM  J  IS IN THE SOLUTION FOUND,
C      = 0 OTHERWISE;
C JUB  = UPPER BOUND ON THE OPTIMAL SOLUTION VALUE (TO EVALUATE  Z
C        WHEN JFO=0).
C
C ARRAYS IA1, IA2, IA3, IA4 AND RA ARE DUMMY.
C
C ALL THE PARAMETERS BUT RA ARE INTEGER. ON RETURN OF MT2 ALL THE
C INPUT PARAMETERS ARE UNCHANGED.
C
      INTEGER P(JDIM),W(JDIM),X(JDIM),C,Z
      INTEGER IA1(JDIM),IA2(JDIM),IA3(JDIM),IA4(JDIM)
      REAL    RA(JDIM)
      INTEGER FN1,FN0
      Z = 0
      JUB = 0
      IF ( JCK .EQ. 1 ) CALL CHMT2(N,P,W,C,JFS,Z,JDIM)
      IF ( Z .LT. 0 ) RETURN
C
C STEP 1 (INITIALIZATION AND DEFINITION OF THE CORE PROBLEM).
C
      JDS = JDIM/3
      NP1 = N + 1
C ON RETURN OF CORE OR CORES, ASSUMING THAT N1 (N0) IS THE SET OF ITEMS
C TEMPORARILY INCLUDED IN (EXCLUDED FROM) THE SOLUTION, WE HAVE:
C IZ1 = TOTAL PROFIT OF ITEMS IN SET N1;
C ICW = C - TOTAL WEIGHT OF ITEMS IN SET N1;
C MINW0 = MINIMUM WEIGHT OF ITEMS IN SET N0;
C IA2(1) TO IA2(NNF) = ITEMS IN THE CORE PROBLEM (FREE ITEMS);
C IA1(I) = SUCCESSOR OF ITEM I IN SETS N0 AND N1 (= N + 1 IF LAST);
C FN1 = POINTER TO THE FIRST ITEM IN SET N1 (= N + 1 IF EMPTY);
C FN0 = POINTER TO THE FIRST ITEM IN SET N0 (= N + 1 IF EMPTY);
C RA(I) = P(I)/W(I).
      IF ( JFS .EQ. 0 ) GO TO 10
      CALL CORES(N,P,W,C,JFO,IZ1,ICW,MINW0,JDIM,IA2,NNF,IA1,FN1,FN0,RA)
      GO TO 20
   10 CALL CORE(N,P,W,C,JFO,IZ1,ICW,MINW0,JDIM,IA2,NNF,IA1,FN1,FN0,RA)
   20 IF ( NNF .EQ. N ) GO TO 130
      IF ( NNF .GT. JDS ) GO TO 130
C
C STEP 2 (HEURISTIC SOLUTION THROUGH THE CORE PROBLEM).
C
C SORT THE ITEMS IN THE CORE PROBLEM.
      CALL SORTR(NNF,RA,IA2,JDIM)
C SOLVE THE CORE PROBLEM, THROUGH KP01M, WITH:
C   IA3(1)  TO  IA3(NNF) = SORTED PROFITS OF FREE ITEMS;
C   IA3(NNF+2)  TO  IA3(2*NNF+1) = DUMMY FOR PS;
C   IA4(1)  TO  IA4(NNF) = SORTED WEIGHTS OF FREE ITEMS;
C   IA4(NNF+2)  TO  IA4(2*NNF+1) = DUMMY FOR WS;
C   IA2(1)  TO  IA2(NNF) = POINTERS TO FREE ITEMS; ON INPUT THEY ARE
C                          .GT. 0 ; ON OUTPUT THEY ARE .GT. 0 (.LT. 0)
C                          IF THE CORRESPONDING ITEM IS SET TO 1 (TO 0);
C   IA2(NNF+2)  TO  IA2(2*NNF+1) = DUMMY FOR ZS;
C   IA2(2*NNF+2)  TO  IA2(3*NNF+1) = DUMMY FOR XX;
C   RA(1)  TO  RA(NNF) = DUMMY FOR MINW;
C   IZC (ON OUTPUT) = SOLUTION VALUE OF THE CORE PROBLEM;
C   IUBF0 (ON OUTPUT) = UPPER BOUND ON THE PROBLEM DEFINED BY THE FREE
C                       ITEMS AND THE ITEMS IN N0.
      NP2 = NNF + 2
      N2P2 = 2*NNF + 2
      N2P1 = 2*NNF + 1
      DO 30 I=1,NNF
        J = IA2(I)
        IA3(I) = P(J)
        IA4(I) = W(J)
   30 CONTINUE
      IZC = 0
      IUBF0 = 0
      CALL KP01M(NNF,IA3,IA4,ICW,MINW0,IZC,IA2,JFO,IUBF0,NNF+1,
     1           IA2(N2P2),IA3(NP2),IA4(NP2),IA2(NP2),RA)
C JUB = UPPER BOUND ON THE ORIGINAL PROBLEM.
      JUB = IZ1 + IUBF0
C DEFINE THE HEURISTIC SOLUTION OF VALUE IZ1 + IZC.
      Z = IZ1 + IZC
      ICWR = ICW
      JVAL = 1
      J = FN1
   40 IF ( J .GT. N ) GO TO 50
      X(J) = JVAL
      J = IA1(J)
      GO TO 40
   50 IF ( JVAL .EQ. 0 ) GO TO 60
      JVAL = 0
      J = FN0
      GO TO 40
   60 DO 80 I=1,NNF
        J = IA2(I)
        IF ( J .LT. 0 ) GO TO 70
        X(J) = 1
        ICWR = ICWR - W(J)
        GO TO 80
   70   MJ = - J
        X(MJ) = 0
        IA2(I) = MJ
   80 CONTINUE
      IF ( ICWR .LT. MINW0 ) GO TO 110
C THE SOLUTION IS NOT MAXIMAL.
      J = FN0
   90 IF ( J .GT. N ) GO TO 110
      IF ( W(J) .GT. ICWR ) GO TO 100
      X(J) = 1
      Z = Z + P(J)
      ICWR = ICWR - W(J)
      IF ( ICWR .LT. MINW0 ) GO TO 110
  100 J = IA1(J)
      GO TO 90
C HALTING TEST.
  110 IF ( Z .EQ. JUB ) GO TO 280
      IF ( JFO .EQ. 0 ) GO TO 280
C
C STEP 3 (REDUCTION, WITHOUT SORTING, OF THE ITEMS NOT IN CORE).
C
      NNFO = NNF
C ON RETURN OF REDNS, ASSUMING THAT N1 (N0) IS THE SET OF ITEMS INCLUDED
C IN (EXCLUDED FROM) THE SOLUTION, THE MEANING OF THE PARAMETERS IS THE
C SAME GIVEN AT STEP 1.
      CALL REDNS(N,P,W,IZC,IZ1,ICW,IA2,NNFO,NNF,IA1,FN1,FN0)
C HALTING TEST.
      IF ( NNF .EQ. NNFO ) GO TO 280
C ITEMS PREVIOUSLY FIXED ARE NOW FREE. SORT THE FREE ITEMS.
      IZN = Z - IZ1 - 1
      DO 120 I=1,NNF
        J = IA2(I)
        RA(J) = FLOAT(P(J))/FLOAT(W(J))
        X(J) = 2
  120 CONTINUE
      CALL SORTR(NNF,RA,IA2,JDIM)
      GO TO 160
C
C STEP 4 (REDUCTION, WITH PRELIMINAR SORTING, OF THE ORIGINAL PROBLEM).
C
C SORT THE ITEMS IN THE ORIGINAL PROBLEM AND DEFINE:
C   IA1(1) TO IA1(N) = POINTERS TO THE ORIGINAL ITEMS;
C   IA3(1) TO IA3(N) = SORTED PROFITS;
C   IA4(1) TO IA4(N) = SORTED WEIGHTS.
  130 DO 140 J=1,N
        IA1(J) = J
  140 CONTINUE
      IF ( JFS .EQ. 0 ) CALL SORTR(N,RA,IA1,JDIM)
      DO 150 I=1,N
        J = IA1(I)
        IA3(I) = P(J)
        IA4(I) = W(J)
  150 CONTINUE
C ON RETURN OF REDS:
C X(J) = 0, 1 OR 2 FOR ITEM J FIXED TO 0, TO 1 OR FREE;
C IA2(1) TO IA2(NNF) = POINTERS TO FREE ITEMS;
C IZH = HEURISTIC SOLUTION;
C ICW = C - TOTAL WEIGHT OF ITEMS FIXED TO 1;
C IZ1 = TOTAL PROFIT OF ITEMS FIXED TO 1.
      CALL REDS(N,IA3,IA4,P,W,C,IA1,NP1,NNF,X,IZ1,IZH,ICW,IA2)
      IZN = IZH - IZ1 - 1
C HALTING TEST.
      IF ( NNF .GT. 0 ) GO TO 160
      Z = IZH
      GO TO 280
  160 IF ( NNF .GT. JDS ) GO TO 200
C
C STEP 5 (EXACT SOLUTION OF THE REDUCED PROBLEM IF NNF IS SMALL).
C
C SOLVE THE REDUCED PROBLEM, THROUGH KP01M, WITH:
C   IA3(1) TO IA3(NNF) = SORTED PROFITS OF FREE ITEMS;
C   IA3(NNF+2) TO IA3(2*NNF+1) = DUMMY FOR PS;
C   IA4(1) TO IA4(NNF) = SORTED WEIGHTS OF FREE ITEMS;
C   IA4(NNF+2) TO IA4(2*NNF+1) = DUMMY FOR WS;
C   IA2(1) TO IA2(NNF) = POINTERS TO FREE ITEMS; ON INPUT THEY ARE
C                        .GT. 0 ; ON OUTPUT THEY ARE .GT. 0 (.LT. 0)
C                        IF THE CORRESPONDING ITEM IS SET TO 1 (TO 0);
C   IA2(NNF+2) TO IA2(2*NNF+1) = DUMMY FOR ZS;
C   IA2(2*NNF+2) TO IA2(3*NNF+1) = DUMMY FOR XX;
C   RA(1) TO RA(NNF) = DUMMY FOR MINW;
C   IZN = SOLUTION OF THE REDUCED PROBLEM.
      NP2 = NNF + 2
      N2P2 = 2*NNF + 2
      N2P1 = 2*NNF + 1
      DO 170 I=1,NNF
        J = IA2(I)
        IA3(I) = P(J)
        IA4(I) = W(J)
        IA2(I) = J
  170 CONTINUE
      IUBF0 = - 1
      CALL KP01M(NNF,IA3,IA4,ICW,C+1,IZN,IA2,JFO,IUBF0,NNF+1,
     1           IA2(N2P2),IA3(NP2),IA4(NP2),IA2(NP2),RA)
C DEFINE THE OPTIMAL SOLUTION.
      Z = IZ1 + IZN
      DO 190 I=1,NNF
        J = IA2(I)
        IF ( J .LT. 0 ) GO TO 180
        X(J) = 1
        GO TO 190
  180   MJ = - J
        X(MJ) = 0
  190 CONTINUE
      GO TO 280
C
C STEP 6 (EXACT SOLUTION OF THE REDUCED PROBLEM IF NNF IS LARGE).
C
C SOLVE THE REDUCED PROBLEM, THROUGH KP01M, WITH:
C   IA3(1) TO IA3(NNF) = SORTED PROFITS OF FREE ITEMS;
C   IA3(NNF+2) TO IA3(NP1) = PROFITS OF FIXED ITEMS;
C   IA4(1) TO IA4(NNF) = SORTED WEIGHTS OF FREE ITEMS;
C   IA4(NNF+2) TO IA4(NP1) = WEIGHTS OF FIXED ITEMS;
C   IA2 = POINTERS TO THE ORIGINAL ITEMS;
C         IA2(1) TO IA2(NNF) CORRESPOND TO FREE ITEMS; ON INPUT THEY
C                              ARE .GT. 0 ; ON OUTPUT THEY ARE .GT. 0
C                              (.LT. 0) IF THE CORRESPONDING ITEM IS
C                              SET TO 1 (TO 0);
C         IA2(NNF+1) TO IA2(N) CORRESPOND TO FIXED ITEMS; THEY ARE
C                               .GT. 0 (.LT. 0) IF THE CORRESPONDING
C                               ITEM IS FIXED TO 1 (TO 0);
C   P(1) TO P(NNF) = DUMMY FOR XX;
C   W(1) TO W(NNF) = DUMMY FOR PS;
C   X(1) TO X(NNF) = DUMMY FOR WS;
C   IA1(1) TO IA1(NNF) = DUMMY FOR ZS;
C   RA(1) TO RA(NNF) = DUMMY FOR MINW;
C   IZN = SOLUTION OF THE REDUCED PROBLEM.
  200 JL = N + 2
      DO 220 J=1,N
        IF ( X(J) .EQ. 2 ) GO TO 220
        JL = JL - 1
        IA3(JL) = P(J)
        IA4(JL) = W(J)
        IF ( X(J) .EQ. 0 ) GO TO 210
        IA2(JL-1) = J
        GO TO 220
  210   IA2(JL-1) = - J
  220 CONTINUE
      DO 230 I=1,NNF
        J = IA2(I)
        IA3(I) = P(J)
        IA4(I) = W(J)
  230 CONTINUE
      IUBF0 = - 1
      CALL KP01M(NNF,IA3,IA4,ICW,ICW+1,IZN,IA2,JFO,IUBF0,NNF+1,
     1           P,W,X,IA1,RA)
C RESET THE ORIGINAL PROBLEM AND DEFINE THE OPTIMAL SOLUTION.
      Z = IZ1 + IZN
      DO 270 I=1,N
        J = IA2(I)
        IF ( J .LT. 0 ) GO TO 240
        X(J) = 1
        GO TO 250
  240   J = - J
        X(J) = 0
  250   IF ( I .GT. NNF ) GO TO 260
        P(J) = IA3(I)
        W(J) = IA4(I)
        GO TO 270
  260   P(J) = IA3(I+1)
        W(J) = IA4(I+1)
  270 CONTINUE
  280 IF ( JUB .EQ. 0 ) JUB = Z
      RETURN
      END
      SUBROUTINE CHMT2(N,P,W,C,JFS,Z,JDIM)
C
C CHECK THE INPUT DATA.
C
      INTEGER P(JDIM),W(JDIM),C,Z
      IF ( N .GE. 2 .AND. N .LE. JDIM - 3 ) GO TO 10
      Z = - 1
      RETURN
   10 IF ( C .GT. 0 ) GO TO 30
   20 Z = - 2
      RETURN
   30 JSW = 0
      DO 40 J=1,N
        IF ( P(J) .LE. 0 ) GO TO 20
        IF ( W(J) .LE. 0 ) GO TO 20
        JSW = JSW + W(J)
        IF ( W(J) .LE. C ) GO TO 40
        Z = - 3
        RETURN
   40 CONTINUE
      IF ( JSW .GT. C ) GO TO 50
      Z = - 4
      RETURN
   50 IF ( JFS .EQ. 0 ) RETURN
      RR = FLOAT(P(1))/FLOAT(W(1))
      DO 60 J=2,N
        R = RR
        RR = FLOAT(P(J))/FLOAT(W(J))
        IF ( RR .GT. R ) GO TO 70
   60 CONTINUE
      RETURN
   70 Z = - 5
      RETURN
      END
      SUBROUTINE CORE(N,P,W,C,JFO,IZ1,ICW,MINW0,JDIM,FF,NNF,NF,FN1,
     1                FN0,PW)
C
C DETERMINE THE CORE PROBLEM.
C
C NF(J) = SUCCESSOR OF ITEM J IN THE CORRESPONDING SET;
C FN1 (LN1) = POINTER TO THE FIRST (LAST) ITEM IN SET N1 (FIXED TO 1);
C FN0 (LN0) = POINTER TO THE FIRST (LAST) ITEM IN SET N0 (FIXED TO 0);
C FNF (LNF) = POINTER TO THE FIRST (LAST) ITEM IN SET NF (FREE ITEMS);
C THE SET OF FREE ITEMS IS PARTITIONED INTO 3 SETS:
C    NG = ITEMS WITH RATIO .GT. LAMBDA,
C    NL = ITEMS WITH RATIO .LT. LAMBDA,
C    NE = ITEMS WITH RATIO .EQ. LAMBDA;
C NF(N+1) (NF(N+2),NF(N+3)) = POINTER TO THE FIRST ITEM IN NG (NL,NE);
C LNG (LNL,LNE) = POINTER TO THE LAST ITEM IN NG (NL,NE);
C MODNG (MODNL,MODNE) = NUMBER OF ITEMS IN NG (NL,NE).
C
      INTEGER P(JDIM),W(JDIM),FF(JDIM),NF(JDIM)
      INTEGER C,S1,S2,FNF,FN0,FN1,SE,TETA,SS
      REAL    PW(JDIM),LAMBDA
C
C STEP 0 (INITIALIZE).
C
      DO 10 J=1,N
        PW(J) = FLOAT(P(J))/FLOAT(W(J))
   10 CONTINUE
      IF ( JFO .EQ. 0 ) GO TO 20
      IF ( N .GE. 200 ) GO TO 20
      NNF = N
      RETURN
   20 TETA = 2.0*SQRT(FLOAT(N))
      IF ( JFO .EQ. 0 ) TETA = 5
      KETA = 20
      IF ( JFO .EQ. 0 ) KETA = 200
      ALPHA = 0.2
      IF ( JFO .EQ. 0 ) ALPHA = 0.
      BETA = 1.
      FN0 = N + 1
      LN0 = N + 1
      FN1 = N + 1
      LN1 = N + 1
      FNF = 1
      LNF = N
      DO 30 J=1,N
        NF(J) = J + 1
   30 CONTINUE
      NP1 = N + 1
      MODNF = N
      LN0P = N + 1
      LN1P = N + 1
      LS1 = 0
      KLAM = 0
C
C STEP 1.
C
C CHOOSE LAMBDA (MEDIAN OF THE RATIOS OF 3 ITEMS IN NF).
   40 LAMBDA = FMED(NF,JDIM,PW,N,FNF,LNF)
C
C STEP 2.
C
C DEFINE NG, NL, NE;
C COMPUTE S1 = SUM OF WEIGHTS OF ITEMS IN (N1 UNION NG);
C COMPUTE S2 = SUM OF WEIGHTS OF ITEMS IN (N1 UNION NG UNION NE).
   50 KLAM = KLAM + 1
      IF ( KLAM .GT. KETA ) GO TO 440
      S1 = LS1
      SE = 0
      LNG = N + 1
      NF(LNG) = NP1
      LNL = N + 2
      NF(LNL) = NP1
      LNE = N + 3
      NF(LNE) = NP1
      NF(LNF) = NP1
      J = FNF
      MODNG = 0
      MODNE = 0
   60 IF ( PW(J) - LAMBDA ) 70,80,90
   70 NF(LNL) = J
      LNL = J
      GO TO 100
   80 NF(LNE) = J
      LNE = J
      SE = SE + W(J)
      MODNE = MODNE + 1
      GO TO 100
   90 NF(LNG) = J
      LNG = J
      S1 = S1 + W(J)
      MODNG = MODNG + 1
  100 J = NF(J)
      IF ( J .LE. N ) GO TO 60
      S2 = S1 + SE
      IF ( S1 .GT. C ) GO TO 310
      IF ( S2 .LE. C ) GO TO 370
C LAMBDA HAS BEEN FOUND.
      IF ( MODNE .GE. TETA ) GO TO 110
      IF ( MODNG + MODNE .LT. TETA ) GO TO 470
      GO TO 460
C
C STEP 2.1 (ADD NG TO N1 AND NL TO N0).
C
  110 IF ( FN1 .LE. N ) GO TO 120
      FN1 = NF(N+1)
      GO TO 130
  120 NF(LN1) = NF(N+1)
  130 IF ( FN0 .LE. N ) GO TO 140
      FN0 = NF(N+2)
      GO TO 150
  140 NF(LN0) = NF(N+2)
  150 JFPR = 2
      IF ( NF(N+1) .GT. N ) GO TO 160
      LN1P = LN1
      LN1 = LNG
  160 IF ( NF(N+2) .GT. N ) GO TO 170
      LN0P = LN0
      LN0 = LNL
  170 NF(LN1) = NP1
      NF(LN0) = NP1
      NOUT = 0
      IF ( MODNE .EQ. TETA ) GO TO 260
C
C STEP 2.2 (SUBSET-SUM TYPE PROBLEMS).
C
      J = NF(N+3)
      SS = S1
      NESS = 0
  180 SS = SS + W(J)
      IF ( SS .GT. C ) GO TO 190
      NESS = NESS + 1
      J = NF(J)
      GO TO 180
  190 IF ( NESS .LE. TETA/2 ) GO TO 260
      IF ( MODNE - NESS .LE. TETA/2 ) GO TO 200
      NOUT = NESS - TETA/2
      GO TO 210
  200 NOUT = MODNE - TETA
C INSERT IN SET N1 THE FIRST NOUT ELEMENTS OF SET NE.
  210 J = NF(N+3)
      NESS = 0
  220 NESS = NESS + 1
      IF ( NESS .EQ. NOUT ) GO TO 230
      J = NF(J)
      GO TO 220
  230 IF ( FN1 .LE. N ) GO TO 240
      FN1 = NF(N+3)
      GO TO 250
  240 NF(LN1) = NF(N+3)
  250 LN1 = J
      NF(N+3) = NF(J)
      NF(LN1) = NP1
C DEFINE THE CORE PROBLEM.
  260 K = 0
      NF(LNE) = NP1
      J = NF(N+3)
  270 K = K + 1
      FF(K) = J
      IF ( K .EQ. TETA ) GO TO 280
      J = NF(J)
      GO TO 270
  280 IF ( MODNE .EQ. NOUT + TETA ) GO TO 530
C INSERT IN SET N0 THE LAST ELEMENTS OF SET NE.
      J = NF(J)
      IF ( FN0 .LE. N ) GO TO 290
      FN0 = J
      GO TO 300
  290 NF(LN0) = J
  300 LN0 = LNE
      GO TO 530
C
C STEP 3 (LAMBDA IS TOO SMALL).
C
  310 IF ( FLOAT(MODNG) .LT. FLOAT(TETA)*(1. - ALPHA) ) GO TO 470
C SET NF EQUAL TO NG.
      FNF = NF(N+1)
      LNF = LNG
      MODNF = MODNG
C ADD (NL UNION NE) TO N0.
      LN0P = LN0
      JFPR = 0
      IF ( FN0 .LE. N ) GO TO 330
      IF ( NF(N+2) .LE. N ) GO TO 320
      FN0 = NF(N+3)
      GO TO 360
  320 FN0 = NF(N+2)
      GO TO 350
  330 IF ( NF(N+2) .LE. N ) GO TO 340
      NF(LN0) = NF(N+3)
      GO TO 360
  340 NF(LN0) = NF(N+2)
  350 NF(LNL) = NF(N+3)
  360 LN0 = LNE
      GO TO 430
C
C STEP 4 (LAMBDA IS TOO LARGE).
C
  370 MODNL = MODNF - MODNG - MODNE
      IF ( FLOAT(MODNL) .LT. FLOAT(TETA)*(1. - ALPHA) ) GO TO 460
C SET NF EQUAL TO NL.
      FNF = NF(N+2)
      LNF = LNL
      MODNF = MODNL
C ADD (NG UNION NE) TO N1.
      LN1P = LN1
      JFPR = 1
      IF ( FN1 .LE. N ) GO TO 390
      IF ( NF(N+1) .LE. N ) GO TO 380
      FN1 = NF(N+3)
      GO TO 420
  380 FN1 = NF(N+1)
      GO TO 410
  390 IF ( NF(N+1) .LE. N ) GO TO 400
      NF(LN1) = NF(N+3)
      GO TO 420
  400 NF(LN1) = NF(N+1)
  410 NF(LNG) = NF(N+3)
  420 LN1 = LNE
      LS1 = S2
  430 IF ( FLOAT(MODNF) .GT. (1. + BETA)*FLOAT(TETA) ) GO TO 40
C
C STEP 5 (SET THE CORE PROBLEM EQUAL TO NF).
C
  440 J = FNF
      NF(LNF) = NP1
      K = 0
  450 K = K + 1
      FF(K) = J
      J = NF(J)
      IF ( J .LE. N ) GO TO 450
      GO TO 530
C
C STEP 6 (UPDATE LAMBDA).
C
C LAMBDA = MEDIAN OF THE RATIOS OF 3 ITEMS IN NG.
  460 JFRET = 1
      IF ( MODNG .LT. 3 ) GO TO 480
      JFRET = 0
      LAMBDA = FMED(NF,JDIM,PW,N,NF(N+1),LNG)
      GO TO 480
C LAMBDA = MEDIAN OF THE RATIOS OF 3 ITEMS IN NL.
  470 JFRET = 1
      MODNL = MODNF - MODNG - MODNE
      IF ( MODNL .LT. 3 ) GO TO 480
      JFRET = 0
      LAMBDA = FMED(NF,JDIM,PW,N,NF(N+2),LNL)
C RE-DEFINE THE PREVIOUS SET OF FREE ITEMS ( = NE UNION NG UNION NL ).
  480 FNF = NF(N+3)
      IF ( NF(N+1) .LE. N ) GO TO 500
      IF ( NF(N+2) .LE. N ) GO TO 490
      LNF = LNE
      GO TO 520
  490 NF(LNE) = NF(N+2)
      LNF = LNL
      GO TO 520
  500 NF(LNE) = NF(N+1)
      IF ( NF(N+2) .LE. N ) GO TO 510
      LNF = LNG
      GO TO 520
  510 NF(LNG) = NF(N+2)
      LNF = LNL
  520 IF ( JFRET .EQ. 1 ) GO TO 440
      GO TO 50
C
C STEP 7.
C
  530 NNF = K
C COMPUTE IZ1 AND ICW.
      IZ1 = 0
      ICW = C
      IF ( FN1 .GT. N ) GO TO 550
      J = FN1
      NF(LN1) = NP1
  540 IF ( J .GT. N ) GO TO 550
      IZ1 = IZ1 + P(J)
      ICW = ICW - W(J)
      J = NF(J)
      GO TO 540
C COMPUTE MINW0.
  550 MINW0 = 10*C
      IF ( FN0 .GT. N ) GO TO 570
      J = FN0
      NF(LN0) = NP1
  560 IF ( J .GT. N ) GO TO 570
      IF ( W(J) .LT. MINW0 ) MINW0 = W(J)
      J = NF(J)
      GO TO 560
C ADD ITEMS TO THE CORE PROBLEM UNTIL THE MAXIMUM WEIGHT
C IN CORE IS .LE. ICW .
  570 MAXWC = 0
      DO 580 K=1,NNF
        J = FF(K)
        IF ( W(J) .GT. MAXWC ) MAXWC = W(J)
  580 CONTINUE
  590 IF ( MAXWC .LE. ICW ) GO TO 640
      J = FN1
      PWMIN = PW(J)
      JMIN = J
      JMINP = 0
  600 IF ( J .GT. N ) GO TO 620
      IF ( PW(J) .GE. PWMIN ) GO TO 610
      PWMIN = PW(J)
      JMIN = J
      JMINP = JP
  610 JP = J
      J = NF(J)
      GO TO 600
  620 NNF = NNF + 1
      FF(NNF) = JMIN
      IZ1 = IZ1 - P(JMIN)
      ICW = ICW + W(JMIN)
      IF ( JMINP .NE. 0 ) GO TO 630
      FN1 = NF(JMIN)
      GO TO 590
  630 NF(JMINP) = NF(JMIN)
      GO TO 590
  640 RETURN
      END
      SUBROUTINE CORES(N,P,W,C,JFO,IZ1,ICW,MINW0,JDIM,FF,NNF,NF,FN1,
     1                 FN0,PW)
C
C DETERMINE THE CORE PROBLEM WHEN THE ITEMS ARE ALREADY SORTED
C ACCORDING TO DECREASING PROFIT PER UNIT WEIGHT.
C
C NF(J) = SUCCESSOR OF ITEM J IN THE CORRESPONDING SET;
C FN1 = POINTER TO THE FIRST ITEM IN SET N1 (FIXED TO 1);
C FN0 = POINTER TO THE FIRST ITEM IN SET N0 (FIXED TO 0);
C
      INTEGER P(JDIM),W(JDIM),FF(JDIM),NF(JDIM)
      INTEGER C,FN0,FN1,TETA
      REAL PW(JDIM)
C FIND THE CRITICAL ITEM.
      IF ( JFO .EQ. 0 ) GO TO 10
      IF ( N .GE. 200 ) GO TO 10
      NNF = N
      RETURN
   10 TETA = 2.0*SQRT(FLOAT(N))
      IF ( JFO .EQ. 0 ) TETA = 8
      IF ( TETA .GT. N ) TETA = N
      IC = C
      DO 20 J=1,N
        IF ( W(J) .GT. IC ) GO TO 30
        IC = IC - W(J)
   20 CONTINUE
C NO EXIT THIS WAY.
   30 LL = J - 1
C DEFINE THE CORE PROBLEM.
      LL1 = LL - TETA/2
      IF ( LL1 .LT. 1 ) LL1 = 1
      LL2 = LL1 + TETA - 1
      IF ( LL2 .LE. N ) GO TO 40
      LL2 = N
      LL1 = N - TETA + 1
   40 NNF = TETA
      IZ1 = 0
      ICW = C
      LL1M1 = LL1 - 1
      IF ( LL1M1 .GT. 0 ) GO TO 50
      FN1 = N + 1
      GO TO 70
   50 FN1 = 1
      DO 60 J=1,LL1M1
        IZ1 = IZ1 + P(J)
        ICW = ICW - W(J)
        NF(J) = J + 1
   60 CONTINUE
      NF(LL1M1) = N + 1
   70 MAXWC = 0
      DO 80 J=1,TETA
        JJ = J + LL1M1
        FF(J) = JJ
        IF ( W(JJ) .GT. MAXWC ) MAXWC = W(JJ)
   80 CONTINUE
      IF ( MAXWC .LE. ICW ) GO TO 110
      J = LL1
   90 J = J - 1
        IZ1 = IZ1 - P(J)
        ICW = ICW + W(J)
        NNF = NNF + 1
        FF(NNF) = J
      IF ( MAXWC .GT. ICW ) GO TO 90
      IF ( J .GT. 1 ) GO TO 100
      FN1 = N + 1
      GO TO 110
  100 NF(J-1) = N + 1
  110 MINW0 = 10*C
      LL2P1 = LL2 + 1
      IF ( LL2P1 .LE. N ) GO TO 120
      FN0 = N + 1
      GO TO 140
  120 FN0 = LL2P1
      DO 130 J=LL2P1,N
        IF ( W(J) .LT. MINW0 ) MINW0 = W(J)
        NF(J) = J + 1
  130 CONTINUE
  140 DO 150 J=1,NNF
        JJ = FF(J)
        PW(JJ) = FLOAT(P(JJ))/FLOAT(W(JJ))
  150 CONTINUE
      RETURN
      END
      FUNCTION FMED(NF,JDIM,PW,N,I1,ILAST)
C
C DEFINE FMED = MEDIAN OF THE RATIOS OF THE FIRST 2 AND THE
C               LAST ITEM.
C
      INTEGER NF(JDIM)
      REAL    PW(N)
      I2 = NF(I1)
      I3 = ILAST
      IF ( PW(I1) .LE. PW(I2) ) GO TO 10
      IF ( PW(I1) .LE. PW(I3) ) GO TO 20
      FMED = PW(I2)
      IF ( FMED .LT. PW(I3) ) FMED = PW(I3)
      RETURN
   10 IF ( PW(I1) .GE. PW(I3) ) GO TO 20
      FMED = PW(I2)
      IF ( FMED .GT. PW(I3) ) FMED = PW(I3)
      RETURN
   20 FMED = PW(I1)
      RETURN
      END
      SUBROUTINE KP01M(N,P,W,C,MINW0,Z,X,JFO,IUBF0,NP1,
     1                 XX,PS,WS,ZS,MINW)
C
C SOLVE, THROUGH BRANCH-AND-BOUND, A 0-1 SINGLE KNAPSACK PROBLEM.
C IT IS ASSUMED THAT THE ITEMS ARE SORTED ACCORDING TO DECREASING
C PROFIT PER UNIT WEIGHT.
C
      INTEGER P(NP1),W(NP1),X(N),C,Z
      INTEGER XX(N),PS(N),WS(N),ZS(N)
      REAL    MINW(N)
      INTEGER CWF,CWS,DIFF,VAL,R,T,UB,UB1
C INITIALIZE.
      CWF = C
      IP = 0
      CWS = C
      DO 10 LL=1,N
        IF ( W(LL) .GT. CWS ) GO TO 20
        IP = IP + P(LL)
        CWS = CWS - W(LL)
   10 CONTINUE
   20 LL = LL - 1
      UB = 0
      IF ( CWS .EQ. 0 ) GO TO 80
      P(NP1) = 0
      W(NP1) = C + 1
      NTEST = 2*N
      IF ( IUBF0 .LT. 0 ) GO TO 30
      KNP1 = (C + 1)/W(N) + 1
      P(NP1) = P(N)*KNP1
      W(NP1) = W(N)*KNP1
      NTEST = N
C NTEST IS USED TO AVOID INTEGER OVERFLOWS WHEN IUBF0 = 0 .
   30 IF ( LL + 2 .GT. NTEST ) GO TO 40
      UB = IP + CWS*P(LL+2)/W(LL+2)
      GO TO 50
   40 UB = IP + CWS*P(N)/W(N)
   50 A = W(LL+1) - CWS
      UB1 = FLOAT(IP + P(LL+1)) - A*FLOAT(P(LL))/FLOAT(W(LL))
      IF ( UB1 .GT. UB ) UB = UB1
      IF ( IUBF0 .EQ. 0 ) IUBF0 = UB
      IF ( JFO .EQ. 0 .AND. N .GT. 10 ) GO TO 80
      MINK = C + 1
      MINW(N) = MINK
      DO 60 J=2,N
        KK = N + 2 - J
        IF ( W(KK) .LT. MINK ) MINK = W(KK)
        MINW(KK-1) = MINK
   60 CONTINUE
      DO 70 J=1,N
        XX(J) = 0
   70 CONTINUE
      VAL = 0
      LOLD = N
      II = 1
      IF ( IUBF0 .LT. 0 ) GO TO 220
      IUBF0 = 0
      FW1 = W(1)
      FP1 = P(1)
      FPN1 = FLOAT(MINW0)*FLOAT(P(N))/FLOAT(W(N))
      A = W(LL+1) - CWS
      IB = FLOAT(IP + P(LL+1)) - A*FP1/FW1
      IF ( IB .GT. IUBF0 ) IUBF0 = IB
      GO TO 220
   80 Z = IP
      NN = LL + 1
      DO 90 J=NN,N
        X(J) = - X(J)
   90 CONTINUE
      GO TO 580
C TRY TO INSERT THE II-TH ITEM INTO THE CURRENT SOLUTION.
  100 IF ( W(II) .LE. C ) GO TO 120
      IF ( IUBF0 .LT. 0 ) GO TO 110
      A = W(II) - C
      IB = FLOAT(VAL + P(II)) - A*FP1/FW1
      IF ( IB .GT. IUBF0 ) IUBF0 = IB
  110 II1 = II + 1
      IF ( Z .GE. C*P(II1)/W(II1) + VAL ) GO TO 360
      II = II1
      GO TO 100
C BUILD A NEW CURRENT SOLUTION.
  120 IP = PS(II)
      CWS = C - WS(II)
      IN = ZS(II)
      DO 130 LL=IN,N
        IF ( W(LL) .GT. CWS ) GO TO 200
        IP = IP + P(LL)
        CWS = CWS - W(LL)
  130 CONTINUE
      LL = N
      IF ( IUBF0 .LT. 0 ) GO TO 140
      CALL NEWB(CWS,IP+VAL,MINW0,P(N),W(N),FP1,FPN1,FW1,IUBF0)
      CALL NEWB(CWS+W(N),IP+VAL-P(N),MINW0,P(N),W(N),FP1,FPN1,
     1          FW1,IUBF0)
  140 IF ( Z .GE. IP + VAL ) GO TO 360
      Z = IP + VAL
      NN = II - 1
      DO 160 J=1,NN
        IF ( XX(J) .EQ. 0 ) GO TO 150
        X(J) = IABS(X(J))
        GO TO 160
  150   X(J) = - IABS(X(J))
  160 CONTINUE
      DO 170 J=II,LL
        X(J) = IABS(X(J))
  170 CONTINUE
      IF ( LL .EQ. N ) GO TO 190
      NN = LL + 1
      DO 180 J=NN,N
        X(J) = - IABS(X(J))
  180 CONTINUE
  190 IF ( Z .NE. UB ) GO TO 360
      C = CWF
      GO TO 580
  200 LL = LL - 1
      IF ( IUBF0 .LT. 0 ) GO TO 210
      A = W(LL+1) - CWS
      IB = FLOAT(VAL + IP + P(LL+1)) - A*FP1/FW1
      IF ( IB .GT. IUBF0 ) IUBF0 = IB
  210 IF ( CWS .EQ. 0 ) GO TO 140
      IF ( Z .GE. VAL + IP + CWS*P(LL+1)/W(LL+1) ) GO TO 360
C SAVE THE CURRENT SOLUTION.
  220 WS(II) = C - CWS
      PS(II) = IP
      ZS(II) = LL + 1
      XX(II) = 1
      NN = LL - 1
      IF ( NN .LT. II ) GO TO 240
      DO 230 J=II,NN
        JP1 = J + 1
        WS(JP1) = WS(J) - W(J)
        PS(JP1) = PS(J) - P(J)
        ZS(JP1) = LL + 1
        XX(JP1) = 1
  230 CONTINUE
  240 J1 = LL + 1
      DO 250 J=J1,LOLD
        WS(J) = 0
        PS(J) = 0
        ZS(J) = J
  250 CONTINUE
      LOLD = LL
      C = CWS
      VAL = VAL + IP
      IF ( LL - (N - 2) ) 300, 270, 260
  260 II = N
      IF ( IUBF0 .GE. 0 ) GO TO 290
      GO TO 320
  270 IF ( C .LT. W(N) ) GO TO 280
      C = C - W(N)
      VAL = VAL + P(N)
      XX(N) = 1
      IF ( IUBF0 .LT. 0 ) GO TO 280
      II = N - 1
      CALL NEWB(C+W(N),VAL-P(N),MINW0,P(N),W(N),FP1,FPN1,FW1,
     1          IUBF0)
      GO TO 290
  280 II = N - 1
      IF ( IUBF0 .LT. 0 ) GO TO 320
      A = W(N) - C
      IB = FLOAT(VAL + P(N)) - A*FP1/FW1
      IF ( IB .GT. IUBF0 ) IUBF0 = IB
  290 CALL NEWB(C,VAL,MINW0,P(N),W(N),FP1,FPN1,FW1,IUBF0)
      GO TO 320
  300 II = LL + 2
      IF ( C .GE. INT(MINW(II-1)) ) GO TO 100
      IF ( IUBF0 .LT. 0 ) GO TO 320
C COMPUTE THE BOUND CORRESPONDING TO THE INSERTION OF ITEMS
C II,...,N.
      DO 310 J=II,N
        IF ( Z .GE. VAL + C*P(J)/W(J) ) GO TO 320
        A = W(J) - C
        IB = FLOAT(VAL + P(J)) - A*FP1/FW1
        IF ( IB .GT. IUBF0 ) IUBF0 = IB
  310 CONTINUE
      CALL NEWB(C,VAL,MINW0,P(N),W(N),FP1,FPN1,FW1,IUBF0)
C SAVE THE CURRENT OPTIMAL SOLUTION.
  320 IF ( Z .GE. VAL ) GO TO 350
      Z = VAL
      DO 340 J=1,N
        IF ( XX(J) .EQ. 0 ) GO TO 330
        X(J) = IABS(X(J))
        GO TO 340
  330   X(J) = - IABS(X(J))
  340 CONTINUE
      IF ( Z .NE. UB ) GO TO 350
      C = CWF
      GO TO 580
  350 IF ( XX(N) .EQ. 0 ) GO TO 360
      XX(N) = 0
      C = C + W(N)
      VAL = VAL - P(N)
C BACKTRACK.
  360 NN = II - 1
      IF ( NN .EQ. 0 ) GO TO 380
      DO 370 J=1,NN
        KK = II - J
        IF ( XX(KK) .EQ. 1 ) GO TO 390
  370 CONTINUE
  380 C = CWF
      GO TO 580
  390 R = C
      C = C + W(KK)
      VAL = VAL - P(KK)
      XX(KK) = 0
      IF ( R .LT. INT(MINW(KK)) ) GO TO 400
      II = KK + 1
      GO TO 100
  400 NN = KK + 1
      II = KK
C TRY TO SUBSTITUTE THE NN-TH ITEM FOR THE KK-TH.
  410 IF ( NN .GT. N ) GO TO 360
      IF ( Z .GE. VAL + C*P(NN)/W(NN) ) GO TO 360
      IF ( IUBF0 .LT. 0 ) GO TO 420
      IF ( NN .EQ. N ) CALL NEWB(C,VAL,MINW0,P(N),W(N),FP1,FPN1,
     1                           FW1,IUBF0)
  420 DIFF = W(NN) - W(KK)
      IF ( DIFF ) 530, 430, 440
  430 NN = NN + 1
      GO TO 410
  440 IF ( DIFF .LE. R ) GO TO 450
      IF ( IUBF0 .LT. 0 ) GO TO 430
      A = W(NN) - C
      IB = FLOAT(VAL + P(NN)) - A*FP1/FW1
      IF ( IB .GT. IUBF0 ) IUBF0 = IB
      GO TO 430
  450 IF ( IUBF0 .LT. 0 ) GO TO 480
C COMPUTE THE BOUND CORRESPONDING TO THE INSERTION OF ITEMS
C NN+1,...,N.
      NPRO = VAL + P(NN)
      NCW = C - W(NN)
      IF ( NN .EQ. N ) GO TO 470
      NN1 = NN + 1
      DO 460 J=NN1,N
        A = W(J) - NCW
        IB = FLOAT(NPRO + P(J)) - A*FP1/FW1
        IF ( IB .GT. IUBF0 ) IUBF0 = IB
  460 CONTINUE
  470 CALL NEWB(NCW,NPRO,MINW0,P(N),W(N),FP1,FPN1,FW1,IUBF0)
  480 IF ( Z .GE. VAL + P(NN) ) GO TO 430
      Z = VAL + P(NN)
      DO 500 J=1,KK
        IF ( XX(J) .EQ. 0 ) GO TO 490
        X(J) = IABS(X(J))
        GO TO 500
  490   X(J) = - IABS(X(J))
  500 CONTINUE
      JJ = KK + 1
      DO 510 J=JJ,N
        X(J) = - IABS(X(J))
  510 CONTINUE
      X(NN) = IABS(X(NN))
      IF ( Z .NE. UB ) GO TO 520
      C = CWF
      GO TO 580
  520 R = R - DIFF
      KK = NN
      NN = NN + 1
      GO TO 410
  530 T = R - DIFF
      IF ( T .GE. INT(MINW(NN)) ) GO TO 560
      IF ( IUBF0 .LT. 0 ) GO TO 430
C COMPUTE THE BOUND CORRESPONDING TO THE INSERTION OF ITEMS
C NN+1,...,N.
      NPRO = VAL + P(NN)
      NCW = C - W(NN)
      IF ( NN .EQ. N ) GO TO 550
      NN1 = NN + 1
      DO 540 J=NN1,N
        IF ( Z .GE. NPRO + NCW*P(J)/W(J) ) GO TO 430
        A = W(J) - NCW
        IB = FLOAT(NPRO + P(J)) - A*FP1/FW1
        IF ( IB .GT. IUBF0 ) IUBF0 = IB
  540 CONTINUE
  550 CALL NEWB(NCW,NPRO,MINW0,P(N),W(N),FP1,FPN1,FW1,IUBF0)
      GO TO 430
  560 IF ( Z .GE. VAL + P(NN) + T*P(NN+1)/W(NN+1) ) GO TO 360
      C = C - W(NN)
      VAL = VAL + P(NN)
      XX(NN) = 1
      II = NN + 1
      WS(NN) = W(NN)
      PS(NN) = P(NN)
      ZS(NN) = II
      N1 = NN + 1
      DO 570 J=N1,LOLD
        WS(J) = 0
        PS(J) = 0
        ZS(J) = J
  570 CONTINUE
      LOLD = NN
      GO TO 100
  580 IF ( IUBF0 .LT. 0 ) IUBF0 = UB
      IF ( IUBF0 .LT. Z ) IUBF0 = Z
      RETURN
      END
      SUBROUTINE NEWB(C,VAL,MINW0,IPN,IWN,FP1,FPN1,FW1,IUBF0)
C
C IMPROVE ON THE CURRENT UPPER BOUND IUBF0 BY TAKING INTO ACCOUNT THE
C ITEMS FOLLOWING THE CORE PROBLEM.
C
      INTEGER C,VAL
      IF ( C .LT. MINW0 ) GO TO 10
      IB = VAL + C*IPN/IWN
      GO TO 20
   10 A = MINW0 - C
      IB = FLOAT(VAL) + FPN1 - A*FP1/FW1
   20 IF ( IB .GT. IUBF0 ) IUBF0 = IB
      RETURN
      END
      SUBROUTINE REDNS(N,P,W,IZC,IZ1,ICW,FF,NNFO,NNF,NF,FN1,FN0)
C
C REDUCE, WITHOUT SORTING,THE ITEMS NOT IN CORE.
C
C FF(1) TO FF(NNFO) = FREE ITEMS (INPUT);
C FF(1) TO FF(NNF) = FREE ITEMS (OUTPUT).
C
      INTEGER P(N),W(N),FF(N),NF(N),FN1,FN0
C DETERMINE THE BREAK ITEM JBR.
      IC = ICW
      IP = 0
      DO 10 I=1,NNFO
        J = FF(I)
        IF ( W(J) .GT. IC ) GO TO 20
        IC = IC - W(J)
        IP = IP + P(J)
   10 CONTINUE
   20 JBR = I
C REDUCE THE ITEMS IN N1.
      IR = FN1
   30 IF ( IR .GT. N ) GO TO 70
      ICR = IC + W(IR)
      IPR = IP - P(IR)
      DO 40 I=JBR,NNFO
        J = FF(I)
        IF ( W(J) .GT. ICR ) GO TO 50
        ICR = ICR - W(J)
        IPR = IPR + P(J)
   40 CONTINUE
      J = FF(NNFO)
   50 IUB = IPR + ICR*P(J)/W(J)
      IF ( IUB .LE. IZC ) GO TO 60
C ITEM IR CANNOT BE FIXED TO 1.
      NNF = NNF + 1
      FF(NNF) = IR
      IZ1 = IZ1 - P(IR)
      ICW = ICW + W(IR)
   60 IR = NF(IR)
      GO TO 30
C REDUCE THE ITEMS IN N0.
   70 IR = FN0
   80 IF ( IR .GT. N ) RETURN
      IF ( W(IR) .GT. ICW ) GO TO 120
      ICR = IC - W(IR)
      IPR = IP + P(IR)
      I = JBR
      J = FF(I)
   90 IF ( ICR .GE. 0 ) GO TO 110
      I = I - 1
      IF ( I .EQ. 0 ) GO TO 100
      J = FF(I)
      ICR = ICR + W(J)
      IPR = IPR - P(J)
      GO TO 90
  100 J = FF(1)
  110 IUB = IPR + ICR*P(J)/W(J)
      IF ( IUB .LE. IZC ) GO TO 120
C ITEM IR CANNOT BE FIXED TO 0.
      NNF = NNF + 1
      FF(NNF) = IR
  120 IR = NF(IR)
      GO TO 80
      END
      SUBROUTINE REDS(N,PS,WS,P,W,C,IA1,NP1,NNF,X,IZ1,IZH,ICW,IA2)
C
C REDUCE THE ORIGINAL PROBLEM. IT IS ASSUMED THAT THE ITEMS ARE SORTED
C ACCORDING TO DECREASING PROFIT PER UNIT WEIGHT.
C
      INTEGER PS(NP1),WS(NP1),P(N),W(N),X(N),IA1(N),IA2(N)
      INTEGER C,CWR,AP,AR,FR,R,SW
C INITIALIZE.
      NNF = 0
      IP = 0
      CWR = C
      J = 1
   10 IF ( CWR .LE. WS(J) ) GO TO 20
      CWR = CWR - WS(J)
      IP = IP + PS(J)
      J = J + 1
      GO TO 10
   20 R = J
      JBR = IA1(R)
      IF ( CWR .LT. WS(J) ) GO TO 50
      IZH = IP + PS(R)
      DO 30 J=1,R
        JJ = IA1(J)
        X(JJ) = 1
   30 CONTINUE
      IR = R + 1
      DO 40 J=IR,N
        JJ = IA1(J)
        X(JJ) = 0
   40 CONTINUE
      RETURN
C COMPUTE THE UPPER BOUND.
   50 LL = J - 1
      PS(NP1) = 0
      WS(NP1) = C + 1
      IUB = IP + CWR*PS(LL+2)/WS(LL+2)
      A = WS(LL+1) - CWR
      IUB1 = FLOAT(IP + PS(LL+1)) - A*FLOAT(PS(LL))/FLOAT(WS(LL))
      IF ( IUB1 .GT. IUB ) IUB = IUB1
C GREEDY SOLUTION.
      IZG = IP
      ICWR = CWR
      IG = LL
      DO 60 I=J,N
        IF ( WS(I) .GT. ICWR ) GO TO 60
        IG = I
        IZG = IZG + PS(I)
        ICWR = ICWR - WS(I)
   60 CONTINUE
      IF ( IZG .LT. IUB ) GO TO 120
C DETERMINE THE GREEDY SOLUTION (OPTIMAL).
      IZH = IZG
      DO 70 I=1,LL
        II = IA1(I)
        X(II) = 1
   70 CONTINUE
      ICW = CWR
      IF ( J .GT. IG ) GO TO 100
      DO 90 I=J,IG
        II = IA1(I)
        IF ( WS(I) .GT. ICW ) GO TO 80
        X(II) = 1
        ICW = ICW - WS(I)
        GO TO 90
   80   X(II) = 0
   90 CONTINUE
  100 IF ( IG .EQ. N ) RETURN
      IG1 = IG + 1
      DO 110 I=IG1,N
        II = IA1(I)
        X(II) = 0
  110 CONTINUE
      RETURN
  120 IZH = IZG
C COMPUTE F(I) FOR I .LE. R .
      JFLAG = 0
      DO 150 I=1,R
        AR = CWR + WS(I)
        AP = IP - PS(I)
        K = R
  130   IF ( AR .LT. WS(K) ) GO TO 140
        AR = AR - WS(K)
        AP = AP + PS(K)
        K = K + 1
        GO TO 130
  140   II = IA1(I)
        X(II) = AP + AR*PS(K)/WS(K)
        IF ( IZH .GE. AP ) GO TO 150
        IZH = AP
        IF ( IZH .EQ. IUB ) GO TO 280
  150 CONTINUE
      FR = X(JBR)
C COMPUTE F(I) FOR I .GE. R .
      JFLAG = 1
      DO 180 I=R,N
        AR = CWR - WS(I)
        AP = IP + PS(I)
        K = R
  160   IF ( AR .GE. 0 ) GO TO 170
        K = K - 1
        AR = AR + WS(K)
        AP = AP - PS(K)
        GO TO 160
  170   II = IA1(I)
        X(II) = AP + AR*PS(K)/WS(K)
        IF ( IZH .GE. AP ) GO TO 180
        IZH = AP
        IF ( IZH .EQ. IUB ) GO TO 280
  180 CONTINUE
C TRY TO INSERT ITEMS IN THE SOLUTION.
      JR = R - 1
      ICW = C
      SW = 0
      IZ1 = 0
      DO 200 I=1,JR
        II = IA1(I)
        IF ( X(II) .LT. IZH ) GO TO 190
        X(II) = 2
        NNF = NNF + 1
        IA2(NNF) = II
        SW = SW + WS(I)
        GO TO 200
  190   X(II) = 1
        ICW = ICW - WS(I)
        IZ1 = IZ1 + PS(I)
  200 CONTINUE
      IR = R
      IF ( FR .GE. IZH ) GO TO 210
      X(JBR) = 1
      ICW = ICW - WS(R)
      IZ1 = IZ1 + PS(R)
      IR = R + 1
      IF ( IR .GT. N ) GO TO 240
  210 DO 230 I=IR,N
        II = IA1(I)
        IF ( X(II) .LT. IZH ) GO TO 220
        X(II) = 2
        NNF = NNF + 1
        IA2(NNF) = II
        SW = SW + WS(I)
        GO TO 230
  220   X(II) = 0
  230 CONTINUE
  240 IF ( NNF .EQ. 0 ) RETURN
      NNN = NNF
      NNF = 0
      DO 260 I=1,NNN
        II = IA2(I)
        IF ( W(II) .LE. ICW ) GO TO 250
        X(II) = 0
        SW = SW - W(II)
        GO TO 260
  250   NNF = NNF + 1
        IA2(NNF) = II
  260 CONTINUE
      IF ( ICW .LT. SW ) RETURN
      IF ( NNF .EQ. 0 ) RETURN
      DO 270 I=1,NNF
        II = IA2(I)
        X(II) = 1
        IZ1 = IZ1 + P(II)
  270 CONTINUE
      IZH = IZ1
      NNF = 0
      RETURN
C DETERMINE THE HEURISTIC SOLUTION (OPTIMAL).
  280 IF ( K .EQ. 1 ) GO TO 300
      K1 = K - 1
      DO 290 J=1,K1
        JJ = IA1(J)
        X(JJ) = 1
  290 CONTINUE
  300 IF ( K .GT. N ) GO TO 320
      DO 310 J=K,N
        JJ = IA1(J)
        X(JJ) = 0
  310 CONTINUE
  320 X(II) = JFLAG
      RETURN
      END
      SUBROUTINE SORTR(N,A,V,JDA)
C
C SORT THE REAL ARRAY A BY DECREASING VALUES (DERIVED FROM SUBROUTINE
C SORTZV OF THE C.E.R.N. LIBRARY).
C
C JDA           = LENGTH OF ARRAY A;
C N             = NUMBER OF ELEMENTS OF A TO BE SORTED;
C V(I) (INPUT)  = POINTER TO THE I-TH ELEMENT TO BE SORTED;
C V(I) (OUTPUT) = POINTER TO THE I-TH ELEMENT OF THE SORTED ARRAY.
C
C ON RETURN, ARRAY A IS UNCHANGED.
C
      INTEGER V(N),IU(20),IL(20)
      REAL    A(JDA)
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
