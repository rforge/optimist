      SUBROUTINE MTP(N,W,C,Z,XSTAR,JDIM,BACK,JCK,LB,
     1               WR,XSTARR,DUM,RES,REL,X,R,WA,WB,KFIX,FIXIT,
     2               XRED,LS,LSB,XHEU)
C
C THIS SUBROUTINE SOLVES THE BIN PACKING PROBLEM
C
C MINIMIZE Z = Y(1) + ... + Y(N)
C
C SUBJECT TO:
C
C        W(1)*X(I,1) + ... + W(N)*X(I,N) .LE. C*Y(I)  FOR I=1,...,N,
C        X(1,J) + ... + X(M,J) = 1                    FOR J=1,...,N,
C        Y(I) = 0 OR 1                                FOR I=1,...,N,
C        X(I,J) = 0 OR 1                   FOR I=1,...,N, J=1,...,N
C
C (I.E., MINIMIZE THE NUMBER OF BINS OF CAPACITY  C  NEEDED TO ALLOCATE
C  N  ITEMS OF SIZE  W(1),...,W(N) ).
C
C THE PROGRAM IS INCLUDED IN THE VOLUME
C   S. MARTELLO, P. TOTH, "KNAPSACK PROBLEMS: ALGORITHMS
C   AND COMPUTER IMPLEMENTATIONS", JOHN WILEY, 1990
C AND IMPLEMENTS THE BRANCH-AND-BOUND ALGORITHM DESCRIBED
C IN SECTION 8.5 .
C
C THE INPUT PROBLEM MUST SATISFY THE CONDITIONS
C
C   1) 2 .LE. N .LE. JDIM;
C   2) W(J) AND C POSITIVE INTEGERS;
C   3) W(J) .LE. C FOR J=1,..., N;
C   4) W(J) .GE. W(J+1) FOR J=1,...,N-1.
C
C IN THE OUTPUT SOLUTION (SEE BELOW) THE  Z  LOWEST INDEXED BINS ARE
C USED.
C
C MTP CALLS 14 PROCEDURES: CHMTP, ENUMER, FFDLS, FIXRED, HBFDS,
C                          INSERT, LCL2, L2, L3, MWFDS, RESTOR,
C                          SEARCH, SORTI2 AND UPDATE.
C
C THE PROGRAM IS COMPLETELY SELF-CONTAINED AND COMMUNICATION TO IT IS
C ACHIEVED SOLELY THROUGH THE PARAMETER LIST OF MTP.
C NO MACHINE-DEPENDENT CONSTANT IS USED .
C THE PROGRAM IS WRITTEN IN 1967 AMERICAN NATIONAL STANDARD FORTRAN
C AND IS ACCEPTED BY THE PFORT VERIFIER (PFORT IS THE PORTABLE
C SUBSET OF ANSI DEFINED BY THE ASSOCIATION FOR COMPUTING MACHINERY).
C THE PROGRAM HAS BEEN TESTED ON A DIGITAL VAX 11/780 AND AN H.P.
C 9000/840.
C
C MTP NEEDS 17 ARRAYS ( W ,  XSTAR ,  WR,  XSTARR ,  DUM ,  RES ,
C                       REL ,  X ,  R ,  WA ,  WB ,  KFIX ,  FIXIT ,
C                       XRED ,  LS ,  LSB  AND  XHEU ) OF LENGTH
C                       AT LEAST  JDIM .
C
C MEANING OF THE INPUT PARAMETERS:
C N        = NUMBER OF ITEMS;
C W(J)     = WEIGHT OF ITEM J;
C C        = CAPACITY OF THE BINS;
C JDIM     = DIMENSION OF THE 17 ARRAYS;
C BACK     = - 1 IF EXACT SOLUTION IS REQUIRED,
C          = MAXIMUM NUMBER OF BACKTRACKINGS TO BE PERFORMED,
C            IF HEURISTIC SOLUTION IS REQUIRED;
C JCK      = 1 IF CHECK ON THE INPUT DATA IS DESIRED,
C          = 0 OTHERWISE.
C
C MEANING OF THE OUTPUT PARAMETERS:
C Z        = VALUE OF THE SOLUTION FOUND IF Z .GT. 0 ,
C          = ERROR IN THE INPUT DATA (WHEN JCK=1) IF Z .LT. 0 : CONDI-
C            TION  - Z  IS VIOLATED;
C XSTAR(J) = BIN WHERE ITEM J IS INSERTED IN THE SOLUTION FOUND;
C LB       = LOWER BOUND ON THE OPTIMAL SOLUTION VALUE.
C
C ALL THE ARRAYS EXCEPT W AND XSTAR ARE DUMMY.
C
C ALL THE PARAMETERS ARE INTEGER. ON RETURN OF MTP ALL THE INPUT
C PARAMETERS ARE UNCHANGED EXCEPT  BACK , WHICH GIVES THE NUMBER OF
C BACKTRACKINGS PERFORMED.
C
      INTEGER W(JDIM),XSTAR(JDIM),C,Z,BACK
      INTEGER WR(JDIM),XSTARR(JDIM),DUM(JDIM),VSTAR
      INTEGER RES(JDIM),REL(JDIM),X(JDIM),R(JDIM),WA(JDIM),WB(JDIM),
     1        KFIX(JDIM),FIXIT(JDIM),XRED(JDIM),LS(JDIM),LSB(JDIM),
     2        XHEU(JDIM)
      Z = 0
      IF ( JCK .EQ. 1 ) CALL CHMTP(N,W,C,JDIM,Z)
      IF ( Z .LT. 0 ) RETURN
      LBSTAR = 0
      NN = 9*N/10
      IF ( W(1) + W(NN) .GE. C ) GO TO 40
C
C TRY A QUICK SOLUTION.
C
C HEURISTIC.
      CALL FFDLS(N,W,C,Z,R,XSTAR,LS,LSB,JDIM)
C LOWER BOUND L1.
      ISUMR = 0
      DO 10 J=1,Z
        ISUMR = ISUMR + R(J)
   10 CONTINUE
      ISUM = Z*C - ISUMR
      LBSTAR = (ISUM - 1)/C + 1
      IF ( LBSTAR .EQ. Z ) GO TO 70
C IMPROVED LOWER BOUND.
      ISS = 0
      DO 20 I=1,N
        IF ( W(I) + W(N) .LE. C ) GO TO 30
        ISS = ISS + W(I)
   20 CONTINUE
      GO TO 70
   30 ISS = ISUM - ISS
      LBSTAR = I - 1 + (ISS - 1)/C + 1
      IF ( LBSTAR .EQ. Z ) GO TO 70
C LOWER BOUND L2.
      CALL L2(N,W,C,LBAL,JDIM)
      IF ( LBAL .LE. LBSTAR ) GO TO 60
      LBSTAR = LBAL
      IF ( LBSTAR .EQ. Z ) GO TO 70
      GO TO 60
C
C REGULAR SOLUTION.
C
C LOWER BOUND L3 AND REDUCTION.
   40 ISUM = 0
      DO 50 I=1,N
        ISUM = ISUM + W(I)
   50 CONTINUE
   60 Z = 0
      CALL L3(N,W,C,0,M,DUM,XSTAR,NF,LB3,N+1,XSTARR,ISUM,Z,
     1        RES,REL,JDIM)
      IF ( LB3 .GT. LBSTAR ) LBSTAR = LB3
      IF ( NF .GT. 0 ) GO TO 80
   70 BACK = 0
      LB = LBSTAR
      RETURN
   80 IF ( NF .EQ. N ) GO TO 100
C DEFINE THE REDUCED PROBLEM.
      II = 0
      DO 90 I=1,N
        IF ( XSTAR(I) .GT. 0 ) GO TO 90
        II = II + 1
        WR(II) = W(I)
        XSTARR(II) = XSTARR(I) - M
   90 CONTINUE
      GO TO 120
  100 DO 110 I=1,N
        WR(I) = W(I)
  110 CONTINUE
  120 VSTAR = Z - M
      LB = LBSTAR - M
C BRANCH-AND-BOUND.
      CALL ENUMER(NF,WR,C,XSTARR,VSTAR,LB,BACK,
     1            X,R,WA,WB,KFIX,FIXIT,XRED,LS,LSB,DUM,XHEU,
     2            RES,REL,JDIM)
C RE-BUILD THE SOLUTION.
      Z = VSTAR + M
      LB = LB + M
      II = 0
      DO 130 I=1,N
        IF ( XSTAR(I) .GT. 0 ) GO TO 130
        II = II + 1
        XSTAR(I) = XSTARR(II) + M
  130 CONTINUE
      RETURN
      END
      SUBROUTINE CHMTP(N,W,C,JDIM,Z)
C
C CHECK THE INPUT DATA.
C
      INTEGER W(JDIM),C,Z
      IF ( N .GE. 2 .AND. N .LE. JDIM ) GO TO 10
      Z = - 1
      RETURN
   10 IF ( C .GT. 0 ) GO TO 30
   20 Z = - 2
      RETURN
   30 JWP = W(1)
      DO 60 J=1,N
        IF ( W(J) .LE. 0 ) GO TO 20
        IF ( W(J) .LE. C ) GO TO 40
        Z = - 3
        RETURN
   40   IF ( W(J) .LE. JWP ) GO TO 50
        Z = - 4
        RETURN
   50   JWP = W(J)
   60 CONTINUE
      RETURN
      END
      SUBROUTINE ENUMER(N,W,C,XSTAR,Z,LBSTAR,BACK,
     1                  X,R,WA,WB,KFIX,FIXIT,XRED,LS,LSB,LOCAL,XHEU,
     2                  RES,REL,JDIM)
C
C PERFORM A BRANCH-AND-BOUND SEARCH.
C COMPUTATION OF LOWER BOUNDS L2 AND L3 AT ALL NODES.
C REDUCTION AT ALL NODES (BUT INITIAL REDUCTION BEFORE CALLING).
C HEURISTIC ALGORITHMS AT ALL NODES.
C DOMINANCE TESTS AT ALL NODES.
C
C KFIX(K)  = 0 IF NO ITEM WAS FIXED AT LEVEL  K ,
C          = POINTER TO THE FIRST ITEM FIXED AT LEVEL  K , OTHERWISE.
C FIXIT(J) = 0  IF ITEM  J  IS NOT FIXED BY REDUCTION,
C          = POINTER TO THE NEXT ITEM FIXED BY REDUCTION AT THE SAME
C            LEVEL, OTHERWISE ( = - 1  FOR THE LAST ITEM FIXED BY
C            REDUCTION AT A LEVEL).
C LS(I)    = POINTER TO THE LAST ITEM INSERTED IN BIN  I  IF
C            LSB(I) = N + 1 ,
C          = POINTER TO THE LAST ITEM WHICH CAN BE (AND WAS) INSERTED
C            WITH ITEM  LSB(I)  AS LAST BUT ONE, IF LSB(I) .LE. N .
C
      INTEGER W(JDIM),C,XSTAR(JDIM),Z,BACK
      INTEGER X(JDIM),R(JDIM),WA(JDIM),WB(JDIM),KFIX(JDIM),FIXIT(JDIM),
     1        XRED(JDIM),LS(JDIM),LSB(JDIM),LOCAL(JDIM),XHEU(JDIM),
     2        RES(JDIM),REL(JDIM)
      INTEGER ZZ
C
C INITIALIZATION.
C
      MBACK = BACK
      BACK = 0
      JDIRK = 1
C FIRST HEURISTIC.
      CALL FFDLS(N,W,C,ZZ,R,X,LS,LSB,JDIM)
      IF ( ZZ .GE. Z ) GO TO 20
      DO 10 I=1,N
        XSTAR(I) = X(I)
   10 CONTINUE
      Z = ZZ
   20 DO 30 I=1,N
        FIXIT(I) = 0
        KFIX(I) = 0
   30 CONTINUE
      LASTW = W(N)
C LOWER BOUND L1.
      ISUMR = 0
      DO 40 J=1,ZZ
        ISUMR = ISUMR + R(J)
   40 CONTINUE
      ISUM = ZZ*C - ISUMR
      LB1 = (ISUM - 1)/C + 1
      IF ( LB1 .GT. LBSTAR ) LBSTAR = LB1
      IF ( LBSTAR .EQ. Z ) RETURN
C IMPROVED LOWER BOUND.
      ISS = 0
      DO 50 I=1,N
        IF ( W(I) + W(N) .LE. C ) GO TO 60
        ISS = ISS + W(I)
   50 CONTINUE
C NO EXIT THIS WAY.
   60 ISS = ISUM - ISS
      LB1I = I - 1 + (ISS - 1)/C + 1
      IF ( LB1I .GT. LBSTAR ) LBSTAR = LB1I
      IF ( LBSTAR .EQ. Z ) RETURN
C LOWER BOUND L2.
      CALL L2(N,W,C,LB2,JDIM)
      IF ( LB2 .GT. LBSTAR ) LBSTAR = LB2
      IF ( LBSTAR .EQ. Z ) RETURN
      KZFFD = ZZ
C SECOND HEURISTC.
      MW = 1
      MB = 1
      WA(1) = C
      CALL HBFDS(N,W,C,MB,WA,LOCAL,WB,JDIM)
      IF ( MB .GE. Z ) GO TO 80
      DO 70 I=1,N
        XSTAR(I) = WB(I)
   70 CONTINUE
      Z = MB
      IF ( Z .EQ. LBSTAR ) RETURN
C THIRD HEURISTIC.
   80 CALL MWFDS(N,W,C,MW,WA,LOCAL,LBSTAR,WB,JDIM)
      IF ( MW .GE. Z ) GO TO 100
      DO 90 I=1,N
        XSTAR(I) = WB(I)
   90 CONTINUE
      Z = MW
      IF ( Z .EQ. LBSTAR ) RETURN
C
C ITERATIVE PART.
C
C BACKTRACKING WHEN THE OPTIMAL SOLUTION HAS BEEN UPDATED.
C FIND THE FIRST ITEM  K  INSERTED IN BIN  Z .
  100 K = N
      ZZ = KZFFD - 1
  110 IF ( FIXIT(K) .NE. 0 ) GO TO 130
        J = X(K)
        R(J) = R(J) + W(K)
        IF ( KFIX(K) .EQ. 0 ) GO TO 120
        CALL RESTOR(K,ZZ,N,C,KFIX,FIXIT,W,X,R,LASTW,JDIM)
  120   IF ( R(KZFFD) .EQ. C ) GO TO 140
  130   K = K - 1
        JDIRK = - 1
      GO TO 110
C FIND THE NEXT ITEM  K  NOT INSERTED IN BIN  Z - 1 .
  140 K = K - 1
        JDIRK = - 1
        IF ( FIXIT(K) .NE. 0 ) GO TO 140
        J = X(K)
        IF ( J .LT. ZZ ) GO TO 150
        R(J) = R(J) + W(K)
        IF ( KFIX(K) .EQ. 0 ) GO TO 140
        CALL RESTOR(K,ZZ,N,C,KFIX,FIXIT,W,X,R,LASTW,JDIM)
      GO TO 140
  150 IF ( R(ZZ) .EQ. C ) ZZ = ZZ - 1
C BACKTRACKING ON ITEM  K.
  160 IF ( K .EQ. 1 ) RETURN
        IF ( FIXIT(K) .NE. 0 ) GO TO 180
        IF ( BACK .EQ. MBACK ) RETURN
        BACK = BACK + 1
        J = X(K)
        R(J) = R(J) + W(K)
        IF ( KFIX(K) .EQ. 0 ) GO TO 170
        CALL RESTOR(K,ZZ,N,C,KFIX,FIXIT,W,X,R,LASTW,JDIM)
        IF ( R(J) .LT. C ) GO TO 190
        GO TO 180
  170   IF ( R(ZZ) .LT. C ) GO TO 190
        LSB(ZZ) = N + 1
        ZZ = ZZ - 1
  180   K = K - 1
        JDIRK = - 1
        GO TO 160
C FIND THE FIRST BIN FOLLOWING BIN  J  WHERE ITEM  K  CAN BE INSERTED.
  190   IF ( J .LT. Z - 1 ) GO TO 200
        K = K - 1
        JDIRK = - 1
      GO TO 160
  200 J = J + 1
      IF ( R(J) .LT. W(K) ) GO TO 190
C DOMINANCE TESTS.
      IF ( LSB(J) .GT. N ) GO TO 250
      IF ( J .GT. ZZ ) GO TO 250
      LSJ = LS(J)
      IF ( W(K) .GT. W(LSJ) ) GO TO 210
      IF ( W(K) + LASTW .LE. R(J) ) GO TO 250
      GO TO 190
  210 LSJ = LSB(J)
      IF ( W(K) .GT. W(LSJ) ) GO TO 250
      NEXT = LS(J) - 1
      LSJ = LS(J)
  220 IF ( NEXT .LE. K ) GO TO 190
        IF ( FIXIT(NEXT) .GT. 0 ) GO TO 230
        IF ( W(NEXT) .GT. W(LSJ) ) GO TO 240
  230   NEXT = NEXT - 1
      GO TO 220
  240 IF ( W(K) + W(NEXT) .GT. R(J) ) GO TO 190
  250 X(K) = J
      R(J) = R(J) - W(K)
      IF ( R(J) .LT. LASTW ) GO TO 260
      LS(J) = K
      LSB(J) = N + 1
      GO TO 270
  260 LSB(J) = LS(J)
      LS(J) = K
  270 IF ( J .LE. ZZ ) GO TO 280
      ZZ = ZZ + 1
C FORWARD STEP.
  280 IF ( K .EQ. N ) GO TO 420
C
C COMPUTATION OF A LOCAL LOWER BOUND.
C ON OUTPUT FROM LCL2, LLB  CONTAINS THE LOWER BOUND, WHILE  NA ,
C WA  AND  WB  DEFINE THE PROBLEM TO BE REDUCED:
C NA  ITEMS WITH WEIGHTS IN  WA . FROM  WA(1)  TO  WA(ZZ)  WEIGHTS
C CORRESPOND TO BINS PARTIALLY FILLED, FROM  WA(ZZ+1)  TO  WA(NA)
C CORRESPOND TO FREE ITEMS. FROM  WB(1)  TO  WB(ZZ)  POINTERS TO THE
C BINS, FROM  WB(ZZ+1)  TO  WB(NA)  POINTERS TO THE ITEMS.
C
      CALL LCL2(N,W,C,ISUM,R,FIXIT,ZZ,Z,K,NA,WA,WB,LLB,JDIM)
      IF ( LLB .GE. Z ) GO TO 160
C
C REDUCTION.
C ON RETURN FROM  L3 : NBIN  BINS (IN TOTAL, OLD + NEW), XRED
C GIVES THE CORRESPONDING PARTIAL SOLUTION, NFREE THE NUMBER OF
C REMAINING FREE ITEMS. IF NFREE .LT. 0 , L3  TRIED TO MATCH
C PAIRS OF BINS.
C
      CALL L3(NA,WA,C,ZZ,NBIN,LOCAL,XRED,NFREE,LBR,Z,XHEU,ISUM,MR,
     1        RES,REL,JDIM)
      IF ( NFREE .LT. 0 ) GO TO 160
      IF ( NFREE .EQ. 0 ) GO TO 330
      IF ( LBR .GE. Z ) GO TO 160
      IF ( MR .GE. Z ) GO TO 320
      Z = MR
      DO 290 I=1,N
        XSTAR(I) = X(I)
  290 CONTINUE
      JZ1 = ZZ + 1
      DO 310 II=JZ1,NA
        I = WB(II)
        J = XHEU(II)
        IF ( J .LE. ZZ ) GO TO 300
        XSTAR(I) = J
        GO TO 310
  300   XSTAR(I) = WB(J)
  310 CONTINUE
      IF ( Z .EQ. LBSTAR ) RETURN
      IF ( LBR .GE. Z ) GO TO 160
  320 CALL FIXRED(W,R,NA,WA,WB,ZZ,X,NBIN,XRED,K,KFIX,FIXIT,JDIM)
      LASTW = WA(NA)
      GO TO 370
  330 IF ( NBIN .GE. Z ) GO TO 160
      Z = NBIN
      DO 340 I=1,N
        XSTAR(I) = X(I)
  340 CONTINUE
      JZ1 = ZZ + 1
      DO 360 II=JZ1,NA
        I = WB(II)
        J = XRED(II)
        IF ( J .LE. ZZ ) GO TO 350
        XSTAR(I) = J
        GO TO 360
  350   XSTAR(I) = WB(J)
  360 CONTINUE
      IF ( Z .EQ. LBSTAR ) RETURN
      GO TO 160
C LOCAL HEURISTICS.
C ON OUTPUT FROM FIXRED, THE  NA  FREE WEIGHTS ARE IN  WA , THE
C CORRESPONDING POINTERS IN  XRED .
  370 MB = ZZ
      CALL HBFDS(NA,WA,C,MB,R,LOCAL,WB,JDIM)
      IF ( MB .GE. Z ) GO TO 380
      CALL UPDATE(N,Z,XSTAR,NA,MB,X,WB,XRED,JDIM)
      IF ( Z .EQ. LBSTAR ) RETURN
      IF ( LLB .GE. Z ) GO TO 160
  380 MW = ZZ
      CALL MWFDS(NA,WA,C,MW,R,LOCAL,LLB,WB,JDIM)
      IF ( MW .GE. Z ) GO TO 390
      CALL UPDATE(N,Z,XSTAR,NA,MW,X,WB,XRED,JDIM)
      IF ( Z .EQ. LBSTAR ) RETURN
      IF ( LLB .GE. Z ) GO TO 160
  390 K = K + 1
      IF ( JDIRK .EQ. 1 ) GO TO 410
      DO 400 II=1,ZZ
        IF ( LS(II) .LT. K ) GO TO 400
        LS(II) = LSB(II)
        LSB(II) = N + 1
  400 CONTINUE
  410 IF ( FIXIT(K) .NE. 0 ) GO TO 370
      J = 0
      GO TO 190
  420 Z = ZZ
      DO 430 I=1,N
        XSTAR(I) = X(I)
  430 CONTINUE
      KZFFD = Z
      IF ( Z .GT. LBSTAR ) GO TO 100
      RETURN
      END
      SUBROUTINE FFDLS(N,W,C,M,K,X,LS,LSB,JDIM)
C
C PERFORM A FIRST-FIT DECREASING HEURISTIC AND INITIALIZE  LS
C AND  LSB .
C TIME COMPLEXITY  O(N**2) .
C
      INTEGER W(JDIM),C,K(JDIM),X(JDIM),LS(JDIM),LSB(JDIM)
      M = 1
      K(1) = C - W(1)
      X(1) = 1
      LS(1) = 1
      LSB(1) = N + 1
      DO 40 I=2,N
C INSERT THE NEXT ITEM.
        IWI = W(I)
        DO 10 J=1,M
          IF ( IWI .LE. K(J) ) GO TO 20
   10   CONTINUE
C INITIALIZE A NEW BIN.
        M = M + 1
        K(M) = C - IWI
        X(I) = M
        LS(M) = I
        LSB(M) = N + 1
        GO TO 40
C INSERT THE ITEM INTO AN OLD BIN.
   20   K(J) = K(J) - IWI
        X(I) = J
        IF ( K(J) .LT. W(N) ) GO TO 30
        LS(J) = I
        GO TO 40
   30   LSB(J) = LS(J)
        LS(J) = I
   40 CONTINUE
      RETURN
      END
      SUBROUTINE FIXRED(W,R,NA,WA,WB,ZZ,X,NBIN,XRED,K,KFIX,FIXIT,
     1                  JDIM)
C
C FIX THE VARIABLES AFTER A LOCAL REDUCTION.
C CURRENT SOLUTION: ZZ ,  X .
C CURRENT LEVEL: K .
C ON OUTPUT, WA  CONTAINS THE WEIGHTS OF THE  NA  FREE ITEMS,
C XRED  THE CORRESPONDING POINTERS.
C
      INTEGER W(JDIM),R(JDIM),WA(JDIM),WB(JDIM),X(JDIM),XRED(JDIM),
     1        KFIX(JDIM),FIXIT(JDIM)
      INTEGER ZZ
C FIND THE FIRST ITEM FIXED.
      IZ1 = ZZ + 1
      IA = 0
      DO 10 II=IZ1,NA
        IF ( XRED(II) .GT. 0 ) GO TO 20
        IA = IA + 1
        WA(IA) = WA(II)
        XRED(IA) = WB(II)
   10 CONTINUE
      NA = IA
      RETURN
   20 I = WB(II)
      KFIX(K) = I
   30 LAST = I
      J = XRED(II)
      IF ( J .LE. ZZ ) GO TO 40
      X(I) = J
      R(J) = R(J) - W(I)
      GO TO 50
   40 X(I) = WB(J)
      IWBJ = WB(J)
      R(IWBJ) = R(IWBJ) - W(I)
   50 IF ( II .EQ. NA ) GO TO 70
      II = II + 1
      IF ( XRED(II) .GT. 0 ) GO TO 60
      IA = IA + 1
      WA(IA) = WA(II)
      XRED(IA) = WB(II)
      GO TO 50
   60 I = WB(II)
      FIXIT(LAST) = I
      GO TO 30
   70 FIXIT(LAST) = - 1
      ZZ = NBIN
      NA = IA
      RETURN
      END
      SUBROUTINE HBFDS(N,W,C,M,KK,K,X,JDIM)
C
C PERFORM A BEST-FIT DECREASING HEURISTIC.
C FOR LOCAL USE WITH CURRENT SOLUTION GIVEN.
C TIME COMPLEXITY  O(N**2) .
C
      INTEGER W(JDIM),C,K(JDIM),X(JDIM),KK(JDIM)
      DO 10 J=1,M
        K(J) = KK(J)
   10 CONTINUE
      DO 40 I=1,N
C INSERT THE NEXT ITEM.
        IWI = W(I)
        MINRES = C
        DO 20 J=1,M
          KRES = K(J) - IWI
          IF ( KRES .LT. 0 ) GO TO 20
          IF ( KRES .GE. MINRES ) GO TO 20
          MINRES = KRES
          JM = J
   20   CONTINUE
        IF ( MINRES . LT. C ) GO TO 30
C INITIALIZE A NEW BIN.
        M = M + 1
        K(M) = C - IWI
        X(I) = M
        GO TO 40
C INSERT THE ITEM  INTO AN OLD BIN.
   30   K(JM) = K(JM) - IWI
        X(I) = JM
   40 CONTINUE
      RETURN
      END
      SUBROUTINE INSERT(I,M,FS,X,IFP,K,JDIM)
C
C INSERT ITEM  I  IN BIN  M  AND UPDATE  FS ,  X ,  IFP ,  K .
C
      INTEGER X(JDIM),K(JDIM),FS
      IS = X(I)
      IP = K(I)
      IF ( IS .GT. 0 ) GO TO 10
      IFP = IP
      GO TO 20
   10 K(IS) = IP
   20 IF ( IP .GT. 0 ) GO TO 30
      FS = IS
      GO TO 40
   30 X(IP) = IS
   40 X(I) = M
      RETURN
      END
      SUBROUTINE LCL2(N,W,C,ISUM,R,FIXIT,ZZ,Z,K,NA,WA,WB,LLB,JDIM)
C
C COMPUTE A LOCAL LOWER BOUND AND EXECUTE A PREPROCESSING
C FOR REDUCTION.
C
      INTEGER W(JDIM),R(JDIM),WA(JDIM),WB(JDIM),FIXIT(JDIM),C,ZZ,Z
      I = N
   10 IF ( FIXIT(I) .EQ. 0 ) GO TO 20
      I = I - 1
      GO TO 10
   20 JWN = W(I)
      ISS = ISUM
      KCB = 0
      DO 30 J=1,ZZ
        WB(J) = J
        IF ( R(J) .GE. JWN ) GO TO 30
        ISS = ISS - (C - R(J))
        KCB = KCB + 1
   30 CONTINUE
      LLB = KCB + (ISS - 1)/C + 1
      IF ( LLB .GE. Z ) RETURN
      CALL SORTI2(ZZ,R,WB,JDIM)
      DO 40 I=1,ZZ
        IWBI = WB(I)
        WA(I) = C - R(IWBI)
   40 CONTINUE
      K1 = K + 1
      NA = ZZ
      DO 50 I=K1,N
        IF ( FIXIT(I) .NE. 0 ) GO TO 50
        NA = NA + 1
        WA(NA) = W(I)
        WB(NA) = I
   50 CONTINUE
      CALL L2(NA,WA,C,LBA,JDIM)
      IF ( LBA .GT. LLB ) LLB = LBA
      RETURN
      END
      SUBROUTINE L2(N,W,C,LB,JDIM)
C
C COMPUTE LOWER BOUND  L2 .
C
      INTEGER W(JDIM),C
      IF ( W(1) .GT. C/2 ) GO TO 40
C
C CASE 1: ALL ITEMS ARE .LE. C/2 .
C
      LB = 0
      IRAT = C/W(1)
      II = 2
   10 DO 20 I=II,N
        IR = C/W(I)
        IF ( IR .GT. IRAT ) GO TO 30
   20 CONTINUE
      NLB = (N - 1)/IRAT + 1
      IF ( LB .LT. NLB ) LB = NLB
      RETURN
   30 NLB = (I - 2)/IRAT + 1
      IF ( LB .LT. NLB ) LB = NLB
      IF ( (N - 1)/IR + 1 .LE. LB ) RETURN
      II = I
      IRAT = IR
      GO TO 10
C
C CASE 2: THERE EXIST ITEMS  .GT. C/2 .
C
   40 CALL SEARCH(N,W,FLOAT(C)/2.,N12,JDIM)
C N12 = N1 + N2 .
      LB = N12
      IF ( N12 .EQ. N ) RETURN
      N12P1 = N12 + 1
      NMN12 = N - N12
      JBW = C - W(N12)
C I2 = NEXT ITEM TO BE CONSIDERED FOR POSSIBLE INCLUSION IN  N2 .
C I3 = NEXT ITEM TO BE CONSIDERED FOR POSSIBLE INCLUSION IN  N3 .
      I2 = N12
      I3 = N12P1
   50 JWW = W(I3)
   60 I3 = I3 + 1
      IF ( I3 .GT. N ) GO TO 70
      IF ( W(I3) .GE. JWW ) GO TO 60
   70 JWST = W(I3-1)
      JBWST = C - JWST
      N3 = I3 - N12P1
   80 IF ( W(I2) .GT. JBWST ) GO TO 100
      IF ( I2 .EQ. 1 ) GO TO 90
      I2 = I2 - 1
      GO TO 80
   90 N2 = N12
      GO TO 110
  100 N2 = N12 - I2
  110 NN = N2*(JBW/JWST)
      JADD = 0
      IF ( N3 .GT. NN ) JADD = (N3 - NN - 1)/(C/JWST) + 1
      NLB = N12 + JADD
      IF ( LB .LT. NLB ) LB = NLB
      NADD = (NMN12 - NN - 1)/(C/JWST) + 1
      IF ( N12 + NADD .LE. LB ) RETURN
      IF ( I3 .GT. N ) RETURN
      GO TO 50
      END
      SUBROUTINE L3(N,W,C,ZZ,MRED,K,X,NFREER,LB,UB,XHEU,ISUM,NUB,
     1              RES,REL,JDIM)
C
C REDUCE THE CURRENT PROBLEM, AND COMPUTE LOWER BOUND  L3  AND A NEW
C UPPER BOUND  NUB .
C ON OUTPUT:
C NFREER    = NUMBER OF UNASSIGNED VARIABLES;
C MRED      = NUMBER OF BINS USED IN THE REDUCTION WITHOUT RELAXATION;
C X(I)      = BIN IN WHICH ITEM I HAS BEEN INSERTED BY THE REDUCTION
C             PROCEDURE WITHOUT RELAXATION,
C           = 0 IF ITEM I HAS NOT BEEN ASSIGNED
C             (DURING EXECUTION,  X(I)  IS USED AS A POINTER TO THE
C             NEXT FREE ITEM FOLLOWING  I  OR AS THE BIN IN WHICH
C             ITEM  I  HAS BEEN INSERTED WITH RELAXATION);
C XHEU(I)   = BIN IN WHICH ITEM  I  HAS BEEN INSERTED IN THE HEURISTIC
C             SOLUTION CORRESPONDING TO NUB
C             (DURING EXECUTION,  XHEU(I)  CONTAINS THE LEVEL OF
C             RELAXATION AT WHICH ITEM  I  HAS BEEN ASSIGNED;
C             IF  XHEU(I) .GT. N , ITEM  I  HAS BEEN RELAXED);
C RES(L)    = RESIDUAL CAPACITY OF BIN  L (L=1,...,M);
C NREL      = NUMBER OF RELAXED ITEMS;
C REL(IREL) = IREL-TH RELAXED ITEM (IREL=1,...,NREL);
C KREL      = LEVEL OF RELAXATION;
C MREL      = MAXIMUM LEVEL OF RELAXATIONS FOR WHICH THE REDUCTION IS
C             VALID;
C KINF      = N + 1  IF NO PAIR OF BINS HAS BEEN MATCHED,
C           = LEVEL AT WHICH THE FIRST PAIR OF BINS HAS BEEN MATCHED,
C             OTHERWISE.
C
C THE ITEMS FROM  1  TO  ZZ  REPRESENT BINS PARTIALLY FILLED,
C THOSE FROM  ZZ + 1  TO  N  ARE REAL ITEMS.
C THE SUBROUTINE TERMINATES AND RETURNS  NFREE = - 1  IF IT
C ATTEMPTS TO MATCH PAIRS OF BINS.
C
      INTEGER W(JDIM),C,K(JDIM),X(JDIM),XHEU(JDIM),ZZ,UB
      INTEGER RES(JDIM),REL(JDIM),FS,TP,WW,QP,WWW
C
C INITIALIZATION.
C
      LB = 0
      ISUMR = ISUM
      MP = 0
      KINF = N + 1
      IF ( ZZ .EQ. 0 ) GO TO 20
      DO 10 I=1,ZZ
        RES(I) = 0
   10 CONTINUE
   20 DO 30 I=1,N
        X(I) = I + 1
        K(I) = I - 1
        XHEU(I) = 0
   30 CONTINUE
      KREL = 0
      MREL = - 1
      MRED = ZZ
      NFREER = N
      NREL = 0
      X(N) = 0
      FS = 1
      M = ZZ
      NFREE = N
      I = 0
      IFP = N
      IF ( N .LE. 2 ) GO TO 470
   40 I = FS
      JJ = IFP
C TRY TO ASSIGN  W(I) .
   50   JWB = C - W(I)
        II = X(I)
C FIRST TEST.
        IF ( W(IFP) .LE. JWB ) GO TO 90
C INSERT ONLY ITEM  I .
        MM = I
        IF ( I .LE. ZZ ) GO TO 60
        M = M + 1
        MM = M
   60   IS = X(I)
        IP = K(I)
        K(IS) = IP
        IF ( IP .GT. 0 ) GO TO 70
        FS = IS
        GO TO 80
   70   X(IP) = IS
   80   X(I) = MM
        XHEU(I) = KREL
        MP = MP + 1
        RES(MM) = C - W(I)
        ISUMR = ISUMR - W(I)
        NFREE = NFREE - 1
        GO TO 450
C FIND THE LARGEST FREE  W(JJ)  WHICH FITS INTO A BIN TOGETHER
C WITH W(I)  ( JJ .GT. I ) .
   90   J = JJ
  100   IF ( J .EQ. I ) GO TO 110
        IF ( W(J) .GT. JWB ) GO TO 120
        J = K(J)
        GO TO 100
  110   IF ( I .EQ. FS ) GO TO 120
        IJ = K(I)
        IF ( W(IJ) .GT. JWB ) GO TO 120
        JJ = X(J)
        GO TO 460
C SECOND TEST.
  120   JJ = X(J)
        IF ( W(JJ) .EQ. JWB ) GO TO 360
C FOURTH TEST.
        ISP = K(IFP)
        IF ( W(IFP) + W(ISP) .GT. JWB ) GO TO 360
        IF ( I .EQ. ISP ) GO TO 360
C FIFTH TEST.
        LLP = IFP
        TP = K(ISP)
        IF ( TP .EQ. I ) GO TO 130
        IF ( W(IFP) + W(TP) .LE. JWB ) GO TO 240
C ONLY PAIR  (ISP,IFP)  CAN BE INSERTED WITH ITEM  I .
        IF ( W(JJ) .GE. W(IFP) + W(ISP) ) GO TO 360
        IF ( W(JJ) .NE. W(ISP) ) GO TO 460
C INSERT ITEMS  I ,  JJ  AND  LLP .
  130   IF ( JJ .GT. ZZ ) GO TO 140
        IF ( KREL .LT. KINF ) KINF = KREL
  140   MM = I
        IF ( I .LE. ZZ ) GO TO 150
        M = M + 1
        MM = M
  150   IS = X(I)
        IP = K(I)
        K(IS) = IP
        IF ( IP .GT. 0 ) GO TO 160
        FS = IS
        GO TO 170
  160   X(IP) = IS
  170   X(I) = MM
        XHEU(I) = KREL
        MP = MP + 1
        IW = W(I) + W(JJ) + W(LLP)
        RES(MM) = C - IW
        ISUMR = ISUMR - IW
        IF ( JJ .EQ. II ) II = X(II)
        IF ( II .EQ. 0 ) II = IFP
        JS = X(JJ)
        JP = K(JJ)
        IF ( JP .GT. 0 ) GO TO 180
        FS = JS
        GO TO 190
  180   X(JP) = JS
  190   K(JS) = JP
        X(JJ) = MM
        XHEU(JJ) = KREL
        JJ = JS
        IF ( LLP .EQ. II ) II = X(II)
        IF ( II .EQ. 0 ) II = IFP
        LS = X(LLP)
        LP = K(LLP)
        IF ( LP .GT. 0 ) GO TO 200
        FS = LS
        GO TO 210
  200   X(LP) = LS
  210   IF ( LS .GT. 0 ) GO TO 220
        IFP = LP
        GO TO 230
  220   K(LS) = LP
  230   X(LLP) = MM
        XHEU(LLP) = KREL
        IF ( JJ .EQ. LLP ) JJ = LS
        IF ( LS .EQ. 0 ) JJ = IFP
        NFREE = NFREE - 3
        GO TO 450
C SIXTH TEST.
  240   WWW = 0
        IF ( W(IFP) + W(ISP) + W(TP) .GT. JWB ) GO TO 280
        QP = K(TP)
        IF ( QP .EQ. I ) GO TO 250
        IF ( W(IFP) + W(ISP) + W(QP) .LE. JWB ) GO TO 550
C ONLY TRIPLET  (TP,ISP,IFP)  CAN BE INSERTED WITH ITEM  I .
        IF ( W(JJ) .EQ. W(TP) ) GO TO 250
        WWW = W(IFP) + W(ISP) + W(TP)
        GO TO 280
C INSERT ITEMS  I ,  TP ,  ISP ,  IFP .
  250   IF ( TP .GT. ZZ ) GO TO 260
        IF ( KREL .LT. KINF ) KINF = KREL
  260   MM = I
        IF ( I .LE. ZZ ) GO TO 270
        M = M + 1
        MM = M
  270   CALL INSERT(I,MM,FS,X,IFP,K,JDIM)
        CALL INSERT(TP,MM,FS,X,IFP,K,JDIM)
        CALL INSERT(ISP,MM,FS,X,IFP,K,JDIM)
        LIFP = IFP
        CALL INSERT(LIFP,MM,FS,X,IFP,K,JDIM)
        XHEU(I) = KREL
        XHEU(TP) = KREL
        XHEU(ISP) = KREL
        XHEU(LIFP) = KREL
        MP = MP + 1
        IW = W(I) + W(TP) + W(ISP) + W(LIFP)
        RES(MM) = C - IW
        ISUMR = ISUMR - IW
        IF ( JJ .EQ. II) II = IFP
        JJ = IFP
        NFREE = NFREE - 4
        GO TO 450
C FIND THE BEST PAIR  (KKP,LLP)  OF TOTAL WEIGHT  WW  WHICH CAN BE
C INSERTED INTO A BIN TOGETHER WITH ITEM  I .
  280   KK = K(JJ)
        LL = IFP
        WW = 0
  290   LL1 = K(LL)
        IF ( KK. GE. LL1 ) GO TO 340
        J = KK
        JWL = JWB - W(LL)
  300   J = X(J)
        IF ( J .GE. LL ) GO TO 340
        IF ( W(J) .GT. JWL ) GO TO 300
        KK = J
        J = LL
        JWK = JWB - W(KK)
  310   J = K(J)
        IF ( J .LE. KK ) GO TO 320
        IF ( W(J) .LE. JWK ) GO TO 310
  320   LL = X(J)
        IF ( W(KK) + W(LL) .LE. WW ) GO TO 330
        WW = W(KK) + W(LL)
        KKP = KK
        LLP = LL
        IF ( W(JJ) .GE. WW ) GO TO 330
        IF ( W(JJ) .NE. W(KKP) ) GO TO 460
  330   LL = K(LL)
        GO TO 290
  340   IF ( (W(JJ) .GE. WW) .AND. (W(JJ) .GE. WWW) ) GO TO 360
        IF ( W(JJ) .NE. W(KKP) ) GO TO 460
C ITEM  JJ  IS IN THE BEST PAIR (JJ=KKP).
C CHECK WHETHER PAIR  (KKP,LLP)  DOMINATES ALL THE FEASIBLE PAIRS.
        JLL = K(LLP)
        JJLL = K(JLL)
        IF ( JJLL .LE. JJ ) GO TO 350
        IF ( W(JLL) + W(JJLL) .LE. JWB ) GO TO 460
  350   IF ( WWW .EQ. 0 ) GO TO 130
C CHECK WHETHER  PAIR  (KKP,LLP)  DOMINATES TRIPLET  (TP,ISP,IFP) .
        IF ( W(LLP) .LT. W(ISP) ) GO TO 460
        IF ( W(KKP) + W(LLP) .GE. WWW ) GO TO 130
        GO TO 460
C INSERT ITEMS  I  AND  JJ .
  360   IF ( JJ .GT. ZZ ) GO TO 370
        IF ( KREL .LT. KINF ) KINF = KREL
  370   MM = I
        IF ( I .LE. ZZ ) GO TO 380
        M = M + 1
        MM = M
  380   IS = X(I)
        IP = K(I)
        K(IS) = IP
        IF ( IP .GT. 0 ) GO TO 390
        FS = IS
        GO TO 400
  390   X(IP) = IS
  400   X(I) = MM
        XHEU(I) = KREL
        MP = MP + 1
        IW = W(I) + W(JJ)
        RES(MM) = C - IW
        ISUMR = ISUMR - IW
        IF ( JJ .EQ. II ) II = X(II)
        IF ( II .EQ. 0 ) II = IFP
        JS = X(JJ)
        JP = K(JJ)
        IF ( JP .GT. 0 ) GO TO 410
        FS = JS
        GO TO 420
  410   X(JP) = JS
  420   IF ( JS .GT. 0 ) GO TO 430
        IFP = JP
        GO TO 440
  430   K(JS) = JP
  440   X(JJ) = MM
        XHEU(JJ) = KREL
        JJ = JS
        IF ( JS .EQ. 0 ) JJ = IFP
        NFREE = NFREE - 2
C STOPPING TEST.
  450   IF ( NFREE .LE. 2 ) GO TO 470
  460   I = II
      IF ( I .LT. IFP ) GO TO 50
      GO TO 670
C OPTIMAL SOLUTION.
  470 IF ( NFREE .EQ. 0 ) GO TO 670
      I1 = 0
      I = FS
  480   II = X(I)
        IF ( I1 .GT. 0 ) GO TO 500
        I1 = I
        MM = I
        IF ( I .LE. ZZ ) GO TO 490
        M = M + 1
        MM = M
  490   X(I1) = MM
        XHEU(I1) = KREL
        MP = MP + 1
        RES(MM) = C - W(I1)
        ISUMR = ISUMR - W(I1)
        IF ( NFREE .EQ. 1 ) GO TO 540
        GO TO 530
  500   IF ( W(I) .GT. C - W(I1) ) GO TO 510
        IF ( I .LE. ZZ ) GO TO 510
        X(I) = MM
        XHEU(I) = KREL
        RES(MM) = RES(MM) - W(I)
        ISUMR = ISUMR - W(I)
        GO TO 540
  510   MM = I
        IF ( I .LE. ZZ ) GO TO 520
        M = M + 1
        MM = M
  520   X(I) = MM
        XHEU(I) = KREL
        MP = MP + 1
        RES(MM) = C - W(I)
        ISUMR = ISUMR - W(I)
        GO TO 540
  530   I = II
      GO TO 480
  540 NFREE = 0
      FS = 0
      GO TO 670
C ONLY SECOND TEST FOR THE REMAINING ITEMS.
  550 IF ( KREL .GT. 0 ) GO TO 670
      J = K(IFP)
      JB2 = (C + 1)/2
  560 I = II
      IF ( I .EQ. IFP ) GO TO 670
      IF ( W(I) .LT. JB2 ) GO TO 670
      II = X(I)
      JWB = C - W(I)
  570 IF ( J .EQ. I ) GO TO 670
      IF ( W(J) .GT. JWB ) GO TO 560
      IF ( W(J) .EQ. JWB ) GO TO 580
      J = K(J)
      GO TO 570
C INSERT ITEMS  I  AND  J .
  580 IF ( J .GT. ZZ ) GO TO 590
      IF ( KREL .LT. KINF ) KINF = KREL
  590 MM = I
      IF ( I .LE. ZZ ) GO TO 600
      M = M + 1
      MM = M
  600 IS = X(I)
      IP = K(I)
      K(IS) = IP
      IF ( IP .GT. 0 ) GO TO 610
      FS = IS
      GO TO 620
  610 X(IP) = IS
  620 X(I) = MM
      XHEU(I) = KREL
      MP = MP + 1
      IW = W(I) + W(J)
      RES(MM) = C - IW
      ISUMR = ISUMR - IW
      IF ( J .EQ. II ) II = X(II)
      IF ( II .EQ. 0 ) II = IFP
      JS = X(J)
      JP = K(J)
      IF ( JP .GT. 0 ) GO TO 630
      FS = JS
      GO TO 640
  630 X(JP) = JS
  640 IF ( JS .GT. 0 ) GO TO 650
      IFP = JP
      GO TO 660
  650 K(JS) = JP
  660 X(J) = MM
      XHEU(J) = KREL
      J = JP
      NFREE = NFREE - 2
      IF ( NFREE .LE. 2 ) GO TO 470
      GO TO 560
  670 IF ( NFREE .EQ. 0 ) GO TO 680
      IF ( W(FS) + W(IFP) .GT. C ) GO TO 40
C CHECK WHETHER ALL THE PREVIOUSLY RELAXED ITEMS CAN BE INSERTED
C INTO THE BINS CURRENTLY USED BY THE REDUCED ITEMS.
  680 IF ( KINF .LE. N ) GO TO 710
      IF ( NREL .EQ. 0 ) GO TO 700
      IX = REL(NREL)
C CHECK WHETHER THE LAST RELAXED ITEM  IX  CAN BE INSERTED INTO
C A FEASIBLE BIN.
      IW = W(IX)
      MIN = C + 1
      DO 690 L=1,M
        IF ( IW .GT. RES(L) ) GO TO 690
        IF ( RES(L) .GE. MIN ) GO TO 690
        MIN = RES(L)
        LMIN = L
  690 CONTINUE
C IF  MIN .GT. C  THEN ITEM  IX  CANNOT BE INSERTED IN ANY CURRENT BIN.
      IF ( MIN .GT. C ) GO TO 710
C INSERT ITEM  IX  IN BIN  LMIN  AND CONSIDER THE NEXT RELAXED ITEM.
      X(IX) = LMIN
      XHEU(IX) = KREL
      RES(LMIN) = RES(LMIN) - IW
      NREL = NREL - 1
      GO TO 680
C ALL THE RELAXED ITEMS HAVE BEEN INSERTED IN THE CURRENT BINS.
  700 MREL = KREL
      MRED = M
      NFREER = NFREE
C COMPUTE THE CURRENT LOWER BOUND  LBC .
  710 LBC = MP + ( ISUMR + C - 1 )/C
      IF ( KREL .EQ. 0 .AND. KINF .LE. N ) GO TO 780
      IF ( LB .LT. LBC ) LB = LBC
      IF ( LB .GE. UB ) RETURN
      IF ( NFREE .EQ. 0 ) GO TO 720
C RELAX THE LAST ITEM  IFP .
      KREL = KREL + 1
      XHEU(IFP) = N + 1
      ISUMR = ISUMR - W(IFP)
      NREL = NREL + 1
      REL(NREL) = IFP
      IFP = K(IFP)
      X(IFP) = 0
      NFREE = NFREE - 1
      GO TO 40
C DEFINE  LB ,  X(I) ,  XHEU(I)  AND  NUB .
  720 IF ( LB .LT. MP ) LB = MP
      NUB = M
      DO 770 I=1,N
        IF ( XHEU(I) .GT. MREL ) GO TO 730
        XHEU(I) = X(I)
        GO TO 770
C ITEM  I  HAS NOT BEEN REDUCED IN A VALID WAY.
  730   IF ( XHEU(I) .GT. N ) GO TO 740
        XHEU(I) = X(I)
        X(I) = 0
        GO TO 770
C TRY TO INSERT THE RELAXED ITEM  I  INTO AN EXISTING BIN
C OR INITIALIZE A NEW BIN.
  740   IW = W(I)
        MIN = C + 1
        DO 750 L=1,NUB
          IF ( IW .GT. RES(L) ) GO TO 750
          IF ( MIN .LE. RES(L) ) GO TO 750
          MIN = RES(L)
          LMIN = L
  750   CONTINUE
        IF ( MIN .GT. C ) GO TO 760
        XHEU(I) = LMIN
        X(I) = 0
        RES(LMIN) = RES(LMIN) - W(I)
        GO TO 770
  760   NUB = NUB + 1
        XHEU(I) = NUB
        X(I) = 0
        RES(NUB) = C - W(I)
  770 CONTINUE
      RETURN
  780 NFREER = - 1
      RETURN
      END
      SUBROUTINE MWFDS(N,W,C,M,KK,K,LLB,X,JDIM)
C
C PERFORM A MODIFIED WORST-FIT DECREASING HEURISTIC.
C FOR LOCAL USE WITH CURRENT SOLUTION GIVEN.
C TIME COMPLEXITY  O(N**2) .
C
      INTEGER W(JDIM),C,K(JDIM),X(JDIM),KK(JDIM)
      DO 10 J=1,M
        K(J) = KK(J)
   10 CONTINUE
      I1 = 1
      IF ( M .GE. LLB ) GO TO 30
      M1 = M + 1
      I = 0
      DO 20 J=M1,LLB
        I = I + 1
        K(J) = C - W(I)
        X(I) = J
   20 CONTINUE
      M = LLB
      I1 = I + 1
      IF ( I1 .GT. N ) RETURN
   30 DO 60 I=I1,N
C INSERT THE NEXT ITEM.
        IWI = W(I)
        MAXRES = - 1
        DO 40 J=1,M
          KRES = K(J) - IWI
          IF ( KRES .LT. 0 ) GO TO 40
          IF ( KRES .LE. MAXRES ) GO TO 40
          MAXRES = KRES
          JM = J
   40   CONTINUE
        IF ( MAXRES .GE. 0 ) GO TO 50
C INITIALIZE A NEW BIN.
        M = M + 1
        K(M) = C - IWI
        X(I) = M
        GO TO 60
C INSERT THE ITEM INTO AN OLD BIN.
   50   K(JM) = K(JM) - IWI
        X(I) = JM
   60 CONTINUE
      RETURN
      END
      SUBROUTINE RESTOR(K,ZZ,N,C,KFIX,FIXIT,W,X,R,LASTW,JDIM)
C
C RE-STORE THE SITUATION PRECEDING THE REDUCTION OF LEVEL  K
C AND UPDATE  LASTW .
C
      INTEGER KFIX(JDIM),FIXIT(JDIM),W(JDIM),X(JDIM),R(JDIM)
      INTEGER ZZ,C
      I = KFIX(K)
      KFIX(K) = 0
   10 NEXT = FIXIT(I)
      FIXIT(I) = 0
      J = X(I)
      R(J) = R(J) + W(I)
      IF ( NEXT .EQ. (- 1) ) GO TO 20
      I = NEXT
      GO TO 10
   20 IF ( R(ZZ) .LT. C ) GO TO 30
      ZZ = ZZ - 1
      GO TO 20
   30 DO 40 II=1,N
        I = N - II + 1
        IF ( FIXIT(I) .EQ. 0 ) GO TO 50
   40 CONTINUE
C NO EXIT THIS WAY.
   50 LASTW = W(I)
      RETURN
      END
      SUBROUTINE SEARCH(N,W,R,NL,JDIM)
C
C GIVEN  W(1),...,W(N) , SORTED BY DECREASING VALUES, FIND, THROUGH
C BINARY SEARCH, THE LARGEST  NL  SUCH THAT  W(NL) .GT. R .
C
      INTEGER W(JDIM)
      N1 = 1
      N2 = N
   10 IF ( N2 - N1 .LE. 2 ) GO TO 30
      NL = (N1 + N2)/2
      IF ( FLOAT(W(NL)) .GT. R ) GO TO 20
      N2 = NL - 1
      GO TO 10
   20 N1 = NL + 1
      GO TO 10
   30 DO 40 I=N1,N2
        IF ( FLOAT(W(I)) .LE. R ) GO TO 50
   40 CONTINUE
      NL = N2
      RETURN
   50 NL = I - 1
      RETURN
      END
      SUBROUTINE SORTI2(N,A,V,JDA)
C
C SORT THE INTEGER ARRAY A BY INCREASING VALUES (DERIVED FROM
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
      IF ( A(KI) .LE. T ) GO TO 30
      V(IJ) = KI
      V(I) = IV
      IV = V(IJ)
      T = A(IV)
   30 L = J
      KI = V(J)
      IF ( A(KI) .GE. T ) GO TO 50
      V(IJ) = KI
      V(J) = IV
      IV = V(IJ)
      T = A(IV)
      KI = V(I)
      IF ( A(KI) .LE. T ) GO TO 50
      V(IJ) = KI
      V(I) = IV
      IV = V(IJ)
      T = A(IV)
      GO TO 50
   40 V(L) = V(K)
      V(K) = IVT
   50 L = L - 1
      KI = V(L)
      IF ( A(KI) .GT. T ) GO TO 50
      IVT = KI
   60 K = K + 1
      KI = V(K)
      IF ( A(KI) .LT. T ) GO TO 60
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
      IF ( A(KI) .LE. T ) GO TO 100
      K = I
  110 V(K+1) = V(K)
      K = K - 1
      KI = V(K)
      IF ( T .LT. A(KI) ) GO TO 110
      V(K+1) = IV
      GO TO 100
      END
      SUBROUTINE UPDATE(N,Z,XSTAR,NA,M,X,WB,XRED,JDIM)
C
C UPDATE THE OPTIMAL SOLUTION AFTER A LOCAL HEURISTIC.
C
      INTEGER XSTAR(JDIM),X(JDIM),WB(JDIM),XRED(JDIM),Z
      DO 10 I=1,N
        XSTAR(I) = X(I)
   10 CONTINUE
      DO 20 II=1,NA
        IXRII = XRED(II)
        XSTAR(IXRII) = WB(II)
   20 CONTINUE
      Z = M
      RETURN
      END
