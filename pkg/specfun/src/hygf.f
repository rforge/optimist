
        SUBROUTINE HYGFX(A,B,C,X,HF)
C
C       ====================================================
C       Purpose: Compute hypergeometric function F(a,b,c,x)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter, c <> 0,-1,-2,...
C                x --- Argument   ( x < 1 )
C       Output:  HF --- F(a,b,c,x)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ====================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        LOGICAL L0,L1,L2,L3,L4,L5
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        L0=C.EQ.INT(C).AND.C.LT.0.0
        L1=1.0D0-X.LT.1.0D-15.AND.C-A-B.LE.0.0
        L2=A.EQ.INT(A).AND.A.LT.0.0
        L3=B.EQ.INT(B).AND.B.LT.0.0
        L4=C-A.EQ.INT(C-A).AND.C-A.LE.0.0
        L5=C-B.EQ.INT(C-B).AND.C-B.LE.0.0
        IF (L0.OR.L1) THEN
           RETURN
        ENDIF
        EPS=1.0D-15
        IF (X.GT.0.95) EPS=1.0D-8
        IF (X.EQ.0.0.OR.A.EQ.0.0.OR.B.EQ.0.0) THEN
           HF=1.0D0
           RETURN
        ELSE IF (1.0D0-X.EQ.EPS.AND.C-A-B.GT.0.0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           HF=GC*GCAB/(GCA*GCB)
           RETURN
        ELSE IF (1.0D0+X.LE.EPS.AND.DABS(C-A+B-1.0).LE.EPS) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0-B,G2)
           CALL GAMMA(0.5D0+0.5*A,G3)
           HF=G0*G1/(G2*G3)
           RETURN
        ELSE IF (L2.OR.L3) THEN
           IF (L2) NM=INT(ABS(A))
           IF (L3) NM=INT(ABS(B))
           HF=1.0D0
           R=1.0D0
           DO 10 K=1,NM
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
10            HF=HF+R
           RETURN
        ELSE IF (L4.OR.L5) THEN
           IF (L4) NM=INT(ABS(C-A))
           IF (L5) NM=INT(ABS(C-B))
           HF=1.0D0
           R=1.0D0
           DO 15 K=1,NM
              R=R*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*X
15            HF=HF+R
           HF=(1.0D0-X)**(C-A-B)*HF
           RETURN
        ENDIF
        AA=A
        BB=B
        X1=X
        IF (X.LT.0.0D0) THEN
           X=X/(X-1.0D0)
           IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN
              A=BB
              B=AA
           ENDIF
           B=C-B
        ENDIF
        IF (X.GE.0.75D0) THEN
           GM=0.0D0
           IF (DABS(C-A-B-INT(C-A-B)).LT.1.0D-15) THEN
              M=INT(C-A-B)
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A+M,GAM)
              CALL GAMMA(B+M,GBM)
              CALL PSI(A,PA)
              CALL PSI(B,PB)
              IF (M.NE.0) GM=1.0D0
              DO 30 J=1,ABS(M)-1
30               GM=GM*J
              RM=1.0D0
              DO 35 J=1,ABS(M)
35               RM=RM*J
              F0=1.0D0
              R0=1.0D0
              R1=1.0D0
              SP0=0.D0
              SP=0.0D0
              IF (M.GE.0) THEN
                 C0=GM*GC/(GAM*GBM)
                 C1=-GC*(X-1.0D0)**M/(GA*GB*RM)
                 DO 40 K=1,M-1
                    R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(K-M))*(1.0-X)
40                  F0=F0+R0
                 DO 45 K=1,M
45                  SP0=SP0+1.0D0/(A+K-1.0)+1.0/(B+K-1.0)-1.0/K
                 F1=PA+PB+SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 55 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 50 J=1,M
50                     SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0))+1.0/
     &                    (B+J+K-1.0)
                    RP=PA+PB+2.0D0*EL+SP+SM+DLOG(1.0D0-X)
                    R1=R1*(A+M+K-1.0D0)*(B+M+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 60
55                  HW=F1
60               HF=F0*C0+F1*C1
              ELSE IF (M.LT.0) THEN
                 M=-M
                 C0=GM*GC/(GA*GB*(1.0D0-X)**M)
                 C1=-(-1)**M*GC/(GAM*GBM*RM)
                 DO 65 K=1,M-1
                    R0=R0*(A-M+K-1.0D0)*(B-M+K-1.0)/(K*(K-M))*(1.0-X)
65                  F0=F0+R0
                 DO 70 K=1,M
70                  SP0=SP0+1.0D0/K
                 F1=PA+PB-SP0+2.0D0*EL+DLOG(1.0D0-X)
                 DO 80 K=1,250
                    SP=SP+(1.0D0-A)/(K*(A+K-1.0))+(1.0-B)/(K*(B+K-1.0))
                    SM=0.0D0
                    DO 75 J=1,M
75                     SM=SM+1.0D0/(J+K)
                    RP=PA+PB+2.0D0*EL+SP-SM+DLOG(1.0D0-X)
                    R1=R1*(A+K-1.0D0)*(B+K-1.0)/(K*(M+K))*(1.0-X)
                    F1=F1+R1*RP
                    IF (DABS(F1-HW).LT.DABS(F1)*EPS) GO TO 85
80                  HW=F1
85               HF=F0*C0+F1*C1
              ENDIF
           ELSE
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              CALL GAMMA(C-A-B,GCAB)
              CALL GAMMA(A+B-C,GABC)
              C0=GC*GCAB/(GCA*GCB)
              C1=GC*GABC/(GA*GB)*(1.0D0-X)**(C-A-B)
              HF=0.0D0
              R0=C0
              R1=C1
              DO 90 K=1,250
                 R0=R0*(A+K-1.0D0)*(B+K-1.0)/(K*(A+B-C+K))*(1.0-X)
                 R1=R1*(C-A+K-1.0D0)*(C-B+K-1.0)/(K*(C-A-B+K))
     &              *(1.0-X)
                 HF=HF+R0+R1
                 IF (DABS(HF-HW).LT.DABS(HF)*EPS) GO TO 95
90               HW=HF
95            HF=HF+C0+C1
           ENDIF
        ELSE
           A0=1.0D0
           IF (C.GT.A.AND.C.LT.2.0D0*A.AND.
     &         C.GT.B.AND.C.LT.2.0D0*B) THEN
              A0=(1.0D0-X)**(C-A-B)
              A=C-A
              B=C-B
           ENDIF
           HF=1.0D0
           R=1.0D0
           DO 100 K=1,250
              R=R*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*X
              HF=HF+R
              IF (DABS(HF-HW).LE.DABS(HF)*EPS) GO TO 105
100           HW=HF
105        HF=A0*HF
        ENDIF
        IF (X1.LT.0.0D0) THEN
           X=X1
           C0=1.0D0/(1.0D0-X)**AA
           HF=C0*HF
        ENDIF
        A=AA
        B=BB
C        IF (K.GT.120) WRITE(*,115)
115     FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END


        SUBROUTINE HYGFZ(A,B,C,Z,ZHF)
C
C       ======================================================
C       Purpose: Compute the hypergeometric function for a 
C                complex argument, F(a,b,c,z)
C       Input :  a --- Parameter
C                b --- Parameter
C                c --- Parameter,  c <> 0,-1,-2,...
C                z --- Complex argument
C       Output:  ZHF --- F(a,b,c,z)
C       Routines called:
C            (1) GAMMA for computing gamma function
C            (2) PSI for computing psi function
C       ======================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Y)
        IMPLICIT COMPLEX *16 (Z)
        LOGICAL L0,L1,L2,L3,L4,L5,L6
        X=REAL(Z)
        Y=DIMAG(Z)
        EPS=1.0D-15
        L0=C.EQ.INT(C).AND.C.LT.0.0D0
        L1=DABS(1.0D0-X).LT.EPS.AND.Y.EQ.0.0D0.AND.C-A-B.LE.0.0D0
        L2=CDABS(Z+1.0D0).LT.EPS.AND.DABS(C-A+B-1.0D0).LT.EPS
        L3=A.EQ.INT(A).AND.A.LT.0.0D0
        L4=B.EQ.INT(B).AND.B.LT.0.0D0
        L5=C-A.EQ.INT(C-A).AND.C-A.LE.0.0D0
        L6=C-B.EQ.INT(C-B).AND.C-B.LE.0.0D0
        AA=A
        BB=B
        A0=CDABS(Z)
        IF (A0.GT.0.95D0) EPS=1.0D-8
        PI=3.141592653589793D0
        EL=.5772156649015329D0
        IF (L0.OR.L1) THEN
C          WRITE(*,*)'The hypergeometric series is divergent'
           RETURN
        ENDIF
        IF (A0.EQ.0.0D0.OR.A.EQ.0.0D0.OR.B.EQ.0.0D0) THEN
           ZHF=(1.0D0,0.0D0)
        ELSE IF (Z.EQ.1.0D0.AND.C-A-B.GT.0.0D0) THEN
           CALL GAMMA(C,GC)
           CALL GAMMA(C-A-B,GCAB)
           CALL GAMMA(C-A,GCA)
           CALL GAMMA(C-B,GCB)
           ZHF=GC*GCAB/(GCA*GCB)
        ELSE IF (L2) THEN
           G0=DSQRT(PI)*2.0D0**(-A)
           CALL GAMMA(C,G1)
           CALL GAMMA(1.0D0+A/2.0D0-B,G2)
           CALL GAMMA(0.5D0+0.5D0*A,G3)
           ZHF=G0*G1/(G2*G3)
        ELSE IF (L3.OR.L4) THEN
           IF (L3) NM=INT(ABS(A))
           IF (L4) NM=INT(ABS(B))
           ZHF=(1.0D0,0.0D0)
           ZR=(1.0D0,0.0D0)
           DO 10 K=1,NM
              ZR=ZR*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*Z
10            ZHF=ZHF+ZR
        ELSE IF (L5.OR.L6) THEN
           IF (L5) NM=INT(ABS(C-A))
           IF (L6) NM=INT(ABS(C-B))
           ZHF=(1.0D0,0.0D0)
           ZR=(1.0D0,0.0D0)
           DO 15 K=1,NM
              ZR=ZR*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*Z
15            ZHF=ZHF+ZR
           ZHF=(1.0D0-Z)**(C-A-B)*ZHF
        ELSE IF (A0.LE.1.0D0) THEN
           IF (X.LT.0.0D0) THEN
              Z1=Z/(Z-1.0D0)
              IF (C.GT.A.AND.B.LT.A.AND.B.GT.0.0) THEN  
                 A=BB
                 B=AA
              ENDIF
              ZC0=1.0D0/((1.0D0-Z)**A)
              ZHF=(1.0D0,0.0D0)
              ZR0=(1.0D0,0.0D0)
              DO 20 K=1,500
                 ZR0=ZR0*(A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C+K-1.0D0))*Z1
                 ZHF=ZHF+ZR0
                 IF (CDABS(ZHF-ZW).LT.CDABS(ZHF)*EPS) GO TO 25
20               ZW=ZHF
25            ZHF=ZC0*ZHF
           ELSE IF (A0.GE.0.90D0) THEN
              GM=0.0D0
              MCAB=INT(C-A-B+EPS*DSIGN(1.0D0,C-A-B))
              IF (DABS(C-A-B-MCAB).LT.EPS) THEN
                 M=INT(C-A-B)
                 CALL GAMMA(A,GA)
                 CALL GAMMA(B,GB)
                 CALL GAMMA(C,GC)
                 CALL GAMMA(A+M,GAM)
                 CALL GAMMA(B+M,GBM)
                 CALL PSI(A,PA)
                 CALL PSI(B,PB)
                 IF (M.NE.0) GM=1.0D0
                 DO 30 J=1,ABS(M)-1
30                  GM=GM*J
                 RM=1.0D0
                 DO 35 J=1,ABS(M)
35                  RM=RM*J
                 ZF0=(1.0D0,0.0D0)
                 ZR0=(1.0D0,0.0D0)
                 ZR1=(1.0D0,0.0D0)
                 SP0=0.D0
                 SP=0.0D0
                 IF (M.GE.0) THEN
                    ZC0=GM*GC/(GAM*GBM)
                    ZC1=-GC*(Z-1.0D0)**M/(GA*GB*RM)
                    DO 40 K=1,M-1
                       ZR0=ZR0*(A+K-1.D0)*(B+K-1.D0)/(K*(K-M))*(1.D0-Z)
40                     ZF0=ZF0+ZR0
                    DO 45 K=1,M
45                     SP0=SP0+1.0D0/(A+K-1.0D0)+1.0/(B+K-1.0D0)-1.D0/K
                    ZF1=PA+PB+SP0+2.0D0*EL+CDLOG(1.0D0-Z)
                    DO 55 K=1,500
                       SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+(1.0D0-B)/
     &                    (K*(B+K-1.0D0))
                       SM=0.0D0
                       DO 50 J=1,M
                          SM=SM+(1.0D0-A)/((J+K)*(A+J+K-1.0D0))
     &                       +1.0D0/(B+J+K-1.0D0)
50                     CONTINUE
                       ZP=PA+PB+2.0D0*EL+SP+SM+CDLOG(1.0D0-Z)
                       ZR1=ZR1*(A+M+K-1.0D0)*(B+M+K-1.0D0)/(K*(M+K))
     &                     *(1.0D0-Z)
                       ZF1=ZF1+ZR1*ZP
                       IF (CDABS(ZF1-ZW).LT.CDABS(ZF1)*EPS) GO TO 60
55                     ZW=ZF1
60                  ZHF=ZF0*ZC0+ZF1*ZC1
                 ELSE IF (M.LT.0) THEN
                    M=-M
                    ZC0=GM*GC/(GA*GB*(1.0D0-Z)**M)
                    ZC1=-(-1)**M*GC/(GAM*GBM*RM)
                    DO 65 K=1,M-1
                       ZR0=ZR0*(A-M+K-1.0D0)*(B-M+K-1.0D0)/(K*(K-M))
     &                     *(1.0D0-Z)
65                     ZF0=ZF0+ZR0
                    DO 70 K=1,M
70                     SP0=SP0+1.0D0/K
                    ZF1=PA+PB-SP0+2.0D0*EL+CDLOG(1.0D0-Z)
                    DO 80 K=1,500
                       SP=SP+(1.0D0-A)/(K*(A+K-1.0D0))+(1.0D0-B)/(K*
     &                    (B+K-1.0D0))
                       SM=0.0D0
                       DO 75 J=1,M
75                        SM=SM+1.0D0/(J+K)
                       ZP=PA+PB+2.0D0*EL+SP-SM+CDLOG(1.0D0-Z)
                       ZR1=ZR1*(A+K-1.D0)*(B+K-1.D0)/(K*(M+K))*(1.D0-Z)
                       ZF1=ZF1+ZR1*ZP
                       IF (CDABS(ZF1-ZW).LT.CDABS(ZF1)*EPS) GO TO 85
80                     ZW=ZF1
85                  ZHF=ZF0*ZC0+ZF1*ZC1
                 ENDIF
              ELSE
                 CALL GAMMA(A,GA)
                 CALL GAMMA(B,GB)
                 CALL GAMMA(C,GC)
                 CALL GAMMA(C-A,GCA)
                 CALL GAMMA(C-B,GCB)
                 CALL GAMMA(C-A-B,GCAB)
                 CALL GAMMA(A+B-C,GABC)
                 ZC0=GC*GCAB/(GCA*GCB)
                 ZC1=GC*GABC/(GA*GB)*(1.0D0-Z)**(C-A-B)
                 ZHF=(0.0D0,0.0D0)
                 ZR0=ZC0
                 ZR1=ZC1
                 DO 90 K=1,500
                    ZR0=ZR0*(A+K-1.D0)*(B+K-1.D0)/(K*(A+B-C+K))*(1.D0-Z)
                    ZR1=ZR1*(C-A+K-1.0D0)*(C-B+K-1.0D0)/(K*(C-A-B+K))
     &                  *(1.0D0-Z)
                    ZHF=ZHF+ZR0+ZR1
                    IF (CDABS(ZHF-ZW).LT.CDABS(ZHF)*EPS) GO TO 95
90                  ZW=ZHF
95               ZHF=ZHF+ZC0+ZC1
              ENDIF
           ELSE
              Z00=(1.0D0,0.0D0)
              IF (C-A.LT.A.AND.C-B.LT.B) THEN
                  Z00=(1.0D0-Z)**(C-A-B)
                  A=C-A
                  B=C-B
              ENDIF
              ZHF=(1.0D0,0.D0)
              ZR=(1.0D0,0.0D0)
              DO 100 K=1,1500
                 ZR=ZR*(A+K-1.0D0)*(B+K-1.0D0)/(K*(C+K-1.0D0))*Z
                 ZHF=ZHF+ZR
                 IF (CDABS(ZHF-ZW).LE.CDABS(ZHF)*EPS) GO TO 105
100              ZW=ZHF
105           ZHF=Z00*ZHF
           ENDIF
        ELSE IF (A0.GT.1.0D0) THEN
           MAB=INT(A-B+EPS*DSIGN(1.0D0,A-B))
           IF (DABS(A-B-MAB).LT.EPS.AND.A0.LE.1.1D0) B=B+EPS
           IF (DABS(A-B-MAB).GT.EPS) THEN
              CALL GAMMA(A,GA)
              CALL GAMMA(B,GB)
              CALL GAMMA(C,GC)
              CALL GAMMA(A-B,GAB)
              CALL GAMMA(B-A,GBA)
              CALL GAMMA(C-A,GCA)
              CALL GAMMA(C-B,GCB)
              ZC0=GC*GBA/(GCA*GB*(-Z)**A)
              ZC1=GC*GAB/(GCB*GA*(-Z)**B)
              ZR0=ZC0
              ZR1=ZC1
              ZHF=(0.0D0,0.0D0)
              DO 110 K=1,500
                 ZR0=ZR0*(A+K-1.0D0)*(A-C+K)/((A-B+K)*K*Z)
                 ZR1=ZR1*(B+K-1.0D0)*(B-C+K)/((B-A+K)*K*Z)
                 ZHF=ZHF+ZR0+ZR1
                 IF (CDABS((ZHF-ZW)/ZHF).LE.EPS) GO TO 115
110              ZW=ZHF
115           ZHF=ZHF+ZC0+ZC1
           ELSE
              IF (A-B.LT.0.0D0) THEN
                 A=BB
                 B=AA
              ENDIF
              CA=C-A
              CB=C-B
              NCA=INT(CA+EPS*DSIGN(1.0D0,CA))
              NCB=INT(CB+EPS*DSIGN(1.0D0,CB))
              IF (DABS(CA-NCA).LT.EPS.OR.DABS(CB-NCB).LT.EPS) C=C+EPS
              CALL GAMMA(A,GA)
              CALL GAMMA(C,GC)
              CALL GAMMA(C-B,GCB)
              CALL PSI(A,PA)
              CALL PSI(C-A,PCA)
              CALL PSI(A-C,PAC)
              MAB=INT(A-B+EPS)
              ZC0=GC/(GA*(-Z)**B)
              CALL GAMMA(A-B,GM)
              ZF0=GM/GCB*ZC0
              ZR=ZC0
              DO 120 K=1,MAB-1
                 ZR=ZR*(B+K-1.0D0)/(K*Z)
                 T0=A-B-K
                 CALL GAMMA(T0,G0)
                 CALL GAMMA(C-B-K,GCBK)
120              ZF0=ZF0+ZR*G0/GCBK
              IF (MAB.EQ.0) ZF0=(0.0D0,0.0D0)
              ZC1=GC/(GA*GCB*(-Z)**A)
              SP=-2.0D0*EL-PA-PCA
              DO 125 J=1,MAB
125              SP=SP+1.0D0/J
              ZP0=SP+CDLOG(-Z)
              SQ=1.0D0
              DO 130 J=1,MAB
130              SQ=SQ*(B+J-1.0D0)*(B-C+J)/J
              ZF1=(SQ*ZP0)*ZC1
              ZR=ZC1
              RK1=1.0D0
              SJ1=0.0D0
              DO 145 K=1,10000
                 ZR=ZR/Z
                 RK1=RK1*(B+K-1.0D0)*(B-C+K)/(K*K)
                 RK2=RK1
                 DO 135 J=K+1,K+MAB
135                 RK2=RK2*(B+J-1.0D0)*(B-C+J)/J
                 SJ1=SJ1+(A-1.0D0)/(K*(A+K-1.0D0))+(A-C-1.0D0)/
     &               (K*(A-C+K-1.0D0))
                 SJ2=SJ1
                 DO 140 J=K+1,K+MAB
140                 SJ2=SJ2+1.0D0/J
                 ZP=-2.0D0*EL-PA-PAC+SJ2-1.0D0/(K+A-C)
     &              -PI/DTAN(PI*(K+A-C))+CDLOG(-Z)
                 ZF1=ZF1+RK2*ZR*ZP
                 WS=CDABS(ZF1)
                 IF (DABS((WS-W0)/WS).LT.EPS) GO TO 150
145              W0=WS
150           ZHF=ZF0+ZF1
           ENDIF
        ENDIF
155     A=AA
        B=BB
C       IF (K.GT.150) WRITE(*,160)
C160    FORMAT(1X,'Warning! You should check the accuracy')
        RETURN
        END
