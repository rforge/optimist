        SUBROUTINE BETA(P,Q,BT)
C
C       ==========================================
C       Purpose: Compute the beta function B(p,q)
C       Input :  p  --- Parameter  ( p > 0 )
C                q  --- Parameter  ( q > 0 )
C       Output:  BT --- B(p,q)
C       Routine called: GAMMA for computing G(x)
C       ==========================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CALL GAMMA(P,GP)
        CALL GAMMA(Q,GQ)
        PPQ=P+Q
        CALL GAMMA(PPQ,GPQ)
        BT=GP*GQ/GPQ
        RETURN
        END


        SUBROUTINE INCOB(A,B,X,BIX)
C
C       ========================================================
C       Purpose: Compute the incomplete beta function Ix(a,b)
C       Input :  a --- Parameter
C                b --- Parameter
C                x --- Argument ( 0 < x < 1 )
C       Output:  BIX --- Ix(a,b)
C       Routine called: BETA for computing beta function B(p,q)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION DK(51),FK(51)
        S0=(A+1.0D0)/(A+B+2.0D0)
        CALL BETA(A,B,BT)
        IF (X.LE.S0) THEN
           DO 10 K=1,20
10            DK(2*K)=K*(B-K)*X/(A+2.0D0*K-1.0D0)/(A+2.0D0*K)
           DO 15 K=0,20
15            DK(2*K+1)=-(A+K)*(A+B+K)*X/(A+2.D0*K)/(A+2.0*K+1.0)
           T1=0.0D0
           DO 20 K=20,1,-1
20            T1=DK(K)/(1.0D0+T1)
           TA=1.0D0/(1.0D0+T1)
           BIX=X**A*(1.0D0-X)**B/(A*BT)*TA
        ELSE
           DO 25 K=1,20
25            FK(2*K)=K*(A-K)*(1.0D0-X)/(B+2.*K-1.0)/(B+2.0*K)
           DO 30 K=0,20
30            FK(2*K+1)=-(B+K)*(A+B+K)*(1.D0-X)/
     &                   (B+2.D0*K)/(B+2.D0*K+1.D0)
           T2=0.0D0
           DO 35 K=20,1,-1
35            T2=FK(K)/(1.0D0+T2)
           TB=1.0D0/(1.0D0+T2)
           BIX=1.0D0-X**A*(1.0D0-X)**B/(B*BT)*TB
        ENDIF
        RETURN
        END
