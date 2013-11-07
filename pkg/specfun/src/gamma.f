
        SUBROUTINE CGAMMA(X,Y,KF,GR,GI)
C
C       =========================================================
C       Purpose: Compute the gamma function G(z) or ln[G(z)]
C                for a complex argument
C       Input :  x  --- Real part of z
C                y  --- Imaginary part of z
C                KF --- Function code
C                       KF=0 for ln[G(z)]
C                       KF=1 for G(z)
C       Output:  GR --- Real part of ln[G(z)] or G(z)
C                GI --- Imaginary part of ln[G(z)] or G(z)
C       ========================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        PI=3.141592653589793D0
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        IF (Y.EQ.0.0D0.AND.X.EQ.INT(X).AND.X.LE.0.0D0) THEN
           GR=1.0D+300
           GI=0.0D0
           RETURN
        ELSE IF (X.LT.0.0D0) THEN
           X1=X
           Y1=Y
           X=-X
           Y=-Y
        ENDIF
        X0=X
        IF (X.LE.7.0) THEN
           NA=INT(7-X)
           X0=X+NA
        ENDIF
        Z1=DSQRT(X0*X0+Y*Y)
        TH=DATAN(Y/X0)
        GR=(X0-.5D0)*DLOG(Z1)-TH*Y-X0+0.5D0*DLOG(2.0D0*PI)
        GI=TH*(X0-0.5D0)+Y*DLOG(Z1)-Y
        DO 10 K=1,10
           T=Z1**(1-2*K)
           GR=GR+A(K)*T*DCOS((2.0D0*K-1.0D0)*TH)
10         GI=GI-A(K)*T*DSIN((2.0D0*K-1.0D0)*TH)
        IF (X.LE.7.0) THEN
           GR1=0.0D0
           GI1=0.0D0
           DO 15 J=0,NA-1
              GR1=GR1+.5D0*DLOG((X+J)**2+Y*Y)
15            GI1=GI1+DATAN(Y/(X+J))
           GR=GR-GR1
           GI=GI-GI1
        ENDIF
        IF (X1.LT.0.0D0) THEN
           Z1=DSQRT(X*X+Y*Y)
           TH1=DATAN(Y/X)
           SR=-DSIN(PI*X)*DCOSH(PI*Y)
           SI=-DCOS(PI*X)*DSINH(PI*Y)
           Z2=DSQRT(SR*SR+SI*SI)
           TH2=DATAN(SI/SR)
           IF (SR.LT.0.0D0) TH2=PI+TH2
           GR=DLOG(PI/(Z1*Z2))-GR
           GI=-TH1-TH2-GI
           X=X1
           Y=Y1
        ENDIF
        IF (KF.EQ.1) THEN
           G0=DEXP(GR)
           GR=G0*DCOS(GI)
           GI=G0*DSIN(GI)
        ENDIF
        RETURN
        END


        SUBROUTINE GAMMA(X,GA)
C
C       ==================================================
C       Purpose: Compute the gamma function G(x)
C       Input :  x  --- Argument of G(x)
C                       ( x is not equal to 0,-1,-2,... )
C       Output:  GA --- G(x)
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION G(26)
        PI=3.141592653589793D0
        IF (X.EQ.INT(X)) THEN
           IF (X.GT.0.0D0) THEN
              GA=1.0D0
              M1=X-1
              DO 10 K=2,M1
10               GA=GA*K
           ELSE
              GA=1.0D+300
           ENDIF
        ELSE
           IF (DABS(X).GT.1.0D0) THEN
              Z=DABS(X)
              M=INT(Z)
              R=1.0D0
              DO 15 K=1,M
15               R=R*(Z-K)
              Z=Z-M
           ELSE
              Z=X
           ENDIF
           DATA G/1.0D0,0.5772156649015329D0,
     &          -0.6558780715202538D0, -0.420026350340952D-1,
     &          0.1665386113822915D0,-.421977345555443D-1,
     &          -.96219715278770D-2, .72189432466630D-2,
     &          -.11651675918591D-2, -.2152416741149D-3,
     &          .1280502823882D-3, -.201348547807D-4,
     &          -.12504934821D-5, .11330272320D-5,
     &          -.2056338417D-6, .61160950D-8,
     &          .50020075D-8, -.11812746D-8,
     &          .1043427D-9, .77823D-11,
     &          -.36968D-11, .51D-12,
     &          -.206D-13, -.54D-14, .14D-14, .1D-15/
           GR=G(26)
           DO 20 K=25,1,-1
20            GR=GR*Z+G(K)
           GA=1.0D0/(GR*Z)
           IF (DABS(X).GT.1.0D0) THEN
              GA=GA*R
              IF (X.LT.0.0D0) GA=-PI/(X*GA*DSIN(PI*X))
           ENDIF
        ENDIF
        RETURN
        END


        SUBROUTINE LGAMMA(KF,X,GL)
C
C       ==================================================
C       Purpose: Compute gamma function G(x) or ln[G(x)]
C       Input:   x  --- Argument of G(x) ( x > 0 )
C                KF --- Function code
C                       KF=1 for G(x); KF=0 for ln[G(x)]
C       Output:  GL --- G(x) or ln[G(x)]
C       ==================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(10)
        DATA A/8.333333333333333D-02,-2.777777777777778D-03,
     &         7.936507936507937D-04,-5.952380952380952D-04,
     &         8.417508417508418D-04,-1.917526917526918D-03,
     &         6.410256410256410D-03,-2.955065359477124D-02,
     &         1.796443723688307D-01,-1.39243221690590D+00/
        X0=X
        IF (X.EQ.1.0.OR.X.EQ.2.0) THEN
           GL=0.0D0
           GO TO 20
        ELSE IF (X.LE.7.0) THEN
           N=INT(7-X)
           X0=X+N
        ENDIF
        X2=1.0D0/(X0*X0)
        XP=6.283185307179586477D0
        GL0=A(10)
        DO 10 K=9,1,-1
10         GL0=GL0*X2+A(K)
        GL=GL0/X0+0.5D0*DLOG(XP)+(X0-.5D0)*DLOG(X0)-X0
        IF (X.LE.7.0) THEN
           DO 15 K=1,N
              GL=GL-DLOG(X0-1.0D0)
15            X0=X0-1.0D0
        ENDIF
20      IF (KF.EQ.1) GL=DEXP(GL)
        RETURN
        END


        SUBROUTINE INCOG(A,X,GIN,GIM,GIP)
C
C       ===================================================
C       Purpose: Compute the incomplete gamma function
C                g(a,x), G(a,x) and P(a,x)
C       Input :  a   --- Parameter ( a > 170 )
C                x   --- Argument 
C       Output:  GIN --- g(a,x)
C                GIM --- G(a,x)
C                GIP --- P(a,x)
C       Routine called: GAMMA for computing G(x)
C       ===================================================
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        XAM=-X+A*DLOG(X)
        IF (X.EQ.0.0) THEN
           GIN=0.0
           CALL GAMMA(A,GA)
           GIM=GA
           GIP=0.0
        ELSE IF (X.LE.1.0+A) THEN
           S=1.0D0/A
           R=S
           DO 10 K=1,60
              R=R*X/(A+K)
              S=S+R
              IF (DABS(R/S).LT.1.0D-15) GO TO 15
10         CONTINUE
15         GIN=DEXP(XAM)*S
           CALL GAMMA(A,GA)
           GIP=GIN/GA
           GIM=GA-GIN
        ELSE IF (X.GT.1.0+A) THEN
           T0=0.0D0
           DO 20 K=60,1,-1
              T0=(K-A)/(1.0D0+K/(X+T0))
20         CONTINUE
           GIM=DEXP(XAM)/(X+T0)
           CALL GAMMA(A,GA)
           GIN=GA-GIM
           GIP=1.0D0-GIM/GA
        ENDIF
        END
