C     ///////////////////////////////////////////////////////////////
C STRE    (OPT2) (NUMB)
      SUBROUTINE VECT1(NAD,K1,K2,V1,XA,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      W=ZERO
      DO 5  I=1,3
      V1(I)=XA(K1,I)-XA(K2,I)
 5    W=W+V1(I)*V1(I)
      W=DSQRT(W)
      C=ONE/W
      DO 10  I=1,3
 10   V1(I)=V1(I)*C
      RETURN
      END
