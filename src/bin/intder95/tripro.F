C     //////////////////////////////////////////////////////////////
      SUBROUTINE TRIPRO(PROD)
C        PROD(I,J,K)=(EI X EJ) * EK
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PROD(3,3,3),RMAT(3,3),VECT(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      DO 10 K=1,3
      DO 5 N=1,3
5     VECT(N)=ZERO
      VECT(K)=ONE
      CALL MAT1(RMAT,VECT)
      DO 10 J=1,3
      DO 10 I=1,3
10    PROD(I,J,K)=RMAT(I,J)
      RETURN
      END
