C     //////////////////////////////////////////////////////////////
      SUBROUTINE MAT1(EM,V)
C       EM(I,J)=(EI X EJ)*V
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION EM(3,3),V(3)
      PARAMETER(ZERO=0.0D0)
      EM(2,1)=-V(3)
      EM(3,1)=V(2)
      EM(3,2)=-V(1)
      EM(1,2)=-EM(2,1)
      EM(1,3)=-EM(3,1)
      EM(2,3)=-EM(3,2)
      DO 5  I=1,3
 5    EM(I,I)=ZERO
      RETURN
      END
