C     //////////////////////////////////////////////////////////////
      SUBROUTINE MAT2(AA,VK)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AA(3,3),VK(3)
      PARAMETER(ZERO=0.0D0)
      DO 10  I=1,3
  10  AA(I,I)=ZERO
      AA(1,2)=VK(1)
      AA(2,1)=-VK(1)
      AA(1,3)=VK(2)
      AA(3,1)=-VK(2)
      AA(2,3)=VK(3)
      AA(3,2)=-VK(3)
      END
