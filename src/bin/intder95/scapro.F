C     //////////////////////////////////////////////////////////////
      SUBROUTINE SCAPRO(U,V,D)
      REAL*8 D,U(3),V(3),ZERO
      PARAMETER(ZERO=0.0D0)
      D=ZERO
      DO 5  I=1,3
 5    D=D+U(I)*V(I)
      RETURN
      END
