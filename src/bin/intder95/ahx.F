C     ////////////////////////////////////////////////////////////
      SUBROUTINE AHX(NC,NSX,SR,A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NC,NC),X(NC,NC),SR(NC,NC)
      PARAMETER(ZERO=0.0D0)
      DO 5  L=1,NSX
      DO 5  M=1,L
         X(L,M)=ZERO
         DO 5  I=1,NC
         DO 5  J=1,NC
 5       X(L,M)=X(L,M)+SR(I,J)*A(I,L)*A(J,M)
      DO 10  L=2,NSX
      DO 10  M=1,L-1
 10      X(M,L)=X(L,M)
      RETURN
      END
