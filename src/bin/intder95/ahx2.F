C     ////////////////////////////////////////////////////////////
      SUBROUTINE AHX2(NC,NSX,L1,L2,H11,A,X)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NC,NC),X(NC,NC),H11(3,3)
      DO 10  N=1,NSX
      DO 10  M=1,N
          DO 15  I=1,3
          DO 15  J=1,3
          W1=(A(L1+I,M)-A(L2+I,M))*(A(L1+J,N)-A(L2+J,N))
  15      X(M,N)=X(M,N)+W1*H11(I,J)
  10  CONTINUE
      DO 20  N=1,NSX
      DO 20  M=1,N
  20  X(N,M)=X(M,N)
      RETURN
      END
