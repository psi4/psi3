C     /////////////////////////////////////////////////////////////
      SUBROUTINE AHY2(NC,NSX,L1,L2,H111,A,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,NSX
      DIMENSION A(NC,NC),Y(NC,NC,NC)
      DIMENSION H111(3,3,3)
      DO 10  P=1,NSX
      DO 10  N=1,P
      DO 10  M=1,N
         DO 12  I=1,3
         DO 12  J=1,3
         DO 12  K=1,3
      W1=A(L1+J,N)*(A(L1+K,P)-A(L2+K,P))-A(L2+J,N)*(A(L1+K,P)-A(L2+K,P))
      W1=(A(L1+I,M)-A(L2+I,M))*W1
 12      Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)
 10   CONTINUE
      RETURN
      END
