C     /////////////////////////////////////////////////////////////////
      SUBROUTINE HSRY2(NC,L1,L2,H111,SR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SR(NC,NC,NC)
      DIMENSION H111(3,3,3)
      DO 10  K=1,3
      DO 10  J=1,3
      DO 10  I=1,3
      Z=H111(I,J,K)
      SR(L1+I,L1+J,L1+K)=Z
      SR(L1+I,L1+J,L2+K)=-Z
      SR(L1+I,L2+J,L1+K)=-Z
      SR(L1+I,L2+J,L2+K)=Z
      SR(L2+I,L1+J,L1+K)=-Z
      SR(L2+I,L1+J,L2+K)=Z
      SR(L2+I,L2+J,L1+K)=Z
 10   SR(L2+I,L2+J,L2+K)=-Z
      RETURN
      END
