C     ///////////////////////////////////////////////////////////////
      SUBROUTINE FILL3A(NX,NY,F3)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,N,P
      DIMENSION F3(NX,NX,NX)
      DO 5  P=1,NY
      DO 6  N=1,P-1
      DO 7  M=1,N-1
      F3(N,M,P)=F3(M,N,P)
      F3(N,P,M)=F3(M,N,P)
      F3(M,P,N)=F3(M,N,P)
      F3(P,M,N)=F3(M,N,P)
   7  F3(P,N,M)=F3(M,N,P)
      F3(N,P,N)=F3(N,N,P)
   6  F3(P,N,N)=F3(N,N,P)
      DO 8  M=1,P-1
      F3(P,M,P)=F3(M,P,P)
   8  F3(P,P,M)=F3(M,P,P)
   5  CONTINUE
      RETURN
      END
