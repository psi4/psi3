C     ///////////////////////////////////////////////////////////////
      SUBROUTINE FILL3B(NX,NY,F3)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,N,P
      DIMENSION F3(NX,NX,NX)
      DO 5   M=1,NY
      DO 10  N=1,M-1
      DO 15  P=1,N-1
      F3(N,M,P)=F3(M,N,P)
      F3(N,P,M)=F3(M,N,P)
      F3(M,P,N)=F3(M,N,P)
      F3(P,M,N)=F3(M,N,P)
  15  F3(P,N,M)=F3(M,N,P)
      F3(N,M,N)=F3(M,N,N)
  10  F3(N,N,M)=F3(M,N,N)
      DO 20  P=1,M-1
      F3(M,P,M)=F3(M,M,P)
  20  F3(P,M,M)=F3(M,M,P)
   5  CONTINUE
      RETURN
      END
