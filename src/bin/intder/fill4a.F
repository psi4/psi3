C     ///////////////////////////////////////////////////////////////
      SUBROUTINE FILL4A(NX,NY,F4)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,N,P,Q
      DIMENSION F4(NX,NX,NX,NX)
      DO 5  Q=1,NY
      DO 5  P=1,Q
      DO 5  N=1,P
      DO 5  M=1,N
      F4(N,M,P,Q)=F4(M,N,P,Q)
      F4(N,P,M,Q)=F4(M,N,P,Q)
      F4(N,P,Q,M)=F4(M,N,P,Q)
      F4(M,P,N,Q)=F4(M,N,P,Q)
      F4(P,M,N,Q)=F4(M,N,P,Q)
      F4(P,N,M,Q)=F4(M,N,P,Q)
      F4(P,N,Q,M)=F4(M,N,P,Q)
      F4(M,P,Q,N)=F4(M,N,P,Q)
      F4(P,M,Q,N)=F4(M,N,P,Q)
      F4(P,Q,M,N)=F4(M,N,P,Q)
      F4(P,Q,N,M)=F4(M,N,P,Q)
      F4(M,N,Q,P)=F4(M,N,P,Q)
      F4(N,M,Q,P)=F4(M,N,P,Q)
      F4(N,Q,M,P)=F4(M,N,P,Q)
      F4(N,Q,P,M)=F4(M,N,P,Q)
      F4(M,Q,N,P)=F4(M,N,P,Q)
      F4(Q,M,N,P)=F4(M,N,P,Q)
      F4(Q,N,M,P)=F4(M,N,P,Q)
      F4(Q,N,P,M)=F4(M,N,P,Q)
      F4(M,Q,P,N)=F4(M,N,P,Q)
      F4(Q,M,P,N)=F4(M,N,P,Q)
      F4(Q,P,M,N)=F4(M,N,P,Q)
      F4(Q,P,N,M)=F4(M,N,P,Q)
   5  CONTINUE
      RETURN
      END
