C     ////////////////////////////////////////////////////////////
C     /////////////////////////////////////////////////////////////////
      SUBROUTINE H4TH1(NAD,K1,K2,XA,H1111)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),H11(3,3),H111(3,3,3),H1111(3,3,3,3)
      PARAMETER(ONE=1.0D0)
      CALL VECT1(NAD,K1,K2,V1,XA,T21)
      CALL HIJS1(NAD,K1,K2,XA,H11)
      CALL HIJKS1(NAD,K1,K2,XA,H111)
      DO 5  L=1,3
      DO 5  K=1,L
      DO 5  J=1,K
      DO 5  I=1,J
      F=H11(I,L)*H11(K,J)+H11(J,L)*H11(K,I)+H11(I,J)*H11(K,L)
      F=F+V1(I)*H111(J,K,L)+V1(J)*H111(I,K,L)+V1(K)*H111(I,J,L)
      F=F+V1(L)*H111(I,J,K)
  5   H1111(I,J,K,L)=-F/T21
      CALL FILL4A(3,3,H1111)
      RETURN
      END
