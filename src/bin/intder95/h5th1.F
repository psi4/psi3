C     /////////////////////////////////////////////////////////////////
      SUBROUTINE H5TH1(NAD,K1,K2,XA,H11111)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),H11(3,3),H111(3,3,3),H1111(3,3,3,3)
      DIMENSION H11111(3,3,3,3,3)
      PARAMETER(ONE=1.0D0)
      CALL VECT1(NAD,K1,K2,V1,XA,T21)
      CALL HIJS1(NAD,K1,K2,XA,H11)
      CALL HIJKS1(NAD,K1,K2,XA,H111)
      CALL H4TH1(NAD,K1,K2,XA,H1111)
      DO 5  M=1,3
      DO 5  L=1,3
      DO 5  K=1,L
      DO 5  J=1,K
      DO 5  I=1,J
      A=H11(I,L)*H111(K,J,M)+H11(J,L)*H111(K,I,M)+H11(I,J)*H111(K,L,M)
      B=H111(I,L,M)*H11(K,J)+H111(J,L,M)*H11(K,I)+H111(I,J,M)*H11(K,L)
      C=H11(I,M)*H111(K,J,L)+H11(J,M)*H111(K,I,L)
      D=H11(K,M)*H111(I,J,L)+H11(L,M)*H111(K,I,J)
      E=V1(I)*H1111(J,K,L,M)+V1(J)*H1111(I,K,L,M)+V1(K)*H1111(I,J,L,M)
      F=V1(L)*H1111(I,J,K,M)+V1(M)*H1111(I,J,K,L)
  5   H11111(I,J,K,L,M)=-(A+B+C+D+E+F)/T21
      DO 6 M=1,3
      CALL FILL4A(3,3,H11111(1,1,1,1,M))
  6   CONTINUE
      RETURN
      END
