C     //////////////////////////////////////////////////////////////
C STRE  (OPT2)   (NUMB)
      SUBROUTINE HIJKS1(NAD,K1,K2,XA,H111)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),H11(3,3),H111(3,3,3)
      PARAMETER(ONE=1.0D0)
      CALL VECT1(NAD,K1,K2,V1,XA,T21)
      CALL HIJS1(NAD,K1,K2,XA,H11)
      W1=ONE/T21
      DO 5  K=1,3
      DO 5  J=K,3
      DO 5  I=J,3
  5   H111(I,J,K)=-(V1(I)*H11(K,J)+V1(J)*H11(K,I)+V1(K)*H11(J,I))*W1
      CALL FILL3B(3,3,H111)
      RETURN
      END
