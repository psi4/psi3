C     ////////////////////////////////////////////////////////////
C STRE    (OPT2)  (NUMB)
      SUBROUTINE HIJS1(NAD,K1,K2,XA,H11)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),H11(3,3)
      PARAMETER(ONE=1.0D0)
      CALL VECT1(NAD,K1,K2,V1,XA,T21)
      DO 5  J=1,3
      DO 5  I=J,3
  5   H11(I,J)=-V1(I)*V1(J)
      DO 10  I=1,3
  10  H11(I,I)=H11(I,I)+ONE
      W=ONE/T21
      DO 15  J=1,3
      DO 15  I=J,3
  15  H11(I,J)=H11(I,J)*W
      DO 20  J=1,2
      DO 20  I=J+1,3
  20  H11(J,I)=H11(I,J)
      RETURN
      END
