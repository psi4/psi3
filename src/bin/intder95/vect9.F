C     ////////////////////////////////////////////////////////////
C LINY    (OPT2) (NUMB)
      SUBROUTINE VECT9(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3)
      DIMENSION E4(3),E3(3),E2(3),E1(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,
     $PI=3.14159265358979323846264D0)
      CALL VECT5(NAD,K4,K3,K2,K1,E4,E3,E2,E1,XA,TOUT)
      W=-DSIN(TOUT)
      COSY=DCOS(TOUT)
      DO 5  I=1,3
      V2(I)=-COSY*E2(I)
      V3(I)=-COSY*E3(I)
      V4(I)=-COSY*E4(I)
 5    V1(I)=-COSY*E1(I)
      RETURN
      END
