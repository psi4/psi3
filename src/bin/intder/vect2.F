C     //////////////////////////////////////////////////////////////
C BEND    (OPT2) (NUMB)
      SUBROUTINE VECT2(NAD,K1,K2,K3,V1,V2,V3,XA,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3)
      PARAMETER(ONE=1.0D0)
      CALL VECT1(NAD,K1,K2,E21,XA,T12)
      CALL VECT1(NAD,K3,K2,E23,XA,T32)
      CALL SCAPRO(E21,E23,W)
      CP=W
      SP=DSQRT(ONE-W*W)
      W=DACOS(W)
      C1=ONE/(T12*SP)
      C2=ONE/(T32*SP)
      DO 5  I=1,3
      V1(I)=(CP*E21(I)-E23(I))*C1
      V3(I)=(CP*E23(I)-E21(I))*C2
  5   V2(I)=-V1(I)-V3(I)
      RETURN
      END
