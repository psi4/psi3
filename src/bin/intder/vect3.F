C     //////////////////////////////////////////////////////////////
C LIN1   (OPT2) (NUMB)
      SUBROUTINE VECT3(NAD,K1,K2,K3,K4,V1,V2,V3,XA,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 V1(3),V2(3),V3(3),XA(NAD,3),E21(3),E23(3)
      DIMENSION E2M(3),EA(3),V4(3),V5(3)
      PARAMETER(ONE=1.0D0)
      CALL VECT1(NAD,K1,K2,E21,XA,T21)
      CALL VECT1(NAD,K3,K2,E23,XA,T23)
      DO 5  J=1,3
  5   EA(J)=XA(K4,J)
      CALL SCAPRO(EA,EA,D)
      D=ONE/DSQRT(D)
      DO 10  I=1,3
 10   EA(I)=D*EA(I)
      CALL VECPRO(E23,E21,E2M)
      CALL SCAPRO(EA,E2M,STHETA)
      W=DASIN(STHETA)
      CTHETA=DCOS(W)
      TTHETA=STHETA/CTHETA
      CALL VECPRO(EA,E23,V4)
      CALL VECPRO(EA,E21,V5)
      C1=ONE/(CTHETA*T21)
      C2=TTHETA/T21
      C3=ONE/(CTHETA*T23)
      C4=TTHETA/T23
      DO 15  I=1,3
      V1(I)=C1*V4(I)-C2*E21(I)
      V3(I)=-(C3*V5(I)+C4*E23(I))
 15   V2(I)=-(V1(I)+V3(I))
      RETURN
      END
