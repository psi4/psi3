C     ///////////////////////////////////////////////////////////////
C OUT     (OPT2)  (NUMB)
      SUBROUTINE VECT5(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V1(3),V2(3),V3(3),V4(3),V5(3),V6(3),V7(3)
      DIMENSION XA(NAD,3),E23(3),E21(3),E24(3),B3P(3),B4P(3)
      PARAMETER(ONE=1.0D0)
      CALL VECT1(NAD,K1,K2,E21,XA,T21)
      CALL VECT1(NAD,K3,K2,E23,XA,T23)
      CALL VECT1(NAD,K4,K2,E24,XA,T24)
      CALL VECPRO(E23,E24,V5)
      CALL VECPRO(E24,E21,V6)
      CALL VECPRO(E21,E23,V7)
      CALL VECT2(NAD,K3,K2,K4,B3P,V2,B4P,XA,PHI)
      SPHI=DSIN(PHI)
      CALL SCAPRO(E21,V5,W)
      W=DASIN(W/SPHI)
      CG=DCOS(W)
      SG=DSIN(W)
      TG=SG/CG
      W1=CG*SPHI
      W2=ONE/(T21*W1)
      W3=TG/T21
      W4=ONE/(T23*W1)
      W5=T24*SG*W4
      W6=ONE/(T24*W1)
      W7=T23*SG*W6
      DO 10  I=1,3
      V1(I)=V5(I)*W2-E21(I)*W3
      V3(I)=V6(I)*W4+B4P(I)*W5
      V4(I)=V7(I)*W6+B3P(I)*W7
 10   V2(I)=-V1(I)-V3(I)-V4(I)
      RETURN
      END
