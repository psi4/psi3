C     /////////////////////////////////////////////////////////////////
C     RCOM
      SUBROUTINE HIJS10(NAD,K1,K2,K3,K4,XA,XMASS,XMA,XMB,H11)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),XMASS(1),RA(3),RB(3),V1(3),H11(3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      XMA=ZERO
      XMB=ZERO
      DO 2  J=1,3
      RA(J)=ZERO
 2    RB(J)=ZERO
      DO 5  K = K1,K2
      XMA=XMA+XMASS(K)
         DO 5  J = 1,3
 5       RA(J)=RA(J)+XMASS(K)*XA(K,J)
      DO 10  K = K3,K4
      XMB=XMB+XMASS(K)
         DO 10  J = 1,3
 10      RB(J)=RB(J)+XMASS(K)*XA(K,J)
      DO 15  J = 1,3
         RA(J)=RA(J)/XMA
         RB(J)=RB(J)/XMB
 15      V1(J)=RA(J)-RB(J)
      CALL SCAPRO(V1,V1,W)
      W=DSQRT(W)
      DO 20  I = 1,3
 20     V1(I)=V1(I)/W
      DO 25  J = 1,3
      DO 25  I = 1,3
 25     H11(I,J)=-V1(I)*V1(J)
      DO 30  I = 1,3
 30     H11(I,I)=H11(I,I)+ONE
      DO 35  J = 1,3
      DO 35  I = 1,3
 35     H11(I,J)=H11(I,J)/W
      RETURN
      END
