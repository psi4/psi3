C     ////////////////////////////////////////////////////////////
C RCOM
      SUBROUTINE VECT10(NAD,K1,K2,K3,K4,V1,XA,XMASS,XMA,XMB,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),XMASS(1),V1(3),RA(3),RB(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      W=ZERO
      XMA=ZERO
      XMB=ZERO
      DO 2  J=1,3
      RA(J)=ZERO
  2   RB(J)=ZERO
      DO 5  K = K1,K2
        XMA=XMA+XMASS(K)
        DO 5  J = 1,3
  5     RA(J)=RA(J)+XA(K,J)*XMASS(K)
      DO 10  K = K3,K4
        XMB=XMB+XMASS(K)
        DO 10  J = 1,3
 10     RB(J)=RB(J)+XA(K,J)*XMASS(K)
      DO 15  J = 1,3
        RA(J)=RA(J)/XMA
        RB(J)=RB(J)/XMB
 15     V1(J)=RA(J)-RB(J)
      CALL SCAPRO(V1,V1,W)
      W=DSQRT(W)
      DO 20  J = 1,3
 20   V1(J)=V1(J)/W
      RETURN
      END
