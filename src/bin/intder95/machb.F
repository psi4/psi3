C     ////////////////////////////////////////////////////////////
      SUBROUTINE MACHB(NAD,NC,NS,XA,XMASS,TYPE,IA,B,S)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IL(5)
      CHARACTER TYPE*5
      DIMENSION TYPE(NS),IA(NS,6),B(NS,NC),XA(NAD,3),S(NS)
      DIMENSION V(3,5),Q(3,5),Z(3,5),XMASS(1)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      DO 100  I=1,NS
         DO 5  J=1,NC
 5       B(I,J)=ZERO
           NIA=0
           K1=IA(I,1)
           K2=IA(I,2)
           K3=IA(I,3)
           K4=IA(I,4)
           K5=IA(I,5)
           IL(1)=3*(K1-1)
           IL(2)=3*(K2-1)
           IL(3)=3*(K3-1)
           IL(4)=3*(K4-1)
           IL(5)=3*(K5-1)
      IF(TYPE(I).NE.' STRE') GO TO 10
         CALL VECT1(NAD,K1,K2,V(1,1),XA,W)
         DO 6  K=1,3
  6      V(K,2)=-V(K,1)
         NIA=2
         GO TO 80
 10   IF(TYPE(I).NE.' BEND') GO TO 15
         CALL VECT2(NAD,K1,K2,K3,V(1,1),V(1,2),V(1,3),XA,W)
         NIA=3
         GO TO 80
 15   IF(TYPE(I).NE.' LIN1') GO TO 25
         CALL VECT3(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),XA,W)
         NIA=3
         GO TO 80
 25   IF(TYPE(I).NE.'  OUT') GO TO 30
         CALL VECT5(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)
         NIA=4
         GO TO 80
 30   IF(TYPE(I).NE.' TORS') GO TO 35
         CALL VECT6(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)
         NIA=4
         GO TO 80
 35   IF(TYPE(I).NE.'  SPF') GO TO 40
         CALL VECT1(NAD,K1,K2,V(1,1),XA,W)
         DO 36  K=1,3
 36      V(K,2)=-V(K,1)
         NIA=2
         FACT=S(I)/(W*W)
         W=ONE-S(I)/W
         DO 38  L=1,2
         DO 38  K=1,3
 38      V(K,L)=V(K,L)*FACT
         GO TO 80
 40   IF(TYPE(I).NE.' LINX') GO TO 45
         CALL VECT8(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)
         NIA=4
         GO TO 80
 45   IF(TYPE(I).NE.' LINY') GO TO 50
         CALL VECT9(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)
         NIA=4
         GO TO 80
 50   IF(TYPE(I).NE.' RCOM') GO TO 55
         CALL VECT10(NAD,K1,K2,K3,K4,V(1,1),XA,XMASS,XMA,XMB,W)
         DO 51  K=K1,K2
           L=3*(K-1)
           DO 51  J=1,3
 51        B(I,L+J)=XMASS(K)*V(J,1)/XMA
         DO 52  K=K3,K4
           L=3*(K-1)
           DO 52  J=1,3
 52        B(I,L+J)=-XMASS(K)*V(J,1)/XMB
         GO TO 90
 55   CONTINUE
 80   DO 85  K=1,NIA
      DO 85  J=1,3
 85   B(I,IL(K)+J)=V(J,K)
 90   S(I)=W
 100  CONTINUE
      RETURN
      END
