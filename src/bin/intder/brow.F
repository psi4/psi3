C     //////////////////////////////////////////////////////////////
      SUBROUTINE BROW(NAD,NC,NS,XA,TYPE,IA,RB,S,IR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER IL(5)
      CHARACTER TYPE*5
      DIMENSION TYPE(NS),IA(NS,6),RB(NC),XA(NAD,3),S(NS)
      DIMENSION V(3,5)
      PARAMETER(ZERO=0.0D0)
         DO 5  J=1,NC
 5       RB(J)=ZERO
           NIA=0
           K1=IA(IR,1)
           K2=IA(IR,2)
           K3=IA(IR,3)
           K4=IA(IR,4)
           K5=IA(IR,5)
           IL(1)=3*(K1-1)
           IL(2)=3*(K2-1)
           IL(3)=3*(K3-1)
           IL(4)=3*(K4-1)
           IL(5)=3*(K5-1)
      IF(TYPE(IR).NE.' STRE') GO TO 10
         CALL VECT1(NAD,K1,K2,V(1,1),XA,W)
         DO 6  K=1,3
  6      V(K,2)=-V(K,1)
         NIA=2
         GO TO 80
 10   IF(TYPE(IR).NE.' BEND') GO TO 15
         CALL VECT2(NAD,K1,K2,K3,V(1,1),V(1,2),V(1,3),XA,W)
         NIA=3
         GO TO 80
 15   IF(TYPE(IR).NE.' LIN1') GO TO 25
           CALL VECT3(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),XA,W)
           NIA=3
           GO TO 80
 25   IF(TYPE(IR).NE.'  OUT') GO TO 30
           CALL VECT5(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)
           NIA=4
           GO TO 80
 30   IF(TYPE(IR).NE.' TORS') GO TO 35
           CALL VECT6(NAD,K1,K2,K3,K4,V(1,1),V(1,2),V(1,3),V(1,4),XA,W)
           NIA=4
           GO TO 80
 35   IF(TYPE(IR).NE.'  SPF') GO TO 40
           CALL VECT1(NAD,K1,K2,V(1,1),XA,W)
           DO 36  K=1,3
 36        V(K,2)=-V(K,1)
           NIA=2
           FACT=S(IR)/(W*W)
           DO 38  L=1,2
           DO 38  K=1,3
 38        V(K,L)=V(K,L)*FACT
           GO TO 80
 40   CONTINUE
 80   DO 85  K=1,NIA
      DO 85  J=1,3
 85   RB(IL(K)+J)=V(J,K)
      RETURN
      END
