C     ///////////////////////////////////////////////////////////////
      SUBROUTINE NUMY(NAD,NC,NS,NSX,XA,TYPE,IA,A,S,Y,SR,XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TYPE*5
      INTEGER R,P,Q
      DIMENSION XA(NAD,3),TYPE(NS),IA(NS,6),XT(NC,NC)
      DIMENSION A(NC,NC),S(NS),SR(NC,NC,NC),Y(NC,NC,NC)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,DELTA=1.0D-4)
      TDELTA=DELTA+DELTA
      DO 10  R=1,NS
              DO 700  K=1,NC
              DO 700  J=1,NC
              DO 700  I=1,NC
              Y(I,J,K)=ZERO
 700          SR(I,J,K)=ZERO
         IF(TYPE(R).EQ.'  SPF') THEN
              CALL VECT1(NAD,IA(R,1),IA(R,2),XT(1,1),XA,T21)
              S(R)=(ONE-S(R))*T21
         END IF
         IF(IA(R,6).EQ.0) GO TO 30
         DO 25  I=1,NC
         IP=(I-1)/3+1
         JP=I-3*(IP-1)
         XA(IP,JP)=XA(IP,JP)+DELTA
         CALL XROW(NAD,NC,NS,XA,TYPE,IA,S,XT,R)
         DO 27  K=1,NC
         DO 27  J=1,NC
 27           SR(I,J,K)=XT(J,K)
         XA(IP,JP)=XA(IP,JP)-TDELTA
         CALL XROW(NAD,NC,NS,XA,TYPE,IA,S,XT,R)
         DO 28  K=1,NC
         DO 28  J=1,NC
 28           SR(I,J,K)=(SR(I,J,K)-XT(J,K))/TDELTA
 25      XA(IP,JP)=XA(IP,JP)+DELTA
 30   CALL YOUT(NC,NC,SR,-R,ISCR7)
      IF(IABS(IA(R,6)).NE.2) GO TO 100
      DO 35  P=1,NSX
      DO 35  J=1,NC
      DO 35  I=1,NC
      DO 35  K=1,NC
 35   Y(I,J,P)=Y(I,J,P)+SR(I,J,K)*A(K,P)
      DO 40  P=1,NSX
      DO 40  N=1,NSX
      DO 40  I=1,NC
      SR(I,N,P)=ZERO
      DO 40  J=1,NC
 40   SR(I,N,P)=SR(I,N,P)+Y(I,J,P)*A(J,N)
      DO 45  P=1,NSX
      DO 45  N=1,NSX
      DO 45  M=1,NSX
      Y(M,N,P)=ZERO
      DO 45  I=1,NC
 45   Y(M,N,P)=Y(M,N,P)+SR(I,N,P)*A(I,M)
 100  CALL YOUT(NC,NSX,Y,R,ISCR7)
         IF(TYPE(R).EQ.'  SPF') THEN
              S(R)=ONE-S(R)/T21
         END IF
 10   CONTINUE
      END
