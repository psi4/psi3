C//////////////////////////////////////////////////////////////////////
      SUBROUTINE NUMX(NAD,NC,NS,NSX,XA,XMASS,TYPE,IA,RB,A,S,X,SR)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TYPE*5
      INTEGER R
      DIMENSION XA(NAD,3),XMASS(1),TYPE(NS),IA(NS,6)
      DIMENSION RB(NC),A(NC,NC),X(NC,NC),S(NS),SR(NC,NC)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,DELTA=1.0D-4)
      TDELTA=DELTA+DELTA
      DO 10  R=1,NS
              DO 700  I=1,NC
              DO 700  J=1,NC
              X(I,J)=ZERO
 700          SR(I,J)=ZERO
         IF(TYPE(R).EQ.'  SPF') THEN
              CALL VECT1(NAD,IA(R,1),IA(R,2),RB,XA,T21)
              S(R)=(ONE-S(R))*T21
         END IF
         IF(IA(R,6).EQ.0) GO TO 30
         DO 25  I=1,NC
         IP=(I-1)/3+1
         JP=I-3*(IP-1)
         XA(IP,JP)=XA(IP,JP)+DELTA
         CALL BROW(NAD,NC,NS,XA,XMASS,TYPE,IA,RB,S,R)
         DO 27  J=1,NC
 27           SR(I,J)=RB(J)
         XA(IP,JP)=XA(IP,JP)-TDELTA
         CALL BROW(NAD,NC,NS,XA,XMASS,TYPE,IA,RB,S,R)
         DO 28  J=1,NC
 28           SR(I,J)=(SR(I,J)-RB(J))/TDELTA
 25      XA(IP,JP)=XA(IP,JP)+DELTA
 30   CALL XOUT(NC,NC,SR,-R,ISCR6)
      IF(IABS(IA(R,6)).NE.1) GO TO 100
      DO 40  N=1,NSX
      DO 40  I=1,NC
      DO 40  J=1,NC
 40   X(I,N)=X(I,N)+SR(I,J)*A(J,N)
      DO 45  N=1,NSX
      DO 45  I=1,NC
      SR(I,N)=X(I,N)
 45   X(I,N)=ZERO
      DO 50  M=1,NSX
      DO 50  N=1,NSX
      DO 50  I=1,NC
 50   X(M,N)=X(M,N)+SR(I,N)*A(I,M)
100   CALL XOUT(NC,NSX,X,R,ISCR6)
         IF(TYPE(R).EQ.'  SPF') THEN
              S(R)=ONE-S(R)/T21
         END IF
 10   CONTINUE
      END
