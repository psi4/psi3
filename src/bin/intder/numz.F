C     /////////////////////////////////////////////////////////////
      SUBROUTINE NUMZ(NAD,NC,NS,NSX,XA,TYPE,IA,A,S,YR,ZR,SR)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TYPE*5
      INTEGER R,P,Q
      DIMENSION XA(NAD,3),TYPE(NS),S(NS),A(NC,NC),IA(NS,6),V1(3)
      DIMENSION YR(NC,NC,NC),SR(NC,NC,NC,NC),ZR(NC,NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,DELTA=1.0D-4)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
      TDELTA=DELTA+DELTA
      DO 10  R=1,NS
              DO 200  L=1,NC
              DO 200  K=1,NC
              DO 200  J=1,NC
              DO 200  I=1,NC
              SR(I,J,K,L)=ZERO
 200          ZR(I,J,K,L)=ZERO
         IF(TYPE(R).EQ.'  SPF') THEN
              CALL VECT1(NAD,IA(R,1),IA(R,2),V1,XA,T21)
              S(R)=(ONE-S(R))*T21
         END IF
         DO 25  I=1,NC
         IP=(I-1)/3+1
         JP=I-3*(IP-1)
         XA(IP,JP)=XA(IP,JP)+DELTA
         CALL YROW(NAD,NC,NS,XA,TYPE,IA,S,YR,R)
            DO 27  L=1,NC
            DO 27  K=1,NC
            DO 27  J=1,NC
  27        SR(I,J,K,L)=YR(J,K,L)
         XA(IP,JP)=XA(IP,JP)-TDELTA
         CALL YROW(NAD,NC,NS,XA,TYPE,IA,S,YR,R)
            DO 28  L=1,NC
            DO 28  K=1,NC
            DO 28  J=1,NC
  28        SR(I,J,K,L)=(SR(I,J,K,L)-YR(J,K,L))/TDELTA
  25     XA(IP,JP)=XA(IP,JP)+DELTA
         CALL ZOUT(NC,NC,SR,-R,ISCR8)
         DO 40  Q=1,NSX
         DO 40  K=1,NC
         DO 40  J=1,NC
         DO 40  I=1,NC
         DO 40  L=1,NC
 40      ZR(I,J,K,Q)=ZR(I,J,K,Q)+SR(I,J,K,L)*A(L,Q)
         DO 45  Q=1,NSX
         DO 45  P=1,NSX
         DO 45  J=1,NC
         DO 45  I=1,NC
         SR(I,J,P,Q)=ZERO
         DO 45  K=1,NC
 45      SR(I,J,P,Q)=SR(I,J,P,Q)+ZR(I,J,K,Q)*A(K,P)
         DO 50  Q=1,NSX
         DO 50  P=1,NSX
         DO 50  N=1,NSX
         DO 50  I=1,NC
         ZR(I,N,P,Q)=ZERO
         DO 50  J=1,NC
 50      ZR(I,N,P,Q)=ZR(I,N,P,Q)+SR(I,J,P,Q)*A(J,N)
         DO 55  Q=1,NSX
         DO 55  P=1,NSX
         DO 55  N=1,NSX
         DO 55  M=1,NSX
         SR(M,N,P,Q)=ZERO
         DO 55  I=1,NC
 55      SR(M,N,P,Q)=SR(M,N,P,Q)+ZR(I,N,P,Q)*A(I,M)
         CALL ZOUT(NC,NS,SR,R,ISCR8)
            IF(TYPE(R).EQ.'  SPF') THEN
                S(R)=ONE-S(R)/T21
            END IF
 10   CONTINUE
      RETURN
      END
