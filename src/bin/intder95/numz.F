C     /////////////////////////////////////////////////////////////
      SUBROUTINE NUMZ(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,YR1,YR2)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TYPE*5
      INTEGER R,P,Q,IOPT(30)
      DIMENSION XA(NAD,3),TYPE(NS),S(NS),A(NC,NC),IA(NS,6),V1(3)
      DIMENSION YR1(NC,NC,NC),YR2(NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,DELTA=1.0D-4)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      TDELTA=DELTA+DELTA
      NSYM=IOPT(3)
      IF(NSYM.NE.0) THEN
         ISCR=ISCR10
      ELSE
         ISCR=ISCR9
      END IF
      DO 10  R=1,NS
         IF(TYPE(R).EQ.'  SPF') THEN
              CALL VECT1(NAD,IA(R,1),IA(R,2),V1,XA,T21)
              S(R)=(ONE-S(R))*T21
         END IF
         DO 25  L=1,NC
         IP=(L-1)/3+1
         JP=L-3*(IP-1)
         XA(IP,JP)=XA(IP,JP)+DELTA
         CALL YROW(NAD,NC,NS,XA,TYPE,IA,S,YR1,R)
         XA(IP,JP)=XA(IP,JP)-TDELTA
         CALL YROW(NAD,NC,NS,XA,TYPE,IA,S,YR2,R)
            DO 28  K=1,NC
            DO 28  J=1,NC
            DO 28  I=1,NC
  28        YR1(I,J,K)=(YR1(I,J,K)-YR2(I,J,K))/TDELTA
         XA(IP,JP)=XA(IP,JP)+DELTA
         CALL YOUT2(NC,NC,YR1,-R,L,ISCR)
  25     CONTINUE
C
         DO 40  L=1,NC
         CALL YIN2(NC,NC,YR2,-R,L,ISCR)
            DO 32  K=1,NC
            DO 32  J=1,NC
            DO 32  M=1,NSX
  32        YR1(M,J,K)=ZERO
         DO 42  K=1,NC
         DO 42  J=1,NC
         DO 42  M=1,NSX
         DO 42  I=1,NC
 42      YR1(M,J,K)=YR1(M,J,K)+YR2(I,J,K)*A(I,M)
 40      CALL YOUT2(NC,NC,YR1,R,L,ISCR)
C
         DO 45  L=1,NC
         CALL YIN2(NC,NC,YR2,R,L,ISCR)
            DO 43  K=1,NC
            DO 43  N=1,NSX
            DO 43  M=1,NSX
 43         YR1(M,N,K)=ZERO
         DO 47  K=1,NC
         DO 47  N=1,NSX
         DO 47  M=1,NSX
         DO 47  J=1,NC
 47      YR1(M,N,K)=YR1(M,N,K)+YR2(M,J,K)*A(J,N)
 45      CALL YOUT2(NC,NC,YR1,R,L,ISCR)
C
         DO 55  L=1,NC
         CALL YIN2(NC,NC,YR2,R,L,ISCR)
              DO 56  P=1,NSX
              DO 56  N=1,NSX
              DO 56  M=1,NSX
 56           YR1(M,N,P)=ZERO
         DO 58  P=1,NSX
         DO 58  N=1,NSX
         DO 58  M=1,NSX
         DO 58  K=1,NC
 58      YR1(M,N,P)=YR1(M,N,P)+YR2(M,N,K)*A(K,P)
         LR=NS+1
 55      CALL YOUT2(NC,NC,YR1,-LR,L,ISCR)
C
         DO 70  Q=1,NSX
            DO 62  P=1,NSX
            DO 62  N=1,NSX
            DO 62  M=1,NSX
 62         YR1(M,N,P)=ZERO
         DO 65  L=1,NC
         LR=NS+1
         CALL YIN2(NC,NC,YR2,-LR,L,ISCR)
            DO 67  P=1,NSX
            DO 67  N=1,NSX
            DO 67  M=1,NSX
 67         YR1(M,N,P)=YR1(M,N,P)+YR2(M,N,P)*A(L,Q)
 65      CONTINUE
 70      CALL YOUT2(NC,NC,YR1,R,Q,ISCR)
C
         IF(TYPE(R).EQ.'  SPF') THEN
              S(R)=ONE-S(R)/T21
         END IF
 10   CONTINUE
      RETURN
      END
