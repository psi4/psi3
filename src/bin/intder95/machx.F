C     ///////////////////////////////////////////////////////////////
      SUBROUTINE MACHX(NAD,NC,NS,NSX,IOPT,XA,XMASS,
     $                  TYPE,IA,A,S,U,IU,X,SR)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TYPE*5
      INTEGER R,RP,NSX,IOPT(30)
      DIMENSION TYPE(NS),IA(NS,6),A(NC,NC),X(NC,NC),SR(NC,NC)
      DIMENSION XA(NAD,3),S(NS),XMASS(1),U(NS,1),IU(NS,0:1)
      DIMENSION H11(3,3),H21(3,3),H22(3,3),H31(3,3),H32(3,3)
      DIMENSION H33(3,3),H41(3,3),H42(3,3),H43(3,3),H44(3,3)
      DIMENSION E21(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(/,1X,'NUMERICAL SR(I,J) AND X(M,N) MATRICES USED FOR',/,
     $  ' SIMPLE INTERNAL COORDINATE',I5)
    2 FORMAT(/,1X,'SR(I,J) AND X(M,N) MATRICES SET TO ZERO FOR',/,
     $  ' SIMPLE INTERNAL COORDINATE',I5)
      NSYM=IOPT(3)
      IF(NSYM.NE.0) THEN
         ISCR=ISCR2
      ELSE
         ISCR=ISCR1
      END IF
C     ***********************************
      DO 500  R=1,NS
         DO 10  N=1,NC
         DO 10  M=1,NC
 10      X(M,N)=ZERO
C  OPTION
         DO 49  J=1,NC
         DO 49  I=1,NC
 49      SR(I,J)=ZERO
C  **************
         IF(IA(R,6).EQ.-1) GO TO 275
         IF(TYPE(R).NE.' STRE') GO TO 25
         K1=IA(R,1)
         K2=IA(R,2)
         L1=3*(K1-1)
         L2=3*(K2-1)
         CALL HIJS1(NAD,K1,K2,XA,H11)
C   OPTION
              DO 28  J=1,3
              DO 28  I=1,3
              SR(L1+I,L1+J)=H11(I,J)
              SR(L2+I,L2+J)=H11(I,J)
              SR(L1+I,L2+J)=-H11(I,J)
 28           SR(L2+I,L1+J)=-H11(I,J)
C  *************
         CALL AHX2(NC,NSX,L1,L2,H11,A,X)
         GO TO 300
 25      IF(TYPE(R).NE.' BEND') GO TO 75
         K1=IA(R,1)
         K2=IA(R,2)
         K3=IA(R,3)
         L1=3*(K1-1)
         L2=3*(K2-1)
         L3=3*(K3-1)
         CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)
C  OPTION
      DO 48  J=1,3
      DO 48  I=1,3
      SR(L1+I,L1+J)=H11(I,J)
      SR(L2+I,L1+J)=H21(I,J)
      SR(L3+I,L1+J)=H31(I,J)
      SR(L1+I,L2+J)=H21(J,I)
      SR(L2+I,L2+J)=H22(I,J)
      SR(L3+I,L2+J)=H32(I,J)
      SR(L1+I,L3+J)=H31(J,I)
      SR(L2+I,L3+J)=H32(J,I)
 48   SR(L3+I,L3+J)=H33(I,J)
C     ***********************************
         CALL AHX3(NC,NSX,L1,L2,L3,H11,H21,H31,H22,H32,H33,A,X)
         GO TO 300
 75      IF(TYPE(R).NE.' LIN1') GO TO 100
         K1=IA(R,1)
         K2=IA(R,2)
         K3=IA(R,3)
         K4=IA(R,4)
         L1=3*(K1-1)
         L2=3*(K2-1)
         L3=3*(K3-1)
         CALL HIJS3(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H22,H32,H33)
C  OPTION
      DO 86  J=1,3
      DO 86  I=1,3
      SR(L1+I,L1+J)=H11(I,J)
      SR(L2+I,L1+J)=H21(I,J)
      SR(L3+I,L1+J)=H31(I,J)
      SR(L1+I,L2+J)=H21(J,I)
      SR(L2+I,L2+J)=H22(I,J)
      SR(L3+I,L2+J)=H32(I,J)
      SR(L1+I,L3+J)=H31(J,I)
      SR(L2+I,L3+J)=H32(J,I)
 86   SR(L3+I,L3+J)=H33(I,J)
C     ***********************************
         CALL AHX3(NC,NSX,L1,L2,L3,H11,H21,H31,H22,H32,H33,A,X)
         GO TO 300
 100     IF(TYPE(R).NE.'  SPF') GO TO 125
              K1=IA(R,1)
              K2=IA(R,2)
              L1=3*(K1-1)
              L2=3*(K2-1)
              CALL VECT1(NAD,K1,K2,E21,XA,T21)
              CALL HIJS1(NAD,K1,K2,XA,H11)
              FACT1=(ONE-S(R))/T21
              FACT2=TWO*FACT1/T21
              DO 102  J=1,3
              DO 104  I=1,3
 104          H11(I,J)=THREE*H11(I,J)*FACT1
 102          H11(J,J)=H11(J,J)-FACT2
C   OPTION
              DO 108  J=1,3
              DO 108  I=1,3
              SR(L1+I,L1+J)=H11(I,J)
              SR(L2+I,L2+J)=H11(I,J)
              SR(L1+I,L2+J)=-H11(I,J)
 108          SR(L2+I,L1+J)=-H11(I,J)
C  *************
         CALL AHX2(NC,NSX,L1,L2,H11,A,X)
         GO TO 300
 125     IF(TYPE(R).NE.' TORS') GO TO 150
         K1=IA(R,1)
         K2=IA(R,2)
         K3=IA(R,3)
         K4=IA(R,4)
         L1=3*(K1-1)
         L2=3*(K2-1)
         L3=3*(K3-1)
         L4=3*(K4-1)
         CALL HIJS6(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,
     $               H33,H43,H44)
C  OPTION
      DO 148  J=1,3
      DO 148  I=1,3
      SR(L1+I,L1+J)=H11(I,J)
      SR(L2+I,L1+J)=H21(I,J)
      SR(L3+I,L1+J)=H31(I,J)
      SR(L4+I,L1+J)=H41(I,J)
      SR(L1+I,L2+J)=H21(J,I)
      SR(L2+I,L2+J)=H22(I,J)
      SR(L3+I,L2+J)=H32(I,J)
      SR(L4+I,L2+J)=H42(I,J)
      SR(L1+I,L3+J)=H31(J,I)
      SR(L2+I,L3+J)=H32(J,I)
      SR(L3+I,L3+J)=H33(I,J)
      SR(L4+I,L3+J)=H43(I,J)
      SR(L1+I,L4+J)=H41(J,I)
      SR(L2+I,L4+J)=H42(J,I)
      SR(L3+I,L4+J)=H43(J,I)
 148  SR(L4+I,L4+J)=H44(I,J)
C     ***********************************
         CALL AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,H22,H32,H42,
     $             H33,H43,H44,A,X)
         GO TO 300
 150     IF(TYPE(R).NE.'  OUT') GO TO 175
         K1=IA(R,1)
         K2=IA(R,2)
         K3=IA(R,3)
         K4=IA(R,4)
         L1=3*(K1-1)
         L2=3*(K2-1)
         L3=3*(K3-1)
         L4=3*(K4-1)
         CALL HIJS7(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,H22,H32,H42,
     $               H33,H43,H44)
C  OPTION
      DO 158  J=1,3
      DO 158  I=1,3
      SR(L1+I,L1+J)=H11(I,J)
      SR(L2+I,L1+J)=H21(I,J)
      SR(L3+I,L1+J)=H31(I,J)
      SR(L4+I,L1+J)=H41(I,J)
      SR(L1+I,L2+J)=H21(J,I)
      SR(L2+I,L2+J)=H22(I,J)
      SR(L3+I,L2+J)=H32(I,J)
      SR(L4+I,L2+J)=H42(I,J)
      SR(L1+I,L3+J)=H31(J,I)
      SR(L2+I,L3+J)=H32(J,I)
      SR(L3+I,L3+J)=H33(I,J)
      SR(L4+I,L3+J)=H43(I,J)
      SR(L1+I,L4+J)=H41(J,I)
      SR(L2+I,L4+J)=H42(J,I)
      SR(L3+I,L4+J)=H43(J,I)
 158  SR(L4+I,L4+J)=H44(I,J)
C     ***********************************
         CALL AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,H22,H32,H42,
     $             H33,H43,H44,A,X)
         GO TO 300
 175     IF(TYPE(R).NE.' LINX'.AND.TYPE(R).NE.' LINY') GO TO 200
         K1=IA(R,1)
         K2=IA(R,2)
         K3=IA(R,3)
         K4=IA(R,4)
         L1=3*(K1-1)
         L2=3*(K2-1)
         L3=3*(K3-1)
         L4=3*(K4-1)
         IF(TYPE(R).EQ.' LINX') CALL HIJS8(NAD,K1,K2,K3,K4,XA,H11,
     $       H21,H31,H41,H22,H32,H42,H33,H43,H44)
         IF(TYPE(R).EQ.' LINY') CALL HIJS9(NAD,K1,K2,K3,K4,XA,H11,
     $       H21,H31,H41,H22,H32,H42,H33,H43,H44)
C  OPTION
      DO 178  J=1,3
      DO 178  I=1,3
      SR(L1+I,L1+J)=H11(I,J)
      SR(L2+I,L1+J)=H21(I,J)
      SR(L3+I,L1+J)=H31(I,J)
      SR(L4+I,L1+J)=H41(I,J)
      SR(L1+I,L2+J)=H21(J,I)
      SR(L2+I,L2+J)=H22(I,J)
      SR(L3+I,L2+J)=H32(I,J)
      SR(L4+I,L2+J)=H42(I,J)
      SR(L1+I,L3+J)=H31(J,I)
      SR(L2+I,L3+J)=H32(J,I)
      SR(L3+I,L3+J)=H33(I,J)
      SR(L4+I,L3+J)=H43(I,J)
      SR(L1+I,L4+J)=H41(J,I)
      SR(L2+I,L4+J)=H42(J,I)
      SR(L3+I,L4+J)=H43(J,I)
 178  SR(L4+I,L4+J)=H44(I,J)
C     ***********************************
         CALL AHX4(NC,NSX,L1,L2,L3,L4,H11,H21,H31,H41,H22,H32,H42,
     $             H33,H43,H44,A,X)
         GO TO 300
 200     IF(TYPE(R).NE.' RCOM') GO TO 225
         K1=IA(R,1)
         K2=IA(R,2)
         K3=IA(R,3)
         K4=IA(R,4)
         CALL HIJS10(NAD,K1,K2,K3,K4,XA,XMASS,XMA,XMB,H11)
         DO 202  K = K1,K2
         DO 202  L = K1,K2
         FF=XMASS(K)*XMASS(L)/(XMA*XMA)
            DO 202  J=1,3
            DO 202  I=1,3
 202     SR(3*(K-1)+I,3*(L-1)+J)=H11(I,J)*FF
         DO 204  K = K3,K4
         DO 204  L = K3,K4
         FF=XMASS(K)*XMASS(L)/(XMB*XMB)
            DO 204  J=1,3
            DO 204  I=1,3
 204     SR(3*(K-1)+I,3*(L-1)+J)=H11(I,J)*FF
         DO 206  K = K1,K2
         DO 206  L = K3,K4
         FF=XMASS(K)*XMASS(L)/(XMA*XMB)
            DO 206  J=1,3
            DO 206  I=1,3
         SR(3*(K-1)+I,3*(L-1)+J)=-H11(I,J)*FF
 206     SR(3*(L-1)+I,3*(K-1)+J)=-H11(I,J)*FF
         CALL AHX(NC,NSX,SR,A,X)
         GO TO 300
C
 225     IF(IABS(IA(R,6)).NE.1) THEN
             WRITE(IOUT,2) R
             GO TO 300
         END IF
 275     CALL XIN(NC,NC,SR,-R,ISCR6)
         CALL XIN(NC,NSX,X,R,ISCR6)
         WRITE(IOUT,1) R
 300     CONTINUE
      CALL XOUT(NC,NC,SR,-R,ISCR)
      CALL XOUT(NC,NSX,X,R,ISCR)
 500  CONTINUE
 600  IF(NSYM.EQ.0) RETURN
      DO 650  R=1,NSYM
         DO 610  N=1,NSYM
         DO 610  M=1,NSYM
 610     X(M,N)=ZERO
         L=IU(R,0)
         DO 620  I=1,L
         RP=IU(R,I)
         CALL XIN(NC,NSYM,SR,RP,ISCR2)
         DO 630  N=1,NSYM
         DO 630  M=1,NSYM
 630     X(M,N)=X(M,N)+U(R,I)*SR(M,N)
 620     CONTINUE
         CALL XOUT(NC,NSYM,X,R,ISCR1)
 650  CONTINUE
C   OPTION
      DO 750  R=1,NSYM
         DO 710  N=1,NC
         DO 710  M=1,NC
 710     SR(M,N)=ZERO
         L=IU(R,0)
         DO 720  I=1,L
         RP=IU(R,I)
         CALL XIN(NC,NC,X,-RP,ISCR2)
         DO 730  N=1,NC
         DO 730  M=1,NC
 730     SR(M,N)=SR(M,N)+U(R,I)*X(M,N)
 720     CONTINUE
         CALL XOUT(NC,NC,SR,-R,ISCR1)
 750  CONTINUE
      END
