C     /////////////////////////////////////////////////////////////
      SUBROUTINE XROW(NAD,NC,NS,XA,TYPE,IA,S,SR,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TYPE*5
      INTEGER R,RP
      DIMENSION TYPE(NS),IA(NS,6),S(NS),XA(NAD,3),SR(NC,NC)
      DIMENSION H11(3,3),H21(3,3),H22(3,3),H31(3,3),H32(3,3)
      DIMENSION H33(3,3),H41(3,3),H42(3,3),H43(3,3),H44(3,3)
      DIMENSION H11A(3,3),H33A(3,3),EA(3),E21(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
         DO 49  J=1,NC
         DO 49  I=1,NC
 49      SR(I,J)=ZERO
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
         GO TO 300
 100  IF(TYPE(R).NE.'  SPF') GO TO 125
         K1=IA(R,1)
         K2=IA(R,2)
         L1=3*(K1-1)
         L2=3*(K2-1)
         CALL VECT1(NAD,K1,K2,E21,XA,T21)
         CALL HIJS1(NAD,K1,K2,XA,H11)
         FACT1=S(R)/(T21*T21)
         FACT2=TWO*FACT1/T21
         DO 102  J=1,3
         DO 104  I=1,3
 104     H11(I,J)=THREE*H11(I,J)*FACT1
 102     H11(J,J)=H11(J,J)-FACT2
         DO 106  J=1,3
         DO 106  I=1,3
         SR(L1+I,L1+J)=H11(I,J)
         SR(L2+I,L2+J)=H11(I,J)
         SR(L1+I,L2+J)=-H11(I,J)
 106     SR(L2+I,L1+J)=-H11(I,J)
         GO TO 300
 125     IF(TYPE(R).NE.' TORS') GO TO 175
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
         GO TO 300
 175  CONTINUE
 300  CONTINUE
      END
