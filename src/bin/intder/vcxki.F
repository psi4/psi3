C     ////////////////////////////////////////////////////////////
      SUBROUTINE VCXKI(NA,NAD,NC,XMASS,XA,XKI)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER QQ
      DIMENSION XMASS(NA),XA(NAD,3),DIP(3),XKI(NC,3)
      DIMENSION W(3,3),XX(3,6),B(3),RCM(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
      COMMON /DIPOLE/DIP,QQ
      RMT=ZERO
      DO 80  K=1,NA
  80  RMT=RMT+XMASS(K)
      DO 90  I=1,3
      RCM(I)=ZERO
      DO 95  K=1,NA
  95  RCM(I)=RCM(I)+XA(K,I)*XMASS(K)
  90  RCM(I)=RCM(I)/RMT
      DO 100  I=1,3
      DO 100  J=1,3
      XX(I,J+3)=ZERO
 100  W(I,J)=-RMT*RCM(I)*RCM(J)
      DO 110  I=1,3
      DO 110  J=I,3
      DO 110  K=1,NA
 110  W(I,J)=W(I,J)+XA(K,I)*XA(K,J)*XMASS(K)
      XX(1,1)=W(1,3)
      XX(2,1)=-W(2,3)
      XX(3,1)=-(W(1,1)+W(2,2))
      XX(1,2)=-W(1,2)
      XX(2,2)=W(1,1)+W(3,3)
      XX(3,2)=-W(2,3)
      XX(1,3)=-(W(2,2)+W(3,3))
      XX(2,3)=W(1,2)
      XX(3,3)=W(1,3)
      XX(1,4)=ONE
      XX(2,5)=ONE
      XX(3,6)=ONE
      CALL FLIN(XX,3,3,3,DD)
      DO 1000  IX=1,NC
      II=(IX-1)/3+1
      IG=IX-3*(II-1)
      IF(IG.EQ.1) THEN
        B(1)=ZERO
        B(2)=XMASS(II)*(XA(II,3)-RCM(3))
        B(3)=-XMASS(II)*(XA(II,2)-RCM(2))
      ELSE IF(IG.EQ.2) THEN
        B(1)=-XMASS(II)*(XA(II,3)-RCM(3))
        B(2)=ZERO
        B(3)=XMASS(II)*(XA(II,1)-RCM(1))
      ELSE IF(IG.EQ.3) THEN
        B(1)=XMASS(II)*(XA(II,2)-RCM(2))
        B(2)=-XMASS(II)*(XA(II,1)-RCM(1))
        B(3)=ZERO
      END IF
      DO 1010 I=1,3
      XKI(IX,I)=ZERO
      DO 1010 J=1,3
 1010 XKI(IX,I)=XKI(IX,I)+XX(I,J+3)*B(J)
 1000 CONTINUE
      END
