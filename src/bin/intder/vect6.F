C     //////////////////////////////////////////////////////////////
C TORS    (OPT2) (NUMB)
      SUBROUTINE VECT6(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),V5(3),V6(3)
      DIMENSION E21(3),E32(3),E43(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,PI=3.14159265358979D0)
      CALL VECT1(NAD,K1,K2,E21,XA,T21)
      CALL VECT1(NAD,K2,K3,E32,XA,T32)
      CALL VECT1(NAD,K3,K4,E43,XA,T43)
      CALL VECPRO(E21,E32,V5)
      CALL VECPRO(E43,E32,V6)
      CALL SCAPRO(E21,E32,W2)
      CALL SCAPRO(E32,E43,W3)
      CP2=-W2
      CP3=-W3
      SP2=DSQRT(ONE-CP2*CP2)
      SP3=DSQRT(ONE-CP3*CP3)
      CALL SCAPRO(E21,V6,W2)
      CALL SCAPRO(V5,V6,W3)
      W3=-W3
C TRAP ROUNDING ERROR
      DASARG = W2/(SP2*SP3)
      IF (DABS(DASARG).GT.ONE) THEN
         IF ((DABS(DASARG)-ONE).LT.1.0D-12) THEN
            IF (DASARG.GT.ZERO) THEN
               DASARG =  1.0D0
               ELSE
               DASARG = -1.0D0
               ENDIF
            ENDIF
         ENDIF
      W=DASIN(DASARG)
C     W=DASIN(W2/(SP2*SP3))
      IF(W3.LT.ZERO) W=PI-W
      W1=ONE/(T21*SP2*SP2)
      W2=ONE/(T43*SP3*SP3)
      DO 5  I=1,3
      V1(I)=-W1*V5(I)
 5    V4(I)=-W2*V6(I)
      W3=(T32-T21*CP2)*W1/T32
      W4=CP3/(T32*SP3*SP3)
      W5=(T32-T43*CP3)*W2/T32
      W6=CP2/(T32*SP2*SP2)
      DO 10  I=1,3
      V2(I)=W3*V5(I)+W4*V6(I)
 10   V3(I)=W5*V6(I)+W6*V5(I)
      RETURN
      END
