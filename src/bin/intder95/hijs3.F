C     //////////////////////////////////////////////////////////////
C LIN1   (OPT2) (NUMB)
      SUBROUTINE HIJS3(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H22,H32,H33)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3),EA(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)
      DIMENSION H11A(3,3),H33A(3,3),EM(3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      CALL VECT3(NAD,K1,K2,K3,K4,V1,V2,V3,XA,TH)
      CALL VECT1(NAD,K1,K2,E21,XA,T21)
      CALL VECT1(NAD,K3,K2,E23,XA,T23)
      CALL HIJS1(NAD,K1,K2,XA,H11A)
      CALL HIJS1(NAD,K3,K2,XA,H33A)
         DO 5  J=1,3
  5      EA(J)=XA(K4,J)
           CALL SCAPRO(EA,EA,D)
           D=ONE/DSQRT(D)
           DO 10  I=1,3
  10       EA(I)=D*EA(I)
      TANTH=DTAN(TH)
      COSTH=DCOS(TH)
      EM(2,1)=EA(3)
      EM(3,1)=-EA(2)
      EM(3,2)=EA(1)
      EM(1,2)=-EM(2,1)
      EM(1,3)=-EM(3,1)
      EM(2,3)=-EM(3,2)
      DO 15  I=1,3
 15   EM(I,I)=ZERO
      DO 20  J=1,3
      DO 20  I=1,3
      H22(I,J)=ZERO
      DO 20  K=1,3
 20   H22(I,J)=H22(I,J)+EM(I,K)*H33A(K,J)
      W1=ONE/T21
      W2=ONE/T23
         DO 25  J=1,3
         DO 25  I=1,3
         H11(I,J)=(-H11A(I,J)*W1+V1(I)*V1(J))*TANTH-(E21(I)*V1(J)
     $                  +V1(I)*E21(J))*W1
         H31(I,J)=(H22(J,I)/COSTH-V3(I)*E21(J))/T21+V3(I)*V1(J)*TANTH
         H33(I,J)=(-H33A(I,J)*W2+V3(I)*V3(J))*TANTH-(E23(I)*V3(J)
     $                  +V3(I)*E23(J))*W2
         H21(I,J)=-(H11(I,J)+H31(I,J))
         H32(I,J)=-(H31(I,J)+H33(I,J))
 25      CONTINUE
         DO 30  J=1,3
         DO 30  I=1,3
         H22(I,J)=-(H21(J,I)+H32(I,J))
 30      CONTINUE
      RETURN
      END
