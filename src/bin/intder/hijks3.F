C     ///////////////////////////////////////////////////////////////
C LIN1  (OPT1)
      SUBROUTINE HIJKS3(NAD,K1,K2,K3,K4,XA,
     $    H111,H112,H113,H123,H221,H222,H223,H331,H332,H333)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3),EA(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)
      DIMENSION H111(3,3,3),H112(3,3,3),H113(3,3,3),H123(3,3,3)
      DIMENSION H221(3,3,3),H222(3,3,3),H223(3,3,3),H331(3,3,3)
      DIMENSION H332(3,3,3),H333(3,3,3)
      DIMENSION H111A(3,3,3),H333A(3,3,3),H11A(3,3),H33A(3,3)
      PARAMETER(ONE=1.0D0)
         CALL VECT3(NAD,K1,K2,K3,K4,V1,V2,V3,XA,TH)
         CALL VECT1(NAD,K1,K2,E21,XA,T21)
         CALL VECT1(NAD,K3,K2,E23,XA,T23)
         CALL HIJS1(NAD,K1,K2,XA,H11A)
         CALL HIJS1(NAD,K3,K2,XA,H33A)
         DO 4  J=1,3
  4      EA(J)=XA(K4,J)
           CALL SCAPRO(EA,EA,D)
           D=ONE/DSQRT(D)
           DO 6  I=1,3
  6        EA(I)=D*EA(I)
         CALL HIJS3(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H22,H32,H33)
         CALL HIJKS1(NAD,K1,K2,XA,H111A)
         CALL HIJKS1(NAD,K3,K2,XA,H333A)
      TANTH=DTAN(TH)
      COSTH=DCOS(TH)
      W1=ONE/T21
      W2=ONE/T23
      W3=TANTH*W1
      W4=TANTH*W2
      DO 10  K=1,3
      DO 10  J=1,3
      DO 10  I=1,3
      H221(I,J,K)=H11(I,J)*(V1(K)*TANTH-E21(K)/T21)
      H221(I,J,K)=H221(I,J,K)+V1(K)*V1(J)*E21(I)*TANTH/T21
      H221(I,J,K)=H221(I,J,K)-(H11A(I,J)*V1(K))/T21
      H223(I,J,K)=H33(I,J)*(V3(K)*TANTH-E23(K)/T23)
      H223(I,J,K)=H223(I,J,K)+V3(K)*V3(J)*E23(I)*TANTH/T23
 10   H223(I,J,K)=H223(I,J,K)-(H33A(I,J)*V3(K))/T23
      DO 15  K=1,3
      DO 15  J=K,3
      DO 15  I=J,3
      H111(I,J,K)=(H221(I,J,K)+H221(J,K,I)+H221(K,I,J))+
     $    V1(I)*V1(J)*V1(K)-H111A(I,J,K)*W3
      H333(I,J,K)=(H223(I,J,K)+H223(J,K,I)+H223(K,I,J))+
     $    V3(I)*V3(J)*V3(K)-H333A(I,J,K)*W4
 15   CONTINUE
      CALL FILL3B(3,3,H111)
      CALL FILL3B(3,3,H333)
      DO 20  I=1,3
      W5=V1(I)*TANTH-E21(I)*W1
      W6=V3(I)*TANTH-E23(I)*W2
      DO 20  J=1,3
      DO 20  K=1,3
      H221(I,J,K)=W5*H31(K,J)
 20   H223(I,J,K)=W6*H31(J,K)
      W5=ONE/(COSTH*COSTH)
      DO 25  K=1,3
      DO 25  J=1,3
      DO 25  I=1,3
      H113(I,J,K)=V3(K)*(V1(I)*V1(J)-H11A(I,J)*W1)*W5
     $        +H221(I,J,K)+H221(J,I,K)
      H331(I,J,K)=V1(K)*(V3(I)*V3(J)-H33A(I,J)*W2)*W5
     $        +H223(I,J,K)+H223(J,I,K)
 25   CONTINUE
      DO 30  K=1,3
      DO 30  J=1,3
      DO 30  I=1,3
      H123(I,J,K)=-(H331(J,K,I)+H113(I,J,K))
      H112(I,J,K)=-(H111(I,J,K)+H113(I,J,K))
 30   H332(I,J,K)=-(H333(I,J,K)+H331(I,J,K))
      DO 35  K=1,3
      DO 35  J=1,3
      DO 35  I=1,3
      H221(J,K,I)=-(H123(I,J,K)+H112(I,K,J))
 35   H223(J,K,I)=-(H332(I,J,K)+H123(J,K,I))
      DO 40  K=1,3
      DO 40  J=1,3
      DO 40  I=1,3
 40   H222(I,J,K)=-(H223(J,K,I)+H221(J,K,I))
      RETURN
      END
