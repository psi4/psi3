C     //////////////////////////////////////////////////////////////
C BEND    (OPT2)  (NUMB)
      SUBROUTINE HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),E21(3),E23(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)
      DIMENSION H11A(3,3),H33A(3,3)
      PARAMETER(ONE=1.0D0)
         CALL VECT2(NAD,K1,K2,K3,V1,V2,V3,XA,PHI)
         CALL VECT1(NAD,K1,K2,E21,XA,T21)
         CALL VECT1(NAD,K3,K2,E23,XA,T23)
         CALL HIJS1(NAD,K1,K2,XA,H11A)
         CALL HIJS1(NAD,K3,K2,XA,H33A)
      SPHI=DSIN(PHI)
      CTPHI=DCOS(PHI)/SPHI
      W1=CTPHI
      W2=ONE/T21
      W3=W1*W2
      W4=ONE/T23
      W5=W1*W4
         DO 5  J=1,3
         DO 5  I=J,3
         H11(I,J)=H11A(I,J)*W3-V1(I)*V1(J)*W1
     $                        -(E21(I)*V1(J)+V1(I)*E21(J))*W2
         H33(I,J)=H33A(I,J)*W5-V3(I)*V3(J)*W1
     $                        -(E23(I)*V3(J)+V3(I)*E23(J))*W4
  5      CONTINUE
         DO 10  J=1,2
         DO 10  I=J+1,3
         H11(J,I)=H11(I,J)
  10     H33(J,I)=H33(I,J)
         W3=ONE/(T21*SPHI)
         DO 15  J=1,3
         W4=W2*E21(J)+W1*V1(J)
         DO 15  I=1,3
         H31(I,J)=-H33A(I,J)*W3-V3(I)*W4
         H21(I,J)=-(H11(I,J)+H31(I,J))
         H32(I,J)=-(H31(I,J)+H33(I,J))
 15      CONTINUE
         DO 20  J=1,3
         DO 20  I=1,3
         H22(I,J)=-(H21(J,I)+H32(I,J))
 20      CONTINUE
      RETURN
      END
