C     //////////////////////////////////////////////////////////////
C TORS   (OPT2)   (NUMB)
      SUBROUTINE HIJS6(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,
     $               H22,H32,H42,H33,H43,H44)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3),V4(3),E21(3),E23(3),E34(3)
      DIMENSION BP21(3),BP22(3),BP23(3),BP32(3),BP34(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
         CALL VECT6(NAD,K1,K2,K3,K4,V1,V2,V3,V4,XA,W)
         CALL VECT1(NAD,K1,K2,E21,XA,T21)
         CALL VECT1(NAD,K3,K2,E23,XA,T23)
         CALL VECT1(NAD,K4,K3,E34,XA,T34)
         CALL VECT2(NAD,K1,K2,K3,BP21,BP22,BP23,XA,P2)
         CALL VECT2(NAD,K2,K3,K4,BP32,V2,BP34,XA,P3)
      CALL MAT1(H11,E23)
      CALL MAT1(H31,E21)
      CALL MAT1(H44,E23)
      CALL MAT1(H42,E34)
      XX=DSIN(P2)
      XY=DSIN(P3)
      XX=T21*XX*XX
      XY=T34*XY*XY
      W1=ONE/(T21*XX)
      W2=ONE/(T23*XX)
      W3=ONE/(T34*XY)
      W4=ONE/(T23*XY)
         DO 5  J=1,3
         DO 5  I=1,3
         H11(I,J)=-H11(I,J)*W1
         H31(I,J)=H31(I,J)*W2
         H44(I,J)=H44(I,J)*W3
  5      H42(I,J)=-H42(I,J)*W4
      XX=DCOS(P2)/DSIN(P2)
      XY=DCOS(P3)/DSIN(P3)
         DO 10  I=1,3
         W1=TWO*(E21(I)/T21+BP21(I)*XX)
         W2=(E23(I)/T23+TWO*BP23(I)*XX)
         W3=TWO*(E34(I)/T34+BP34(I)*XY)
         W4=(E23(I)/T23-TWO*BP32(I)*XY)
         DO 10  J=1,3
         H11(I,J)=H11(I,J)-W1*V1(J)
         H31(I,J)=H31(I,J)-W2*V1(J)
         H44(I,J)=H44(I,J)-W3*V4(J)
 10      H42(J,I)=H42(J,I)+W4*V4(J)
         DO 15  J=1,3
         DO 15  I=1,3
         H41(I,J)=ZERO
         H21(I,J)=-(H11(I,J)+H31(I,J))
 15      H43(I,J)=-(H44(I,J)+H42(I,J))
      X1=T21/T23
      Y1=T34/T23
      X2=DCOS(P2)
      Y2=DSIN(P2)
      X3=DCOS(P3)
      Y3=DSIN(P3)
      C1=X1*X2-ONE
      C2=-X3*Y1
      C3=-X2/T23
      C4=-X1*Y2
      C5=X1*X2/T23
      C6=Y1*Y3
      C7=-Y1*X3/T23
         DO 20  I=1,3
         W1=C3*E21(I)+C4*BP22(I)+C5*E23(I)
         W2=C6*BP32(I)+C7*E23(I)
         DO 20  J=1,3
 20      H22(I,J)=C1*H21(I,J)+C2*H42(J,I)+W1*V1(J)+W2*V4(J)
         DO 25  J=1,3
         DO 25  I=1,3
 25      H32(I,J)=-(H21(J,I)+H22(I,J)+H42(I,J))
         DO 30  J=1,3
         DO 30  I=1,3
 30      H33(I,J)=-(H31(I,J)+H32(I,J)+H43(J,I))
      RETURN
      END
