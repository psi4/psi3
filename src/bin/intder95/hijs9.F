C     //////////////////////////////////////////////////////////////////
      SUBROUTINE HIJS9(NAD,K1,K2,K3,K4,XA,H11,H21,H31,H41,
     $               H22,H32,H42,H33,H43,H44)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),E1(3),E2(3),E3(3),E4(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H41(3,3),H22(3,3),H32(3,3)
      DIMENSION H42(3,3),H33(3,3),H43(3,3),H44(3,3)
      DIMENSION Q11(3,3),Q12(3,3),Q13(3,3),Q14(3,3),Q22(3,3),Q23(3,3)
      DIMENSION Q24(3,3),Q33(3,3),Q34(3,3),Q44(3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      CALL VECT5(NAD,K4,K3,K2,K1,E4,E3,E2,E1,XA,TOUT)
      W=-DSIN(TOUT)
      COSY=DCOS(TOUT)
      CALL HIJS7(NAD,K4,K3,K2,K1,XA,Q44,Q34,Q24,Q14,
     $           Q33,Q23,Q13,Q22,Q12,Q11)
      DO 1 K=1,3
      DO 1 J=1,3
      H22(J,K)=-W*E2(J)*E2(K)-COSY*Q22(K,J)
      H32(J,K)=-W*E3(J)*E2(K)-COSY*Q23(K,J)
      H42(J,K)=-W*E4(J)*E2(K)-COSY*Q24(K,J)
      H33(J,K)=-W*E3(J)*E3(K)-COSY*Q33(K,J)
      H43(J,K)=-W*E4(J)*E3(K)-COSY*Q34(K,J)
      H44(J,K)=-W*E4(J)*E4(K)-COSY*Q44(K,J)
      H41(J,K)=-W*E4(J)*E1(K)-COSY*Q14(K,J)
      H31(J,K)=-W*E3(J)*E1(K)-COSY*Q13(K,J)
      H21(J,K)=-W*E2(J)*E1(K)-COSY*Q12(K,J)
  1   H11(J,K)=-W*E1(J)*E1(K)-COSY*Q11(K,J)
      RETURN
      END
