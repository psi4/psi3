C     //////////////////////////////////////////////////////////////
      SUBROUTINE H4TH2(NAD,K1,K2,K3,XA,H1111,H1112,H1113,H1122,
     $    H1123,H1133,H1222,H1223,H1233,H1333,H2222,H2223,H2233,
     $    H2333,H3333)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),V1(3),V2(3),V3(3)
      DIMENSION H11(3,3),H21(3,3),H31(3,3),H22(3,3),H32(3,3),H33(3,3)
      DIMENSION H111(3,3,3),H112(3,3,3),H113(3,3,3),H123(3,3,3)
      DIMENSION H221(3,3,3),H222(3,3,3),H223(3,3,3),H331(3,3,3)
      DIMENSION H332(3,3,3),H333(3,3,3)
      DIMENSION H1111(3,3,3,3),H1112(3,3,3,3),H1113(3,3,3,3)
      DIMENSION H1122(3,3,3,3),H1123(3,3,3,3),H1133(3,3,3,3)
      DIMENSION H1222(3,3,3,3),H1223(3,3,3,3),H1233(3,3,3,3)
      DIMENSION H1333(3,3,3,3),H2222(3,3,3,3),H2223(3,3,3,3)
      DIMENSION H2233(3,3,3,3),H2333(3,3,3,3),H3333(3,3,3,3)
      DIMENSION Q11111(3,3,3,3,3),Q33333(3,3,3,3,3)
      DIMENSION Q1111(3,3,3,3),Q3333(3,3,3,3)
      PARAMETER(ONE=1.0D0)
      CALL VECT2(NAD,K1,K2,K3,V1,V2,V3,XA,PHI)
      CALL HIJS2(NAD,K1,K2,K3,XA,H11,H21,H31,H22,H32,H33)
      CALL HIJKS2(NAD,K1,K2,K3,XA,H111,H112,H113,H123,H221,
     $    H222,H223,H331,H332,H333)
      CSCP=ONE/DSIN(PHI)
      COTP=DCOS(PHI)*CSCP
      DO 10  L=1,3
      DO 10  K=1,3
      DO 10  J=1,3
      DO 10  I=1,3
      H1111(I,J,K,L)=V1(I)*V1(J)*V1(K)*V1(L) - H111(J,K,L)*V1(I) -
     $  H111(I,J,K)*V1(L) - H111(I,J,L)*V1(K) - H111(I,K,L)*V1(J) -
     $  H11(I,J)*H11(K,L) - H11(I,L)*H11(K,J) - H11(I,K)*H11(J,L)
      H1111(I,J,K,L)=H1111(I,J,K,L)*COTP +
     $  H11(I,J)*V1(K)*V1(L) + H11(I,K)*V1(J)*V1(L) +
     $  H11(I,L)*V1(J)*V1(K) + H11(J,K)*V1(I)*V1(L) +
     $  H11(J,L)*V1(K)*V1(I) + H11(K,L)*V1(J)*V1(I)
      H1113(I,J,K,L)=V1(I)*V1(J)*V1(K)*V3(L) - H113(J,K,L)*V1(I) -
     $  H111(I,J,K)*V3(L) - H113(I,J,L)*V1(K) - H113(I,K,L)*V1(J) -
     $  H11(I,J)*H31(L,K) - H31(L,I)*H11(K,J) - H11(I,K)*H31(L,J)
      H1113(I,J,K,L)=H1113(I,J,K,L)*COTP +
     $  H11(I,J)*V1(K)*V3(L) + H11(I,K)*V1(J)*V3(L) +
     $  H31(L,I)*V1(J)*V1(K) + H11(J,K)*V1(I)*V3(L) +
     $  H31(L,J)*V1(K)*V1(I) + H31(L,K)*V1(J)*V1(I)
      H1133(I,J,K,L)=V1(I)*V1(J)*V3(K)*V3(L) - H331(L,K,J)*V1(I) -
     $  H113(I,J,K)*V3(L) - H113(I,J,L)*V3(K) - H331(L,K,I)*V1(J) -
     $  H11(I,J)*H33(L,K) - H31(L,I)*H31(K,J) - H31(K,I)*H31(L,J)
      H1133(I,J,K,L)=H1133(I,J,K,L)*COTP +
     $  H11(I,J)*V3(K)*V3(L) + H31(K,I)*V1(J)*V3(L) +
     $  H31(L,I)*V1(J)*V3(K) + H31(K,J)*V1(I)*V3(L) +
     $  H31(L,J)*V3(K)*V1(I) + H33(L,K)*V1(J)*V1(I)
      H1333(I,J,K,L)=V1(I)*V3(J)*V3(K)*V3(L) - H333(L,K,J)*V1(I) -
     $  H331(K,J,I)*V3(L) - H331(L,J,I)*V3(K) - H331(L,K,I)*V3(J) -
     $  H31(J,I)*H33(L,K) - H31(L,I)*H33(K,J) - H31(K,I)*H33(L,J)
      H1333(I,J,K,L)=H1333(I,J,K,L)*COTP +
     $  H31(J,I)*V3(K)*V3(L) + H31(K,I)*V3(J)*V3(L) +
     $  H31(L,I)*V3(J)*V3(K) + H33(K,J)*V1(I)*V3(L) +
     $  H33(L,J)*V3(K)*V1(I) + H33(L,K)*V3(J)*V1(I)
      H3333(I,J,K,L)=V3(I)*V3(J)*V3(K)*V3(L) - H333(L,K,J)*V3(I) -
     $  H333(K,J,I)*V3(L) - H333(L,J,I)*V3(K) - H333(L,K,I)*V3(J) -
     $  H33(J,I)*H33(L,K) - H33(L,I)*H33(K,J) - H33(K,I)*H33(L,J)
      H3333(I,J,K,L)=H3333(I,J,K,L)*COTP +
     $  H33(J,I)*V3(K)*V3(L) + H33(K,I)*V3(J)*V3(L) +
     $  H33(L,I)*V3(J)*V3(K) + H33(K,J)*V3(I)*V3(L) +
     $  H33(L,J)*V3(K)*V3(I) + H33(L,K)*V3(J)*V3(I)
 10   CONTINUE
      CALL VECT1(NAD,K1,K2,V1,XA,R1)
      CALL VECT1(NAD,K3,K2,V3,XA,R3)
      CALL HIJS1(NAD,K1,K2,XA,H11)
      CALL HIJS1(NAD,K3,K2,XA,H33)
      CALL HIJKS1(NAD,K1,K2,XA,H111)
      CALL HIJKS1(NAD,K3,K2,XA,H333)
      CALL H4TH1(NAD,K1,K2,XA,Q1111)
      CALL H4TH1(NAD,K3,K2,XA,Q3333)
      CALL H5TH1(NAD,K1,K2,XA,Q11111)
      CALL H5TH1(NAD,K3,K2,XA,Q33333)
      DO 15  M=1,3
      DO 15  L=1,3
      DO 15  K=1,3
      DO 15  J=1,3
      DO 15  I=1,3
      H1111(I,J,K,L)=H1111(I,J,K,L)-CSCP*Q11111(M,I,J,K,L)*V3(M)
      H1113(I,J,K,L)=H1113(I,J,K,L)-CSCP*Q1111(M,I,J,K)*H33(L,M)
      H1133(I,J,K,L)=H1133(I,J,K,L)-CSCP*H111(M,I,J)*H333(K,L,M)
      H1333(I,J,K,L)=H1333(I,J,K,L)-CSCP*H11(M,I)*Q3333(J,K,L,M)
      H3333(I,J,K,L)=H3333(I,J,K,L)-CSCP*V1(M)*Q33333(I,J,K,L,M)
 15   CONTINUE
      DO 16  L=1,3
      DO 16  K=1,3
      DO 16  J=1,3
      DO 16  I=1,3
      H1112(I,J,K,L)=-H1111(I,J,K,L)-H1113(I,J,K,L)
      H1123(I,J,K,L)=-H1113(I,J,K,L)-H1133(I,J,K,L)
      H1233(I,J,K,L)=-H1133(I,J,K,L)-H1333(I,J,K,L)
 16   H2333(I,J,K,L)=-H1333(I,J,K,L)-H3333(I,J,K,L)
      DO 17  L=1,3
      DO 17  K=1,3
      DO 17  J=1,3
      DO 17  I=1,3
      H1122(I,J,K,L)=-H1112(I,J,K,L)-H1123(I,J,L,K)
      H1223(I,J,K,L)=-H1123(I,J,K,L)-H1233(I,K,J,L)
 17   H2233(I,J,K,L)=-H1233(I,J,K,L)-H2333(J,I,K,L)
      DO 18  L=1,3
      DO 18  K=1,3
      DO 18  J=1,3
      DO 18  I=1,3
      H1222(I,J,K,L)=-H1122(I,J,K,L)-H1223(I,K,L,J)
 18   H2223(I,J,K,L)=-H1223(I,J,K,L)-H2233(J,K,I,L)
      DO 19  L=1,3
      DO 19  K=1,3
      DO 19  J=1,3
      DO 19  I=1,3
 19   H2222(I,J,K,L)=-H1222(I,J,K,L)-H2223(J,K,L,I)
      RETURN
      END
