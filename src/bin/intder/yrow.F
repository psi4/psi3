C////////////////////////////////////////////////////////////////////
      SUBROUTINE YROW(NAD,NC,NS,XA,TYPE,IA,S,SR,R)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER R,P,Q
      CHARACTER TYPE*5
      DIMENSION TYPE(NS),S(NS),IA(NS,6)
      DIMENSION XA(NAD,3),SR(NC,NC,NC)
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)
      DIMENSION H332(3,3,3),H333(3,3,3),H411(3,3,3),H421(3,3,3)
      DIMENSION H422(3,3,3),H431(3,3,3),H432(3,3,3),H433(3,3,3)
      DIMENSION H441(3,3,3),H442(3,3,3),H443(3,3,3),H444(3,3,3)
      DIMENSION E21(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
         DO 10  I=1,NC
         DO 10  J=1,NC
         DO 10  K=1,NC
 10      SR(I,J,K)=ZERO
      IF(TYPE(R).NE.' STRE') GO TO 25
      K1=IA(R,1)
      K2=IA(R,2)
      L1=3*(K1-1)
      L2=3*(K2-1)
      CALL HIJKS1(NAD,K1,K2,XA,H111)
C  OPTION
      DO 17  K=1,3
      DO 17  J=1,3
      DO 17  I=1,3
      Z=H111(I,J,K)
      SR(L1+I,L1+J,L1+K)=Z
      SR(L1+I,L1+J,L2+K)=-Z
      SR(L1+I,L2+J,L1+K)=-Z
      SR(L1+I,L2+J,L2+K)=Z
      SR(L2+I,L1+J,L1+K)=-Z
      SR(L2+I,L1+J,L2+K)=Z
      SR(L2+I,L2+J,L1+K)=Z
 17   SR(L2+I,L2+J,L2+K)=-Z
C  OPTION
      GO TO 300
 25   IF(TYPE(R).NE.' BEND') GO TO 75
      K1=IA(R,1)
      K2=IA(R,2)
      K3=IA(R,3)
      L1=3*(K1-1)
      L2=3*(K2-1)
      L3=3*(K3-1)
      CALL HIJKS2(NAD,K1,K2,K3,XA,H111,H112,H113,H123,H221,
     $    H222,H223,H331,H332,H333)
C  OPTION
      DO 37  K=1,3
      DO 37  J=1,3
      DO 37  I=1,3
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)
 37   SR(L3+I,L3+J,L3+K)=H333(I,J,K)
C  OPTION
         GO TO 300
 75      IF(TYPE(R).NE.' LIN1') GO TO 100
         K1=IA(R,1)
         K2=IA(R,2)
         K3=IA(R,3)
         K4=IA(R,4)
         L1=3*(K1-1)
         L2=3*(K2-1)
         L3=3*(K3-1)
      CALL HIJKS3(NAD,K1,K2,K3,K4,XA,H111,H112,H113,H123,
     $            H221,H222,H223,H331,H332,H333)
C  OPTION
      DO 88  K=1,3
      DO 88  J=1,3
      DO 88  I=1,3
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)
 88   SR(L3+I,L3+J,L3+K)=H333(I,J,K)
C  OPTION
      GO TO 300
 100  IF(TYPE(R).NE.'  SPF') GO TO 125
           K1=IA(R,1)
           K2=IA(R,2)
           L1=3*(K1-1)
           L2=3*(K2-1)
           CALL VECT1(NAD,K1,K2,E21,XA,T21)
           CALL HIJKS1(NAD,K1,K2,XA,H111)
           FACT1=(ONE-S(R))*THREE/T21
           FACT2=TWO/(T21*T21)
           DO 102  K=1,3
           DO 102  J=1,3
           DO 102  I=1,3
 102       H111(I,J,K)=FACT1*(H111(I,J,K)+FACT2*E21(I)*E21(J)*E21(K))
C  OPTION
      DO 107  K=1,3
      DO 107  J=1,3
      DO 107  I=1,3
      Z=H111(I,J,K)
      SR(L1+I,L1+J,L1+K)=Z
      SR(L1+I,L1+J,L2+K)=-Z
      SR(L1+I,L2+J,L1+K)=-Z
      SR(L1+I,L2+J,L2+K)=Z
      SR(L2+I,L1+J,L1+K)=-Z
      SR(L2+I,L1+J,L2+K)=Z
      SR(L2+I,L2+J,L1+K)=Z
 107  SR(L2+I,L2+J,L2+K)=-Z
      GO TO 300
 125  IF(TYPE(R).NE.' TORS') GO TO 175
      K1=IA(R,1)
      K2=IA(R,2)
      K3=IA(R,3)
      K4=IA(R,4)
      L1=3*(K1-1)
      L2=3*(K2-1)
      L3=3*(K3-1)
      L4=3*(K4-1)
      CALL HIJKS6(NAD,K1,K2,K3,K4,XA,H111,H112,H221,H222,H113,H123,
     $      H223,H331,H332,H333,H411,H421,H422,H431,H432,H433,H441,
     $      H442,H443,H444)
C  OPTION
      DO 137  K=1,3
      DO 137  J=1,3
      DO 137  I=1,3
      SR(L1+I,L1+J,L1+K)=H111(I,J,K)
      SR(L1+I,L1+J,L2+K)=H112(I,J,K)
      SR(L1+I,L1+J,L3+K)=H113(I,J,K)
      SR(L1+I,L1+J,L4+K)=H411(K,J,I)
      SR(L1+I,L2+J,L1+K)=H112(I,K,J)
      SR(L1+I,L2+J,L2+K)=H221(J,K,I)
      SR(L1+I,L2+J,L3+K)=H123(I,J,K)
      SR(L1+I,L2+J,L4+K)=H421(K,J,I)
      SR(L1+I,L3+J,L1+K)=H113(I,K,J)
      SR(L1+I,L3+J,L2+K)=H123(I,K,J)
      SR(L1+I,L3+J,L3+K)=H331(J,K,I)
      SR(L1+I,L3+J,L4+K)=H431(K,J,I)
      SR(L1+I,L4+J,L1+K)=H411(J,K,I)
      SR(L1+I,L4+J,L2+K)=H421(J,K,I)
      SR(L1+I,L4+J,L3+K)=H431(J,K,I)
      SR(L1+I,L4+J,L4+K)=H441(J,K,I)
      SR(L2+I,L1+J,L1+K)=H112(J,K,I)
      SR(L2+I,L1+J,L2+K)=H221(I,K,J)
      SR(L2+I,L1+J,L3+K)=H123(J,I,K)
      SR(L2+I,L1+J,L4+K)=H421(K,I,J)
      SR(L2+I,L2+J,L1+K)=H221(I,J,K)
      SR(L2+I,L2+J,L2+K)=H222(I,J,K)
      SR(L2+I,L2+J,L3+K)=H223(I,J,K)
      SR(L2+I,L2+J,L4+K)=H422(K,J,I)
      SR(L2+I,L3+J,L1+K)=H123(K,I,J)
      SR(L2+I,L3+J,L2+K)=H223(I,K,J)
      SR(L2+I,L3+J,L3+K)=H332(J,K,I)
      SR(L2+I,L3+J,L4+K)=H432(K,J,I)
      SR(L2+I,L4+J,L1+K)=H421(J,I,K)
      SR(L2+I,L4+J,L2+K)=H422(J,I,K)
      SR(L2+I,L4+J,L3+K)=H432(J,K,I)
      SR(L2+I,L4+J,L4+K)=H442(J,K,I)
      SR(L3+I,L1+J,L1+K)=H113(J,K,I)
      SR(L3+I,L1+J,L2+K)=H123(J,K,I)
      SR(L3+I,L1+J,L3+K)=H331(I,K,J)
      SR(L3+I,L1+J,L4+K)=H431(K,I,J)
      SR(L3+I,L2+J,L1+K)=H123(K,J,I)
      SR(L3+I,L2+J,L2+K)=H223(J,K,I)
      SR(L3+I,L2+J,L3+K)=H332(I,K,J)
      SR(L3+I,L2+J,L4+K)=H432(K,I,J)
      SR(L3+I,L3+J,L1+K)=H331(I,J,K)
      SR(L3+I,L3+J,L2+K)=H332(I,J,K)
      SR(L3+I,L3+J,L3+K)=H333(I,J,K)
      SR(L3+I,L3+J,L4+K)=H433(K,I,J)
      SR(L3+I,L4+J,L1+K)=H431(J,I,K)
      SR(L3+I,L4+J,L2+K)=H432(J,I,K)
      SR(L3+I,L4+J,L3+K)=H433(J,I,K)
      SR(L3+I,L4+J,L4+K)=H443(J,K,I)
      SR(L4+I,L1+J,L1+K)=H411(I,J,K)
      SR(L4+I,L1+J,L2+K)=H421(I,K,J)
      SR(L4+I,L1+J,L3+K)=H431(I,K,J)
      SR(L4+I,L1+J,L4+K)=H441(I,K,J)
      SR(L4+I,L2+J,L1+K)=H421(I,J,K)
      SR(L4+I,L2+J,L2+K)=H422(I,J,K)
      SR(L4+I,L2+J,L3+K)=H432(I,K,J)
      SR(L4+I,L2+J,L4+K)=H442(I,K,J)
      SR(L4+I,L3+J,L1+K)=H431(I,J,K)
      SR(L4+I,L3+J,L2+K)=H432(I,J,K)
      SR(L4+I,L3+J,L3+K)=H433(I,J,K)
      SR(L4+I,L3+J,L4+K)=H443(I,K,J)
      SR(L4+I,L4+J,L1+K)=H441(I,J,K)
      SR(L4+I,L4+J,L2+K)=H442(I,J,K)
      SR(L4+I,L4+J,L3+K)=H443(I,J,K)
 137  SR(L4+I,L4+J,L4+K)=H444(I,J,K)
 175  CONTINUE
      GO TO 300
 300  CONTINUE
      RETURN
      END
