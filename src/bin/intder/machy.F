C     /////////////////////////////////////////////////////////////
      SUBROUTINE MACHY(NAD,NC,NS,NSX,IOPT,XA,TYPE,IA,A,S,U,IU,Y,SR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER R,P,Q,NSX,RP,IOPT(16)
      CHARACTER TYPE*5
      DIMENSION TYPE(NS),IA(NS,6),S(NS),U(NS,1),IU(NS,0:1)
      DIMENSION XA(NAD,3),A(NC,NC),Y(NC,NC,NC),SR(NC,NC,NC)
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
    1 FORMAT(/,1X,'NUMERICAL SR(I,J,K) AND Y(M,N,P) MATRICES USED FOR',
     $ /, ' SIMPLE INTERNAL COORDINATE',I5)
    2 FORMAT(/,1X,'SR(I,J,K) AND Y(M,N,P) MATRICES SET TO ZERO FOR',
     $ /, ' SIMPLE INTERNAL COORDINATE',I5)
      NSYM=IOPT(3)
      IF(NSYM.NE.0) THEN
         ISCR=ISCR4
      ELSE
         ISCR=ISCR3
      END IF
      DO 500  R=1,NS
         DO 10  P=1,NC
         DO 10  N=1,NC
         DO 10  M=1,NC
 10      Y(M,N,P)=ZERO
         DO 11  I=1,NC
         DO 11  J=1,NC
         DO 11  K=1,NC
 11      SR(I,J,K)=ZERO
      IF(IA(R,6).EQ.-2) GO TO 275
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
      DO 20  P=1,NSX
      DO 20  N=1,P
      DO 20  M=1,N
         DO 22  I=1,3
         DO 22  J=1,3
         DO 22  K=1,3
      W1=A(L1+J,N)*(A(L1+K,P)-A(L2+K,P))-A(L2+J,N)*(A(L1+K,P)-A(L2+K,P))
      W1=(A(L1+I,M)-A(L2+I,M))*W1
 22      Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)
 20   CONTINUE
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
      DO 40  P=1,NSX
      DO 40  N=1,P
      DO 40  M=1,N
         DO 42  I=1,3
         DO 42  J=1,3
         DO 42  K=1,3
         W1=A(L1+I,M)*A(L1+J,N)*A(L1+K,P)
         W2=A(L2+I,M)*A(L2+J,N)*A(L2+K,P)
         W3=A(L3+I,M)*A(L3+J,N)*A(L3+K,P)
         Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)+W2*H222(I,J,K)+W3*H333(I,J,K)
         W1=A(L1+I,M)*(A(L1+J,N)*A(L2+K,P)+A(L2+K,N)*A(L1+J,P))
         W1=W1+A(L2+K,M)*A(L1+I,N)*A(L1+J,P)
         W2=A(L1+I,M)*(A(L1+J,N)*A(L3+K,P)+A(L3+K,N)*A(L1+J,P))
         W2=W2+A(L3+K,M)*A(L1+I,N)*A(L1+J,P)
         W3=A(L3+I,M)*(A(L3+J,N)*A(L2+K,P)+A(L2+K,N)*A(L3+J,P))
         W3=W3+A(L2+K,M)*A(L3+I,N)*A(L3+J,P)
         W4=A(L3+I,M)*(A(L3+J,N)*A(L1+K,P)+A(L1+K,N)*A(L3+J,P))
         W4=W4+A(L1+K,M)*A(L3+I,N)*A(L3+J,P)
         W5=A(L2+I,M)*(A(L2+J,N)*A(L1+K,P)+A(L1+K,N)*A(L2+J,P))
         W5=W5+A(L1+K,M)*A(L2+I,N)*A(L2+J,P)
         W6=A(L2+I,M)*(A(L2+J,N)*A(L3+K,P)+A(L3+K,N)*A(L2+J,P))
         W6=W6+A(L3+K,M)*A(L2+I,N)*A(L2+J,P)
         Y(M,N,P)=Y(M,N,P)+W1*H112(I,J,K)+W2*H113(I,J,K)+W3*H332(I,J,K)
         Y(M,N,P)=Y(M,N,P)+W4*H331(I,J,K)+W5*H221(I,J,K)+W6*H223(I,J,K)
         W1=A(L1+I,M)*(A(L2+J,N)*A(L3+K,P)+A(L3+K,N)*A(L2+J,P))
         W1=W1+A(L2+J,M)*(A(L1+I,N)*A(L3+K,P)+A(L3+K,N)*A(L1+I,P))
         W1=W1+A(L3+K,M)*(A(L1+I,N)*A(L2+J,P)+A(L2+J,N)*A(L1+I,P))
 42      Y(M,N,P)=Y(M,N,P)+W1*H123(I,J,K)
 40   CONTINUE
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
      DO 90  P=1,NSX
      DO 90  N=1,P
      DO 90  M=1,N
         DO 92  I=1,3
         DO 92  J=1,3
         DO 92  K=1,3
         W1=A(L1+I,M)*A(L1+J,N)*A(L1+K,P)
         W2=A(L2+I,M)*A(L2+J,N)*A(L2+K,P)
         W3=A(L3+I,M)*A(L3+J,N)*A(L3+K,P)
         Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)+W2*H222(I,J,K)+W3*H333(I,J,K)
         W1=A(L1+I,M)*(A(L1+J,N)*A(L2+K,P)+A(L2+K,N)*A(L1+J,P))
         W1=W1+A(L2+K,M)*A(L1+I,N)*A(L1+J,P)
         W2=A(L1+I,M)*(A(L1+J,N)*A(L3+K,P)+A(L3+K,N)*A(L1+J,P))
         W2=W2+A(L3+K,M)*A(L1+I,N)*A(L1+J,P)
         W3=A(L3+I,M)*(A(L3+J,N)*A(L2+K,P)+A(L2+K,N)*A(L3+J,P))
         W3=W3+A(L2+K,M)*A(L3+I,N)*A(L3+J,P)
         W4=A(L3+I,M)*(A(L3+J,N)*A(L1+K,P)+A(L1+K,N)*A(L3+J,P))
         W4=W4+A(L1+K,M)*A(L3+I,N)*A(L3+J,P)
         W5=A(L2+I,M)*(A(L2+J,N)*A(L1+K,P)+A(L1+K,N)*A(L2+J,P))
         W5=W5+A(L1+K,M)*A(L2+I,N)*A(L2+J,P)
         W6=A(L2+I,M)*(A(L2+J,N)*A(L3+K,P)+A(L3+K,N)*A(L2+J,P))
         W6=W6+A(L3+K,M)*A(L2+I,N)*A(L2+J,P)
         Y(M,N,P)=Y(M,N,P)+W1*H112(I,J,K)+W2*H113(I,J,K)+W3*H332(I,J,K)
         Y(M,N,P)=Y(M,N,P)+W4*H331(I,J,K)+W5*H221(I,J,K)+W6*H223(I,J,K)
         W1=A(L1+I,M)*(A(L2+J,N)*A(L3+K,P)+A(L3+K,N)*A(L2+J,P))
         W1=W1+A(L2+J,M)*(A(L1+I,N)*A(L3+K,P)+A(L3+K,N)*A(L1+I,P))
         W1=W1+A(L3+K,M)*(A(L1+I,N)*A(L2+J,P)+A(L2+J,N)*A(L1+I,P))
 92      Y(M,N,P)=Y(M,N,P)+W1*H123(I,J,K)
 90   CONTINUE
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
C  OPTION
      DO 110  P=1,NSX
      DO 110  N=1,P
      DO 110  M=1,N
         DO 112  I=1,3
         DO 112  J=1,3
         DO 112  K=1,3
         W=A(L1+J,N)*A(L1+K,P)-A(L1+J,N)*A(L2+K,P)
         W=W-A(L2+J,N)*A(L1+K,P)+A(L2+J,N)*A(L2+K,P)
         W=(A(L1+I,M)-A(L2+I,M))*W
 112     Y(M,N,P)=Y(M,N,P)+W*H111(I,J,K)
 110  CONTINUE
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
C  OPTION
      DO 140  P=1,NSX
      DO 140  N=1,P
      DO 140  M=1,N
         DO 142  I=1,3
         DO 142  J=1,3
         DO 142  K=1,3
         W1=A(L1+I,M)*A(L1+J,N)*A(L1+K,P)
         W2=A(L2+I,M)*A(L2+J,N)*A(L2+K,P)
         W3=A(L3+I,M)*A(L3+J,N)*A(L3+K,P)
         W4=A(L4+I,M)*A(L4+J,N)*A(L4+K,P)
         Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)+W2*H222(I,J,K)
         Y(M,N,P)=Y(M,N,P)+W3*H333(I,J,K)+W4*H444(I,J,K)
         W1=A(L1+I,M)*(A(L1+J,N)*A(L2+K,P)+A(L2+K,N)*A(L1+J,P))
         W1=W1+A(L2+K,M)*A(L1+I,N)*A(L1+J,P)
         W2=A(L1+I,M)*(A(L1+J,N)*A(L3+K,P)+A(L3+K,N)*A(L1+J,P))
         W2=W2+A(L3+K,M)*A(L1+I,N)*A(L1+J,P)
         W3=A(L3+I,M)*(A(L3+J,N)*A(L2+K,P)+A(L2+K,N)*A(L3+J,P))
         W3=W3+A(L2+K,M)*A(L3+I,N)*A(L3+J,P)
         W4=A(L3+I,M)*(A(L3+J,N)*A(L1+K,P)+A(L1+K,N)*A(L3+J,P))
         W4=W4+A(L1+K,M)*A(L3+I,N)*A(L3+J,P)
         W5=A(L2+I,M)*(A(L2+J,N)*A(L1+K,P)+A(L1+K,N)*A(L2+J,P))
         W5=W5+A(L1+K,M)*A(L2+I,N)*A(L2+J,P)
         W6=A(L2+I,M)*(A(L2+J,N)*A(L3+K,P)+A(L3+K,N)*A(L2+J,P))
         W6=W6+A(L3+K,M)*A(L2+I,N)*A(L2+J,P)
         W7=A(L1+I,M)*(A(L1+J,N)*A(L4+K,P)+A(L4+K,N)*A(L1+J,P))
         W7=W7+A(L4+K,M)*A(L1+I,N)*A(L1+J,P)
         W8=A(L2+I,M)*(A(L2+J,N)*A(L4+K,P)+A(L4+K,N)*A(L2+J,P))
         W8=W8+A(L4+K,M)*A(L2+I,N)*A(L2+J,P)
         W9=A(L3+I,M)*(A(L3+J,N)*A(L4+K,P)+A(L4+K,N)*A(L3+J,P))
         W9=W9+A(L4+K,M)*A(L3+I,N)*A(L3+J,P)
         W10=A(L4+I,M)*(A(L4+J,N)*A(L1+K,P)+A(L1+K,N)*A(L4+J,P))
         W10=W10+A(L1+K,M)*A(L4+I,N)*A(L4+J,P)
         W11=A(L4+I,M)*(A(L4+J,N)*A(L2+K,P)+A(L2+K,N)*A(L4+J,P))
         W11=W11+A(L2+K,M)*A(L4+I,N)*A(L4+J,P)
         W12=A(L4+I,M)*(A(L4+J,N)*A(L3+K,P)+A(L3+K,N)*A(L4+J,P))
         W12=W12+A(L3+K,M)*A(L4+I,N)*A(L4+J,P)
         Y(M,N,P)=Y(M,N,P)+W1*H112(I,J,K)+W2*H113(I,J,K)+W3*H332(I,J,K)
         Y(M,N,P)=Y(M,N,P)+W4*H331(I,J,K)+W5*H221(I,J,K)+W6*H223(I,J,K)
         Y(M,N,P)=Y(M,N,P)+W7*H411(K,I,J)+W8*H422(K,I,J)+W9*H433(K,I,J)
      Y(M,N,P)=Y(M,N,P)+W10*H441(I,J,K)+W11*H442(I,J,K)+W12*H443(I,J,K)
         W1=A(L1+I,M)*(A(L2+J,N)*A(L3+K,P)+A(L3+K,N)*A(L2+J,P))
         W1=W1+A(L2+J,M)*(A(L1+I,N)*A(L3+K,P)+A(L3+K,N)*A(L1+I,P))
         W1=W1+A(L3+K,M)*(A(L1+I,N)*A(L2+J,P)+A(L2+J,N)*A(L1+I,P))
         W2=A(L4+I,M)*(A(L2+J,N)*A(L1+K,P)+A(L1+K,N)*A(L2+J,P))
         W2=W2+A(L2+J,M)*(A(L4+I,N)*A(L1+K,P)+A(L1+K,N)*A(L4+I,P))
         W2=W2+A(L1+K,M)*(A(L4+I,N)*A(L2+J,P)+A(L2+J,N)*A(L4+I,P))
         W3=A(L4+I,M)*(A(L3+J,N)*A(L1+K,P)+A(L1+K,N)*A(L3+J,P))
         W3=W3+A(L3+J,M)*(A(L4+I,N)*A(L1+K,P)+A(L1+K,N)*A(L4+I,P))
         W3=W3+A(L1+K,M)*(A(L4+I,N)*A(L3+J,P)+A(L3+J,N)*A(L4+I,P))
         W4=A(L4+I,M)*(A(L3+J,N)*A(L2+K,P)+A(L2+K,N)*A(L3+J,P))
         W4=W4+A(L3+J,M)*(A(L4+I,N)*A(L2+K,P)+A(L2+K,N)*A(L4+I,P))
         W4=W4+A(L2+K,M)*(A(L4+I,N)*A(L3+J,P)+A(L3+J,N)*A(L4+I,P))
         Y(M,N,P)=Y(M,N,P)+W1*H123(I,J,K)+W2*H421(I,J,K)
 142     Y(M,N,P)=Y(M,N,P)+W3*H431(I,J,K)+W4*H432(I,J,K)
 140  CONTINUE
      GO TO 300
 175  IF(IABS(IA(R,6)).NE.2) THEN
         WRITE(IOUT,2) R
         GO TO 300
      END IF
 275  CALL YIN(NC,NC,SR,-R,ISCR7)
      CALL YIN(NC,NSX,Y,R,ISCR7)
      WRITE(IOUT,1) R
 300  CONTINUE
      CALL FILL3A(NC,NSX,Y)
      CALL YOUT(NC,NC,SR,-R,ISCR)
      CALL YOUT(NC,NSX,Y,R,ISCR)
 500  CONTINUE
 600  IF(NSYM.EQ.0) RETURN
      DO 650  R=1,NSYM
         DO 610  P=1,NSYM
         DO 610  N=1,P
         DO 610  M=1,N
 610     Y(M,N,P)=ZERO
         L=IU(R,0)
         DO 620  I=1,L
         RP=IU(R,I)
         CALL YIN(NC,NSYM,SR,RP,ISCR4)
         W1=U(R,I)
         DO 630  P=1,NSYM
         DO 630  N=1,P
         DO 630  M=1,N
 630     Y(M,N,P)=Y(M,N,P)+W1*SR(M,N,P)
 620     CONTINUE
         CALL FILL3A(NC,NSYM,Y)
         CALL YOUT(NC,NSYM,Y,R,ISCR3)
 650  CONTINUE
C   OPTION
      DO 750  R=1,NSYM
         DO 710  P=1,NC
         DO 710  N=1,P
         DO 710  M=1,N
 710     SR(M,N,P)=ZERO
         L=IU(R,0)
         DO 720  I=1,L
         RP=IU(R,I)
         CALL YIN(NC,NC,Y,-RP,ISCR4)
         W1=U(R,I)
         DO 730  P=1,NC
         DO 730  N=1,P
         DO 730  M=1,N
 730     SR(M,N,P)=SR(M,N,P)+W1*Y(M,N,P)
 720     CONTINUE
         CALL FILL3A(NC,NC,SR)
         CALL YOUT(NC,NC,SR,-R,ISCR3)
 750  CONTINUE
      END
