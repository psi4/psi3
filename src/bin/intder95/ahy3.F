C     /////////////////////////////////////////////////////////////////
      SUBROUTINE AHY3(NC,NSX,L1,L2,L3,H111,H112,H113,H123,H221,H222,
     $      H223,H331,H332,H333,A,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,NSX
      DIMENSION A(NC,NC),Y(NC,NC,NC)
      DIMENSION H111(3,3,3),H112(3,3,3),H221(3,3,3),H222(3,3,3)
      DIMENSION H113(3,3,3),H123(3,3,3),H223(3,3,3),H331(3,3,3)
      DIMENSION H332(3,3,3),H333(3,3,3)
      DO 10  P=1,NSX
      DO 10  N=1,P
      DO 10  M=1,N
         DO 12  I=1,3
         DO 12  J=1,3
         V1=A(L1+I,M)*A(L1+J,N)
         V2=A(L2+I,M)*A(L2+J,N)
         V3=A(L3+I,M)*A(L3+J,N)
         V4=A(L1+I,M)*A(L1+J,P)
         V5=A(L1+I,N)*A(L1+J,P)
         V6=A(L2+I,M)*A(L2+J,P)
         V7=A(L2+I,N)*A(L2+J,P)
         V8=A(L3+I,M)*A(L3+J,P)
         V9=A(L3+I,N)*A(L3+J,P)
         V10=A(L1+I,M)*A(L2+J,N)
         V11=A(L1+I,M)*A(L2+J,P)
         V12=A(L1+I,N)*A(L2+J,M)
         V13=A(L1+I,N)*A(L2+J,P)
         V14=A(L1+I,P)*A(L2+J,M)
         V15=A(L1+I,P)*A(L2+J,N)
C
         DO 14  K=1,3
         W1=V1*A(L1+K,P)
         W2=V2*A(L2+K,P)
         W3=V3*A(L3+K,P)
         Y(M,N,P)=Y(M,N,P)+W1*H111(I,J,K)+W2*H222(I,J,K)+W3*H333(I,J,K)
         W1=V1*A(L2+K,P)+V4*A(L2+K,N)+V5*A(L2+K,M)
         W2=V1*A(L3+K,P)+V4*A(L3+K,N)+V5*A(L3+K,M)
         W3=V3*A(L2+K,P)+V8*A(L2+K,N)+V9*A(L2+K,M)
         W4=V3*A(L1+K,P)+V8*A(L1+K,N)+V9*A(L1+K,M)
         W5=V2*A(L1+K,P)+V6*A(L1+K,N)+V7*A(L1+K,M)
         W6=V2*A(L3+K,P)+V6*A(L3+K,N)+V7*A(L3+K,M)
         Y(M,N,P)=Y(M,N,P)+W1*H112(I,J,K)+W2*H113(I,J,K)+W3*H332(I,J,K)
         Y(M,N,P)=Y(M,N,P)+W4*H331(I,J,K)+W5*H221(I,J,K)+W6*H223(I,J,K)
         W1=A(L3+K,P)*(V10+V12)+A(L3+K,N)*(V11+V14)+A(L3+K,M)*(V13+V15)
 14      Y(M,N,P)=Y(M,N,P)+W1*H123(I,J,K)
 12   CONTINUE
 10   CONTINUE
      RETURN
      END
