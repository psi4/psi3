C     ////////////////////////////////////////////////////////////////
      SUBROUTINE XF3(NC,NS,NINV,BS,XR,XS,XT,F3,F4,YR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,Q,R,S
      DIMENSION F3(NC,NC,NC),F4(NC,NC,NC,NC),YR(NC,NC,NC)
      DIMENSION XR(NC,NC),XS(NC,NC),XT(NC,NC),BS(NC,NC)
      PARAMETER(ZERO=0.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
  1   FORMAT(2I5)
  2   FORMAT(3F20.10)
      IF(NINV.GT.0) GO TO 100
      DO 10  R=1,NS
         CALL XIN(NC,NS,XR,R,ISCR1)
         DO 20  Q=1,NS
         DO 20  P=1,Q
         DO 20  N=1,P
         DO 20  M=1,N
         W=F3(M,N,R)*XR(P,Q)+F3(M,P,R)*XR(N,Q)+F3(M,Q,R)*XR(N,P)
         W=W+F3(N,P,R)*XR(M,Q)+F3(N,Q,R)*XR(M,P)+F3(P,Q,R)*XR(M,N)
 20      F4(M,N,P,Q)=F4(M,N,P,Q)-W
 10   CONTINUE
      CALL FILL4A(NC,NS,F4)
      RETURN
C
cbfy
c100  REWIND I21
c     READ(I21,1) M,N
 100  CONTINUE
      NLINES = (NS*(NS+1)*(NS+2)/6 +2)/3
      DO 105 M=1,NLINES
 105       BACKSPACE I21
cbfy
      READ(I21,2) (((YR(M,N,P),P=1,N),N=1,M),M=1,NS)
      CALL FILL3B(NC,NS,YR)
      DO 110  R=1,NS
         CALL XIN(NC,NC,XR,-R,ISCR1)
         DO 120  L=1,NC
         DO 120  N=1,NC
           XT(N,L)=ZERO
           DO 120  P=1,NC
 120       XT(N,L)=XT(N,L)+YR(R,N,P)*BS(P,L)
         DO 125  L=1,NC
         DO 125  K=1,NC
           XS(K,L)=ZERO
           DO 125  N=1,NC
 125       XS(K,L)=XS(K,L)+XT(N,L)*BS(N,K)
         DO 130  L=1,NC
         DO 130  K=1,L
         DO 130  J=1,K
         DO 130  I=1,J
         W=XR(I,J)*XS(K,L)+XR(I,K)*XS(J,L)+XR(J,K)*XS(I,L)
         W=W+XR(K,L)*XS(I,J)+XR(J,L)*XS(I,K)+XR(I,L)*XS(J,K)
 130     F4(I,J,K,L)=F4(I,J,K,L)+W
 110  CONTINUE
      CALL FILL4A(NC,NC,F4)
      RETURN
      END
