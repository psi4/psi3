C     //////////////////////////////////////////////////////////////
      SUBROUTINE XF2(NDER,NC,NS,NINV,BS,F2,F3,F4,V,XR,XS,XT)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,Q,R,S
      DIMENSION F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)
      DIMENSION XR(NC,NC),XS(NC,NC),XT(NC,NC),BS(NC,NC),V(NC,NC)
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
         DO 20  P=1,NS
         DO 20  N=1,P
         DO 20  M=1,N
         W=F2(P,R)*XR(M,N)+F2(N,R)*XR(M,P)+F2(M,R)*XR(N,P)
 20      F3(M,N,P)=F3(M,N,P)-W
 10   CONTINUE
      CALL FILL3A(NC,NS,F3)
      IF(NDER.LE.3) RETURN
      DO 60  R=1,NS
         DO 30  N=1,NS
         DO 30  M=1,NS
 30      XR(M,N)=ZERO
           DO 35  S=1,NS
           CALL XIN(NC,NS,XS,S,ISCR1)
           IF(R.EQ.S) THEN
             DO 40  N=1,NS
             DO 40  M=1,NS
 40          XT(M,N)=XS(M,N)
           END IF
           DO 45  N=1,NS
           DO 45  M=1,NS
 45        XR(M,N)=XR(M,N)+F2(R,S)*XS(M,N)
 35        CONTINUE
         DO 55  Q=1,NS
         DO 55  P=1,Q
         DO 55  N=1,P
         DO 55  M=1,N
         W=XT(M,N)*XR(P,Q)+XT(M,P)*XR(N,Q)+XT(M,Q)*XR(N,P)
 55      F4(M,N,P,Q)=F4(M,N,P,Q)-W
 60      CONTINUE
      CALL FILL4A(NC,NS,F4)
      RETURN
C
cbfy
c100  REWIND I16
c     READ(I16,1) M,N
 100  CONTINUE
      NLINES = (NS*NS +2)/3
      DO 105 M=1,NLINES
 105       BACKSPACE I16
cbfy
      READ(I16,2) ((V(M,N),N=1,NS),M=1,NS)
      DO 110  M=1,NS
      DO 110  J=1,NC
      XS(M,J)=ZERO
      DO 110  N=1,NS
 110  XS(M,J)=XS(M,J)+V(M,N)*BS(N,J)
      DO 120  R=1,NS
         CALL XIN(NC,NC,XR,-R,ISCR1)
         DO 125  K=1,NC
         DO 125  J=1,K
         DO 125  I=1,J
         W=XR(I,J)*XS(R,K)+XR(I,K)*XS(R,J)+XR(J,K)*XS(R,I)
 125     F3(I,J,K)=F3(I,J,K)+W
 120  CONTINUE
      CALL FILL3A(NC,NC,F3)
      IF(NDER.LE.3) RETURN
      DO 130  R=1,NS
         DO 140  J=1,NC
         DO 140  I=1,NC
 140     XR(I,J)=ZERO
         DO 145  S=1,NS
         CALL XIN(NC,NC,XS,-S,ISCR1)
         IF(R.EQ.S) THEN
            DO 150  J=1,NC
            DO 150  I=1,NC
 150        XT(I,J)=XS(I,J)
         END IF
         DO 155  J=1,NC
         DO 155  I=1,NC
 155     XR(I,J)=XR(I,J)+V(R,S)*XS(I,J)
 145     CONTINUE
         DO 160  L=1,NC
         DO 160  K=1,L
         DO 160  J=1,K
         DO 160  I=1,J
         W=XT(I,J)*XR(K,L)+XT(I,K)*XR(J,L)+XT(I,L)*XR(J,K)
 160     F4(I,J,K,L)=F4(I,J,K,L)+W
 130  CONTINUE
      CALL FILL4A(NC,NC,F4)
      RETURN
      END
