C     ////////////////////////////////////////////////////////////////
C     (OPT1) (NUMB)
      SUBROUTINE YF2(NC,NS,NINV,BS,XR,XS,F2,F4,YR)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,Q,R,S
      DIMENSION F2(NC,NC),F4(NC,NC,NC,NC)
      DIMENSION YR(NC,NC,NC),BS(NC,NC),XR(NC,NC),XS(NC,NC)
      PARAMETER(ZERO=0.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
  1   FORMAT(2I5)
  2   FORMAT(3F20.10)
      IF(NINV.GT.0) GO TO 200
      DO 110  R=1,NS
         CALL YIN(NC,NS,YR,R,ISCR3)
         DO 120  Q=1,NS
         DO 120  P=1,Q
         DO 120  N=1,P
         DO 120  M=1,N
         W=F2(M,R)*YR(N,P,Q)+F2(N,R)*YR(M,P,Q)
         W=W+F2(P,R)*YR(M,N,Q)+F2(Q,R)*YR(M,N,P)
 120     F4(M,N,P,Q)=F4(M,N,P,Q)-W
 110  CONTINUE
      CALL FILL4A(NC,NS,F4)
      RETURN
C
cbfy
c200  REWIND I16
c     READ(I16,1) M,N
 200  CONTINUE
      NLINES = (NS*NS +2)/3
      DO 205 M=1,NLINES
 205       BACKSPACE I16
cbfy
      READ(I16,2) ((XS(M,N),N=1,NS),M=1,NS)
      DO 210  M=1,NS
      DO 210  J=1,NC
      XR(M,J)=ZERO
      DO 210  N=1,NS
 210  XR(M,J)=XR(M,J)+XS(M,N)*BS(N,J)
      DO 220  R=1,NS
         CALL YIN(NC,NC,YR,-R,ISCR3)
         DO 230  L=1,NC
         DO 230  K=1,L
         DO 230  J=1,K
         DO 230  I=1,J
         W=YR(I,J,K)*XR(R,L)+YR(I,J,L)*XR(R,K)
         W=W+YR(I,K,L)*XR(R,J)+YR(J,K,L)*XR(R,I)
 230     F4(I,J,K,L)=F4(I,J,K,L)+W
 220  CONTINUE
      CALL FILL4A(NC,NC,F4)
      RETURN
      END
