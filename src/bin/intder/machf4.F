C//////////////////////////////////////////////////////////////////////
      SUBROUTINE MACHF4(NC,NS,NINV,F1,F4,ZR,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,N,P,Q,R
      DIMENSION F1(NC),F4(NC,NC,NC,NC)
      DIMENSION ZR(NC,NC,NC,NC),V(NC)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
   4  FORMAT(I5)
   5  FORMAT(20X,F20.10)
      IF(NINV.GT.0) GO TO 200
      DO 100  R=1,NS
         CALL ZIN(NC,NS,ZR,R,ISCR9)
         DO 110  Q=1,NS
         DO 110  P=1,Q
         DO 110  N=1,P
         DO 110  M=1,N
 110     F4(M,N,P,Q)=F4(M,N,P,Q)-F1(R)*ZR(M,N,P,Q)
 100  CONTINUE
      CALL FILL4A(NC,NS,F4)
      RETURN
 200  REWIND I12
      READ(I12,4) R
      READ(I12,5)(V(R),R=1,NS)
      DO 210  R=1,NS
         CALL ZIN(NC,NC,ZR,-R,ISCR9)
         DO 220  Q=1,NC
         DO 220  P=1,Q
         DO 220  N=1,P
         DO 220  M=1,N
 220     F4(M,N,P,Q)=F4(M,N,P,Q)+V(R)*ZR(M,N,P,Q)
 210  CONTINUE
      CALL FILL4A(NC,NC,F4)
      RETURN
      END
