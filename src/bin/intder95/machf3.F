C     ////////////////////////////////////////////////////////////
C             (NUMB)
      SUBROUTINE MACHF3(NC,NS,NINV,F1,F3,YR,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER M,N,P,R
      DIMENSION F1(NC),F3(NC,NC,NC)
      DIMENSION YR(NC,NC,NC),V(NC)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
   4  FORMAT(I5)
   5  FORMAT(20X,F20.10)
      IF(NINV.GT.0) GO TO 200
      DO 100  R=1,NS
         CALL YIN(NC,NS,YR,R,ISCR3)
         DO 110  P=1,NS
         DO 110  N=1,P
         DO 110  M=1,N
 110     F3(M,N,P)=F3(M,N,P)-F1(R)*YR(M,N,P)
 100  CONTINUE
      CALL FILL3A(NC,NS,F3)
      RETURN
 200  REWIND I12
      READ(I12,4) R
      READ(I12,5)(V(R),R=1,NS)
      DO 210  R=1,NS
         CALL YIN(NC,NC,YR,-R,ISCR3)
         DO 220  P=1,NC
         DO 220  N=1,P
         DO 220  M=1,N
 220     F3(M,N,P)=F3(M,N,P)+V(R)*YR(M,N,P)
 210  CONTINUE
      CALL FILL3A(NC,NC,F3)
      RETURN
      END
