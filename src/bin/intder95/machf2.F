C     ////////////////////////////////////////////////////////////
      SUBROUTINE MACHF2(NC,NS,NINV,F1,F2,X,V)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER R
      DIMENSION F1(NC),F2(NC,NC),X(NC,NC),V(NC)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
   4  FORMAT(I5)
   5  FORMAT(20X,F20.10)
      IF(NINV.GT.0) GO TO 30
      DO 10  R=1,NS
         CALL XIN(NC,NS,X,R,ISCR1)
         DO 20  N=1,NS
         DO 20  M=1,NS
 20      F2(M,N)=F2(M,N)-F1(R)*X(M,N)
 10   CONTINUE
      RETURN
 30   REWIND I12
      READ(I12,4) R
      READ(I12,5)(V(R),R=1,NS)
 40   DO 50  R=1,NS
         CALL XIN(NC,NC,X,-R,ISCR1)
         DO 60  N=1,NC
         DO 60  M=1,NC
 60      F2(M,N)=F2(M,N)+V(R)*X(M,N)
 50   CONTINUE
      RETURN
      END
