
C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PDER2(NC,IOPT,F2,XS,P,P3,F3)
C /////////////////////////////////////////////////////////////////////
C ADD REMAINING TERMS TO V*ijk WHICH CONTAIN Pijk
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NC        INTEGER      NUMBER OF ATOMS * 3
C F2        (NC,NC)      SCRATCH SPACE
C XS        (NC,NC)      SCRATCH SPACE
C P         (NC,NC)      2ND ORDER PROJECTION MATRIX
C P3        (NC,NC,NC)   3RD ORDER PROJECTION MARTIX
C F3        (NC,NC,NC)   V*ijk WITH TERMS CONTAINING Pij ONLY
C ON RETURN:
C F3        (NC,NC,NC)   COMPLETE 3RD ORDER FORCE CONSTANTS
C /////////////////////////////////////////////////////////////////////

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XS(NC,NC),P(NC,NC),F2(NC,NC),P3(NC,NC,NC),F3(NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT(2I5)
  2   FORMAT(3F20.10)
      REWIND I15
      READ(I15,1) M,N
      READ(I15,2) ((F2(M,N),N=1,NC),M=1,NC)
      DO 10  M=2,NC
      DO 10  N=1,M-1
      XX=(F2(M,N)+F2(N,M))/TWO
      F2(M,N)=XX
 10   F2(N,M)=XX

      DO 25  I=1,NC
      DO 25  M=1,NC
      XS(I,M)=ZERO
      DO 25  L=1,NC
  25  XS(I,M)=XS(I,M)+F2(M,L)*P(L,I)
      DO 30  I=1,NC
      DO 30  J=1,NC
      DO 30  K=1,NC
      DO 30  M=1,NC
  30  F3(I,J,K)=F3(I,J,K)+XS(I,M)*P3(M,J,K)+XS(J,M)*P3(M,I,K)
     $          +XS(K,M)*P3(M,I,J)
CWA
      DO 150  I=1,NC
      WRITE(IOUT,*) 'I=',I
        DO 152  J=1,NC
  152   WRITE(IOUT,888) (P3(I,J,K),K=1,NC)
 150  CONTINUE
 888  FORMAT(9F10.6)
      RETURN
      END
