

C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PDER3(NC,IOPT,P2,P3,P4,F2,F3,F4P)
C /////////////////////////////////////////////////////////////////////
C ADD REMAINING TERMS TO V*ijkl WHICH CONTAIN Pijk AND Pijkl
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NC        INTEGER         NUM ATOMS * 3
C P2        (NC,NC)         2ND ORDER PROJECTION MATRIX
C P3        (NC,NC,NC)      3RD   "               "
C P4        (NC,NC,NC,NC)   4TH   "               "
C F2        (NC,NC)         SCRATCH SPACE
C F3        (NC,NC)         SCRATCH SPACE
C F4P       (NC,NC,NC,NC)   V*ijkl WITH Vijkl CONTRIBUTIONS ONLY
C ON EXIT:
C F4P       (NC,NC,NC,NC)   V*ijkl COMPLETE
C /////////////////////////////////////////////////////////////////////

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      INTEGER*4 P
      DIMENSION F2(NC,NC), P2(NC,NC),F3(NC,NC,NC),P3(NC,NC,NC)
      DIMENSION F4P(NC,NC,NC,NC), P4(NC,NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT (2I5)
  2   FORMAT (3F20.10)

C READ IN UNPROJECTED FORCE CONSTANTS
      REWIND I15
      READ (I15,1) M,N
      READ (I15,2) ((F2(I,J), J=1,NC), I=1,NC)

      REWIND I20
      READ (I20,1) L,M
      READ (I20,2) (((F3(I,J,K), K=1,J),  J=1,I),  I=1,NC)

C FILL OUT F3 MATRIX
      DO 5 I = 1, NC
      DO 5 J = 1, I-1
        F3(J,I,J)=F3(J,J,I)
        F3(I,J,J)=F3(J,J,I)
        F3(I,J,I)=F3(I,I,J)
        F3(J,I,I)=F3(I,I,J)
        DO 5 K = 1,J-1
          F3(I,K,J)=F3(I,J,K)
          F3(J,I,K)=F3(I,J,K)
          F3(J,K,I)=F3(I,J,K)
          F3(K,I,J)=F3(I,J,K)
          F3(K,J,I)=F3(I,J,K)
  5   CONTINUE

C ADD THE LAST TWO TERMS INTO V*
      DO 10 I= 1, NC
      DO 10 J= 1, NC
      DO 10 K= 1, NC
      DO 10 L= 1, NC

      DO 10 M= 1, NC
      DO 10 N= 1, NC

C         TERMS FROM 3RD SUMMATION
          XX = F2(M,N)*P2(M,I)*P4(N,J,K,L)
          F4P(I,J,K,L)=F4P(I,J,K,L)+XX
          F4P(J,I,K,L)=F4P(J,I,K,L)+XX
          F4P(K,I,J,L)=F4P(K,I,J,L)+XX
          F4P(L,I,J,K)=F4P(L,I,J,K)+XX

          XX = F2(M,N)*P3(M,I,J)*P3(N,K,L)
          F4P(I,J,K,L)=F4P(I,J,K,L)+XX
          F4P(I,K,J,L)=F4P(I,K,J,L)+XX
          F4P(I,L,J,K)=F4P(I,L,J,K)+XX

      DO 10  P = 1, NC

C         TERMS FROM 2ND SUMMATION
          XX =  F3(M,N,P)*P2(M,I)*P2(N,J)*P3(P,K,L)
          F4P(I,J,K,L)=F4P(I,J,K,L)+XX
          F4P(I,K,J,L)=F4P(I,K,J,L)+XX
          F4P(I,L,J,K)=F4P(I,L,J,K)+XX
          F4P(J,K,I,L)=F4P(J,K,I,L)+XX
          F4P(J,L,I,K)=F4P(J,L,I,K)+XX
          F4P(K,L,I,J)=F4P(K,L,I,J)+XX

   10 CONTINUE
      RETURN
      END
