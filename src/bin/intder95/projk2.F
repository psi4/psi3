C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PROJK2(NAD,NA,NC,IOPT,XA,B2,BB2,DK,A,P,P3,P3A)
C /////////////////////////////////////////////////////////////////////
C COMPUTE SECOND ORDER PROJECTION MATRIX TO REMOVE ROTATIONAL AND
C TRANSLATIONAL COMPONENTS IF THE FORCE CONSTANT MATRICES WERE
C COMPUTED AT A NON-STATIONARY POINT
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NAD      INTEGER      NA+NDUM
C NA       INTEGER      NUMBER OF ATOMS
C NC       INTEGER      NA*3
C XA       (NAD,3)      CARTESIAN COORDINATES (IN ANGSTROMS)
C B2       (6,NC)       B2 DERIVATIVES
C A        (NC,6)       B2 "INVERSE"
C DK       (NC,NC)      SCRATCH SPACE
C P        (NC,NC)      SECOND ORDER PROJECTION MATRIX
C ON RETURN:
C BB2      (3,NC,NC)    B2 2ND DERIVATIVES ROTATIONS ONLY
C P3       (NC,NC,NC)   3RD ORDER PROJECTION MATRIX
C P3A      (NC,NC,NC)   Gijk DERIVATIVES
C /////////////////////////////////////////////////////////////////////

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),B2(6, NC),A(NC,6 ),P(NC,NC),DK(NC,NC)
      DIMENSION P3(NC,NC,NC), P3A(NC,NC,NC), BB2(3,NC,NC)
      DIMENSION W(6,12),IOPT(30),T0(3,6),T1(3,3),V(3),CM(3)
      DIMENSION V1(3),V2(3),V3(3),V4(3),V5(3),V6(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
C COMPUTE CENTER OF MOLECULE
      TM=DBLE(NA)
      DO 40 K=1,3
         CM(K)=ZERO
         DO 42  I=1,NA
  42     CM(K)=CM(K)+XA(I,K)
  40     CM(K)=CM(K)/TM

C FORM AND INVERT THE INERTIA TENSOR
      DO 62  I=1,3
      DO 62  J=1,6
  62     T0(I,J)=ZERO
      DO 63  I=1,3
  63     T0(I,I+3)=ONE
C
      DO 65  K=1,NA
         DO 67  L=1,3
  67     V1(L)=XA(K,L)-CM(L)
         DO 69  L=2,3
         DO 69  M=1,L-1
  69     T0(L,M)=T0(L,M)-V1(L)*V1(M)
         T0(1,1)=T0(1,1)+V1(2)*V1(2)+V1(3)*V1(3)
         T0(2,2)=T0(2,2)+V1(1)*V1(1)+V1(3)*V1(3)
  65     T0(3,3)=T0(3,3)+V1(1)*V1(1)+V1(2)*V1(2)
      T0(1,2)=T0(2,1)
      T0(1,3)=T0(3,1)
      T0(2,3)=T0(3,2)
C
      CALL FLIN(T0,3,3,3,DET)
      DO 68  I=1,3
      DO 68  J=1,3
  68  T1(I,J)=T0(I,J+3)
C
C     FORM Gijk DERIVATIVES
C
      DO 104  I=1,NA
      DO 104  J=1,NA
      LI=3*(I-1)
      LJ=3*(J-1)
      DO 104  IB=1,3
      DO 104  JG=1,3
        LL=LI+IB
        LM=LJ+JG
        DO 108  K=1,3
        V1(K)=B2(K+3,LL)
        V2(K)=B2(K+3,LM)
        V3(K)=XA(I,K)-CM(K)
 108    V4(K)=XA(J,K)-CM(K)
        CALL SCAPRO(V3,V2,D1)
        CALL SCAPRO(V4,V1,D2)
        DO 109  K=1,3
  109   V(K)=-V2(K)*V3(IB)-V1(K)*V4(JG)
        V(IB)=V(IB)+D1
        V(JG)=V(JG)+D2
        DO 111  IA=1,NA
           DO 115  K=1,3
  115      V3(K)=XA(IA,K)-CM(K)
           CALL SCAPRO(V1,V3,D1)
           CALL SCAPRO(V2,V3,D2)
           CALL VECPRO(V3,V2,V4)
           CALL VECPRO(V3,V1,V5)
           DO 113  K=1,3
  113      V(K)=V(K)+(V4(K)*D1+V5(K)*D2)/TWO
  111   CONTINUE
        DO 110  K=1,3
        V1(K)=ZERO
        DO 110  L=1,3
        V1(K)=V1(K)+T1(K,L)*V(L)
  110   BB2(L,I,J)=V1(L)
        DO 112  K=1,NC
        XX=ZERO
        DO 114  N=1,3
  114   XX=XX+A(K,N+3)*V1(N)
  112   P3A(K,LL,LM)=XX
  104   CONTINUE
C
C   FORM Pijk
C
      DO 120  I=1,NC
        DO 122  K=1,NC
        DO 122  L=1,NC
        DK(K,L)=ZERO
        DO 122  M=1,NC
  122   DK(K,L)=DK(K,L)+P3A(I,L,M)*P(M,K)
        DO 124  J=1,NC
        DO 124  K=1,NC
        P3(I,J,K)=ZERO
        DO 124  L=1,NC
  124   P3(I,J,K)=P3(I,J,K)-DK(K,L)*P(L,J)
  120 CONTINUE
      DO  130  K=1,NC
        DO 132  J=1,NC
        DO 132  L=1,NC
        DK(J,L)=ZERO
        DO 132  M=1,NC
  132   DK(J,L)=DK(J,L)+P3A(K,L,M)*P(M,J)
        DO 134  I=1,NC
        DO 134  J=1,NC
        DO 134  L=1,NC
  134   P3(I,J,K)=P3(I,J,K)-DK(J,L)*P(L,I)
  130 CONTINUE
      DO 140  I=1,NC
      DO 140  J=1,NC
      DO 140  K=1,NC
      DO 142  L=1,NC
  142 P3(I,J,K)=P3(I,J,K)-P3A(J,L,K)*P(L,I)
  140 P3(I,J,K)=P3(I,J,K)*BOHR

      WRITE (IOUT,*) 'P3:'
      DO 3000 I = 1, NC
 3001 FORMAT (I5,I5,I5,'    ',F13.10)

      WRITE (IOUT,3001) I,7,7,P3(I,7,7)
      WRITE (IOUT,3001) I,8,8,P3(I,8,8)
      WRITE (IOUT,3001) I,9,9,P3(I,9,9)
      WRITE (IOUT,3001) I,7,8,P3(I,7,8)
      WRITE (IOUT,3001) I,7,9,P3(I,7,9)
 3000 WRITE (IOUT,3001) I,8,9,P3(I,9,8)


      RETURN

      END
