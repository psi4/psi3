C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PROJK3(NAD,NA,NC,IOPT,XA,BETA,M2,B2,A,BB2,P2,P3,
     $                  G3,D,P4)
C /////////////////////////////////////////////////////////////////////
C COMPUTE THIRD  ORDER PROJECTION MATRIX TO REMOVE ROTATIONAL AND
C TRANSLATIONAL COMPONENTS IF THE FORCE CONSTANT MATRICES WERE
C COMPUTED AT A NON-STATIONARY POINT
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C    NAD    INTEGER      NA+NDUM
C    NA     INTEGER      NUMBER OF ATOMS
C    NC     INTEGER      NA*3
C    XA     INTEGER      CARTESIAN COORS OF ATOMS
C    BETA   (3,NC,NC)    SCRATCH ARRAY
C    M2     SCRATCH ARRAY PASS IN SAME ADDRESS AS PASSED FOR BETA
C    B2     (6,NC)       B 1ST DERIVATIVES
C    A      (NC,6)       B2 "INVERSE"
C    BB2    (3,NC,NC)    B 2ND DERIVATIVES ROTATIONS ONLY
C    P2     (NC,NC)      2ND ORDER PROJECTION MATRIX
C    P3     (NC,NC,NC)   3RD ORDER PROJECTION MATRIX
C    G3     (NC,NC,NC)   Gijk
C    D      (NC,NC,NC)   SCRATCH SPACE
C ON RETURN:
C    P4     (NC,NC,NC,NC) 4TH ORDER PROJECTION MATRIX
C /////////////////////////////////////////////////////////////////////

      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION XA(NAD,3),B2( 6,NC),A(NC,6),P2(NC,NC),P4(NC,NC,NC,NC)
      DIMENSION P3(NC,NC,NC), G3(NC,NC,NC),M2(NC,NC), BETA(3,NC,NC)
      DIMENSION BB2(3,NC,NC), W(6,12),IOPT(30), CM(3), D(NC,NC,NC)
      DIMENSION V1(3),V2(3),V3(3),T1(3),T2(3),T3(3), AI(3), AJ(3), AK(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
       DO 100 I = 1, NC
       DO 100 J = 1, NC
       DO 100 K = 1, NC
       DO 100 L = 1, NC
 100     P4(I,J,K,L) = ZERO

C
C 7TH TERM OF P4
C

C FORM BETAi,j'S
      DO 300 I = 1, NC
      DO 300 J = 1, NC
      DO 300 L = 1, 3
          BETA(L,I,J) = ZERO
          DO 300 M = 1, NC
          DO 300 N = 1, NC
 300          BETA(L,I,J) = BETA(L,I,J)+P2(M,I)*P2(N,J)*BB2(L,M,N)

C FORM (B*BT)-1
      DO 305  I=1,6
      DO 305  J=1,12
 305    W(I,J)=ZERO
      DO 307  I=1,6
 307    W(I,I+6)=ONE
      DO 309 I=1,6
      DO 309 J=1,6
      DO 309 K=1,NC
 309    W(I,J)=W(I,J)+B2(I,K)*B2(J,K)
      CALL FLIN(W,6,6,6,DET)

C MOVE IN 7TH TERM OF P4
       DO 310 L = 1, NC
       DO 310 I = 1, NC
       DO 310 J = 1, NC
       DO 310 K = 1, NC

       DO 310 N = 1, 3
       DO 310 M = 1, 3

           XX = -W(N,M+3)*BETA(N,I,L)*BETA(M,J,K)
           P4(L,I,J,K) = P4(L,I,J,K) + XX
           P4(L,J,I,K) = P4(L,J,I,K) + XX
 310       P4(L,K,I,J) = P4(L,K,I,J) + XX

C
C 5TH AND 6TH TERMS OF P4
C

C CALCULTE Dijk'S
       DO 312 I = 1, NC
       DO 312 J = 1, NC
       DO 312 K = 1, NC
           D(I,J,K) = G3(I,J,K)
           DO 312 L = 1, NC
  312          D(I,J,K) = D(I,J,K) + P2(L,I)*G3(K,L,J)

C SUM IN TERMS 5 AND 6 DOING NC**2 M'S AT A TIME
       DO 325 L = 1, NC

       DO 325 LL = 1, NC
       DO 325 MM = 1, NC

       DO 320 II = 1, NC
       DO 320 JJ = 1, NC
 320      M2(II,JJ) = D(L,LL,MM)*D(II,MM,JJ) -
     $                (P3(MM,II,JJ)+D(II,MM,JJ))*G3(L,LL,MM)

       DO 325 I = 1, NC
       DO 325 J = 1, NC
       DO 325 K = 1, NC
C TERM 5
           XX = P2(LL,I)*M2(J,K)
           P4(L,I,J,K) = P4(L,I,J,K) + XX
           P4(L,J,I,K) = P4(L,J,I,K) + XX
           P4(L,K,I,J) = P4(L,K,I,J) + XX
C TERM 6
           XX = -P2(LL,L)*G3(I,MM,K)*G3(J,LL,M)
           P4(L,I,J,K) = P4(L,I,J,K) + XX
           P4(L,I,K,J) = P4(L,I,K,J) + XX

 325   CONTINUE


C
C TERMS 1-4 OF P4
C

C FIND MOLECULE'S CENTER OF MASS
      TM=DBLE(NA)
      DO 332 K=1,3
         CM(K)=ZERO
         DO 333 I=1,NA
 333     CM(K)=CM(K)+XA(I,K)
 332     CM(K)=CM(K)/TM

C FORM AND INVERT THE INERTIA TENSOR

      DO 326  I=1,3
      DO 326  J=1,6
 326     W(I,J)=ZERO
      DO 327  I=1,3
 327     W(I,I+3)=ONE

      DO 328  K=1,NA
         DO 329  L=1,3
 329     V1(L)=XA(K,L)-CM(L)
         DO 330  L=2,3
         DO 330  M=1,L-1
 330     W(L,M)=W(L,M)-V1(L)*V1(M)
         W(1,1)=W(1,1)+V1(2)*V1(2)+V1(3)*V1(3)
         W(2,2)=W(2,2)+V1(1)*V1(1)+V1(3)*V1(3)
 328     W(3,3)=W(3,3)+V1(1)*V1(1)+V1(2)*V1(2)
      W(1,2)=W(2,1)
      W(1,3)=W(3,1)
      W(2,3)=W(3,2)
      CALL FLIN(W,3,3,3,DET)

C  ONE Gl,i,j,k WILL BE FORMED AT A TIME AND THEN APPLIED WHERE REQUIRED
       DO 400 LL = 1, NC
       DO 400 II = 1, NC
       DO 400 JJ = 1, NC
       DO 400 KK = 1, NC

C CALC BBB2 DERIVATIVES IN RESPECT TO II, JJ, KK AND PUT IN V2

C      CARTESIAN COORD
       IIX = 1 +MOD(II-1,3)
       JJX = 1 +MOD(JJ-1,3)
       KKX = 1 +MOD(KK-1,3)
C      WHICH ATOM
       IIN = (II-IIX)/3 + 1
       JJN = (JJ-JJX)/3 + 1
       KKN = (KK-KKX)/3 + 1

       DO 335 K = 1, 3
           AI(K) = XA(IIN,K) - CM(K)
           AJ(K) = XA(JJN,K) - CM(K)
           AK(K) = XA(KKN,K) - CM(K)
           V1(K) = ZERO
           V2(K) = ZERO
 335       V3(K) = ZERO

           V1(IIX) = BB2(IIX,JJ,KK)
           V2(JJX) = BB2(JJX,II,KK)
           V3(KKX) = BB2(KKX,II,JJ)

       CALL SCAPRO(AI, V1, T1)
       CALL SCAPRO(AJ, V2, T2)
       CALL SCAPRO(AK, V3, T3)

       DO 337 I = 1, 3
 337       V1(I) = -T1(I) - T2(I) - T3(I)

       DO 339 I = 1, 3
           V2(I) = ZERO
       DO 339 J = 1, 3
 339       V2(I) = V2(I) + W(I,J+3)*V1(J)


C CALC G(LL,II,JJ,KK)
       G4 = ZERO
       DO 340 I = 1, 3
 340       G4 = G4 + A(LL,I+3)*V2(I)


C SUM IN GIVEN G4(LL,II,JJ,KK) IN TERMS AS REQUIRED:

C TERM 1
       DO 350 I = 1, NC
       DO 350 J = 1, NC
       DO 350 K = 1, NC
 350       P4(LL,I,J,K) = P4(LL,I,J,K) - P2(II,I)*P2(JJ,J)*P2(KK,K)*G4

C TERM 2
       DO 360 L = 1, NC
       DO 360 I = 1, NC
       DO 360 J = 1, NC
 360       P4(L,I,J,LL) = P4(L,I,J,LL) - P2(II,L)*P2(JJ,I)*P2(KK,J)*G4

C TERM 3
       DO 370 L = 1, NC
       DO 370 I = 1, NC
 370       P4(L,I,LL,II) = P4(L,I,LL,II) - P2(JJ,L)*P2(KK,I)*G4

C TERM 4
       DO 380 L = 1, NC
 380       P4(L,LL,II,JJ) = P4(L,LL,II,JJ) - P2(KK,L)*G4


C CYCLE TO NEXT G
 400   CONTINUE


 501   FORMAT(3F20.10)
       WRITE (IOUT,*) '4TH ORDER PROJECTION MATRIX'
       DO 500 I = 1, NC
       DO 500 J = 1, NC
            WRITE (IOUT,*)  'I=',I,'J=',J
       DO 499 K = 1, NC
       DO 499 L = 1, NC
 499        P4(I,J,K,L) =  P4(I,J,K,L)*BOHR
 500        WRITE (IOUT, 501)((P4(I,J,K,L),K=1,NC),L=1,NC)
       RETURN
       END
