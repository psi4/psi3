C /////////////////////////////////////////////////////////////////////
      SUBROUTINE PROJK1(NAD,NA,NC,IOPT,XA,B2,A ,P)
C /////////////////////////////////////////////////////////////////////
C COMPUTES THE PROJECTION MATRIX (P(NC,NC)) WHICH REMOVES TRANSLATIONAL
C AND ROTATIONAL COMPONENTS IF THE FORCE CONSTANT MATRICES WERE
C OBTAINED AT A NON-STATIONARY POINT
C /////////////////////////////////////////////////////////////////////
C ON ENTRY:
C NAD         INTEGER        NA+NDUM
C NA          INTEGER        NUMBER OF ATOMS
C NC          INTEGER        3*NA
C XA          (NAD,3)        CARTESIAN COORDINATES (IN ANGSTROMS)
C ON RETURN:
C B2          ( 6,NC)        Bij DERIVATIVES
C A           (NC,6 )        B2 "INVERSE"
C P           (NC,NC)        PROJECTION MATRIX
C /////////////////////////////////////////////////////////////////////
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),B2( 6,NC), A(NC,6 ),P(NC,NC)
      DIMENSION W(6,12),IOPT(30),T0(3,6),V(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,TWO=2.0D0,THREE=3.0D0,FOUR=4.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
C
    1 FORMAT(//1X,'B2 MATRIX',/1X,'(ROWS 1-3: TRANSLATION, ROWS 4-6:',
     $            ' ROTATIONS)',/)
    2 FORMAT(12F10.6)
    3 FORMAT(/1X,'PROJECTION MATRIX')

C COMPUTE CENTER OF MOLECULE
      CX=ZERO
      CY=ZERO
      CZ=ZERO
      TM=DBLE(NA)
      DO 40 I=1,NA
         CX=CX+XA(I,1)
         CY=CY+XA(I,2)
  40     CZ=CZ+XA(I,3)
      CX=CX/TM
      CY=CY/TM
      CZ=CZ/TM

C TRANSLATIONAL PART OF THE B2 MATRIX
      DO 60  K=1,NA
      DO 60  I=1,3
         L=3*(K-1)+I
  60     B2(I,L)=ONE/TM
C FORM AND INVERT THE INERTIA TENSOR
      DO 62  I=1,3
      DO 62  J=1,6
  62     T0(I,J)=ZERO
      DO 63  I=1,3
  63     T0(I,I+3)=ONE
C
      DO 70  K=1,NA
         T0(1,1)=T0(1,1)+(XA(K,2)-CY)**2+(XA(K,3)-CZ)**2
         T0(2,2)=T0(2,2)+(XA(K,1)-CX)**2+(XA(K,3)-CZ)**2
         T0(3,3)=T0(3,3)+(XA(K,1)-CX)**2+(XA(K,2)-CY)**2
         T0(2,1)=T0(2,1)-(XA(K,1)-CX)*(XA(K,2)-CY)
         T0(3,1)=T0(3,1)-(XA(K,1)-CX)*(XA(K,3)-CZ)
  70     T0(3,2)=T0(3,2)-(XA(K,2)-CY)*(XA(K,3)-CZ)
      T0(1,2)=T0(2,1)
      T0(1,3)=T0(3,1)
      T0(2,3)=T0(3,2)
C
      CALL FLIN(T0,3,3,3,DET)
C
C  FORM THE B2 MATRIX FOR ROTATIONS
      DO 72  I=1,NA
        DO 74  K=1,3
        DO 74  L=1,3
  74    T0(K,L)=ZERO
C
        V(1)=ZERO
        V(2)=(XA(I,3)-CZ)
        V(3)=-(XA(I,2)-CY)
        DO 76  K=1,3
        DO 76  L=1,3
  76    T0(K,1)=T0(K,1)+T0(K,L+3)*V(L)
C
        V(1)=-(XA(I,3)-CZ)
        V(2)=ZERO
        V(3)=(XA(I,1)-CX)
        DO 78  K=1,3
        DO 78  L=1,3
  78    T0(K,2)=T0(K,2)+T0(K,L+3)*V(L)
C
        V(1)=(XA(I,2)-CY)
        V(2)=-(XA(I,1)-CX)
        V(3)=ZERO
        DO 80  K=1,3
        DO 80  L=1,3
  80    T0(K,3)=T0(K,3)+T0(K,L+3)*V(L)
C
        L=3*(I-1)
        DO 82  K=1,3
        B2(K+3,L+1)=T0(K,1)
        B2(K+3,L+2)=T0(K,2)
  82    B2(K+3,L+3)=T0(K,3)
  72  CONTINUE
C
C FORM (BT*B)-1
      DO 90  I=1,6
      DO 90  J=1,12
  90    W(I,J)=ZERO
      DO 92  I=1,6
  92    W(I,I+6)=ONE
      DO 94 I=1,6
      DO 94 J=1,6
      DO 94 K=1,NC
  94    W(I,J)=W(I,J)+B2(I,K)*B2(J,K)
      CALL FLIN(W,6,6,6,DET)
C
      DO 96  I=1,6
      DO 96  J=1,NC
          A(J,I)=ZERO
      DO 96  K=1,6
 96       A(J,I)= A(J,I)+W(I,K+6)*B2(K,J)
      DO 98  I=1,NC
      DO 98  J=1,NC
         P(I,J)=ZERO
      DO 98  K=1,6
 98      P(I,J)=P(I,J)+B2(K,I)* A(J,K)
      DO 100 I=2,NC
      DO 100 J=1,I-1
         P(I,J)=-P(I,J)
 100     P(J,I)=-P(J,I)
      DO 102 I=1,NC
 102     P(I,I)=ONE-P(I,I)
C
      WRITE(IOUT,*)'P2:'
      DO 200 I = 1, NC
      DO 200 J = 1, NC
 200    WRITE(IOUT,201) I,J,P(I,J)
 201    FORMAT (I5,I5,'   ',F12.10)

      RETURN
      END
