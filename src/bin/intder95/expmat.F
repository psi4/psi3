C     //////////////////////////////////////////////////////////////////
      SUBROUTINE EXPMAT(RK,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 K(3,3)
      DIMENSION RK(3),W(3,3),C(3,3),D(3,3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      PARAMETER(XTOL=1.0D-20,MAXIT=50)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
  1   FORMAT(/'CONVERGENCE NOT REACHED IN EXPMAT AFTER',I5,
     $   ' ITERATIONS.'/' CNORM= ',E10.3,' XTOL= ',E10.3)
      DO 10  I=1,3
 10   K(I,I)=ZERO
      K(1,2)=RK(1)
      K(2,1)=-RK(1)
      K(1,3)=RK(2)
      K(3,1)=-RK(2)
      K(2,3)=RK(3)
      K(3,2)=-RK(3)
      DO 20  I=1,3
      DO 30  J=1,3
      C(I,J)=ZERO
 30   W(I,J)=ZERO
      C(I,I)=ONE
 20   W(I,I)=ONE
      ITER=0
C
 100  CNORM=ZERO
      ITER=ITER+1
      DO 50  I=1,3
      DO 50  J=1,3
      D(I,J)=ZERO
      DO 60  L=1,3
 60   D(I,J)=D(I,J)+C(I,L)*K(L,J)
      D(I,J)=D(I,J)/DFLOAT(ITER)
      IF(DABS(D(I,J)).GT.CNORM) CNORM=DABS(D(I,J))
 50   CONTINUE
      DO 70  I=1,3
      DO 70  J=1,3
      W(I,J)=W(I,J)+D(I,J)
 70   C(I,J)=D(I,J)
      IF(ITER.LT.MAXIT.AND.CNORM.GT.XTOL) GO TO 100
      IF(CNORM.GT.XTOL) THEN
      WRITE(IOUT,1) ITER,CNORM,XTOL
      END IF
      RETURN
      END
