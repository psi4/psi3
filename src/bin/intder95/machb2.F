C     ///////////////////////////////////////////////////////////////
      SUBROUTINE MACHB2(NAD,NC,XA,B2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XA(NAD,3),B2(6,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      NA=NC/3
      T=DSQRT(DBLE(NA))
      DO 10  J=1,NC
      DO 10  I=1,6
 10   B2(I,J)=ZERO
      DO 15  K=1,NA
      DO 15  I=1,3
      L=3*(K-1)+I
 15   B2(I,L)=ONE/T
      DO 20  K=1,NA
      L=3*K
      B2(4,L-1)=-XA(K,3)
      B2(4,L)=XA(K,2)
      B2(5,L)=-XA(K,1)
      B2(5,L-2)=XA(K,3)
      B2(6,L-1)=XA(K,1)
 20   B2(6,L-2)=-XA(K,2)
      DO 25  I=4,6
      T=ZERO
      DO 30  J=1,NC
 30   T=T+B2(I,J)*B2(I,J)
      T=DSQRT(T)
      DO 35  J=1,NC
 35   B2(I,J)=B2(I,J)/T
 25   CONTINUE
      RETURN
      END
