C     ///////////////////////////////////////////////////////////////
cbfy  FUNCTION DOTX(A,NA,B,NB,N)
c     IMPLICIT REAL*8 (A-H,O-Z)
c     DIMENSION A(1),B(1)
c     PARAMETER(ZERO=0.0D0)
c     IAPT=1
c     IBPT=1
c     D=ZERO
c     DO 10 I=1,N
c     D=D+A(IAPT)*B(IBPT)
c     IAPT=IAPT+NA
c     IBPT=IBPT+NB
c10   CONTINUE
c     DOTX=D
c     RETURN
cbfy  END
C     ///////////////////////////////////////////////////////////////
      SUBROUTINE VECPRO(U,V,W)
      REAL*8 U(3),V(3),W(3)
      W(1)=U(2)*V(3)-V(2)*U(3)
      W(2)=U(3)*V(1)-U(1)*V(3)
      W(3)=U(1)*V(2)-U(2)*V(1)
      RETURN
      END
