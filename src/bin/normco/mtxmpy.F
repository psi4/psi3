      SUBROUTINE MTXMPY(A,NAD,B,NBD,C,NCD,D,NDD,N,IMPY)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(NAD,NAD),B(NBD,NBD),C(NCD,NCD),D(NDD,NDD)
      DATA A00 / 0.0D+00 /
C
      DO 101 I=1,NCD
      DO 101 J=1,NCD
  101 C(I,J)=A00
      GO TO (201,202,203,204,205,206),IMPY
C     C=A*B
  201 DO 102 I=1,N
      DO 102 J=1,N
      DO 102 K=1,N
  102 C(I,J)=C(I,J)+A(I,K)*B(K,J)
      RETURN
C     C=TA*B
  202 DO 103 I=1,N
      DO 103 J=1,N
      DO 103 K=1,N
  103 C(I,J)=C(I,J)+A(K,I)*B(K,J)
      RETURN
C     C=A*TB
  203 DO 104 I=1,N
      DO 104 J=1,N
      DO 104 K=1,N
  104 C(I,J)=C(I,J)+A(I,K)*B(J,K)
      RETURN
C     C=TA*TB
  204 DO 105 I=1,N
      DO 105 J=1,N
      DO 105 K=1,N
  105 C(I,J)=C(I,J)+A(K,I)*B(J,K)
      RETURN
C     C=TA*B*A
  205 DO 106 I=1,N
      DO 106 J=1,N
      D(I,J)=A00
      DO 106 K=1,N
  106 D(I,J)=D(I,J)+A(K,I)*B(K,J)
      DO 107 I=1,N
      DO 107 J=1,N
      DO 107 K=1,N
  107 C(I,J)=C(I,J)+D(I,K)*A(K,J)
      RETURN
C     C=A*B*TA
  206 DO 108 I=1,N
      DO 108 J=1,N
      D(I,J)=A00
      DO 108 K=1,N
  108 D(I,J)=D(I,J)+A(I,K)*B(K,J)
      DO 109 I=1,N
      DO 109 J=1,N
      DO 109 K=1,N
  109 C(I,J)=C(I,J)+D(I,K)*A(J,K)
      RETURN
      END
