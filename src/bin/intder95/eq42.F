C/////////////////////////////////////////////////////////////////////
      SUBROUTINE EQ42 (NC,NA,XA,P2,P3,P4,XAS)
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      DIMENSION P2(NC,NC),P3(NC,NC,NC),P4(NC,NC,NC,NC)
      DIMENSION XA(NA,3), XAS(NC), X(3)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      PARAMETER (ZERO=0.0D0)
      WRITE (IOUT,1)
   1  FORMAT ('APPROX PROJECTED CART COORD DISP VIA EQ42')

      DO 2000 I = 1,NC
 2000    XAS(I) = ZERO

      XAS(7)=.1091033790
      XAS(8)=.1091033790
      XAS(9)=.1091033790

C     DO 10 I = 1, NA
C     DO 10 J = 1, 3
C          K = (I-1)*3+J
C 10       XAS(K) = XA(I,J)

      DO 100 I = 1, NA
      DO  90 J = 1, 3
        K = (I-1)*3+J
        T1 = ZERO
        T2 = ZERO
        T3 = ZERO

        DO 50 L = 1, NC
          T1 = T1 + P2(K,L)*XAS(L)

          DO 50 M = 1, NC
            T2 = T2 + P3(K,L,M)*XAS(L)*XAS(M)

             DO 50 N = 1, NC
               T3 = T3 + P4(K,L,M,N)*XAS(L)*XAS(M)*XAS(N)

  50    CONTINUE

  90    X(J) = T1 + T2/2 + T3/6

 100  WRITE (IOUT,101) I,X(1),X(2),X(3)
 101  FORMAT (I3,3F20.10)

      RETURN
      END
