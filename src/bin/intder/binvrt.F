C     ////////////////////////////////////////////////////////////
      SUBROUTINE BINVRT(NS,NC,XMASS,B,D,A,IFLAG,TOLINV,ITST)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(NC,NC),A(NC,NC),D(NC,1),XMASS(1)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
 1    FORMAT(//' B*BT MATRIX FOR (SYMMETRY) INTERNAL COORDINATES',/)
 2    FORMAT(8F12.6)
 3    FORMAT(/1X,'DETERMINANT OF B*BT MATRIX=',G12.5)
 4    FORMAT(//,' A MATRIX FOR (SYMMETRY) INTERNAL COORDINATES',/)
 5    FORMAT(8F12.6)
      IFLAG=0
      DO 10  J=1,NS
      DO 10  I=1,NS
      D(I,J)=ZERO
      D(I,J+NS)=ZERO
      DO 15  K=1,NC
      IK=(K-1)/3+1
 15   D(I,J)=D(I,J)+B(I,K)*B(J,K)/XMASS(IK)
 10   CONTINUE
      IF(LPRT(1,NPRT).GE.4) THEN
         WRITE(IOUT,1)
         DO 20  I=1,NS
         WRITE(IOUT,*) 'I=',I
 20      WRITE(IOUT,2) (D(I,J),J=1,NS)
      END IF
      DO 25  K=1,NS
 25   D(K,K+NS)=ONE
      CALL FLIN(D,NC,NS,NS,DET)
      WRITE(IOUT,3) DET
      IF(DET.EQ.ZERO) THEN
         IFLAG=1
         RETURN
      END IF
      DO 35  J=1,NS
      DO 35  I=1,NC
      A(I,J)=ZERO
      DO 40  K=1,NS
 40   A(I,J)=A(I,J)+B(K,I)*D(K,J+NS)
      IK=(I-1)/3+1
      A(I,J)=A(I,J)/XMASS(IK)
 35   CONTINUE
      IF(LPRT(1,NPRT).GE.3) THEN
         WRITE(IOUT,4)
         DO 45  I=1,NC
         WRITE(IOUT,*) 'I=',I
 45      WRITE(IOUT,5) (A(I,J),J=1,NS)
      END IF
      IF(ITST.EQ.0) RETURN
      DO 50  J=1,NS
      DO 50  I=1,NS
      W=ZERO
      DO 55  K=1,NC
 55   W=W+B(I,K)*A(K,J)
      D(I,J)=W
      IF(I.EQ.J) W=W-ONE
      IF(DABS(W).GT.TOLINV) IFLAG=2
 50   CONTINUE
      RETURN
      END
