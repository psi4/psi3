C     /////////////////////////////////////////////////////////////////
      SUBROUTINE BINVRT(NS,NC,XMASS,B,D,A,AV,V,IFLAG,TOLINV,ITST)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(NC,NC),A(NC,NC),D(NC,1),XMASS(1)
      DIMENSION AV(1),V(NC,3)
C  Note that AV is stored in the second part of D.
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(//' B*BT MATRIX FOR (SYMMETRY) INTERNAL COORDINATES',/)
 2    FORMAT(8F12.6)
 3    FORMAT(/1X,'DETERMINANT OF B*BT MATRIX=',G12.5)
 4    FORMAT(//,' A MATRIX FOR (SYMMETRY) INTERNAL COORDINATES',/)
 5    FORMAT(8F12.6)
 6    FORMAT(//,' NORMALIZED OVERLAP MATRIX FOR INTERNAL COORDINATES',/)
 7    FORMAT(8F12.6)
 8    FORMAT(//,' EIGENVECTORS AND EIGENVALUES OF',
     $          ' NORMALIZED OVERLAP MATRIX')
 9    FORMAT(//,' DETERMINANT OF OVERLAP MATRIX   ', E12.5)
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
CSA   IF(ITST.GT.1) GO TO 60
      DO 50  J=1,NS
      DO 50  I=1,NS
      W=ZERO
      DO 55  K=1,NC
 55   W=W+B(I,K)*A(K,J)
      D(I,J)=W
      IF(I.EQ.J) W=W-ONE
      IF(DABS(W).GT.TOLINV) IFLAG=2
 50   CONTINUE
 60   IF(ITST.LT.2) GO TO 80
      DO 62  J=1,NS
      DO 62  I=1,NS
      D(I,J)=ZERO
      DO 63  K=1,NC
      IK=(K-1)/3+1
 63   D(I,J)=D(I,J)+B(I,K)*B(J,K)/XMASS(IK)
 62   CONTINUE
      DO 64  I=1,NS
      DO 64  J=1,NS
      IF(I.EQ.J) GO TO 64
      D(I,J)=D(I,J)/DSQRT(D(I,I)*D(J,J))
 64   CONTINUE
      DO 65  I=1,NS
 65   D(I,I)=ONE
      WRITE(IOUT,6)
      DO 66  I=1,NS
 66   WRITE(IOUT,7) (D(I,J),J=1,NS)
      II=0
      DO 67  I=1,NS
      DO 67  J=1,I
      II=II+1
 67   AV(II)=D(I,J)
      NV=NS*(NS+1)/2
      CALL RSP(NC,NS,NV,AV,V(1,1),1,D,V(1,2),V(1,3))
      WRITE(IOUT,8)
      CALL TABLE5(NC,NS,NS,V(1,1),D)
      DETS=ONE
      DO  68  I=1,NS
 68   DETS=DETS*V(I,1)
      WRITE(IOUT,9) DETS
 80   CONTINUE
      RETURN
      END
