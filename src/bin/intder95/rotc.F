C     //////////////////////////////////////////////////////////////////
      SUBROUTINE ROTC(NA,NAD,XMASS,XA,XD,IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XMASS(NA),XA(NAD,3),XD(NAD,3)
      DIMENSION W(3,3),XX(3,6),B(3),RK(3),RCM(3)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      PARAMETER(MAXIT=40,TOLDISP=1.0D-12)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(/'ROTATION OF COORDINATES TO SATISFY ECKART CONDITIONS'/)
 2    FORMAT('CONVERGENCE NOT REACHED IN ROTC')
 3    FORMAT('XDISP= ',E10.3,' TOLDISP= ',E10.3)
 4    FORMAT('K VECTOR'/3E12.5)
 5    FORMAT('CENTER OF MASS'/3F20.10)
      ITER=0
      IFLAG=0
      RMT=ZERO
      DO 15  I=1,NA
 15   RMT=RMT+XMASS(I)
      DO 25  J=1,3
      RCM(J)=ZERO
      DO 20  I=1,NA
 20   RCM(J)=RCM(J)+XA(I,J)*XMASS(I)
 25   RCM(J)=RCM(J)/RMT
      DO 30  I=1,NA
      DO 30  J=1,3
      XA(I,J)=XA(I,J)-RCM(J)
 30   XD(I,J)=XD(I,J)-RCM(J)
      WRITE(IOUT,1)
      WRITE(IOUT,5) (RCM(I),I=1,3)
C
  50  ITER=ITER+1
      DO 100  I=1,3
      B(I)=ZERO
      DO 100  J=1,3
      XX(I,J+3)=ZERO
 100  W(I,J)=ZERO
      DO 110  I=1,3
      DO 110  J=I,3
      DO 110  K=1,NA
 110  W(I,J)=W(I,J)+XA(K,I)*XD(K,J)*XMASS(K)
      XX(1,1)=W(3,1)
      XX(2,1)=W(3,2)
      XX(3,1)=-(W(1,1)+W(2,2))
      XX(1,2)=-W(2,1)
      XX(2,2)=W(1,1)+W(3,3)
      XX(3,2)=-W(2,3)
      XX(1,3)=-(W(2,2)+W(3,3))
      XX(2,3)=W(1,2)
      XX(3,3)=W(1,3)
      XX(1,4)=ONE
      XX(2,5)=ONE
      XX(3,6)=ONE
      CALL FLIN(XX,3,3,3,DD)
      DO 120  I=1,NA
      B(1)=B(1)+XMASS(I)*(XD(I,2)*XA(I,3)-XD(I,3)*XA(I,2))
      B(2)=B(2)+XMASS(I)*(XD(I,3)*XA(I,1)-XD(I,1)*XA(I,3))
 120  B(3)=B(3)+XMASS(I)*(XD(I,1)*XA(I,2)-XD(I,2)*XA(I,1))
      DO 130  I=1,3
      RK(I)=ZERO
      DO 130  J=1,3
 130  RK(I)=RK(I)+XX(I,J+3)*B(J)
      CALL EXPMAT(RK,W)
      XDISP=ZERO
      DO 200  I=1,NA
      DO 205  J=1,3
      XY=ZERO
      DO 210  K=1,3
  210 XY=XY+W(J,K)*XD(I,K)
      XZ=DABS(XY-XD(I,J))
      IF(XZ.GT.XDISP) XDISP=XZ
  205 B(J)=XY
      DO 215  J=1,3
  215 XD(I,J)=B(J)
  200 CONTINUE
      IF(XDISP.GT.TOLDISP.AND.ITER.LT.MAXIT) GO TO 50
      DO 300 I=1,NA
      DO 300 J=1,3
      XA(I,J)=XA(I,J)+RCM(J)
 300  XD(I,J)=XD(I,J)+RCM(J)
      IF(XDISP.GT.TOLDISP) THEN
         WRITE(IOUT,2)
         WRITE(IOUT,3) XDISP,TOLDISP
         WRITE(IOUT,4) (RK(I), I=1,3)
         IFLAG=1
         RETURN
      END IF
      RETURN
      END
