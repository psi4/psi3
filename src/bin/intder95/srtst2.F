C     ///////////////////////////////////////////////////////////////
      SUBROUTINE SRTST2(NC,NS,NSX,NSYM,IA,SRA,SRN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER R
      DIMENSION IA(NS,6),SRA(NC,NC,NC),SRN(NC,NC,NC)
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,AVCUT=1.0D-5,ERCUT=1.0D-5)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
    1 FORMAT(//1X,'TESTING OF BR(I,J,K) AND CR(M,N,P) MATRICES',/)
    2 FORMAT(/1X,'R=',I3,' BR:  AVRG=',F12.8,'  RMSDIF=',F12.8)
    3 FORMAT(//1X,'BR(I,J,K) MATRIX FOR COORDINATE ',I5/)
    4 FORMAT(6F12.6)
    5 FORMAT(//1X,'BR(I,J,K) ANALYTIC - BR(I,J,K) NUMERICAL'/)
    6 FORMAT(6F12.6)
    7 FORMAT(1X,' I= ',I5)
    8 FORMAT(1X,' ERROR IN H',3I1)
    9 FORMAT(/1X,'R=',I3,' CR:  AVRG=',F12.8,'  RMSDIF=',F12.8)
   10 FORMAT(//1X,'CR(M,N,P) MATRIX FOR COORDINATE ',I5/)
   11 FORMAT(//1X,'CR(M,N,P) ANALYTIC - CR(M,N,P) NUMERICAL'/)
      IF(NSYM.EQ.0) THEN
         ISCR=ISCR3
      ELSE
         ISCR=ISCR4
      END IF
      WRITE(IOUT,1)
      NA=NC/3
      DO 100  R=1,NS
         IF(IA(R,6).NE.2) GO TO 100
         CALL YIN(NC,NC,SRA,-R,ISCR)
         CALL YIN(NC,NC,SRN,-R,ISCR7)
         RMSDIF=ZERO
         AVRG=ZERO
         NNZ=0
         DO 20  I=1,NC
         DO 20  J=1,NC
         DO 20  K=1,NC
         SRN(I,J,K)=SRA(I,J,K)-SRN(I,J,K)
         IF(DABS(SRA(I,J,K)).GT.AVCUT) THEN
            NNZ=NNZ+1
            AVRG=AVRG+DABS(SRA(I,J,K))
            RMSDIF=RMSDIF+SRN(I,J,K)**2
         END IF
   20    CONTINUE
         IF(NNZ.EQ.0) NNZ=1
         AVRG=AVRG/NNZ
         RMSDIF=DSQRT(RMSDIF)/NNZ
         WRITE(IOUT,2) R,AVRG,RMSDIF
         IF(LPRT(2,NPRT).GE.1) THEN
           WRITE(IOUT,3) R
           DO 30  I=1,NC
           WRITE(IOUT,7) I
           DO 30  J=1,NC
   30      WRITE(IOUT,4) (SRA(I,J,K),K=1,NC)
         END IF
         IF(LPRT(2,NPRT).GE.2) THEN
           WRITE(IOUT,5)
           DO 35  I=1,NC
           WRITE(IOUT,7) I
           DO 35  J=1,NC
   35      WRITE(IOUT,6) (SRN(I,J,K),K=1,NC)
         END IF
         DO 40  I=1,NC
         DO 40  J=1,NC
         DO 40  K=1,NC
   40    SRA(I,J,K)=ZERO
         DO 50  I=1,NC
         DO 50  J=1,NC
         DO 50  K=1,NC
         IF(DABS(SRN(I,J,K)).LT.ERCUT) GO TO 50
         INA=(I-1)/3+1
         JNA=(J-1)/3+1
         KNA=(K-1)/3+1
           DO 45  MM=1,5
         IF(INA.EQ.IA(R,MM)) INA=-MM
         IF(JNA.EQ.IA(R,MM)) JNA=-MM
         IF(KNA.EQ.IA(R,MM)) KNA=-MM
   45      CONTINUE
         INA=-INA
         JNA=-JNA
         KNA=-KNA
         SRA(INA,JNA,KNA)=ONE
   50    CONTINUE
         DO 55  I=1,NA
         DO 55  J=1,I
         DO 55  K=1,J
         IF(SRA(I,J,K).GT.ZERO) THEN
         WRITE(IOUT,8) I,J,K
         END IF
   55    CONTINUE
         CALL YIN(NC,NSX,SRA,R,ISCR)
         CALL YIN(NC,NSX,SRN,R,ISCR7)
         RMSDIF=ZERO
         AVRG=ZERO
         NNZ=0
         DO 60  I=1,NSX
         DO 60  J=1,NSX
         DO 60  K=1,NSX
         SRN(I,J,K)=SRA(I,J,K)-SRN(I,J,K)
         IF(DABS(SRA(I,J,K)).GT.AVCUT) THEN
            NNZ=NNZ+1
            AVRG=AVRG+DABS(SRA(I,J,K))
            RMSDIF=RMSDIF+SRN(I,J,K)**2
         END IF
   60    CONTINUE
         IF(NNZ.EQ.0) NNZ=1
         AVRG=AVRG/NNZ
         RMSDIF=DSQRT(RMSDIF)/NNZ
         WRITE(IOUT,9) R,AVRG,RMSDIF
         IF(LPRT(2,NPRT).GE.1) THEN
           WRITE(IOUT,10) R
           DO 65  I=1,NSX
           WRITE(IOUT,7) I
           DO 65  J=1,NSX
   65      WRITE(IOUT,4) (SRA(I,J,K),K=1,NSX)
         END IF
         IF(LPRT(2,NPRT).GE.2) THEN
           WRITE(IOUT,11)
           DO 70  I=1,NSX
           WRITE(IOUT,7) I
           DO 70  J=1,NSX
   70      WRITE(IOUT,6) (SRN(I,J,K),K=1,NSX)
         END IF
  100 CONTINUE
      END
