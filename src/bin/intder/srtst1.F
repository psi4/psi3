C     ///////////////////////////////////////////////////////////////
      SUBROUTINE SRTST1(NC,NS,NSX,NSYM,IA,SRA,SRN)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER R
      DIMENSION IA(NS,6),SRA(NC,NC),SRN(NC,NC)
      PARAMETER(ZERO=0.0D0,AVCUT=1.0D-6)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
    1 FORMAT(//1X,'TESTING OF SR(I,J) AND X(M,N) MATRICES')
    2 FORMAT(/1X,'R=',I3,' SR:  AVRG=',F12.8,'  RMSDIF=',F12.8)
    3 FORMAT(//1X,'SR(I,J) MATRIX FOR COORDINATE ',I5/)
    4 FORMAT(6F12.6)
    5 FORMAT(//1X,'SR(I,J) ANALYTIC - SR(I,J) NUMERICAL'/)
    6 FORMAT(6F12.6)
    7 FORMAT(/1X,'R=',I3,' X:  AVRG=',F12.8,'  RMSDIF=',F12.8)
    8 FORMAT(//1X,'X(M,N) MATRIX FOR COORDINATE ',I5/)
    9 FORMAT(//1X,'X(M,N) ANALYTIC - X(M,N) NUMERICAL'/)
      IF(NSYM.EQ.0) THEN
         ISCR=ISCR1
      ELSE
         ISCR=ISCR2
      END IF
      WRITE(IOUT,1)
      DO 100  R=1,NS
         IF(IA(R,6).NE.1) GO TO 100
         CALL XIN(NC,NC,SRA,-R,ISCR)
         CALL XIN(NC,NC,SRN,-R,ISCR6)
         RMSDIF=ZERO
         AVRG=ZERO
         NNZ=0
         DO 20  I=1,NC
         DO 20  J=1,NC
         SRN(I,J)=SRA(I,J)-SRN(I,J)
         IF(DABS(SRA(I,J)).GT.AVCUT) THEN
            NNZ=NNZ+1
            AVRG=AVRG+DABS(SRA(I,J))
            RMSDIF=RMSDIF+SRN(I,J)**2
         END IF
   20    CONTINUE
         IF(NNZ.EQ.0) NNZ=1
         AVRG=AVRG/NNZ
         RMSDIF=DSQRT(RMSDIF)/NNZ
         WRITE(IOUT,2) R,AVRG,RMSDIF
         IF(LPRT(2,NPRT).GE.1) THEN
           WRITE(IOUT,3) R
           DO 30  I=1,NC
   30      WRITE(IOUT,4) (SRA(I,J),J=1,NC)
         END IF
         IF(LPRT(2,NPRT).GE.2) THEN
           WRITE(IOUT,5)
           DO 35  I=1,NC
   35      WRITE(IOUT,6) (SRN(I,J),J=1,NC)
         END IF
         CALL XIN(NC,NC,SRA,R,ISCR)
         CALL XIN(NC,NC,SRN,R,ISCR6)
         RMSDIF=ZERO
         AVRG=ZERO
         NNZ=0
         DO 40  I=1,NSX
         DO 40  J=1,NSX
         SRN(I,J)=SRA(I,J)-SRN(I,J)
         IF(DABS(SRA(I,J)).GT.AVCUT) THEN
            NNZ=NNZ+1
            AVRG=AVRG+DABS(SRA(I,J))
            RMSDIF=RMSDIF+SRN(I,J)**2
         END IF
   40    CONTINUE
         IF(NNZ.EQ.0) NNZ=1
         AVRG=AVRG/NNZ
         RMSDIF=DSQRT(RMSDIF)/NNZ
         WRITE(IOUT,7) R,AVRG,RMSDIF
         IF(LPRT(2,NPRT).GE.1) THEN
           WRITE(IOUT,8) R
           DO 50  I=1,NSX
   50      WRITE(IOUT,4) (SRA(I,J),J=1,NSX)
         END IF
         IF(LPRT(2,NPRT).GE.2) THEN
           WRITE(IOUT,9)
           DO 55  I=1,NSX
   55      WRITE(IOUT,6) (SRN(I,J),J=1,NSX)
         END IF
  100 CONTINUE
      END
