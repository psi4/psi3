C     ///////////////////////////////////////////////////////////
      SUBROUTINE FCOUT(NC,NS,NSY,NEQ,NDER,NINV,NVEC,NRUN,SS,ENERGY,
     $                 F1,F2,F3,F4)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,Q,QQ
cbfy
croag  the ltext parameter is the maximum number of lines in file12
croag      PARAMETER (LTEXT = 1200)
      PARAMETER (LTEXT = 3000)
      CHARACTER TEXT(LTEXT)*80, LABEL*80, WORD*4
cbfy
      DIMENSION SS(NS)
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)
      DIMENSION DIP(3)
      PARAMETER(CF1=0.121376050883D0,CF2=0.064229421739D0)
      PARAMETER(CF3=0.033988736554D0,CF4=0.017986059681D0)
      PARAMETER(ONE=1.0D0)
      PARAMETER(F1CUT=1.0D-8,F2CUT=1.0D-6,F3CUT=1.0D-5,F4CUT=1.0D-5)
      COMMON /IO/ IIN1,IOUT,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10
      COMMON /DIPOLE/DIP,QQ
 1    FORMAT(/' GRADIENTS')
 2    FORMAT(6F12.8)
 3    FORMAT(/' QUADRATIC FORCE CONSTANTS'/)
 4    FORMAT(6F12.6)
 5    FORMAT(/' CUBIC FORCE CONSTANTS'/)
 6    FORMAT('M=',I5)
 7    FORMAT(/' QUARTIC FORCE CONSTANTS'/)
 8    FORMAT('M=',I5,'  N=',I5)
 9    FORMAT(//' FINAL INTERNAL (SYMMETRY) COORDINATE FORCE CONSTANTS')
 10   FORMAT(7F10.6)
cbfy 11   FORMAT(8F10.7)
 11   FORMAT(8F10.6)
 12   FORMAT(//' FINAL CARTESIAN COORDINATE FORCE CONSTANTS')
 13   FORMAT(3F20.10)
 14   FORMAT(5I5)
 15   FORMAT(20X,3F20.10)
 16   FORMAT(3F20.10)
 17   FORMAT(I5,12X,F20.10)
 18   FORMAT(//' PRINTING OF CUBIC FORCE CONSTANTS',
     $         ' TO OUTPUT FILE SUPPRESSED.'/' LPRT(1,NPRT)=',I3)
 19   FORMAT(//' PRINTING OF QUARTIC FORCE CONSTANTS',
     $         ' TO OUTPUT FILE SUPPRESSED.'/' LPRT(1,NPRT)=',I3)
 20   FORMAT(I5,15X,F20.10)
 21   FORMAT(2I5,10X,F20.10)
 22   FORMAT(3I5,5X,F20.10)
 23   FORMAT(4I5,F20.10)
 24   FORMAT(2I5,3F20.10)
cbfy
 25   FORMAT(A80)
 26   FORMAT(A4)
cbfy
      NA=NC/3
      IF(NINV.LE.0) THEN
         WRITE(IOUT,9)
         IOUT1=I12
         IOUT2=I16
         IOUT3=I21
         IOUT4=I25
         IOUT5=I18
      ELSE
         WRITE(IOUT,12)
         IOUT1=I11
         IOUT2=I15
         IOUT3=I20
         IOUT4=I24
         IOUT5=I17
      END IF
      IF(NEQ.NE.0) THEN
         IF(NINV.LE.0) THEN
             XX=ONE
             IF(NVEC.EQ.1) XX=CF1
         ELSE
             XX=CF1
             IF(NVEC.EQ.1) XX=ONE
         END IF
         DO 200 M=1,NSY
 200     F1(M)=F1(M)*XX
      END IF
      IF(NINV.LE.0) GO TO 220
      IF(NDER.GE.2) THEN
      DO 205  N=1,NSY
      DO 205  M=1,NSY
 205  F2(M,N)=F2(M,N)*CF2
      END IF
      IF(NDER.GE.3) THEN
      DO 210  P=1,NSY
      DO 210  N=1,NSY
      DO 210  M=1,NSY
 210  F3(M,N,P)=F3(M,N,P)*CF3
      END IF
      IF(NDER.GE.4) THEN
      DO 215  Q=1,NSY
      DO 215  P=1,NSY
      DO 215  N=1,NSY
      DO 215  M=1,NSY
 215  F4(M,N,P,Q)=F4(M,N,P,Q)*CF4
      END IF
 220  CONTINUE
      IF(NEQ.NE.0) THEN
           WRITE(IOUT,1)
           WRITE(IOUT,2) (F1(M),M=1,NSY)
           IF(NVEC.EQ.1) THEN
                IF(NRUN.LE.1) WRITE(IOUT5,24) NA,QQ,(DIP(KK),KK=1,3)
                WRITE(IOUT5,16) (F1(M),M=1,NSY)
                GO TO 235
           END IF
           IF(NINV.GT.0) THEN
cbfy
c               First, skip over the gradients
                DO 221  M=1,NA
 221            READ(IOUT1,25) LABEL
c               Save the rest of the file in a text buffer
                N=1
 222            CONTINUE
                IF (N.GT.LTEXT) THEN
             WRITE(IOUT,*)'FCOUT: ERROR: The text buffer is too small,'
             WRITE(IOUT,*)'       increase LTEXT to the number of lines'
             WRITE(IOUT,*)'       in file12a.dat and recompile. The '
             WRITE(IOUT,*)'       current limit is 3000 lines of data.'
             WRITE(IOUT,*)'Subroutine FCOUT: Aborting.'
                   STOP
                   ENDIF
clj030589
clj Currently some FILE11's (those produced by SCF DERIV) do not
clj use 80 columns for all records.  Others (CI CC NGRCPHF) do use
clj 80 columns for all records.  The following patch will still give
clj error messages, but no data in file11 will be lost.
clj             READ(IOUT1,25,END=223) TEXT(N)
                READ(IOUT1,8891,END=223) (TEXT(N)(ihi:ihi),ihi=1,80)
 8891           format(80a1)
                N=N+1
                GO TO 222
 223            DO 224  M=1,N+NA
 224            BACKSPACE IOUT1
c               Write out new gradients
                WRITE(IOUT1,15) (F1(M),M=1,NSY)
c               Copy the text buffer back to File11
                DO 225  M=1,N-1
 225            WRITE(IOUT1,25) TEXT(M)
                DO 226  M=1,N-1
 226            BACKSPACE IOUT1
cbfy
           ELSE
                WRITE(IOUT1,17) NA,ENERGY
                DO 230  M=1,NSY
 230            WRITE(IOUT1,13) SS(M),F1(M)
           END IF
      END IF
 235  IF(NDER.GE.2) THEN
           WRITE(IOUT,3)
           DO 100  M=1,NSY
 100       WRITE(IOUT,4) (F2(M,N),N=1,NSY)
           WRITE (IOUT2,14) NA,NA*6
           WRITE (IOUT2,13) ((F2(M,N),N=1,NSY),M=1,NSY)
           IF(LPRT(4,NPRT).GE.4)
     $             WRITE(ICHECK,10)((F2(M,N),M=N,NSY),N=1,NSY)
           IF(LPRT(4,NPRT).GE.3) THEN
           DO 102  M=1,NSY
 102       WRITE(ICHECK,11)(F2(M,N),N=1,M)
cbfy
c  Write force constants to BMAT file
              IBMAT = 40
cets122690 write force constants to file fconst.dat
              call ffile(ibmat,'fconst',0)
c          Read until STOP
cets 111          READ(IBMAT,26,END=112) WORD
cets              IF (WORD.EQ.'STOP') GO TO 113
cets              GO TO 111
cets 112          WRITE(IOUT,*) 'FCOUT: No STOP found in BMAT file.'
cets              WRITE(IOUT,*) 'FCOUT: Aborting.'
cets              STOP
cets 113          BACKSPACE IBMAT
c          Read rest of text into buffer
cets              K=1
cets 114          READ(IBMAT,25,END=115) TEXT(K)
cets              K=K+1
cets              GO TO 114
cets 115          CONTINUE
              rewind(ibmat)
cets              call bottom(ibmat)
c          Backspace to force constants
cets              DO 116 M=1,K - 1 + ((NSY/8 + 1)*(NSY - (NSY/8)*4))
cets 116          BACKSPACE IBMAT
c          Now write new force constants and rest of text to BMAT
              DO 117  M=1,NSY
 117          WRITE(IBMAT,11)(F2(M,N),N=1,M)
cets              DO 118 M=1,K-1
cets 118          WRITE(IBMAT,25) TEXT(M)
cbfy
           END IF
      END IF
      IF(NDER.GE.3) THEN
           IF(LPRT(1,NPRT).GE.1) THEN
             WRITE(IOUT,5)
             DO 105  M=1,NSY
             WRITE(IOUT,6) M
             DO 105  N=1,NSY
 105         WRITE(IOUT,4) (F3(M,N,P),P=1,NSY)
           ELSE
             WRITE(IOUT,18) LPRT(1,NPRT)
           END IF
           WRITE (IOUT3,14) NA,(NSY*(NSY+1)*(NSY+2))/6
           WRITE (IOUT3,13) (((F3(M,N,P),P=1,N),N=1,M),M=1,NSY)
      END IF
      IF(NDER.GE.4) THEN
           IF(LPRT(1,NPRT).GE.1) THEN
             WRITE(IOUT,7)
             DO 110  M=1,NSY
             DO 110  N=1,NSY
             WRITE(IOUT,8) M,N
             DO 110  P=1,NSY
 110         WRITE(IOUT,4) (F4(M,N,P,Q),Q=1,NSY)
           ELSE
             WRITE(IOUT,19) LPRT(1,NPRT)
           END IF
      WRITE (IOUT4,14) NA,(NSY*(NSY+1)*(NSY+2)*(NSY+3))/24
      WRITE (IOUT4,13) ((((F4(M,N,P,Q),Q=1,P),P=1,N),N=1,M),M=1,NSY)
      END IF
      IF(LPRT(4,NPRT).NE.2) RETURN
      NZ=0
      IF(NEQ.NE.0) THEN
         DO 300  M=1,NSY
       IF(DABS(F1(M)).GT.F1CUT) WRITE(ICHECK,20) M,F1(M)
 300     CONTINUE
         WRITE(ICHECK,20) NZ
      END IF
      IF(NDER.GE.2) THEN
         DO 310  M=1,NSY
         DO 310  N=1,M
       IF(DABS(F2(M,N)).GT.F2CUT) WRITE(ICHECK,21) M,N,F2(M,N)
 310     CONTINUE
         WRITE(ICHECK,20) NZ
      END IF
      IF(NDER.GE.3) THEN
         DO 320  M=1,NSY
         DO 320  N=1,M
         DO 320  P=1,N
       IF(DABS(F3(M,N,P)).GT.F3CUT) WRITE(ICHECK,22) M,N,P,F3(M,N,P)
 320     CONTINUE
         WRITE(ICHECK,20) NZ
      END IF
      IF(NDER.GE.4) THEN
         DO 330  M=1,NSY
         DO 330  N=1,M
         DO 330  P=1,N
         DO 330  Q=1,P
       IF(DABS(F4(M,N,P,Q)).GT.F4CUT)
     $     WRITE(ICHECK,23) M,N,P,Q,F4(M,N,P,Q)
 330     CONTINUE
         WRITE(ICHECK,20) NZ
      END IF
      RETURN
      END
