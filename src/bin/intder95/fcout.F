C     ///////////////////////////////////////////////////////////
      SUBROUTINE FCOUT(NC,NS,NSY,NEQ,NDER,NINV,NVEC,NRUN,SS,ENERGY,
     $                 F1,F2,F3,F4)
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER P,Q,QQ
      DIMENSION SS(NS)
      DIMENSION F1(NC),F2(NC,NC),F3(NC,NC,NC),F4(NC,NC,NC,NC)
      DIMENSION DIP(3)
CSA   PARAMETER(CF1=0.121376050883D0,CF2=0.064229421739D0)
CSA   PARAMETER(CF3=0.033988736554D0,CF4=0.017986059681D0)
      PARAMETER(ONE=1.0D0)
      PARAMETER(F1CUT=1.0D-8,F2CUT=1.0D-6,F3CUT=1.0D-5,F4CUT=1.0D-5)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
      COMMON /PHYCON/ BOHR,DEBYE,HART,WAVE0,CINT
      COMMON /DIPOLE/ DIP,QQ
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
 11   FORMAT(8F10.7)
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
      CF1=BOHR/HART
      CF2=CF1*BOHR
      CF3=CF2*BOHR
      CF4=CF3*BOHR
      NA=NC/3
      NNINV=NINV
CSA
      IF(ABS(NINV).NE.3) GO TO 55
         NNINV=0
         WRITE(IOUT,12)
         IOUT1=I12
         IOUT2=I16
         IOUT3=I21
         IOUT4=I25
      GO TO 60
CSA
 55   IF(NNINV.LE.0) THEN
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
  60  IF(NEQ.NE.0) THEN
         IF(NNINV.LE.0) THEN
             XX=ONE
             IF(NVEC.EQ.1) XX=CF1
         ELSE
             XX=CF1
             IF(NVEC.EQ.1) XX=ONE
         END IF
         DO 200 M=1,NSY
 200     F1(M)=F1(M)*XX
      END IF
      IF(NNINV.LE.0) GO TO 220
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
                GO TO 225
           END IF
           IF(NNINV.GT.0) THEN
                WRITE(IOUT1,15) (F1(M),M=1,NSY)
           ELSE
                WRITE(IOUT1,17) NA,ENERGY
                DO 230  M=1,NSY
 230            WRITE(IOUT1,13) SS(M),F1(M)
           END IF
      END IF
 225  IF(NDER.GE.2) THEN
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
           END IF
      END IF
      IF(NDER.GE.3) THEN
           LLL=LPRT(1,NPRT)
           IF(LLL.GE.1) THEN
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
