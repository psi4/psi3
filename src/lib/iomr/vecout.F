      SUBROUTINE VECOUT(C,EIGVAL,NUMBLK,NUM,NUMSYM,NSYM,IOUT)
C
C***PURPOSE: TO PRINT A SYMMETRY-BLOCKED MATRIX AND ITS ASSOCIATED
C            EIGENVALUES
C
C PAUL SAXE                 1 NOVEMBER 1984               LANL
C
      IMPLICIT INTEGER (A-Z)
C
      REAL*8 C(NUMBLK),EIGVAL(NUM)
      INTEGER NUMSYM(NSYM)
      CHARACTER*1 LINE
c     The number of data per line
      integer ndpl
#include <misc.h>
C
      DATA LINE /'-'/
c
      ndpl = (linewd-7)/11

      if ((ndpl.le.0).or.(ndpl.gt.200)) then
        write(6,*) 'ERROR: bogus ndpl in vecout: check ''LINEWIDTH'''
        call qabort
        endif
c
C
      PTEIG=0
      PTSQR=0
      DO 100 SYM=1,NSYM
         NUMS=NUMSYM(SYM)
         IF (NUMS.LE.0) GO TO 100
         WRITE (IOUT,90) SYM,NUMS
   90    FORMAT (//,' SYMMETRY BLOCK',I2,' OF DIMENSION',I4,/)
         MX=0
    1    CONTINUE
            MN=MX+1
            MX=MIN(MX+6,NUMS)
            WRITE (IOUT,2) (IQ,IQ=MN,MX)
    2       FORMAT (/,1X,10(I6,5X))
            WRITE (IOUT,3) (LINE,IQ=1,(MX-MN+1)*11)
    3       FORMAT (1X,80A1)
            DO 5 I=1,NUMS
               WRITE (IOUT,4) (C(PTSQR+I+IQ),IQ=(MN-1)*NUMS,(MX-1)*NUMS,
     #                                                             NUMS)
    4          FORMAT (1X,7F11.6)
    5       CONTINUE
            WRITE (IOUT,6) (EIGVAL(PTEIG+IQ),IQ=MN,MX)
    6       FORMAT (/,1X,7F11.6)
C
         IF (MX.LT.NUMS) GO TO 1
         PTEIG=PTEIG+NUMS
         PTSQR=PTSQR+NUMS**2
  100 CONTINUE
C
C
      RETURN
      END
