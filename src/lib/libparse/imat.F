      SUBROUTINE imat(A,NAD,NBD,M,N,IOUT)
      IMPLICIT integer (a-z)
      CHARACTER*3 LINE
      DIMENSION A(NAD,NBD)
c     The number of data per line
      integer ndpl
#include "misc.h"
    1 FORMAT(7X,10(1X,I5))
    2 FORMAT(7X,21A3)
    3 FORMAT(2X,I3,2X,20i6)
    4 FORMAT(/)
C                                                                       
      DATA LINE / '---' /

      ndpl = (linewd-7)/6

      if ((ndpl.le.0).or.(ndpl.gt.200)) then
        write(6,*) 'ERROR: bogus ndpl in imat: check ''LINEWIDTH'''
        call qabort
        endif
      II=0
      JJ=0
  200 II=II+1
      JJ=JJ+1
      KK=10*JJ
      NN=N
      IF(N.GT.KK) NN=KK
      LL=2*(NN-II+1)+1
      WRITE (IOUT,1) (I,I=II,NN)
      WRITE (IOUT,2) (LINE,I=1,LL)
      DO 101 I=1,M
      WRITE (IOUT,3) I,(A(I,J),J=II,NN)
  101 CONTINUE
      IF(N.LE.KK) GO TO 201
      WRITE (IOUT,4)
      II=KK
      GO TO 200
  201 RETURN
      END
