C////////////////////////////////////////////////////////////////////
      INTEGER FUNCTION LPRT(K,NPRT)
      DIMENSION N(4)
      NR=NPRT
      N(4)=NR/1000
      NR=NR-N(4)*1000
      N(3)=NR/100
      NR=NR-N(3)*100
      N(2)=NR/10
      NR=NR-N(2)*10
      N(1)=NR
      LPRT=N(K)
      END
