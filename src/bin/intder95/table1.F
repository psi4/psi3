C//////////////////////////////////////////////////////////////////////
      SUBROUTINE TABLE1(NN,N,A)
      REAL*8 A(NN,NN)
      COMMON /IO/ IIN1,IOUT,IIN2,ICHECK,NPRT,
     $   I11,I15,I17,I20,I24,
     $   I12,I16,I18,I21,I25,
     $   I31,I32,I33,I35,I36,I37,
     $   ISCR1,ISCR2,ISCR3,ISCR4,ISCR5,ISCR6,ISCR7,ISCR8,ISCR9,ISCR10,
     $   ISCR11,ISCR12,ISCR13,ISCR14
 1    FORMAT(6(I5,I3,F10.6))
      DO 10 I=1,N
 10   WRITE(IOUT,1) (I,J,A(I,J),J=1,N)
      END
