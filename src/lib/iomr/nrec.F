      INTEGER FUNCTION NREC(NLEN)                                       
      IMPLICIT INTEGER (A-Z)                                            
CC    IADD=(NLEN+8188)/8189                                             
      IADD=(NLEN+255)/256                                               
      NREC=IADD                                                         
      RETURN                                                            
      END                                                               
