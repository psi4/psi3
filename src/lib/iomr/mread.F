      SUBROUTINE Mread(IA,IADD)                                         
      IMPLICIT INTEGER (A-Z)                                            
      DIMENSION IA(1)                                                   
      COMMON/LOCAT/LOCS(1024)                                           
      COMMON/MFSEC/MFILE,NSECT                                          
C                                                                       
      IADDR=LOCS(IADD)                                                  
      NSECTH=NSECT/2                                                    
      NLEN=LOCS(IADD+NSECTH)                                            
      CALL RREAD(MFILE,IA,NLEN,IADDR)                                   
      RETURN                                                            
      END                                                               
