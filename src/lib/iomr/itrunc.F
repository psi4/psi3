clj122490 XLF (AIX) uses itrunc for internal purposes
      INTEGER FUNCTION ITRNC(N)                                        
CC    ITRUNC=N/192*192                                                  
      ITRNC=N/1024*1024                                                
      RETURN                                                            
      END                                                               
