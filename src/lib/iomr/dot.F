      double precision FUNCTION DOT(A,B,N)                                               
C                                                                       
      REAL*8 A(N),B(N),T                                                
C                                                                       
      T=0.0D+00                                                         
      DO 1 I=1,N                                                        
           T=T+A(I)*B(I)                                                
    1 CONTINUE                                                          
      DOT=T                                                             
      RETURN                                                            
      END                                                               
