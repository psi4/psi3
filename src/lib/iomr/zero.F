      SUBROUTINE ZERO(A,N)                                              
C                                                                       
      REAL*8 A(N)                                                       
C                                                                       
      NBYTES = 8 * N                                                    
cps      CALL VXINIT(A,0,NBYTES)                                           
      DO 1 I=1,N                                                        
           A(I)=0.0D+00                                                 
    1 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
