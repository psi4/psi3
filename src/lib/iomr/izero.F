      SUBROUTINE IZERO(A,N)                                             
C                                                                       
      INTEGER A(N)                                                      
C                                                                       
      NBYTES = 4 * N                                                    
cps      CALL VXINIT(A,0,NBYTES)                                           
      DO 1 I=1,N                                                        
           A(I)=0                                                       
    1 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
