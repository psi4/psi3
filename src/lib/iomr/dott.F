      FUNCTION DOTT(A,NA,B,NB,N)                                        
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(1),B(1)                                               
      DATA ZERO / 0.0D+00 /                                             
C                                                                       
      IAPT=1                                                            
      IBPT=1                                                            
      D=ZERO                                                            
      DO 10 I=1,N                                                       
      D=D+A(IAPT)*B(IBPT)                                               
      IAPT=IAPT+NA                                                      
      IBPT=IBPT+NB                                                      
 10   CONTINUE                                                          
      DOTT=D                                                            
      RETURN                                                            
      END                                                               
