      SUBROUTINE EIVOUT(A,B,NAD,NBD,M,N,IOUT)                           
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      CHARACTER*6 LINE(21)                                              
      DIMENSION A(NAD,NBD),B(NAD)                                       
      DATA LINE / 21*'------' /                                         
    1 FORMAT(2X,10(7X,I5))                                              
    2 FORMAT(2X,21A6)                                                   
    3 FORMAT(2X,I2,2X,10F12.7)                                          
    4 FORMAT(/,6X,10F12.7)                                              
    5 FORMAT(/)                                                         
C                                                                       
      II=0                                                              
      JJ=0                                                              
  200 II=II+1                                                           
      JJ=JJ+1                                                           
      KK=10*JJ                                                          
      NN=N                                                              
      IF(N.GT.KK) NN=KK                                                 
      LL=2*(NN-II+1)+1                                                  
      WRITE(IOUT,1) (I,I=II,NN)                                         
      WRITE(IOUT,2) (LINE(I),I=1,LL)                                    
      DO 101 I=1,M                                                      
  101 WRITE(IOUT,3) I,(A(I,J),J=II,NN)                                  
      WRITE(IOUT,4) (B(J),J=II,NN)                                      
      IF(N.LE.KK) GO TO 201                                             
      WRITE(IOUT,5)                                                     
      II=KK                                                             
      GO TO 200                                                         
  201 RETURN                                                            
      END                                                               
