      SUBROUTINE DNODE(A,RT,K)                                          
C                                                                       
C                                                                       
C     RETURNS IN RT(I) THE ITH ROOT OF A POLYNOMIAL OF ORDER K WHOSE    
C     MTH COEFFICIENT IS STORED IN A(M+1). IT IS ASSUMED THAT THE       
C     INITIAL VALUES IN RT BRACKET THE FINAL VALUES.                    
C                                                                       
C                                                                       
C     IMPLICIT REAL*12 (A-H,O-Z)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION A(10),RT(10)                                            
C                                                                       
      DATA A0,TM21,A1S16,A1S4,A3S4 /0.0D0,1.0D-13,6.25D-2,0.25D0,0.75D0/
C                                                                       
      TOL=TM21                                                          
      K1=K+1                                                            
      R2=A0                                                             
      P2=A(1)                                                           
      DO 100 M=1,K                                                      
      R1=R2                                                             
      P1=P2                                                             
      R2=RT(M)                                                          
      P2=A(K1)                                                          
      DO 10 I=1,K                                                       
   10 P2=P2*R2+A(K1-I)                                                  
      PROD=P1*P2                                                        
      IF(PROD.LT.A0) GO TO 20                                           
      WRITE(6,15) M,K                                                   
   15 FORMAT(/12H0ROOT NUMBER,I4,                                       
     1   38H WAS NOT FOUND FOR POLYNOMIAL OF ORDER,I4//)                
      STOP                                                              
   20 R5=R1                                                             
      P5=P1                                                             
      R6=R2                                                             
      P6=P2                                                             
   30 R3=R5                                                             
      P3=P5                                                             
      R4=R6                                                             
      P4=P6                                                             
      R =(R3*P4-R4*P3)/(P4-P3)                                          
      DR=R4-R3                                                          
      DELTA=DR                                                          
      IF(DABS(DELTA).LT.TOL) GO TO 90                                   
      DR=A1S16*DR                                                       
      R5=R-DR                                                           
      IF(R5.LT.R3) R5=R3                                                
      R6=R+DR                                                           
      IF(R6.GT.R4) R6=R4                                                
      P5=A(K1)                                                          
      P6=P5                                                             
      DO 40 I=1,K                                                       
      P5=P5*R5+A(K1-I)                                                  
   40 P6=P6*R6+A(K1-I)                                                  
   45 PROD=P5*P6                                                        
      IF(PROD.LT.A0) GO TO 30                                           
      PROD=P3*P5                                                        
      IF(PROD.GT.A0) GO TO 60                                           
      R5=A1S4*R3+A3S4*R5                                                
      P5=A(K1)                                                          
      DO 50 I=1,K                                                       
   50 P5=P5*R5+A(K1-I)                                                  
      GO TO 45                                                          
   60 R6=A1S4*R4+A3S4*R6                                                
      P6=A(K1)                                                          
      DO 70 I=1,K                                                       
   70 P6=P6*R6+A(K1-I)                                                  
      GO TO 45                                                          
   90 RT(M)=R                                                           
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
