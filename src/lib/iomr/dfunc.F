      SUBROUTINE DFUNC(X,N)                                             
C                                                                       
C                                                                       
C                                                                       
C                                                                       
C                                                                       
C     IMPLICIT REAL*12 (A-H,O-Z)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      COMMON /FFM/    FF(19)                                            
C                                                                       
      DATA TM78, TM29, A1S2, PIE4, A1, A2, A80, A360 /1.0D-78, 1.0D-15, 
     1   0.5D0,0.7853981633974483096156608D0,1.0D0,2.0D0,80.0D0,360.0D0/
C                                                                       
      TOL=TM29                                                          
      XX=X+X                                                            
      FACMIN=XX                                                         
      E=TM78                                                            
      QX=-X                                                             
      IF(FACMIN.LT.A360) E=DEXP(QX)                                     
      IF (E.EQ.0.0D+00) E=TM78                                          
      IF(FACMIN.GT.A80) GO TO 100                                       
      TERM=A1                                                           
      SUM=A1                                                            
      FAC=N                                                             
      FAC=FAC+A1S2                                                      
   10 FAC=FAC+A1                                                        
      TERM=TERM*X/FAC                                                   
      SUM=SUM+TERM                                                      
      IF(FAC.LE.FACMIN) GO TO 10                                        
      T=TERM                                                            
      S=SUM                                                             
      IF(T.GT.S*TOL) GO TO 10                                           
      FAC=N+N+1                                                         
      FF(N+1)=SUM*E/FAC                                                 
      M=N-1                                                             
      FAC=M+M+1                                                         
   20 IF(M.LT.0) RETURN                                                 
      FF(M+1)=(E+XX*FF(M+2))/FAC                                        
      M=M-1                                                             
      FAC=FAC-A2                                                        
      GO TO 20                                                          
C                                                                       
C     USE ASYMPTOTIC EXPANSION FOR LARGE ARGUMENTS.                     
C                                                                       
  100 A=DSQRT(PIE4/X)                                                   
      TMAX=A*TOL/E                                                      
      TERM=A1/XX                                                        
      SUM=TERM                                                          
      FAC=A1                                                            
  110 FAC=FAC-A2                                                        
      TERM=FAC*TERM/XX                                                  
      SUM=TERM+SUM                                                      
      T=TERM                                                            
      IF(DABS(T).GT.TMAX) GO TO 110                                     
      FF(1)=A-E*SUM                                                     
      FAC=-A1                                                           
      M=0                                                               
  120 IF(M.EQ.N) RETURN                                                 
      M=M+1                                                             
      FAC=FAC+A2                                                        
      FF(M+1)=(FAC*FF(M)-E)/XX                                          
      GO TO 120                                                         
      END                                                               
