      SUBROUTINE DSMIT(C,S,N)                                           
C                                                                       
C     RETURNS AN N*N TRIANGULAR MATRIX C SUCH THAT C(TRANSPOSE)SC=I,    
C     WHERE I IS A N*N IDENTITY MATRIX                                  
C                                                                       
C                                                                       
C     IMPLICIT REAL*12 (A-H,O-Z)                                        
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      DIMENSION C(10,10),S(10,10),V(10),Y(10)                           
C                                                                       
      DATA A0, A1 /0.0D0, 1.0D0/                                        
C                                                                       
      DO 10 I=1,N                                                       
      DO 10 J=1,I                                                       
   10 C(I,J)=A0                                                         
      DO 100 J=1,N                                                      
      KMAX=J-1                                                          
      FAC=S(J,J)                                                        
      IF(KMAX.EQ.0) GO TO 60                                            
      DO 20 K=1,KMAX                                                    
      V(K)=A0                                                           
   20 Y(K)=S(K,J)                                                       
      DO 50 K=1,KMAX                                                    
      DOT=A0                                                            
      DO 30 I=1,K                                                       
   30 DOT=C(I,K)*Y(I)+DOT                                               
      DO 40 I=1,K                                                       
   40 V(I)=V(I)-DOT*C(I,K)                                              
   50 FAC=FAC-DOT*DOT                                                   
   60 FAC=A1/DSQRT(FAC)                                                 
      C(J,J)=FAC                                                        
      IF(KMAX.EQ.0) GO TO 100                                           
      DO 70 K=1,KMAX                                                    
   70 C(K,J)=FAC*V(K)                                                   
  100 CONTINUE                                                          
      RETURN                                                            
      END                                                               
