      SUBROUTINE MATINV(A,B,PIVOT,INDEX,M,N,DETERM,NMAX)                
C     MATRIX INVERSION WITH ACCOMPANYING SOLUTION OF LINEAR EQUATIONS   
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      DIMENSION A(NMAX,NMAX),B(NMAX,NMAX),PIVOT(NMAX),INDEX(NMAX)       
      DATA ZERO,ONE / 0.0D+00 , 1.0D+00 /                               
C                                                                       
      DETERM=ONE                                                        
      DO 20 I=1,N                                                       
      PIVOT(I)=ZERO                                                     
   20 CONTINUE                                                          
C     PERFORM SUCCESSIVE PIVOT OPERATIONS (GRAND LOOP)                  
      DO 550 I=1,N                                                      
C     SEARCH FOR PIVOT ELEMENT AND EXTEND DETERMINANT PARTIAL PRODUCT   
      AMAX=ZERO                                                         
      DO 105 J=1,N                                                      
      IF (PIVOT(J).NE.ZERO) GO TO 105                                   
      DO 100 K=1,N                                                      
      IF (PIVOT(K).NE.ZERO) GO TO 100                                   
      AVAL=A(J,K)                                                       
      TEMP=DABS(AVAL)                                                   
      IF (TEMP.LT.AMAX) GO TO 100                                       
      IROW=J                                                            
      ICOLUM=K                                                          
      AMAX=TEMP                                                         
  100 CONTINUE                                                          
  105 CONTINUE                                                          
      INDEX(I)=4096*IROW+ICOLUM                                         
      J=IROW                                                            
      AMAX=A(J,ICOLUM)                                                  
      DETERM=AMAX*DETERM                                                
C     RETURN IF MATRIX IS SINGULAR (ZERO PIVOT) AFTER COLUMN INTERCHANGE
      IF(DETERM.EQ.ZERO) GO TO 600                                      
      PIVOT(ICOLUM)=AMAX                                                
C     INTERCHANGE ROWS TO PUT PIVOT ELEMENT ON DIAGONAL                 
      IF (IROW.EQ.ICOLUM) GO TO 260                                     
      DETERM=-DETERM                                                    
      DO 200 K=1,N                                                      
      SWAP=A(J,K)                                                       
      A(J,K)=A(ICOLUM,K)                                                
  200 A(ICOLUM,K)=SWAP                                                  
      IF (M.LE.0) GO TO 260                                             
      DO 250 K=1,M                                                      
      SWAP=B(J,K)                                                       
      B(J,K)=B(ICOLUM,K)                                                
  250 B(ICOLUM,K)=SWAP                                                  
C     DIVIDE PIVOT ROW BY PIVOT ELEMENT                                 
  260 K=ICOLUM                                                          
      A(ICOLUM,K)=ONE                                                   
      DO 350 K=1,N                                                      
  350 A(ICOLUM,K)=A(ICOLUM,K)/AMAX                                      
      IF (M.LE.0) GO TO 380                                             
      DO 370 K=1,M                                                      
  370 B(ICOLUM,K)=B(ICOLUM,K)/AMAX                                      
C     REDUCE NON-PIVOT ROWS                                             
  380 DO 550 J=1,N                                                      
      IF (J.EQ.ICOLUM) GO TO 550                                        
      T=A( J,ICOLUM)                                                    
      A( J,ICOLUM)=ZERO                                                 
      DO 450 K=1,N                                                      
  450 A( J,K)=A( J,K)-A(ICOLUM,K)*T                                     
      IF (M.LE.0) GO TO 550                                             
      DO 500 K=1,M                                                      
  500 B( J,K)=B( J,K)-B(ICOLUM,K)*T                                     
  550 CONTINUE                                                          
C     INTERCHANGE COLUMNS AFTER ALL PIVOT OPERATIONS HAVE BEEN PERFORMED
  600 DO 710 I=1,N                                                      
      I1=N+1-I                                                          
      K=INDEX(I1)/4096                                                  
      ICOLUM=INDEX(I1)-4096*K                                           
      IF (K.EQ.ICOLUM) GO TO 710                                        
      DO 705 J=1,N                                                      
      SWAP=A(J,K)                                                       
      A(J,K)=A(J,ICOLUM)                                                
  705 A(J,ICOLUM)=SWAP                                                  
  710 CONTINUE                                                          
      RETURN                                                            
      END                                                               
