C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE TQLRAT(N,D,E2,IERR)                                    
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,L,M,N,II,L1,MML,IERR                                  
      REAL*8 D(N),E2(N)                                                 
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                     
      DATA ZERO,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/                       
C     REAL*8 SQRT,ABS,SIGN                                              
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQLRAT,   
C     ALGORITHM 464, COMM. ACM 16, 689(1973) BY REINSCH.                
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES OF A SYMMETRIC              
C     TRIDIAGONAL MATRIX BY THE RATIONAL QL METHOD.                     
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E2 CONTAINS THE SQUARES OF THE SUBDIAGONAL ELEMENTS OF THE     
C          INPUT MATRIX IN ITS LAST N-1 POSITIONS.  E2(1) IS ARBITRARY. 
C                                                                       
C      ON OUTPUT-                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT AND          
C          ORDERED FOR INDICES 1,2,...IERR-1, BUT MAY NOT BE            
C          THE SMALLEST EIGENVALUES,                                    
C                                                                       
C        E2 HAS BEEN DESTROYED,                                         
C                                                                       
C        IERR IS SET TO                                                 
C          ZERO       FOR NORMAL RETURN,                                
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               
C                     DETERMINED AFTER 30 ITERATIONS.                   
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** MACHEP IS A MACHINE DEPENDENT PARAMETER SPECIFYING     
C                THE RELATIVE PRECISION OF FLOATING POINT ARITHMETIC.   
C                                                                       
C                **********                                             
      MACHEP = TWO**(-47)                                               
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      DO 100 I = 2, N                                                   
  100 E2(I-1) = E2(I)                                                   
C                                                                       
      F =ZERO                                                           
      B =ZERO                                                           
      E2(N) =ZERO                                                       
C                                                                       
      DO 290 L = 1, N                                                   
           J = 0                                                        
           H = MACHEP * (DABS(D(L)) + DSQRT(E2(L)))                     
           IF (B .GT. H) GO TO 105                                      
           B = H                                                        
           C = B * B                                                    
C     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT ********** 
  105      DO 110 M = L, N                                              
                IF (E2(M) .LE. C) GO TO 120                             
C     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT              
C                THROUGH THE BOTTOM OF THE LOOP **********              
  110      CONTINUE                                                     
C                                                                       
  120      IF (M .EQ. L) GO TO 210                                      
  130      IF (J .EQ. 30) GO TO 1000                                    
           J = J + 1                                                    
C     ********** FORM SHIFT **********                                  
           L1 = L + 1                                                   
           S = DSQRT(E2(L))                                             
           G = D(L)                                                     
           P = (D(L1) - G) / (TWO * S)                                  
           R = DSQRT(P*P+ONE)                                           
           D(L) = S / (P + DSIGN(R,P))                                  
           H = G - D(L)                                                 
C                                                                       
           DO 140 I = L1, N                                             
  140      D(I) = D(I) - H                                              
C                                                                       
           F = F + H                                                    
C     ********** RATIONAL QL TRANSFORMATION **********                  
           G = D(M)                                                     
           IF (G .EQ.ZERO) G = B                                        
           H = G                                                        
           S =ZERO                                                      
           MML = M - L                                                  
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
           DO 200 II = 1, MML                                           
                I = M - II                                              
                P = G * H                                               
                R = P + E2(I)                                           
                E2(I+1) = S * R                                         
                S = E2(I) / R                                           
                D(I+1) = H + S * (H + D(I))                             
                G = D(I) - E2(I) / G                                    
                IF (G .EQ.ZERO) G = B                                   
                H = G * P / R                                           
  200      CONTINUE                                                     
C                                                                       
           E2(L) = S * G                                                
           D(L) = H                                                     
C     ********** GUARD AGAINST UNDERFLOW IN CONVERGENCE TEST ********** 
           IF (H .EQ.ZERO) GO TO 210                                    
           IF (DABS(E2(L)) .LE. DABS(C/H)) GO TO 210                    
           E2(L) = H * E2(L)                                            
           IF (E2(L) .NE.ZERO) GO TO 130                                
  210      P = D(L) + F                                                 
C     ********** ORDER EIGENVALUES **********                           
           IF (L .EQ. 1) GO TO 250                                      
C     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********               
           DO 230 II = 2, L                                             
                I = L + 2 - II                                          
                IF (P .GE. D(I-1)) GO TO 270                            
                D(I) = D(I-1)                                           
  230      CONTINUE                                                     
C                                                                       
  250      I = 1                                                        
  270      D(I) = P                                                     
  290 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
 1000 IERR = L                                                          
 1001 RETURN                                                            
C     ********** LAST CARD OF TQLRAT **********                         
      END                                                               
