C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,K,L,M,N,II,L1,NM,MML,IERR                             
      REAL*8 D(N),E(N),Z(NM,N)                                          
      REAL*8 B,C,F,G,H,P,R,S,MACHEP                                     
      DATA ZERO,ONE,TWO /0.0D+00,1.0D+00,2.0D+00/                       
C     REAL*8 SQRT,ABS,SIGN                                              
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     
C     WILKINSON.                                                        
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   
C                                                                       
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          
C          DIMENSION STATEMENT,                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX,          
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY,               
C                                                                       
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        
C          THE IDENTITY MATRIX.                                         
C                                                                       
C      ON OUTPUT-                                                       
C                                                                       
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          
C          UNORDERED FOR INDICES 1,2,...,IERR-1,                        
C                                                                       
C        E HAS BEEN DESTROYED,                                          
C                                                                       
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       
C          EIGENVALUES,                                                 
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
      MACHEP =TWO**(-47)                                                
C                                                                       
      IERR = 0                                                          
      IF (N .EQ. 1) GO TO 1001                                          
C                                                                       
      DO 100 I = 2, N                                                   
  100 E(I-1) = E(I)                                                     
C                                                                       
      F =ZERO                                                           
      B =ZERO                                                           
      E(N) = ZERO                                                       
C                                                                       
      DO 240 L = 1, N                                                   
           J = 0                                                        
           H = MACHEP * (DABS(D(L)) + DABS(E(L)))                       
           IF (B .LT. H) B = H                                          
C     ********** LOOK FOR SMALL SUB-DIAGONAL ELEMENT **********         
           DO 110 M = L, N                                              
                IF (DABS(E(M)) .LE. B) GO TO 120                        
C     ********** E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               
C                THROUGH THE BOTTOM OF THE LOOP **********              
  110      CONTINUE                                                     
C                                                                       
  120      IF (M .EQ. L) GO TO 220                                      
  130      IF (J .EQ. 30) GO TO 1000                                    
           J = J + 1                                                    
C     ********** FORM SHIFT **********                                  
           L1 = L + 1                                                   
           G = D(L)                                                     
           P = (D(L1) - G) / (TWO * E(L))                               
           R = DSQRT(P*P+ONE)                                           
           D(L) = E(L) / (P + DSIGN(R,P))                               
           H = G - D(L)                                                 
C                                                                       
           DO 140 I = L1, N                                             
  140      D(I) = D(I) - H                                              
C                                                                       
           F = F + H                                                    
C     ********** QL TRANSFORMATION **********                           
           P = D(M)                                                     
           C = ONE                                                      
           S = ZERO                                                     
           MML = M - L                                                  
C     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********             
           DO 200 II = 1, MML                                           
                I = M - II                                              
                G = C * E(I)                                            
                H = C * P                                               
                IF (DABS(P) .LT. DABS(E(I))) GO TO 150                  
                C = E(I) / P                                            
                R = DSQRT(C*C+ONE)                                      
                E(I+1) = S * P * R                                      
                S = C / R                                               
                C = ONE / R                                             
                GO TO 160                                               
  150           C = P / E(I)                                            
                R = DSQRT(C*C+ONE)                                      
                E(I+1) = S * E(I) * R                                   
                S = ONE / R                                             
                C = C * S                                               
  160           P = C * D(I) - S * G                                    
                D(I+1) = H + S * (C * G + S * D(I))                     
C     ********** FORM VECTOR **********                                 
                DO 180 K = 1, N                                         
                     H = Z(K,I+1)                                       
                     Z(K,I+1) = S * Z(K,I) + C * H                      
                     Z(K,I) = C * Z(K,I) - S * H                        
  180           CONTINUE                                                
C                                                                       
  200      CONTINUE                                                     
C                                                                       
           E(L) = S * P                                                 
           D(L) = C * P                                                 
           IF (DABS(E(L)) .GT. B) GO TO 130                             
  220      D(L) = D(L) + F                                              
  240 CONTINUE                                                          
C     ********** ORDER EIGENVALUES AND EIGENVECTORS **********          
      DO 300 II = 2, N                                                  
           I = II - 1                                                   
           K = I                                                        
           P = D(I)                                                     
C                                                                       
           DO 260 J = II, N                                             
                IF (D(J) .GE. P) GO TO 260                              
                K = J                                                   
                P = D(J)                                                
  260      CONTINUE                                                     
C                                                                       
           IF (K .EQ. I) GO TO 300                                      
           D(K) = D(I)                                                  
           D(I) = P                                                     
C                                                                       
           DO 280 J = 1, N                                              
                P = Z(J,I)                                              
                Z(J,I) = Z(J,K)                                         
                Z(J,K) = P                                              
  280      CONTINUE                                                     
C                                                                       
  300 CONTINUE                                                          
C                                                                       
      GO TO 1001                                                        
C     ********** SET ERROR -- NO CONVERGENCE TO AN                      
C                EIGENVALUE AFTER 30 ITERATIONS **********              
 1000 IERR = L                                                          
 1001 RETURN                                                            
C     ********** LAST CARD OF TQL2 **********                           
      END                                                               
