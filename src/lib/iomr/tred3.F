C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      SUBROUTINE TRED3(N,NV,A,D,E,E2)                                   
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,K,L,N,II,IZ,JK,NV                                     
      REAL*8 A(NV),D(N),E(N),E2(N)                                      
      REAL*8 F,G,H,HH,SCALE                                             
      DATA ZERO /0.0D+00/                                               
C     REAL*8 SQRT,ABS,SIGN                                              
C                                                                       
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED3,    
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   
C                                                                       
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX, STORED AS        
C     A ONE-DIMENSIONAL ARRAY, TO A SYMMETRIC TRIDIAGONAL MATRIX        
C     USING ORTHOGONAL SIMILARITY TRANSFORMATIONS.                      
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        N IS THE ORDER OF THE MATRIX,                                  
C                                                                       
C        NV MUST BE SET TO THE DIMENSION OF THE ARRAY PARAMETER A       
C          AS DECLARED IN THE CALLING PROGRAM DIMENSION STATEMENT,      
C                                                                       
C        A CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC            
C          INPUT MATRIX, STORED ROW-WISE AS A ONE-DIMENSIONAL           
C          ARRAY, IN ITS FIRST N*(N+1)/2 POSITIONS.                     
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        A CONTAINS INFORMATION ABOUT THE ORTHOGONAL                    
C          TRANSFORMATIONS USED IN THE REDUCTION,                       
C                                                                       
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX,    
C                                                                       
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO,      
C                                                                       
C        E2 CONTAINS THE SQUARES OF THE CORRESPONDING ELEMENTS OF E.    
C          E2 MAY COINCIDE WITH E IF THE SQUARES ARE NOT NEEDED.        
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
C     ********** FOR I=N STEP -1 UNTIL 1 DO -- **********               
      DO  300 II = 1, N                                                 
           I = N + 1 - II                                               
           L = I - 1                                                    
           IZ = (I * L) / 2                                             
           H =ZERO                                                      
           SCALE = ZERO                                                 
           IF (L .LT. 1) GO TO 130                                      
C     ********** SCALE ROW (ALGOL TOL THEN NOT NEEDED) **********       
           DO 120 K = 1, L                                              
                IZ = IZ + 1                                             
                D(K) = A(IZ)                                            
                SCALE = SCALE + DABS(D(K))                              
  120      CONTINUE                                                     
C                                                                       
           IF (SCALE .NE.ZERO) GO TO 140                                
  130      E(I) =ZERO                                                   
           E2(I) = ZERO                                                 
           GO TO 290                                                    
C                                                                       
  140      DO 150 K = 1, L                                              
                D(K) = D(K) / SCALE                                     
                H = H + D(K) * D(K)                                     
  150      CONTINUE                                                     
C                                                                       
           E2(I) = SCALE * SCALE * H                                    
           F = D(L)                                                     
           G = -DSIGN(DSQRT(H),F)                                       
           E(I) = SCALE * G                                             
           H = H - F * G                                                
           D(L) = F - G                                                 
           A(IZ) = SCALE * D(L)                                         
           IF (L .EQ. 1) GO TO 290                                      
           F =ZERO                                                      
C                                                                       
           DO 240 J = 1, L                                              
                G =ZERO                                                 
                JK = (J * (J-1)) / 2                                    
C     ********** FORM ELEMENT OF A*U **********                         
                DO 180 K = 1, L                                         
                     JK = JK + 1                                        
                     IF (K .GT. J) JK = JK + K - 2                      
                     G = G + A(JK) * D(K)                               
  180           CONTINUE                                                
C     ********** FORM ELEMENT OF P **********                           
                E(J) = G / H                                            
                F = F + E(J) * D(J)                                     
  240      CONTINUE                                                     
C                                                                       
           HH = F / (H + H)                                             
           JK = 0                                                       
C     ********** FORM REDUCED A **********                              
           DO 260 J = 1, L                                              
                F = D(J)                                                
                G = E(J) - HH * F                                       
                E(J) = G                                                
C                                                                       
                DO 260 K = 1, J                                         
                     JK = JK + 1                                        
                     A(JK) = A(JK) - F * E(K) - G * D(K)                
  260      CONTINUE                                                     
C                                                                       
  290      D(I) = A(IZ+1)                                               
           A(IZ+1) = SCALE * DSQRT(H)                                   
  300 CONTINUE                                                          
C                                                                       
      RETURN                                                            
C     ********** LAST CARD OF TRED3 **********                          
      END                                                               
