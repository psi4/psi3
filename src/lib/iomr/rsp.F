      SUBROUTINE RSP(NM,N,NV,A,W,MATZ,Z,FV1,FV2)                        
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
C                                                                       
      INTEGER I,J,N,NM,NV,IERR,MATZ                                     
      REAL*8 A(NV),W(N),Z(NM,N),FV1(N),FV2(N)                           
      DATA ZERO,ONE / 0.0D+00,1.0D+00 /                                 
    1 FORMAT(//,2X,' IERR = ',I5//)                                     
C                                                                       
C     THIS SUBROUTINE CALLS THE RECOMMENDED SEQUENCE OF                 
C     SUBROUTINES FROM THE EIGENSYSTEM SUBROUTINE PACKAGE (EISPACK)     
C     TO FIND THE EIGENVALUES AND EIGENVECTORS (IF DESIRED)             
C     OF A REAL SYMMETRIC PACKED MATRIX.                                
C                                                                       
C     ON INPUT-                                                         
C                                                                       
C        NM  MUST BE SET TO THE ROW DIMENSION OF THE TWO-DIMENSIONAL    
C        ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM            
C        DIMENSION STATEMENT,                                           
C                                                                       
C        N  IS THE ORDER OF THE MATRIX  A,                              
C                                                                       
C        NV  IS AN INTEGER VARIABLE SET EQUAL TO THE                    
C        DIMENSION OF THE ARRAY  A  AS SPECIFIED FOR                    
C        A  IN THE CALLING PROGRAM.  NV  MUST NOT BE                    
C        LESS THAN  N*(N+1)/2,                                          
C                                                                       
C        A  CONTAINS THE LOWER TRIANGLE OF THE REAL SYMMETRIC           
C        PACKED MATRIX STORED ROW-WISE,                                 
C                                                                       
C        MATZ  IS AN INTEGER VARIABLE SET EQUAL TO ZERO IF              
C        ONLY EIGENVALUES ARE DESIRED,  OTHERWISE IT IS SET TO          
C        ANY NON-ZERO INTEGER FOR BOTH EIGENVALUES AND EIGENVECTORS.    
C                                                                       
C     ON OUTPUT-                                                        
C                                                                       
C        W  CONTAINS THE EIGENVALUES IN ASCENDING ORDER,                
C                                                                       
C        Z  CONTAINS THE EIGENVECTORS IF MATZ IS NOT ZERO,              
C                                                                       
C        IERR  IS AN INTEGER OUTPUT VARIABLE SET EQUAL TO AN            
C        ERROR COMPLETION CODE DESCRIBED IN SECTION 2B OF THE           
C        DOCUMENTATION.  THE NORMAL COMPLETION CODE IS ZERO,            
C                                                                       
C        FV1  AND  FV2  ARE TEMPORARY STORAGE ARRAYS.                   
C                                                                       
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         
C                                                                       
C     ------------------------------------------------------------------
C                                                                       
      IF (N .LE. NM) GO TO 5                                            
      IERR = 10 * N                                                     
      GO TO 50                                                          
    5 IF (NV .GE. (N * (N + 1)) / 2) GO TO 10                           
      IERR = 20 * N                                                     
      GO TO 50                                                          
C                                                                       
   10 CALL  TRED3(N,NV,A,W,FV1,FV2)                                     
      IF (MATZ .NE. 0) GO TO 20                                         
C     ********** FIND EIGENVALUES ONLY **********                       
      CALL  TQLRAT(N,W,FV2,IERR)                                        
      GO TO 50                                                          
C     ********** FIND BOTH EIGENVALUES AND EIGENVECTORS **********      
   20 DO 40 I = 1, N                                                    
C                                                                       
           DO 30 J = 1, N                                               
                Z(J,I) =ZERO                                            
   30      CONTINUE                                                     
C                                                                       
           Z(I,I) =ONE                                                  
   40 CONTINUE                                                          
C                                                                       
      CALL  TQL2(NM,N,W,FV1,Z,IERR)                                     
      IF (IERR .NE. 0) GO TO 50                                         
      CALL  TRBAK3(NM,N,NV,A,N,Z)                                       
   50 IF(IERR.NE.0) GO TO 60                                            
      RETURN                                                            
   60 WRITE(6,1) IERR                                                   
      RETURN                                                            
C     ********** LAST CARD OF RSP **********                            
      END                                                               
