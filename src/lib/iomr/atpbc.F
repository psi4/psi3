      SUBROUTINE ATPBC(A,B,C,NI,NK,NJ)                                  
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      REAL*8 A(NJ,NI),B(NI,NK),C(NK,NJ)                                 
C                                                                       
      CALL MXMB(B,1,NI,C,1,NK,A,NJ,1,NI,NK,NJ)                          
C                                             
cets111490      DO 3 J=1,NJ                    
cets111490         DO 2 K=1,NK                  
cets111490            CKJ=C(K,J)                 
cets111490            IF (CKJ.EQ.0.0D+00) GO TO 2 
cets111490            DO 1 I=1,NI                
cets111490               A(J,I)=A(J,I)+B(I,K)*CKJ 
cets111490    1       CONTINUE  
cets111490    2    CONTINUE 
cets111490    3 CONTINUE 
C                                                                       
      RETURN                                                            
      END                                                               
