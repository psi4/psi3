      SUBROUTINE ATPBCT(A,B,C,NI,NK,NJ)                                 
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      REAL*8 A(NJ,NI),B(NI,NK),C(NJ,NK)                                 
C                                                                       
      CALL MXMB(B,1,NI,C,NJ,1,A,NJ,1,NI,NK,NJ)                          
C                                                                       
cets111490      DO 3 K=1,NK                   
cets111490         DO 2 J=1,NJ                 
cets111490            CJK=C(J,K)                
cets111490            IF (CJK.EQ.0.0D+00) GO TO 2
cets111490            DO 1 I=1,NI                
cets111490               A(J,I)=A(J,I)+B(I,K)*CJK 
cets111490    1       CONTINUE                     
cets111490    2    CONTINUE                         
cets111490    3 CONTINUE                             
C                                                                       
      RETURN                                                            
      END                                                               
