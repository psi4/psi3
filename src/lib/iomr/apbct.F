      SUBROUTINE APBCT(A,B,C,NI,NK,NJ)                                  
C                                                                       
C                                                                       
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      DIMENSION A(NI,NJ),B(NI,NK),C(NJ,NK)                              
C                                                                       
      CALL MXMB(B,1,NI,C,NJ,1,A,1,NI,NI,NK,NJ)                          
C                                                                       
cets111490      DO 3 I=1,NI                  
cets111490           DO 2 J=1,NJ              
cets111490                T=A(I,J)             
cets111490                DO 1 K=1,NK           
cets111490                     T=T+B(I,K)*C(J,K) 
cets111490    1           CONTINUE 
cets111490                A(I,J)=T  
cets111490    2      CONTINUE  
cets111490    3 CONTINUE  
      RETURN                                                            
C     END                                                               
C                                                                       
C     DIMENSION A(1),B(1),C(1)                                          
C                                                                       
C     IJ=1                                                              
C     DO 3 J=1,NJ                                                       
C     DO 2 I=1,NI                                                       
C     IK=I                                                              
C     JK=J                                                              
C     T=A(IJ)                                                           
C     DO 1 K=1,NK                                                       
C     T=T+B(IK)*C(JK)                                                   
C     IK=IK+NI                                                          
C     JK=JK+NJ                                                          
C   1 CONTINUE                                                          
C     A(IJ)=T                                                           
C     IJ=IJ+1                                                           
C   2 CONTINUE                                                          
C   3 CONTINUE                                                          
C     RETURN                                                            
      END                                                               
