      SUBROUTINE APBCTT(A,B,C,NI,NK,NJ)                                 
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      REAL*8 A(NJ,NI),B(NK,NI),C(NJ,NK)                                 
C                                                                       
CTJL  CALL MXMB(B,1,NI,C,NK,1,A,NI,1,NI,NK,NJ)                          
      call mxmb(b,nk,1,c,nj,1,a,nj,1,ni,nk,nj)
C                                                                       
cets111490      DO 3 K=1,NK                    
cets111490         DO 2 J=1,NJ                  
cets111490            CJK=C(J,K)                 
cets111490            IF (CJK.EQ.0.0D+00) GO TO 2 
cets111490            DO 1 I=1,NI                  
cets111490               A(J,I)=A(J,I)+B(K,I)*CJK   
cets111490    1       CONTINUE
cets111490    2    CONTINUE
cets111490    3 CONTINUE 
C                                                                       
      RETURN                                                            
      END                                                               
