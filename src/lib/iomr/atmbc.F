      SUBROUTINE ATMBC(A,B,C,NI,NK,NJ)                                  
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      REAL*8 A(NJ,NI),B(NI,NK),C(NK,NJ)                                 
C                                                                       
      isize=nk*(ni+nj)
      if(isize.lt.80000) then
         call minmul(b,0,c,0,a,1,ni,nk,nj,1)
      else
      DO 3 J=1,NJ                                                       
         DO 2 K=1,NK                                                    
            CKJ=C(K,J)                                                  
            IF (CKJ.EQ.0.0D+00) GO TO 2                                 
            DO 1 I=1,NI                                                 
               A(J,I)=A(J,I)-B(I,K)*CKJ                                 
    1       CONTINUE                                                    
    2    CONTINUE                                                       
    3 CONTINUE                                                          
      endif
C                                                                       
      RETURN                                                            
      END                                                               
