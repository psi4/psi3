      SUBROUTINE ATEBC(A,B,C,NI,NK,NJ)                                  
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      REAL*8 A(NJ,NI),B(NI,NK),C(NK,NJ)                                 
C                                                                       
      isize=nk*(ni+nj)
      if(isize.lt.80000) then
         call mmult(b,ni,0,c,nk,0,a,nj,1,ni,nk,nj,0)
      else
         CALL ZERO(A,NJ*NI)                                                
C                                                                       
         CALL MXMB(B,1,NI,C,1,NK,A,NJ,1,NI,NK,NJ)                          
      endif
      return
C                                                                       
      DO 3 J=1,NJ                                                       
         DO 2 K=1,NK                                                    
            CKJ=C(K,J)                                                  
            IF (CKJ.EQ.0.0D+00) GO TO 2                                 
            DO 1 I=1,NI                                                 
               A(J,I)=A(J,I)+B(I,K)*CKJ                                 
    1       CONTINUE                                                    
    2    CONTINUE                                                       
    3 CONTINUE                                                          
C                                                                       
      RETURN                                                            
      END                                                               
