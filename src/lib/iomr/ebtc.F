      SUBROUTINE EBTC(A,B,C,NI,NK,NJ)                                   
C                                                                       
C                                                                       
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      DIMENSION A(NI,NJ),B(NK,NI),C(NK,NJ)                              
C                                                                       
      isize=nk*(ni+nj)
      if(isize.lt.80000) then
         call mmult(b,nk,1,c,nk,0,a,ni,0,ni,nk,nj,0)
      else
         call zero(a,ni*nj)
         call mxmb(b,nk,1,c,1,nk,a,1,ni,ni,nk,nj)
      endif
      return
cps      CALL ZERO(A,NI*NJ)                                                
cps      CALL MXMB(B,NK,1,C,1,NK,A,1,NI,NI,NK,NJ)                          
C                                                                       
      DO 3 I=1,NI                                                       
           DO 2 J=1,NJ                                                  
                T=0.0D+00                                               
                DO 1 K=1,NK                                             
                     T=T+B(K,I)*C(K,J)                                  
    1           CONTINUE                                                
                A(I,J)=T                                                
    2      CONTINUE                                                     
    3 CONTINUE                                                          
      RETURN                                                            
C     END                                                               
C                                                                       
C     DIMENSION A(1),B(1),C(1)                                          
C                                                                       
C     IJ=1                                                              
C     JA=1                                                              
C     DO 3 J=1,NJ                                                       
C     KI=1                                                              
C     DO 2 I=1,NI                                                       
C     KJ=JA                                                             
C     T=0.0D+00                                                         
C     DO 1 K=1,NK                                                       
C     T=T+B(KI)*C(KJ)                                                   
C     KI=KI+1                                                           
C     KJ=KJ+1                                                           
C   1 CONTINUE                                                          
C     A(IJ)=T                                                           
C     IJ=IJ+1                                                           
C   2 CONTINUE                                                          
C     JA=JA+NK                                                          
C   3 CONTINUE                                                          
C     RETURN                                                            
      END                                                               
