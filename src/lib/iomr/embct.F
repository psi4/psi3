      SUBROUTINE EMBCT(A,B,C,NI,NK,NJ)                                  
C                                                                       
C                                                                       
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      DIMENSION A(NI,NJ),B(NI,NK),C(NJ,NK)                              
C                                                                       
      isize=nk*(ni+nj)
      if(isize.lt.80000) then
         call minmul(b,0,c,1,a,0,ni,nk,nj,0)
      else
      DO 3 I=1,NI                                                       
           DO 2 J=1,NJ                                                  
                T=0.0D+00                                               
                DO 1 K=1,NK                                             
                     T=T-B(I,K)*C(J,K)                                  
    1           CONTINUE                                                
                A(I,J)=T                                                
    2      CONTINUE                                                     
    3 CONTINUE                                                          
      endif
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
C     T=0.0D+00                                                         
C     DO 1 K=1,NK                                                       
C     T=T-B(IK)*C(JK)                                                   
C     IK=IK+NI                                                          
C     JK=JK+NJ                                                          
C   1 CONTINUE                                                          
C     A(IJ)=T                                                           
C     IJ=IJ+1                                                           
C   2 CONTINUE                                                          
C   3 CONTINUE                                                          
C     RETURN                                                            
      END                                                               
