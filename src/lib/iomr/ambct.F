      SUBROUTINE AMBCT(A,B,C,NI,NK,NJ)                                  
C                                                                       
C                                                                       
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      DIMENSION A(NI,NJ),B(NI,NK),C(NJ,NK)                              
C                                                                       
#if 1
      isize=nk*(ni+nj)
      if(isize.lt.80000) then
         call minmul(b,0,c,1,a,0,ni,nk,nj,1)
      else
         DO 3 I=1,NI                                                       
              DO 2 J=1,NJ                                                  
                   T=A(I,J)                                                
                   DO 1 K=1,NK                                             
                        T=T-B(I,K)*C(J,K)                                  
    1              CONTINUE                                                
                   A(I,J)=T                                                
    2         CONTINUE                                                     
    3    CONTINUE                                                          
      endif
#else
      call dgemm('n','t',ni,nj,nk,-1.0,b,ni,c,nj,1.0,a,ni)
#endif
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
