      SUBROUTINE AMBC(A,B,C,NI,NK,NJ)                                   
C                                                                       
C                                                                       
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
C                                                                       
      DIMENSION A(NI,NJ),B(NI,NK),C(NK,NJ)                              
C                                                                       
      isize=nk*(ni+nj)

#if 1
      if(isize.lt.80000) then
         call minmul(b,0,c,0,a,0,ni,nk,nj,1)
      else
         DO 3 I=1,NI                                                       
              DO 2 J=1,NJ                                                  
                   T=A(I,J)                                                
                   DO 1 K=1,NK                                             
                        T=T-B(I,K)*C(K,J)                                  
    1              CONTINUE                                                
                   A(I,J)=T                                                
    2         CONTINUE                                                     
    3    CONTINUE                                                          
      endif
#else
      call dgemm('n','n',ni,nj,nk,-1.0,b,ni,c,nk,1.0,a,ni)
#endif
      RETURN                                                            
C     END                                                               
C                                                                       
C     DIMENSION A(1),B(1),C(1)                                          
C                                                                       
C     IJ=1                                                              
C     JA=1                                                              
C     DO 3 J=1,NJ                                                       
C     DO 2 I=1,NI                                                       
C     IK=I                                                              
C     KJ=JA                                                             
C     T=A(IJ)                                                           
C     DO 1 K=1,NK                                                       
C     T=T-B(IK)*C(KJ)                                                   
C     IK=IK+NI                                                          
C     KJ=KJ+1                                                           
C   1 CONTINUE                                                          
C     A(IJ)=T                                                           
C     IJ=IJ+1                                                           
C   2 CONTINUE                                                          
C     JA=JA+NK                                                          
C   3 CONTINUE                                                          
C     RETURN                                                            
      END                                                               
