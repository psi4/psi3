      SUBROUTINE NCDLBL(ITAPE,NA,NB,NC,ND,NE,RF,RG,RH)                  
      REAL*8 RF,RG,RH,RN                                                
      DIMENSION NA(26),NB(26),N(112),RN(10)                             
      EQUIVALENCE(N(57),RN(1))                                          
      L = 26                                                            
      DO 10 I=1,26                                                      
           L = L + 1                                                    
           N(I) = NA(I)                                                 
 10   N(L) = NB(I)                                                      
      N(53) = NC                                                        
      N(54) = ND                                                        
      N(55) = NE                                                        
      RN(1) = RF                                                        
      RN(2) = RG                                                        
      RN(3) = RH                                                        
      CALL SWRIT(ITAPE,N,112)                                           
      RETURN                                                            
      end
C                                                                       
cets081591      ENTRY DCDLBL(ITAPE,NA,NB,NC,ND,NE,RF,RG,RH)
      subroutine DCDLBL(ITAPE,NA,NB,NC,ND,NE,RF,RG,RH)                       
      REAL*8 RF,RG,RH,RN                                                
      DIMENSION NA(26),NB(26),N(112),RN(10)                             
      EQUIVALENCE(N(57),RN(1))                                          

      CALL SREAD(ITAPE,N,112)                                           
      L = 26                                                            
      DO 20 I=1,26                                                      
           L = L + 1                                                    
           NA(I) = N(I)                                                 
 20   NB(I) = N(L)                                                      
      NC = N(53)                                                        
      ND = N(54)                                                        
      NE = N(55)                                                        
      RF = RN(1)                                                        
      RG = RN(2)                                                        
      RH = RN(3)                                                        
      RETURN                                                            
      END                                                               
