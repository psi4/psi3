      SUBROUTINE RSETSA(ITAPE,IADR)                                     
C                                                                       
C THIS ROUTINE SETS THE POINTER OF A UTILITIES FILE TO THE DESIRED      
C SECTOR ADDRESS.                                                       
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      common /pointr/ wptr(128),tptr(128)
      COMMON /SECT/ SECTOR                                              
C                                                                       
      IPOS = SEC2I(IADR-1) + 1                                          
      wptr(itape)=ipos
C                                                                       
      RETURN                                                            
      END                                                               
