      SUBROUTINE SREW(ITAPE)                                            
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      common /pointr/ wptr(128),tptr(128)
C                                                                       
C THIS ROUTINE REWINDS FILES USING THE FORTRAN UTILITIES DIRECT ACCESS  
C                                                                       
C PLACE THE POINTER OF ITAPE AT THE BEGINNING                           
C                                                                       
      wptr(itape)=1
C                                                                       
      RETURN                                                            
      END                                                               
