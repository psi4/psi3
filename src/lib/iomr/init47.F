      SUBROUTINE INIT47                                                 
C                                                                       
C   THIS PROGRAM INTIALIZES USEAGE OF TAPE 47                           
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      COMMON/LOC47/LOCS(192)                                            
      COMMON/T47SEC/ITAP47                                              
C                                                                       
C   DEFINE CONSTANTS                                                    
      ITAP47=47                                                         
C                                                                       
C   FIND LOCATIONS OF PARAMETERS AND MATRICES IN THE 47 FILE            
      CALL RFILE(ITAP47)                                                
      CALL WREADW(ITAP47,LOCS,192,1,JUNK)                               
CC    WRITE (*,*) LOCS                                                  
C                                                                       
  201 CONTINUE                                                          
      RETURN                                                            
      END                                                               
