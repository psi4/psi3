      SUBROUTINE INITMF(ISTART)                                         
C   THIS PROGRAM INTIALIZES USEAGE OF THE MASTER FILE                   
C**********************************************************             
C***LAST UPDATED ON FEBRUARY 06, 1985 BY YUKIO YAMAGUCHI***             
C**********************************************************             
C                                                                       
      IMPLICIT REAL*8 (A-H,O-Z)                                         
      COMMON/CALIF/LPARA(1024),APARA(1024)                              
      COMMON/LIMIT/MAXBAS,MAXBUF,MAXBF2,MAXBF4                          
      COMMON/LOCAT/LOCS(1024)                                           
      COMMON/MFSEC/MFILE,NSECT                                          
C                                                                       
C   DEFINE CONSTANTS                                                    
      MFILE=40                                                          
      NSECT=1024                                                        
      MAXBAS=52                                                         
      MAXBUF=4096                                                       
      MAXBF2=MAXBUF*2                                                   
      MAXBF4=MAXBUF*4                                                   
      IF(ISTART.EQ.0) GO TO 201                                         
C                                                                       
C   FIND LOCATIONS OF PARAMETERS AND MATRICES IN THE MASTER FILE        
      CALL RFILE(MFILE)                                                 
      CALL RREAD(MFILE,LOCS,NSECT,1)                                    
      CALL MREAD(LPARA,2)                                               
      CALL MREAD(APARA,3)                                               
C                                                                       
  201 CONTINUE                                                          
      RETURN                                                            
      END                                                               
