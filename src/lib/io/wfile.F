      SUBROUTINE WFILE(ITAPE,ISECT)                                     
C                                                                       
C THIS ROUTINE ZEROES OUT A FILE ON DISK                                
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
      INTEGER BUF(224)                                                  
C                                                                       
      CALL IZERO(BUF,224)                                               
      NSEC = ISECT / 2                                                  
      IF(2*NSEC.NE.ISECT) NSEC = NSEC + 1                               
      DO 1  I=1,NSEC                                                    
   1  CALL SWRIT(ITAPE,BUF,224)                                         
C                                                                       
      RETURN                                                            
      END                                                               
