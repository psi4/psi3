      SUBROUTINE ZFILE(ITAPE,ISECT)                                     
C                                                                       
C THIS ROUTINE ZEROES OUT A FILE ON DISK                                
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
      COMMON /IOBUFF/ BUFF(1024)                                        
CTJL  INTEGER BUF(1024)                                                 
C         IN BLOCK DATA IOBUF IS 1024 INTEGER WORD ZEROES               
CTJL  CALL IZERO(BUF,1024)                                              
C     NSEC = ISECT * 2                                                  
C     IF(2*NSEC.NE.ISECT) NSEC = NSEC + 1                               
      DO 1  I=1,ISECT                                                   
   1  CALL SWRIT(ITAPE,BUFF,1024)                                       
C                                                                       
      RETURN                                                            
      END                                                               
