      SUBROUTINE LOCATE(INPUT,TOKEN,IERROR)                             
C                                                                       
C     ----- SEARCH THROUGH INPUT FILE FOR TOKEN BEGINNING               
C           WITH # TO LOCATE INPUT FOR PROGRAM.  IERROR IS              
C           SET TO 0 IF NO ERRORS, 1 IF ANY ERROR OCCURS.               
C                                                                       
C                                                                       
      CHARACTER*10 TOKEN,LINE                                           
C                                                                       
      REWIND (UNIT=INPUT,ERR=99)                                        
C                                                                       
    1 CONTINUE                                                          
      READ (UNIT=INPUT,FMT='(A10)',END=99,ERR=99) LINE                  
clj   write(6,*) 'locate:token,line:',token,',',line
clj   write(*,*) 'locate:token,line:',token,',',line
      IF (LINE.NE.TOKEN) GO TO 1                                        
C                                                                       
      IERROR=0                                                          
      RETURN                                                            
C                                                                       
C                                                                       
   99 CONTINUE                                                          
      IERROR=1                                                          
      RETURN                                                            
      END                                                               
