      SUBROUTINE RGETSA(ITAPE,IADR)                                     
C                                                                       
C THIS ROUTINE RETURNS THE CURRENT SECTOR ADDRESS FOR FILES USING       
C THE SREAD/SWRIT AND/OR RREAD/RWRIT IO ROUTINES.                       
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      common /pointr/ wptr(128),tptr(128)
      COMMON /SECT/ SECTOR                                              
C                                                                       
      ipos=wptr(itape)
      IREC = (IPOS-1)/SECTOR + 1                                        
      TEST = SECTOR*(IREC-1) + 1                                        
      IF(IPOS.NE.TEST) THEN                                
      WRITE(*,*) ' ERROR ENCOUNTERED IN RGETSA FOR FILE ',ITAPE         
      WRITE(*,*) ' IPOS,TEST = ',IPOS,TEST                    
      STOP                                                              
      END IF                                                            
C                                                                       
      IADR=IREC                                                         
C                                                                       
      RETURN                                                            
      END                                                               
