      SUBROUTINE RWRIT(ITAPE,ARRAY,NLEN,IREC)                           
C                                                                       
C THIS ROUTINE WRITES NLEN INTEGER WORDS FROM ARRAY TO FILE ITAPE       
C STARTING AT SECTOR IREC.                                              
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      COMMON /SECT/ SECTOR                                              
      COMMON /IOBUFF/ BUFF(1024)                                        
C                                                                       
      DIMENSION ARRAY(NLEN)                                             
C                                                                       
      IPOS = (IREC-1)*SECTOR + 1                                        
      call wwritw(itape,array,nlen,ipos,junk)
C                                                                       
      RETURN                                                            
      END                                                               
