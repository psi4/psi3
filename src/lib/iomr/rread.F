      SUBROUTINE RREAD(ITAPE,ARRAY,NLEN,IREC)                           
C                                                                       
C THIS ROUTINE READS NLEN INTEGER WORDS INTO ARRAY STARTING AT          
C SECTOR IREC.                                                          
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      COMMON /SECT/ SECTOR                                              
      COMMON /IOBUFF/ BUFF(1024)                                        
C                                                                       
      DIMENSION ARRAY(NLEN)                                             
C                                                                       
      IPOS = (IREC-1)*SECTOR + 1                                        
      call wreadw(itape,array,nlen,ipos,junk)
C                                                                       
      RETURN                                                            
      END                                                               
