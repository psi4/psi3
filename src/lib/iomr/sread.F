      SUBROUTINE SREAD(ITAPE,ARRAY,NLEN)                                
C                                                                       
C THIS ROUTINE READS NLEN INTEGER WORDS FROM ITAPE INTO ARRAY           
C STARTING AT THE CURRENT POINTER LOCATION.                             
C                                                                       
      IMPLICIT INTEGER (A-Z)                                            
C                                                                       
      common /pointr/ wptr(128),tptr(128)
      COMMON /IOBUFF/ BUFF(1024)                                        
      COMMON /SECT/ SECTOR                                              
C                                                                       
      DIMENSION ARRAY(NLEN)                                             
C                                                                       
CTJL  WRITE(*,*) ' IN SREAD,  ITAPE,NLEN = ',ITAPE,NLEN                 
      start=wptr(itape)
      call wreadw(itape,array,nlen,start,end)
c
      wptr(itape)=(wptr(itape)-1+1023)/1024*1024+1
C                                                                       
      RETURN                                                            
      END                                                               
