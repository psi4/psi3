*Deck Pak
*Deck UnPak
      SUBROUTINE PAK                                                    
C                                                                       
      implicit integer (a-z)
C                                                                       
      COMMON /INOUT/  IOUT,ICHECK,ITAP30                                
      COMMON /PACKQ/  I1,i2,i3,i4,I5                                       
C                                                                       
      i5=or(or(or(i1,lshift(i2,8)),lshift(i3,16)),lshift(i4,24))
C                                                                       
      RETURN                                                            
C                                                                       
C----UNPAK-----------------------------------------------UNPAK          
C                                                                       
      ENTRY UNPAK                                                       
C                                                                       
      i1=and(i5,255)
      i2=and(rshift(i5,8),255)
      i3=and(rshift(i5,16),255)
      i4=and(rshift(i5,24),255)
C                                                                       
      RETURN                                                            
      END                                                               
