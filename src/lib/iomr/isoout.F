*Deck IsoOut
*Deck IsoIn
      SUBROUTINE ISOOUT                                                 
C                                                                       
C                                                                       
C                                                                       
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               
C                                                                       
      COMMON /INOUT/  IOUT,ICHECK,ITAP30                                
      COMMON /SYMTRY/ NT                                                
      COMMON /ISOPAC/ INDIN(120),INDOUT(30)                            
      integer rshift
C                                                                       
C***  INTEGER*1 INDOUT                                                  
C     INTEGER*2 INDOUT                                                  
cps      LOGICAL*4 TRANS(120)                                              
cps      LOGICAL*1 INDOUT                                                  
cps      EQUIVALENCE (INDIN(1),TRANS(1))                                   
C                                                                       
      iword=0
      do 20 j=1,nt
         if (mod(j-1,4).eq.0) then
            iword=iword+1
            indin(j)=and(indout(iword),255)
         else
            indin(j)=and(rshift(indout(iword),8*mod(j-1,4)),255)
         end if
   20 continue
c
      RETURN                                                            
C                                                                       
C--ISOIN-------------------------------------------------ISOIN          
C                                                                       
      ENTRY ISOIN                                                       
C                                                                       
      iword=0
      DO 40 J=1,NT                                                      
         if (mod(j-1,4).eq.0) then
            iword=iword+1
            indout(iword)=indin(j)
         else
            indout(iword)=or(indout(iword),lshift(indin(j),
     #                       8*mod(j-1,4)))
         end if
   40 CONTINUE                                                          
      RETURN                                                            
      END                                                               
