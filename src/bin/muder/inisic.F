      subroutine INISIC(SIC,ns,npar)                        
c/////////////////////////////////////////////////////////////////////  
c                                                                   //  
c  A subroutine to initialize the SIC array to 1, then read in      //  
c  the non-symmetric SIC's and set their value in SIC to 0.         //  
c                                                                   //  
c/////////////////////////////////////////////////////////////////////  
c                  Variable Definitions                             //  
c                                                                   //  
c ns = number of symmetrized internal coordinates                   //  
c SIC(ns) = an integer array indicating if the SIC ns is totally    //  
c           symmetric (1) or not (0)                                //  
c                                                                   //  
c/////////////////////////////////////////////////////////////////////  
      implicit real*8(a-h,o-z)                                          
      integer ns, SIC(ns),unit5,npar,unit6                              
      integer unit18,unit28,unit11,unit12,intco
#include <error.h>
      integer errcod, frdc, frdi, fcnt, frdia1, frdi2
      integer prdup, prcntr, ip
      character*8 opentp
      logical oldinp
      common/croag/ oldinp
      common/unit/ intco,unit5,unit6,unit11,unit12,unit18,unit28
c
      do 10 i = 1,ns                                                    
        SIC(i) = 1                                                      
 10   continue                                                          
c                                                                       
c.......read in the SIC's which are not totally symmetric and set       
c.......their value in SIC to 0                                         
c                                                                       
      if (oldinp) then
        do 15 i = 1,npar                                                  
           read (unit5,2) nsic                                            
 2         format(i5)                                                     
           SIC(nsic) = 0                                                  
           read (unit5,*)                                                 
 15     continue                                                          
      else
	do 600 i=1,npar
	  errcod = frdi2('DIP_PAR',nsic,i,1)
	  SIC(nsic) = 0
600	continue
      endif
      ip=0
      ip = prcntr('IS_ON VERBOSE')
      if (ip.eq.1) then
        write (6,*) 
	write (6,*) ' The symmetries of each IntCoord. '
        do 33 i = 1,ns                                                    
          write (6,*) sic(i)                                                
 33     continue                                                          
      endif
      return                                                            
      end                                                               
c                                                                       
