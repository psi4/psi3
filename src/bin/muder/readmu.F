      subroutine READMU(MU,nblock,natom)                   
c/////////////////////////////////////////////////////////////////////  
c                                                                   //  
c  A subroutine to read in the dipole moments from file11           //  
c                                                                   //  
c/////////////////////////////////////////////////////////////////////  
c                  Variable Definitions                             //  
c                                                                   //  
c See subroutine MUDER                                              //  
c                                                                   //  
c///////////////////////////////////////////////////////////////////// 
      implicit real*8(a-h,o-z)                                          
      integer nblock                                       
      integer unit11,unit6,unit5,unit12,unit28,unit30,intco            
      integer unit18
      real*8 mu(3,nblock)                                               
#include <error.h>
      integer errcod, frdc, frdi, fcnt, frdia1
      integer prdup, prcntr, ip
      character*8 opentp, cityp
      logical oldinp
      common/croag/ oldinp
      common/unit/ intco,unit5,unit6,unit11,unit12,unit18,unit28
Ceeb
      common/type/ cityp
c
      ip = 0
      ip = prcntr('IS_ON DEBUG')
      if (ip.eq.1) then
        write(unit6,*)
        if (oldinp) then
          write(unit6,*) ' oldinp = true '
        else
          write(unit6,*) ' oldinp = false '
        endif
      endif
      if (oldinp) then
	write(unit6,*) ' If old style input is used, '
	write(unit6,*) ' this code can only handle cisd cases. '
	cityp = 'CI'
      else
	errcod = frdc('WFN',cityp)
	if(errcod.ne.EOK) then
          write(iout,*) ' no WFN input found in muder'
          call fioflu(unit6)
          call qabort
        endif
      endif

Ceeb
      if(cityp.eq.'SCF     ') then
        call ffile(unit11,'file11',1)
        rewind unit11  
        do 19 j = 1,nblock
           read (unit11,*)
           mu(1,j) = 0.0D+00
           mu(2,j) = 0.0D+00
           mu(3,j) = 0.0D+00
           if(j.ne.nblock) then
              do 21 k = 1,2*natom
                 read (unit11,*)
 21           continue
           end if
 19     continue
Ceeb
      else if(cityp.eq.'CI') then
        call ffile(unit11,'file11',1)
        rewind unit11                                                     
        do 10 j = 1,nblock                                                
           read (unit11,*)                                                
           read (unit11,1) (mu(i,j), i = 1,3)                             
 1         format(25x,f15.10,5x,f15.10,5x,f15.10)                         
           if(j.ne.nblock) then                                           
              do 20 k = 1,2*natom                                         
                 read (unit11,*)                                          
 20           continue                                                    
           end if                                                         
 10     continue                                                          
      else if (cityp.eq.'CCSD'.or.cityp.eq.'CCSD_T') then
        call ffile(unit28,'file28',1)
	rewind unit28
	do 15 j=1,nblock
	  read (unit28,2) (mu(i,j), i=1,3)
2	  format(f20.10,f20.10,f20.10)
15	continue
      else
	write (unit6,*) ' WFN is not SCF, CI, CCSD, or CCSD_T. '
	write (unit6,*) ' Sorry, but this code just can not handle it.'
        call fioflu(unit6)
        call qabort
      endif
      return                                                            
      end                                                               
c                                                                       
