      subroutine READR(RVAL,PAIRS,dprty,      
     & NCOOR,DISP,MU,SIC,REFR,MUDER,COUNT,ns,nblock,icharg,natom,       
     & DER,ENER,npar)                                              
c/////////////////////////////////////////////////////////////////////  
c                                                                   //  
c  A subroutine to read the values of the internal coordinates      //  
c  and calculate the displacement relative to the reference geometry//  
c                                                                   //  
c/////////////////////////////////////////////////////////////////////  
c                  Variable Definitions                             //  
c                                                                   //  
c See subroutine MUDER                                              //  
c                                                                   //  
croag  dprty(j,i) = parity of ith (x,y,z) dipole moment component   //
croag            for jth positive displacement under the symmetry   //
croag            operation which generates the neg. displacemnt     //
c                                                                   //
c/////////////////////////////////////////////////////////////////////  
      implicit real*8(a-h,o-z)                                          
      integer ns, unit12, check, unit6, iblock, nblock, unit18,icharg,nc
      integer unit5 , dprty(npar,3), intco
      integer unit28, unit11
      integer NCOOR(nblock),COUNT(ns),PAIRS(ns,2),SIC(ns)               
      real*8 MU(3,nblock),MUDER(3,ns),DER(nblock,ns),ENER(nblock)       
      dimension REFR(ns),RVAL(nblock,ns),DISP(nblock,ns)                
#include <error.h>
      integer errcod, frdc,frdi,fcnt,frdia1,frdi2,frdlpn,frdint
      integer prdup, prcntr, ip
      integer un12a
      character*8 opentp
      logical oldinp
      common/croag/ oldinp
      common/unit/ intco,unit5,unit6,unit11,unit12,unit18,unit28
      common/un12a/ un12a
c                                                                       
c.......Read in the reference values for the SIC's assuming they        
c.......appear first in file12                                          
c                                                                       
      rewind unit12                                                     
      rewind unit18                                                     
      iblock = 1                                                        
      ncoor(iblock) = 0                                                 
      DISP(1,0) = 0.0d+00                                               
      read (unit12,11) ener(1)                                          
 11   format(17x,f20.10)                                                
      write (unit6,*)                                                       
      write (unit6,*) ' Energy(ref) = ',ener(1)                              
      ip = 0
      ip = prcntr('IS_ON VERBOSE')
      if (ip.eq.1) then
        write (6,*)                                                       
        write (6,*) ' Reference Value of the IntCoords '                   
        write (6,*)                                                       
      endif
      do 10 i = 1,ns                                                    
	 RVAL(iblock,i) = 0.0d+00
	 DER(iblock,i) = 0.0d+00
         read (unit12,1) RVAL(iblock,i),DER(iblock,i)                   
 1       format(2f20.10)                                                
 10   continue                                                          
      do 20 i = 1,ns                                                    
         REFR(i) = RVAL(iblock,i)                                       
         if (ip.eq.1) then
           write (6,*) refr(i)                                            
         endif
 20   continue                                                          
c........initialize PAIRS(i,j) array to 1                              
      do 21 i = 1,ns                                                    
      do 22 j = 1,2                                                     
         PAIRS(i,j) = 1                                                 
 22   continue                                                          
 21   continue                                                          
c........initialize COUNT(i) to 0                                       
      do 30 i = 1,ns                                                    
         count(i) = 0                                                   
 30   continue                                                          
c                                                                       
c.......Read in the next set of SIC values, calculate the displacement, 
c.......and determine if it corresponds to a totally symmetric          
c.......displacement.                                                   
c                                                                       
         write (6,*)                                                    
         write (6,*) ' Calculating displacements'                        
         write (6,*)                                                    
      do 40 l = 2,nblock                                               
         iblock = iblock + 1                                           
         read (unit12,11,end=999) ener(iblock)                          
	 if (ip.eq.1) then
           write (unit6,*) ' ener',iblock, ener(iblock)                    
           call fioflu(unit6)
	 endif
         check = 0                                                      
         do 50 i = 1,ns                                                 
	    RVAL(iblock,i) = 0.0d+00
	    DER(iblock,i) = 0.0d+00
            read (unit12,1) RVAL(iblock,i), DER(iblock,i)               
            DISP(iblock,i) = 0.0d+00
            DISP(iblock,i) = RVAL(iblock,i) - REFR(i)                   
            if (dabs(DISP(iblock,i)).gt.0.0000002) then                
               check = check + 1                                        
               ncoor(iblock) = i                                        
               count(i) = count(i) + 1                                  
            end if                                                      
            if (check.gt.1) then                                        
               write(unit6,*) 'Error!!!'                                
               write(unit6,*)                                           
               write (unit6,5) iblock                                   
               write (unit6,6)
               write (unit6,7)
               write (unit6,8)
 5    format('READR:  More than one displacement at geometry #',i3,     
     & ' is nonzero.')                                                  
 6     format('        Make sure reference geometry is first in')       
 7     format(' file11 {or 28} and file12 and that the displacements')   
 8     format(' have been generated properly.')                  
               call fioflu(unit6)
               call qabort                                                     
            end if                                                      
            if(dabs(DISP(iblock,i)).gt.0.00000001.and.SIC(i).eq.1.and.  
     &         count(i).gt.2) then                                      
               write(unit6,*) 'Error!!!'                                
               write(unit6,*)                                           
               write(unit6,*) 'READR:  More than two of the displaced ge
     &ometries'                                                         
              write(unit6,*) '        belong to symmetric coordinate ',i
            stop                                                        
            end if                                                      
            if(dabs(DISP(iblock,i)).gt.0.00000001.and.SIC(i).eq.0.and.  
     &      count(i).gt.1) then                                         
               write(unit6,*) 'Error!!!'                                
               write(unit6,*)                                           
               write(unit6,*) 'READR: More than one of the displaced geo
     &metries'                                                          
              write(unit6,*) '       belong to non-symmetric coor. ',i  
            stop                                                        
            end if                                                      
         if(dabs(DISP(iblock,i)).gt.0.00000001) PAIRS(i,count(i))=iblock
 50      continue                                                       
 40   continue                                                          
c........Check to see that there are two displaced geometries for       
c........a totally symmetric SIC and only one for a non-symmetric one   
      if (ip.eq.1) then
        write (6,*)                                                       
        write (6,*) ' Checking for correct number of displacements'        
      endif
      do 55 i = 1,ns                                                    
         if (SIC(i).eq.1.and.count(i).ne.2) then                        
            write(unit6,*) 'Error!!!'                                   
            write(unit6,*)                                              
            write(unit6,*) 'READR:  Require 2 displacements for SIC ',i 
            stop                                                        
         end if                                                         
         if (SIC(i).eq.0.and.count(i).ne.1) then                        
            write(unit6,*) 'Error!!!'                                   
            write(unit6,*)                                              
            write(unit6,*) 'READR:  Require 1 displacement for SIC ', i 
            stop                                                        
         end if                                                         
 55   continue                                                          
c........Calculate the dipole moment derivative wrt the SIC's and       
c........write this information to file18                               
      write (6,*)                                                       
      write (6,*) ' Calculating the dipole derivatives'                  
      write(unit18,2) natom,icharg,(mu(i,1),i=1,3)                      
 2    format(2i5,3f20.10)                                               
croag      call LOCATE(unit5,'# MUDER ##',ierr)                              
croag      if(ierr.ne.0) then                                                
croag         write(unit6,*) 'Error!!!'                                      
croag         write(unit6,*)                                                 
croag         write (unit6,*) 'READR:  No MUDER input found'                 
croag         stop                                                           
croag      end if                                                            
croag      do 95 i = 1,npar                                                  
croag         read (5,14) nc,axis2                                           
croag 14      format(i5,a5)                                                  
croag         axis(nc) = axis2                                               
croag         write (6,*)                                                    
croag         write (6,*) i,axis2                                            
croag         write (6,*)                                                    
croag 95   continue                                                          
      do 600 i=1,npar
        errcod = frdi2('DIP_PAR',ndp,i,1)
        errcod = frdlpn()
        if(errcod.ne.EOK) then
	  write(6,*) 'trouble parsing DIP_PAR vector'
          call fioflu(unit6)
	  call pabort
	  call qabort
	endif
        do 601 j=1,3
	  errcod = frdint(dprty(ndp,j))
	  if(errcod.ne.EOK) then
	    write(6,*) 'trouble parsing DIP_PAR vector'
            call fioflu(unit6)
	    call pabort
	    call qabort
	  endif
	  if ((dprty(ndp,j).ne.1).and.(dprty(ndp,j).ne.-1)) then
	    write(unit6,*) ' DIP_PAR accepts parities of 1 or -1. '
            call fioflu(unit6)
	    call qabort
	  endif
	  write(6,*) ' dprty(', ndp,',',j,') = ',dprty(ndp,j)
601     continue
600   continue
      write(6,*)
croag for displacements which are only positive, and if the change to a neg.
croag   displacement does not change the sign of the dipole component, then
croag   the value delmu = 0.0
croag all other cases delmu = mu1 - mu2
      do 60 i = 1,3                                                     
      delr = 0.0d+00
      delmu = 0.0d+00
      do 70 j = 1,ns                                                    
         delr = DISP(PAIRS(j,1),j) - DISP(PAIRS(j,2),j)                 
         if (sic(j).eq.0) then                                          
            if (dprty(j,i).eq.1) then 
	       delmu = 0.0d+00
            else
               delmu = MU(i,PAIRS(j,1)) - MU(i,PAIRS(j,2))
	    endif
         else
	    delmu = MU(i,PAIRS(j,1)) - MU(i,PAIRS(j,2))
	 endif
	 write(6,*) ' mu(',i,',pairs(', j,',1)) = ',MU(i,PAIRS(j,1))
	 write(6,*) ' mu(',i,',pairs(', j,',2)) = ',MU(i,PAIRS(j,2))
	 write(6,*)
	 write(6,*) ' delr = ', delr
	 write(6,*) ' delmu = ', delmu
	 write(6,*)
         MUDER(i,j) = delmu/delr                                        
 	 write(6,*) ' MUDER(', i,',', j, ') = ', MUDER(i,j)
	 write(6,*)
 70   continue                                                          
 60   continue                                                          
      write(unit18,3)((MUDER(i,j),j=1,ns),i=1,3)                        
 3    format(3f20.10)                                                   
c                                                                       
c    write a modified file12 to file12a for user friendliness           
c                                                                       
      rewind un12a                                                     
      do 80 i = 1,nblock                                                
         Write (un12a,4) ncoor(i),disp(i,ncoor(i)),ener(i)             
 4       format (i5,f12.8,f20.10)                                       
         do 90 j = 1,ns                                                 
            write (un12a,9) RVAL(i,j),DER(i,j)                         
 9          format (2f20.10)                                            
 90      continue                                                       
 80   continue                                                          
      call fioflu(unit6)
      return                                                           
 999  write(unit6,*) 'Error!!!'                                         
      write(unit6,*)                                                    
      write(unit6,12)                                                   
 12   format('READR:  End of file encountered reading file12')          
      write(unit6,13) nblock                                            
 13   format('        Expecting ',i3,' geometry points.')               
      call fioflu(unit6)
      stop                                                              
      end                                                               
