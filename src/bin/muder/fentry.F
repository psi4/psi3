       subroutine fentry(ra,ia,maxcor)
c      program MUDER                                                   
c     Subroutine MUDER(unit5,unit6,unit11,unit12,unit18,natom,multi)    
c
c  put up for rs6000's on 8-31-91 by croag
c   but it only works shitty as of now,  and needs alot of work to 
c    handle ccsd ccsd(t) and the input style of psi2.0
c
C**********************************************************************
C*   NOTICE OF PROGRAM MODIFICATION                                   *
C**********************************************************************
C   Last updated on Febuary 13, 1995 by Evan Bolton                   *
C   Reason:  Added calls to fioflu throughout program.  Why?  You'd   *
C              get an error byt never see it.  Also changed the       *
C              tolerance for changes in internal coordinates due to   *
C              differences between BMAT and INTDER when computing the *
C              internal coordinates of OUT of plane bend.  (six zeroes*
C              now okay rather than the stringent nine)  Put in the   *
C              WFN type SCF.  Won't read dipoles from file11.dat if   *
C              WFN = SCF.  Man page reflects these changes.           *
C              See comments labeled Ceeb.                             *
C**********************************************************************
croag is currently adding psi2.0 input parsing to this beast 2-13
croag only gets around to this once a month do to fucked courses
c
croag it works fine for most cases, but seems to have difficulty
croag  with cases that have degenerate modes like hcch.
c
croag will work on it as he has time.
croag  date is now 8-31-92 // god time flies.  over a year.
c/////////////////////////////////////////////////////////////////////  
c                                                                   //  
c  Written:  Brenda S. Colegrove   (b.t-c.)                         //  
c            December 14,1988                                       //  
c                                                                   //  
c  A Subroutine to calculate the dipole moment derivative from      //  
c  finite displacement points and put them into file18              //  
c                                                                   //
croag      also this modifies file12.dat to format for file12a.dat  //
c                                                                   //
c                                                                   //  
c/////////////////////////////////////////////////////////////////////  
c                  Variable Definitions                             //  
c                                                                   //  
c AXIS(i) = character array which indicated along which axes        //  
c           the dipole moment changes sign for displacement i       //  
c CHAR(maxchr) = array to store character data                      //  
c charge = total molecular charge (real value)                      //  
c check = keeps track of number of non-zero displacements           //  
c count(i) = keeps track of number of displacement points per SIC   //  
c DISP(i,j) = displacement along SIC i for geometry j               //  
c iblock = keeps track of current block of file12 you are at        //  
c icharg = total molecular charge (integral value)                  //  
c maxcor = maximum number of real objects which can be stored       //  
c mu(i,j) = dipole component i at geometry point j (i =x,y,z)       //  
c muder(i,j) = deriv of dipole component i wrt SIC j                //  
c natom = number of atoms                                           //  
c nblock = number of gradient points                                //  
c nc = serial number for SIC                                        // 
c ncoor(i) = the number of the SIC which was displaced in the       //  
c            ith block of file12(i.e. geometry point i)             //  
c npar = number of non-totally symmetric SIC's                      //  
c ns = number of symmetrized internal coordinates (SIC's)           //  
c PAIRS(i,j) = for SIC i, PAIRS(i,1) and PAIRS(i,2) are the numbers //  
c              the geometry points corresponding to the + and -     //  
c              displacements (for symmetric SIC's)                  //  
c pcend = pointer to end of character data                          //  
c pcount = pointer to beginning of COUNT array                      //  
c pend = pointer to end of real data                                //  
c pmu = pointer to beginning of MU array                            //  
c prefr = pointer to beginning of REFR array                        //  
c prval = pointer to beginning of RVAL array                        //  
c psic = pointer to beginning of SIC array                          //  
c RA(maxcor) = array to store real data                             //  
c RVAL(i,j) = value of SIC i at geometry point j                    //  
c SIC(i) = integer array indicating if sic(i) is totally            //  
c           symmetric (1) or non-totally symmetric (0)              //  
c intco =          intco.dat                                        //
c unit5 =          input.dat                                        //  
c unit6 =          output.dat                                       //  
c unit11 =          file11.dat                                      //  
c unit12 =          file12.dat                                      //  
c un12a =          file12a.dat                                    //  
c unit18 =          file18.dat                                      //  
c unit28 =          file28.dat                                      //
c unit30 =          file30.dat                                      //
c                                                                   //  
c/////////////////////////////////////////////////////////////////////  
      implicit real*8 (a-h,o-z)                                         
c      parameter(maxcor=100000)                                          
      parameter(maxchr=1000)                                            
      dimension RA(maxcor),ia(*)
      dimension i30(200)                                                
      real*8 poschg,ZAT(200),charge,negchg                              
      integer unit11, unit12, unit5, unit6, nblock, ns, intco,           
     & npar, pmu, pend, psym,unit18,icharg,unit30,unit28                       
      integer CSHELL(20),OSHELL(20),point,point2,pcshel,pdpar         
      character CHAR(maxchr)                                            
#include <libparse/error.h>
      integer errcod, frdc, frdi, fcnt, frdia1, fwkadd, fwkclr
      integer prdup, prcntr, ip, prup
      integer un12a
      character*8 opentp, citype
      logical oldinp
      common/croag/ oldinp
      common/un12a/ un12a
      common/unit/ intco,unit5,unit6,unit11,unit12,unit18,unit28
C     S.Brown: Added in order for the program to read occupations from file30
      integer ITEMP(200),OCCTEMP(8)
      integer NSYMHF,opoint
      character*4 LAB(8),HFLAB(8)
Ceeb
      common/type/ citype
c
      call psinit('APPEND')
      call zero(ra,maxcor)
      call zero(zat,200)
      call izero(cshell,20)
      call izero(oshell,20)
      mpoint = 200                                                      
      mconst=200                                                        
      mcalcs=100                                                        
      unit5=5                                                           
      unit6=6                                                           
      unit11=11                                                         
      unit12=12                                                         
      un12a=14                                                         
      unit18=18                                                        
      unit28=28                                                        
      unit30=30                                                        
      intco=7
      call tstart(unit6)
c
c  print control options
c
c have three modes for now, DEBUG, VERBOSE, and BRIEF
c       BRIEF  prints only the barest of information to output.dat
c	VERBOSE  prints everything that b.t-c. had plus some
c       DEBUG  will continue to change...prints everything I can find
c
      errcod = prcntr('ON BRIEF')
c        write(unit6,*)
c    	write(unit6,*) ' after on brief, returned is: '
c	write(unit6,*) ' brief ',prcntr('IS_ON BRIEF')
      errcod = prcntr('OFF VERBOSE DEBUG')
c        write(unit6,*)
c    	write(unit6,*) ' after off verb debug, returned is: '
c	write(unit6,*) ' verbose ',prcntr('IS_ON VERBOSE')
c	write(unit6,*) ' debug ',prcntr('IS_ON DEBUG')
      errcod = prdup()
c        write(unit6,*)
c        write(unit6,*) ' after off read input, returned is: '
c        write(unit6,*) ' brief ',prcntr('IS_ON BRIEF')
c        write(unit6,*) ' verbose ',prcntr('IS_ON VERBOSE')
c        write(unit6,*) ' debug ',prcntr('IS_ON DEBUG')
c
      if(prcntr('IS_ON DEBUG').eq.1) errcod = prcntr('ON VERBOSE')
      if(prcntr('IS_ON VERBOSE').eq.1) errcod = prcntr('OFF BRIEF')
c        write(unit6,*)
c        write(unit6,*) ' after off reset, returned is: '
c        write(unit6,*) ' brief ',prcntr('IS_ON BRIEF')
c        write(unit6,*) ' verbose ',prcntr('IS_ON VERBOSE')
c        write(unit6,*) ' debug ',prcntr('IS_ON DEBUG')
c
      call ffile(unit18,'file18',0)
      call ffile(unit12,'file12',1)
      call ffile(un12a,'file12a',0)
      call rfile(unit30)                                                
      ip = 0
      ip = prcntr('IS_ON DEBUG')
      if(ip.eq.1) then
        write(unit6,*)
        write(unit6,*) ' finished through rfile of unit30 '
      endif
      call izero(i30,mconst)
c                                                                       
c   Read in a bunch of stuff to locate the nuclear charges and the      
c   orbital occupancies in file30                                       
c                                                                       
      call wreadw(unit30,i30,mconst,101,junk)                              
      if(ip.eq.1) then
        write(unit6,*)
        write(unit6,*) ' finished through read of i30'
      endif
croag  -total length of file30
      itotal=i30(1)                                                     
croag  -number of calc. sections??
      ncalcs=i30(5)                                                     
croag  -number of atoms
      natom=i30(19)                                                     
croag  -number of M.O.'s
      msfru=i30(18)                                                     
croag  -number of irreps in SCF (after elimin. of irreps w/ no so//mo)
      nsymhf=i30(41)                                                    
croag  -number of m.o. coeff. (after all elim. due to symm.)
      mxcoef=i30(42)                                                    
croag  -number of open shells
      iopen=i30(43)                                                     
croag  -????????????????????????
      nt=i30(29)                                                        
croag  -maximum number of primitive basisfn.
      mxpr=i30(32)                                                      
croag  -number of shells
      nshell=i30(27)                                                    
croag  -total number of irreps.
      nirred=i30(28)                                                    
croag  -????????????????????????
      mxcf=i30(35)                                                      
croag  -????????????????????????
      mxt=i30(36)     
      
      write(unit6,*) 'nirred = ',nirred,'nsymhf = ',nsymhf,'mxcoef = ',mxcoef                    
c                                                                       
      npsym = (nirred + 3)/4                                            
      niso = (nt + 3)/4                                                 
c                                                                       
      if (itotal.eq.0) then                                             
        write(unit6,*)
        write (unit6,*) 'MUDER:  Nothing in file30'                      
        stop                                                           
      end if                                                            
c    Calculate total nuclear charge                                     
c                                                                       
c     ichg = pointer to position of atomic charges in file30            
c                                                                       
      ichg=mcalcs + mpoint + mconst +101                                
      call wreadw(unit30,ZAT,intowp(natom),ichg,junk)                   
      call rclose(unit30,3)
      if (ip.eq.1) then
        write(unit6,*)
        write(unit6,*) ' finished through read of ZAT'
      endif
      poschg = 0.0d0                                                    
      do 3 i = 1,natom                                                  
        poschg = poschg + ZAT(i)                                       
 3    continue                                                          
c                                                                       
c     Calculate number of electrons                                     
c                                                                       
croag
croag   we can just find cshell on oshell in the input, so let's parse it
croag   there.  thus we do not need to worry about funky pointers in file30
croag   for each different case of whether there are f & g and 5d vs. 6d, etc.
croag
croag   the following series of commented code is b.t-c. original.
croag   following it is standard psi2.0 parsing.
croag
croag point = 161 + mconst + mpoint + mcalcs + intowp(natom) + nt*natom 
croag& + 11*natom + intowp(mxpr) + 3*(intowp(mxpr)) + 7*nshell          
croag& + niso*nshell + npsym*nshell*2 + 4*nirred + intowp(mxcf)         
croag& + mxcf + intowp(9*nt) + nt + intowp(3*3*nt) + intowp(6*6*nt)     
croag call wreadw(unit30,point2,1,point,junk)                           
croag pcshel = point2 + intowp(mxcoef) + intowp(msfru) + 2*nsymhf       
ccroag                                                                  
ccroagpcshel = pointer to where the orbital occupancies are in file30.  
ccroag                                                                  
croag call wreadw(unit30,CSHELL,nsymhf,pcshel,junk)                     
croag    call wreadw(unit30,OSHELL,nsymhf,junk,junk)                    
croag
croag we'll check for old style input, but this code will not support it.
croag
      if (ip.eq.1) then
        write(unit6,*)
	write(unit6,*) ' we are at the beginning of input parsing. '
      endif
      oldinp = .true.
      CALL LOCATE(unit5,'# MUDER ##',IERR)
      IF(IERR .ne. 0) then
        oldinp = .false.
      endif
C
      ip = 0
      ip = prcntr('IS_ON VERBOSE')
      if (ip.eq.1) then
        write(unit6,*)
        if (oldinp) then
          write(unit6,*) ' oldinp = true '
        else
          write(unit6,*) ' oldinp = false '
        endif
      endif
      if(oldinp) then
	write(unit6,*) ' This muder code is only psi2.0 compatible, '
        write(unit6,*) ' and will not support psi1.0 input format. '
        call fioflu(unit6)
	call qabort
      else
C     S.Brown: Get those occupations from File30 
C        errcod = fcnt('DOCC',nclshl)
C	errcod = frdia1('DOCC',cshell,nclshl)
C	if (iopen.ne.0) then
C          errcod = fcnt('SOCC',nopshl)
C          errcod = frdia1('SOCC',oshell,nopshl)      
         DO 7000 I=1,MPOINT
            ITEMP(I)=0
 7000    continue
         call rfile(unit30)
         call wreadw (unit30,itemp,mpoint,301,junk)
         call wreadw (unit30,lab,nirred,itemp(16),junk)
         call wreadw (unit30,itemp,mcalcs,501,junk)
         call wreadw (unit30,opoint,1,itemp(1)+60,junk)
         opoint = opoint+(2*mxcoef)+(2*msfru)
         call wreadw (unit30,hflab,nsymhf,opoint,opoint)
         call wreadw (unit30,cshell,nsymhf,opoint+nsymhf,opoint)
         if(iopen.ne.0) call wreadw (unit30,oshell,nsymhf,opoint,opoint)
         M=0
         if(nirred.ne.nsymhf) THEN
            do 7100 I=1,nirred
               if(lab(i+m).ne.hflab(i)) THEN
                  do 7110 k=nirred,i+1,-1
                     cshell(k)=cshell(k-1)
                     oshell(k)=oshell(k-1)
 7110             continue
                  cshell(i)=0
                  oshell(i)=0
                  m=m+1
               endif
 7100       continue
         endif
      call rclose(unit30,3)  
C     endif
      endif
      if (ip.eq.1) then
        write(unit6,*)
c        write(unit6,*)' DOCC = ',cshell(1),cshell(2),cshell(3),cshell(4)
	write(unit6,*) 'DOCC = ',(cshell(j), j=1,nirred)
        if (iopen.ne.0) then
          write(unit6,*)
	  write(unit6,*) 'SOCC = ',(oshell(j), j=1,nirred)
c          write(unit6,*)' SOCC = ',
c     &                    oshell(1),oshell(2),oshell(3),oshell(4)
        endif
      endif
      negchg = 0                                                        
      do 4 i = 1,nirred                                                 
         negchg = negchg + 2*CSHELL(i)                                  
 4    continue                                                          
      if (iopen.ne.0) then                                              
         do 5 i = i,nirred                                              
            negchg = negchg + oshell(i)                                 
 5       continue                                                       
      end if                                                            
c                                                                       
c     total molecular charge                                            
c                                                                       
c     Since the total molecular charge may be non-integer, we convert   
c     to an integer value and abort for non-integer values.             
c                                                                       
c	write(unit6,*) 'poschg = ', poschg
c	write(unit6,*) 'negchg = ', negchg
      charge = poschg - negchg                                          
      if((charge-int(charge)).ne.0.0) then                              
        write(6,*)'MUDER:  Non-integer value for total positive charge.'
      end if                                                            
      icharg = int(charge)                                              
croag  more of the new style input parsing.  i'll leave in the old stuff,
croag    simply b/c i may need it for testing.  nobody need know it's here.
      ns = 0
      npar = 0
      nas = 0
      oldinp = .true.
      CALL LOCATE(unit5,'# INTDIF ##',IERR)
      IF(IERR .ne. 0) then
        oldinp = .false.
      endif
C
      if (ip.eq.1) then
        write(unit6,*)
        if (oldinp) then
          write(unit6,*) ' oldinp = true '
        else
          write(unit6,*) ' oldinp = false '
        endif
      endif
      if(oldinp) then
        write(unit6,*) ' Getting input from # INTDIF ##. '
	read (unit5,*)
	read (unit5,2) ns,npar
croag
croag  i'm not sure that this will work for diatomics.  in fact i'm fairly 
croag  sure that it won't.  sorry.  suffer.
croag
      else
        call ffile(intco,'intco',0)
	errcod = fwkadd(unit5,'INTCO')
	errcod = fwkadd(intco,'INTCO')
c	
	errcod = fcnt('SYMM',ns)
	errcod = fcnt('ASYMM',nas)
	ns=ns+nas
c
	errcod = fwkclr()
        errcod = fwkadd(unit5,'MUDER')
        errcod = fwkadd(unit5,'DEFAULT')
	errcod = frdi('N_INTCO',ns)
        ip = 0
        ip = prcntr('IS_ON DEBUG')
        if (ip.eq.1) then
          write(unit6,*) ' frdi NINTCO '
          if (oldinp) then
            write(unit6,*) ' oldinp = true '
          else
            write(unit6,*) ' oldinp = false '
          endif
        endif
	errcod = fcnt('DIP_PAR',npar)
        if (ip.eq.1) then
          write(unit6,*) ' fcnt DIPPAR '
          if (oldinp) then
            write(unit6,*) ' oldinp = true '
          else
            write(unit6,*) ' oldinp = false '
          endif
        endif
      end if
c
      ip = 0
      ip = prcntr('IS_ON VERBOSE')
      if (ip.eq.1) then
        write(unit6,*)
        write (6,*) ' ns =',ns,' npar =',npar                                 
      endif
 2    format(i5,10x,i5)                                                 
croag
croag  this next line is botched for cases with degenerate modes, 
croag   for such cases nblock = nblock - 2*n.degenerate.pairs
croag
      nblock = 1 + 2*(ns-npar) + npar                                   
c     if(nblock.ne.multi) return                                        
c........Set up pointers to arrays                                      
      kdpar = 1                                                         
      kcend = kaxis + npar*3                                              
      kmu = 1                                                           
      ksic = kmu + 3*nblock                                             
      krval = ksic + ns                                                 
      kpair = krval + ns*nblock                                         
      kncoor = kpair + 2*ns                                             
      kdisp = kncoor + nblock                                           
      kmuder = kdisp + ns*nblock                                        
      krefr = kmuder + 3*ns                                             
      kcount = krefr + ns                                               
      kder = kcount + ns                                                
      kener = kder + ns*nblock                                          
      kend = kener + nblock                                             
      call zero(ra,kend+1)
c                                                                       
c.......test to see if amount of data is larger than allocated space    
c                                                                       
      if(kend.gt.maxcor) then                                           
         write(unit6,*) 'Error!!!'                                      
         write(unit6,*)                                                 
         write (unit6,*) 'MUDER:  You bozo!  Storage capacity exceded.' 
         write (unit6,*) '        Increase maxcor to ',pend             
         stop                                                           
      end if                                                            
      if(kcend.gt.maxchr) then                                          
         write(unit6,*) 'Error!!!'                                      
         write(unit6,*)                                                 
         write (unit6,*) 'MUDER:  You bozo!  Storage capacity exceded.' 
         write (unit6,*) '        Increase maxchr to ',pcend            
         stop                                                           
      end if                                                            
      if (ip.eq.1) then
        write(unit6,*)
        write(unit6,*) ' just finished setting pointers '
        write(unit6,*)
        if (oldinp) then
          write(unit6,*) ' oldinp = true '
        else
          write(unit6,*) ' oldinp = false '
        endif
      endif
      ip = 0
      ip = prcntr('IS_ON DEBUG')
      if (ip.eq.1) then
        write (6,*)                                                       
        write (6,*) ' Entering READMU'                                     
      endif                                                             
      call READMU(RA(kmu),nblock,natom)                    
      if (ip.eq.1) then
        write (6,*)                                                       
        write (6,*) ' Exited READMU'                                       
        if (oldinp) then
          write(unit6,*) ' oldinp = true '
        else
          write(unit6,*) ' oldinp = false '
        endif
        write (6,*) ' Entering INISIC'                                     
      endif
      call INISIC(RA(ksic),ns,npar)                         
      if (ip.eq.1) then
        write (6,*)                                                       
        write (6,*) ' Exited INISIC'                                       
        write (6,*) ' Entering READR'                                      
      endif
      call READR(RA(krval),RA(kpair),char(kdpar),
     & RA(kncoor),RA(kdisp),RA(kmu),RA(ksic),RA(krefr),RA(kmuder),      
     & RA(kcount),ns,nblock,icharg,natom,RA(kder),RA(kener),npar)       
      if (ip.eq.1) then
        write (6,*)                                                      
        write (6,*) ' Exited READR'                                        
      endif
Ceeb
      if (citype.eq.'SCF     ') then
        write(6,*) '  '
        write(6,*) ' CITYP = SCF so delete file18'
        write(6,*) '  '
        call fioflu(unit6)
        close (unit=unit18,status='DELETE',iostat=ierr)
      endif
Ceeb
      call tstop(unit6)
      call psdone()
      stop                                                              
      end                                                               
c                                                                       
