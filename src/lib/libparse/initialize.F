c This is the master initialization routine.  It calls any other needed
c initialization routines.
c  status = 'WRITE' will delete the output file.
c  status = 'APPEND' will append to the output file.
c Error recovery method: print a message and continue
      subroutine psinit(status)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) status
#include "error.h"
#include <misc.h>

      integer errcod
      character*80 newstt

      integer fiocls,fioopn,frdi,chrlen

      integer first
      logical inited
      save first
      save inited
      data first/0/
      data inited/.false./

#ifdef sun
      call ieee_init()
#endif

c Set up the primary input and output.
      if (inited) errcod = fiocls(5, 'KEEP')
      errcod = fioopn(5,'input.dat','USER','READ')
      if (errcod .ne. EOK) then
         call prnerr()
         call qabort
      endif

c Set the working keyword list to the default list.
      call fwkdef()

c Set up the print control info to the default values.
c This is needed because the library routines may need to use print
c control
      call prdup()

c Try to read a nonexistant keyword to parse thru the entire input file
c and search for errors.
      errcod = frdc('YYYNO_SUCH_KEYWORDXXXX',newstt)
      if ((first .eq. 0).and.(errcod .ne. EOK)
     &    .and.(errcod .ne. EEOD)) then
        first = 1
        write(6,*) 'ERROR: There has been an error in a test parse'
        write(6,*) '    of your input file.'
        write(6,*) 'NOTE: If you do not intend on using the new free'
        write(6,*) '    format input style, then make sure that the'
        write(6,*) '    first character on the first line of your input'
        write(6,*) '    file contains a ''#''.'
        call pabort
        call qabort
        endif

c See if the output status cannot be overridden
      if (status(chrlen(status):chrlen(status)).eq.'!') then
        newstt = status(1:chrlen(status)-1)
      else
c See if we want to override the output status.
        errcod = frdc('OUTPUT',newstt)
        if (errcod .ne. EOK) then
          newstt = status
        else if (     (newstt .ne. 'TERMINAL')
     &           .and.(newstt .ne. 'WRITE')
     &           .and.(newstt .ne. 'APPEND') ) then
          write(6,*) 'ERROR: an attempt was made to set ''OUTPUT'''
          write(6,*) '      to something other than ''TERMINAL'','
          write(6,*) '      ''WRITE'', or ''APPEND''.'
          write(6,*) '      The request has been ignored.'
          newstt = status
          endif
        endif

c Don't open 6 if status is terminal.  This will hopefully write the
c output to the screen.
      if (newstt.ne.'TERMINAL') then
        if (inited) errcod = fiocls(6, 'KEEP')
        errcod = fioopn(6,'output.dat','USER',newstt)
        if (errcod .ne. EOK) call prnerr()
        endif

c Open the check file.
c Don't open 3 if status is terminal.  This will keep AIX 4.1 from screwing up
c I hope. (TDC)
      if (newstt.ne.'TERMINAL') then
        if (inited) errcod = fiocls(3, 'KEEP')
        errcod = fioopn(3,'tape3.dat','USER','WRITE')
        if (errcod .ne. EOK) call prnerr()
      endif

c Determine the linewidth of the screen.
      errcod = frdi('LINEWIDTH',linewd)
      if (errcod .eq. EEOD) then
        linewd = 113
      else if (errcod .ne. EOK) then
        linewd = 113
        write(6,*)'''LINEWIDTH'' could not be set as requested(psinit):'
        call prnerr()
        endif

      inited = .true.
      return
      end
c--------------------------------------------------------------------------
c This routine is used to clean up when we are done.
      subroutine psdone()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif

#ifdef sun
      call ieee_done()
#endif

      call fioflu(6)

crawdad Trying to fix the errors in multiple psinit() calls
	errcod = fiocls(6, 'KEEP')
c      close(6)
	errcod = fiocls(5, 'KEEP')
c      close(5)

      return
      end
