clj The find core routine that uses the new style input format to look
clj for the MEMORY keyword.  This allows units to be specified with this
clj keyword.
      subroutine fndcor(maxcrr,maxcri,sizerl)
      integer maxcrr,maxcri,sizerl
#include <libparse/error.h>
      character*10 type
      integer errcod,count
      integer defcri
      parameter (defcri=100000)

      real*8 size

      integer fcnt,isizof,frdf1,frdi1,frdi,frdc1

c The user can call psinit again if output to the terminal is not
c desirable
      call psinit('TERMINAL')

      sizerl = isizof('REAL_IN_BYTES')

      errcod = fcnt('MEMORY',count)
      if (errcod.eq.EEOD) then
        maxcrr = defcri
      else if (errcod.ne.EOK) then
        goto 901
      else if (count.eq.-1) then
        errcod = frdi('MEMORY',maxcrr)
        if (errcod.ne.EOK) goto 901
      else if (count.eq.0) then
        maxcrr = defcri
      else if (count.eq.1) then
        errcod = frdi1('MEMORY',maxcrr,1)
        if (errcod.ne.EOK) goto 901
      else if (count.eq.2) then
        errcod = frdf1('MEMORY',size,1)
        if (errcod.ne.EOK) goto 901
        errcod = frdc1('MEMORY',type,2)
        if (errcod.ne.EOK) goto 901
        call touppr(type)
        if ((type.eq.'REAL').or.(type.eq.'R')) then
          maxcrr = size
        else if ((type.eq.'INTEGER').or.(type.eq.'I')) then
          maxcrr = size*isizof('INT')/isizof('REAL')
        else if ((type.eq.'BYTES').or.(type.eq.'B')) then
          maxcrr = size/isizof('REAL_IN_BYTES')
        else if ((type.eq.'KBYTES').or.(type.eq.'KB')) then
          maxcrr = 1000*size/isizof('REAL_IN_BYTES')
        else if ((type.eq.'MBYTES').or.(type.eq.'MB')) then
          maxcrr = 1000000*size/isizof('REAL_IN_BYTES')
        else
          write(6,*) 'bad data type, specify one of:'
          write(6,*) 'REAL, INTEGER, BYTES, KBYTES, or MBYTES'
          goto 901
          endif
        endif

clj problems with multiple psinit calls should be handled in
clj psinit by examining the "inited" variable
clj calling psdone now will, in some cases, reset the output to fort.6
clj      call psdone()

      maxcri = maxcrr*isizof('REAL')/isizof('INT')
      return
  901 continue
      write(6,*) 'FATAL ERROR: problem reading ''MEMORY'''
      write(6,*) 'NOTE: If you do not intend on using the new free'
      write(6,*) '    format input style, then make sure that the'
      write(6,*) '    first character on the first line of your input'
      write(6,*) '    file contains a ''#''.  Otherwise, there has been'
      write(6,*) '    some error in parsing your input.  Hopefully, the'
      write(6,*) '    following information will be useful in finding'
      write(6,*) '    the problem.'
      call pabort
      call qabort
      return
      end
