      subroutine fioflu(unit)
      implicit integer (a-z)
#if !defined(AIX)
      call flush(unit)
#else
#include <error.h>
      integer errcod, frdc, fioopn
      character*8 term

      term = 'APPEND  '
      errcod = frdc('OUTPUT',term)

      if(term.eq.'TERMINAL') then
        return
      else
        errcod = fiocls(6,'KEEP')
        errcod = fioopn(6,'output.dat','USER','APPEND')
      endif
#endif
      return
      end
