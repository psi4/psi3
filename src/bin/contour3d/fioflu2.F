      subroutine fioflu2(unit)
      implicit integer (a-z)
#if !defined(AIX)
      call flush(unit)
#else
#include <libparse/error.h>
      integer errcod, frdc, fioopn
      character*8 term

      errcod = fiocls(7,'KEEP')
      errcod = fioopn(7,'contour.ps','USER','APPEND')
#endif
      return
      end

