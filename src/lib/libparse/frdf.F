c This reads a integer datum from the file unit and places it in result.
c An error code is returned.
      integer function frdf(keywrd,result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      real*8 result

      integer idx(1)
      real*8 res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdfm

      frdf = frdfm(keywrd,1,res,0,idx)
      if (frdf .ne. EOK) return
      result = res(1)

      return
      end
