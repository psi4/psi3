c This reads a integer datum from the file unit and places it in result.
c An error code is returned.
      integer function frdi(keywrd,result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer result

      integer idx(1)
      integer res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdim

      frdi = frdim(keywrd,1,res,0,idx)
      if (frdi .ne. EOK) return
      result = res(1)

      return
      end
