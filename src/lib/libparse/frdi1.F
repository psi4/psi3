c----------------------------------------------------------------------------
c This reads a integer datum from the file unit and places it in result.
c An error code is returned. The datum is the i1th element of a vector
      integer function frdi1(keywrd,result,i1)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1
      character*(*) keywrd
      integer result

      integer idx(1)
      integer res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdim

      idx(1) = i1

      frdi1 = frdim(keywrd,1,res,1,idx)
      result = res(1)

      return
      end
