c----------------------------------------------------------------------------
c This reads a integer datum from the file unit and places it in result.
c An error code is returned. The datum is the i2th element of a vector found
c in the i1th position of the top level vector
      integer function frdf2(keywrd,result,i1,i2)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2
      character*(*) keywrd
      real*8 result

      integer idx(2)
      real*8 res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdfm

      idx(1) = i1
      idx(2) = i2

      frdf2 = frdfm(keywrd,1,res,2,idx)
      result = res(1)

      return
      end
