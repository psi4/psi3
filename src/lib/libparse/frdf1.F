c----------------------------------------------------------------------------
c This reads a integer datum from the file unit and places it in result.
c An error code is returned. The datum is the i1th element of a vector
      integer function frdf1(keywrd,result,i1)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1
      character*(*) keywrd
      real*8 result

      integer idx(1)
      real*8 res(1)

      integer frdfm

      idx(1) = i1

      frdf1 = frdfm(keywrd,1,res,1,idx)
      result = res(1)

      return
      end
