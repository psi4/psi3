c----------------------------------------------------------------------------
c This reads a integer datum from the file unit and places it in result.
c An error code is returned. The datum is the i2th element of a vector found
c in the i1th position of the top level vector
      integer function frdi2(keywrd,result,i1,i2)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2
      character*(*) keywrd
      integer result

      integer idx(2)
      integer res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdim

      idx(1) = i1
      idx(2) = i2

      frdi2 = frdim(keywrd,1,res,2,idx)
      result = res(1)

      return
      end
