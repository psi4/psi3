c----------------------------------------------------------------------------
c See frdf2
      integer function frdf4(keywrd,result,i1,i2,i3,i4)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2,i3,i4
      character*(*) keywrd
      real*8 result

      integer idx(4)
      real*8 res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdfm

      idx(1) = i1
      idx(2) = i2
      idx(3) = i3
      idx(4) = i4

      frdf4 = frdfm(keywrd,1,res,4,idx)
      result = res(1)

      return
      end
