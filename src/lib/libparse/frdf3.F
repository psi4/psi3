c----------------------------------------------------------------------------
c See frdf2
      integer function frdf3(keywrd,result,i1,i2,i3)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2,i3
      character*(*) keywrd
      real*8 result

      integer idx(3)
      real*8 res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdfm

      idx(1) = i1
      idx(2) = i2
      idx(3) = i3

      frdf3 = frdfm(keywrd,1,res,3,idx)
      result = res(1)

      return
      end
