c----------------------------------------------------------------------------
c See frdi2
      integer function frdi4(keywrd,result,i1,i2,i3,i4)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2,i3,i4
      character*(*) keywrd
      integer result

      integer idx(4)
      integer res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdim

      idx(1) = i1
      idx(2) = i2
      idx(3) = i3
      idx(4) = i4

      frdi4 = frdim(keywrd,1,res,4,idx)
      result = res(1)

      return
      end
