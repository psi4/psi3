c----------------------------------------------------------------------------
c See frdi2
      integer function frdi3(keywrd,result,i1,i2,i3)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2,i3
      character*(*) keywrd
      integer result

      integer idx(3)
      integer res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdim

      idx(1) = i1
      idx(2) = i2
      idx(3) = i3

      frdi3 = frdim(keywrd,1,res,3,idx)
      result = res(1)

      return
      end
