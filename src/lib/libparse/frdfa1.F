c----------------------------------------------------------------------------
c  The formatted read integer array of dimension 1 routine.
c Read a one dimensional character array of length n1 into the result array.
      integer function frdfa1(keywrd,result,n1)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1
      character*(*) keywrd
      real*8 result(n1)

      integer i
      real*8 tmp

      integer frdpos,frdlpn,frdflt,frdrpn

      frdfa1 = frdpos(keywrd)
      if (frdfa1 .ne. EOK) return

      frdfa1 = frdlpn()
      if (frdfa1 .ne. EOK) return

      do 101 i = 1, n1
        frdfa1 = frdflt(tmp)
        if (frdfa1.ne.EOK) return
        result(i) = tmp
 101    continue

      frdfa1 = frdrpn()
      end
