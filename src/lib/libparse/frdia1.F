c----------------------------------------------------------------------------
c  The formatted read integer array of dimension 1 routine.
c Read a one dimensional character array of length n1 into the result array.
      integer function frdia1(keywrd,result,n1)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1
      character*(*) keywrd
      integer result(n1)

      integer i, tmp

      integer frdpos,frdlpn,frdint,frdrpn

      frdia1 = frdpos(keywrd)
      if (frdia1 .ne. EOK) return

      frdia1 = frdlpn()
      if (frdia1 .ne. EOK) return

      do 101 i = 1, n1
        frdia1 = frdint(tmp)
        if (frdia1.ne.EOK) return
        result(i) = tmp
 101    continue

      frdia1 = frdrpn()
      end
