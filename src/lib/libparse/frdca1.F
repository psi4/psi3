c----------------------------------------------------------------------------
c  The formatted read character array of dimension 1 routine.
c Read a one dimensional character array of length n1 into the result array.
      integer function frdca1(keywrd,result,n1)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1
      character*(*) keywrd, result(n1)

      integer i
      character*(lentok) tmp

      integer frdpos,frdlpn,frdchr,frdrpn

      frdca1 = frdpos(keywrd)
      if (frdca1 .ne. EOK) return

      frdca1 = frdlpn()
      if (frdca1 .ne. EOK) return

      do 101 i = 1, n1
        frdca1 = frdchr(tmp)
        if (frdca1.ne.EOK) return
        result(i) = tmp
 101    continue

      frdca1 = frdrpn()
      end
