c----------------------------------------------------------------------------
c  The formatted read character array of dimension 2 routine.
      integer function frdca2(keywrd,result,n1,n2)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1,n2
      character*(*) keywrd, result(n1,n2)

      integer i,j
      character*(lentok) tmp

      integer frdpos,frdlpn,frdchr,frdrpn

      frdca2 = frdpos(keywrd)
      if (frdca2 .ne. EOK) return

      frdca2 = frdlpn()
      if (frdca2 .ne. EOK) return

      do 101 i = 1, n2
        frdca2 = frdlpn()
        if (frdca2 .ne. EOK) return
        do 201 j = 1, n1
          frdca2 = frdchr(tmp)
          if (frdca2.ne.EOK) return
          result(j,i) = tmp
 201      continue
        frdca2 = frdrpn()
        if (frdca2 .ne. EOK) return
 101    continue

      frdca2 = frdrpn()
      end
