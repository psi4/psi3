c----------------------------------------------------------------------------
c  The formatted read character triangular array of dimension 2 routine.
      integer function frdct2(keywrd,result,n1,ntri)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1,ntri
      character*(*) keywrd, result(ntri)

      integer i,j,ij
      character*(lentok) tmp

      integer frdpos,frdlpn,frdchr,frdrpn

      frdct2 = frdpos(keywrd)
      if (frdct2 .ne. EOK) return

      frdct2 = frdlpn()
      if (frdct2 .ne. EOK) return

      ij = 0
      do 101 i = 1, n1
        frdct2 = frdlpn()
        if (frdct2 .ne. EOK) return
        do 201 j = 1, i
          ij = ij + 1
          frdct2 = frdchr(tmp)
          if (frdct2.ne.EOK) return
          result(ij) = tmp
 201      continue
        frdct2 = frdrpn()
        if (frdct2 .ne. EOK) return
 101    continue

      frdct2 = frdrpn()
      end
