c----------------------------------------------------------------------------
c  The formatted read integer triangular array of dimension 2 routine.
      integer function frdft2(keywrd,result,n1,ntri)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1,ntri
      character*(*) keywrd
      real*8 result(ntri)

      integer i,j,ij
      real*8 tmp

      integer frdpos,frdlpn,frdflt,frdrpn

      frdft2 = frdpos(keywrd)
      if (frdft2 .ne. EOK) return

      frdft2 = frdlpn()
      if (frdft2 .ne. EOK) return

      ij = 0
      do 101 i = 1, n1
        frdft2 = frdlpn()
        if (frdft2 .ne. EOK) return
        do 201 j = 1, i
          ij = ij + 1
          frdft2 = frdflt(tmp)
          if (frdft2.ne.EOK) return
          result(ij) = tmp
 201      continue
        frdft2 = frdrpn()
        if (frdft2 .ne. EOK) return
 101    continue

      frdft2 = frdrpn()
      end
