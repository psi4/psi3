c----------------------------------------------------------------------------
c  The formatted read integer array of dimension 2 routine.
      integer function frdfa2(keywrd,result,n1,n2)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1,n2
      character*(*) keywrd
      real*8 result(n1,n2)

      integer i,j
      real*8 tmp

      integer frdpos,frdlpn,frdflt,frdrpn

      frdfa2 = frdpos(keywrd)
      if (frdfa2 .ne. EOK) return

      frdfa2 = frdlpn()
      if (frdfa2 .ne. EOK) return

      do 101 i = 1, n2
        frdfa2 = frdlpn()
        if (frdfa2 .ne. EOK) return
        do 201 j = 1, n1
          frdfa2 = frdflt(tmp)
          if (frdfa2.ne.EOK) return
          result(j,i) = tmp
 201      continue
        frdfa2 = frdrpn()
        if (frdfa2 .ne. EOK) return
 101    continue

      frdfa2 = frdrpn()
      end
