c----------------------------------------------------------------------------
c  The formatted read integer array of dimension 2 routine.
      integer function frdia2(keywrd,result,n1,n2)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer n1,n2
      character*(*) keywrd
      integer result(n1,n2)

      integer i,j,tmp

      integer frdpos,frdlpn,frdint,frdrpn

      frdia2 = frdpos(keywrd)
      if (frdia2 .ne. EOK) return

      frdia2 = frdlpn()
      if (frdia2 .ne. EOK) return

      do 101 i = 1, n2
        frdia2 = frdlpn()
        if (frdia2 .ne. EOK) return
        do 201 j = 1, n1
          frdia2 = frdint(tmp)
          if (frdia2.ne.EOK) return
          result(j,i) = tmp
 201      continue
        frdia2 = frdrpn()
        if (frdia2 .ne. EOK) return
 101    continue

      frdia2 = frdrpn()
      end
