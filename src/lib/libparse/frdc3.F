c----------------------------------------------------------------------------
c See frdc2
      integer function frdc3(keywrd,result,i1,i2,i3)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1,i2,i3
      character*(*) keywrd, result

      integer idx(3)
      character*(lentok) res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdcm

      idx(1) = i1
      idx(2) = i2
      idx(3) = i3

      frdc3 = frdcm(keywrd,1,res,3,idx)
      if (frdc3 .ne. EOK) return
      result = res(1)

      if (result .ne. res(1)) then
        errcod = Etrunc
        frdc3 = errcod
        emess = 'Character data for '''
        call scat(emess,keywrd)
        call scat(emess,''' truncated (frdc3)')
        endif

      return
      end
