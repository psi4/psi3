c----------------------------------------------------------------------------
c This reads a character datum from the file unit and places it in result.
c An error code is returned. The datum is the i1th element of a vector
      integer function frdc1(keywrd,result,i1)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer i1
      character*(*) keywrd, result

      integer idx(1)
      character*(lentok) res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdcm

      idx(1) = i1

      frdc1 = frdcm(keywrd,1,res,1,idx)
      if (frdc1 .ne. EOK) return
      result = res(1)

      if (result .ne. res(1)) then
        errcod = Etrunc
        frdc1 = errcod
        emess = 'Character data for '''
        call scat(emess,keywrd)
        call scat(emess,''' truncated (frdc1)')
        endif

      return
      end
