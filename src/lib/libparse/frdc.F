c This reads a character datum from the file unit and places it in result.
c An error code is returned.
      integer function frdc(keywrd,result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd, result

      integer idx(1)
      character*(lentok) res(1)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdcm

      frdc = frdcm(keywrd,1,res,0,idx)
      if (frdc .ne. EOK) return
      result = res(1)

      if (result .ne. res(1)) then
        errcod = Etrunc
        frdc = errcod
        emess = 'Character data for '''
        call scat(emess,keywrd)
        call scat(emess,''' truncated (frdc)')
        endif

      return
      end
