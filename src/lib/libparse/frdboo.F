c This reads a boolean datum from the file unit and places it in result.
c An error code is returned.  It can be given as YES NO TRUE FALSE 1 or 0.
      integer function frdboo(keywrd,result)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd
      integer result

      character*80 cres
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdc

      frdboo = frdc(keywrd,cres)
      if (frdboo .ne. EOK) return

      if (    (cres .eq. 'YES')
     &    .or.(cres .eq. 'TRUE')
     &    .or.(cres .eq. '1')) then
        result = 1
        return
        endif

      if (    (cres .eq. 'NO')
     &    .or.(cres .eq. 'FALSE')
     &    .or.(cres .eq. '0')) then
        result = 0
        return
        endif

      frdboo = Etype
      emess = 'Expected a boolean type but got, '''
      call scat(emess,cres)
      call scat(emess,''' (frdboo)')

      return
      end
