c This contains the error code processing utilities
c------------------------------------------------------------------------
c  Print out an error message
      subroutine prnerr()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      write(6,'(A80)') emess

      return
      end
