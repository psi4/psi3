      integer function cosize(c,ptr)
#if defined(sun)
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      integer ptr
      integer c(1)
#include "share.in"

      cosize = c(ptr-3)

      return
      end

