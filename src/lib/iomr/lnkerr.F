*Deck LnkErr
      subroutine lnkerr(error)
c
      character*(*) error
c
#ifndef IBMMVS
clj      call abort(error)
      call mabort
      stop
#else
      stop 876
#endif
      end
