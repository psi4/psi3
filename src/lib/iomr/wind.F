*Deck Wind
      subroutine wind(unit)
c
      implicit integer (a-z)
c
      character*1 temp
c
    1 continue
         read (unit,2,end=3,err=4) temp
      go to 1
c
    2 format (a1)
c
    3 continue
      return
c
c
    4 continue
      call lnkerr('I/O error in WIND')
c
      stop
      end
