c The parser abort sequence
      subroutine pabort
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer lineno, linept, errcod
      character*(lentok) line, where
      character*(80) emess
      character*(lentok) progid
      integer cunit, ncwk, units(maxcwk)
      common /parse/ lineno, linept, cunit, ncwk, units
      common /curlin/ line
      common /lparer/ errcod, emess

      integer pt,i

      if (linept .gt. 1) then
        pt = linept - 1
      else
        pt = 1
        endif

      write(6,*)
      write(6,*) 'pabort: status of parser:'
      write(6,*) ' errcod = ', errcod
      write(6,*) ' on line ',lineno,' of input at position ',pt,':'
      write(6,'(A80)') line

      where = ' '

      do 101 i = 1, pt
        where(i:i) = '-'
 101    continue

      where(pt:pt) = '|'

      write(6,'(A80)') where
      write(6,*)
      write(6,*) ' the error message buffer contains:'
      write(6,*) emess
      write(6,*)
      write(6,*) ' the program identifier for the aborting program:'
      call gprgid(progid)
      write(6,*) progid
      write(6,*)

      return
      end
