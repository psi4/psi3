c  The Track Unit routine.
c This keeps information about what units are in use and what they are
c used for.
c  unit is the unit number 0 means a new unit number is to be allocated.
c  name is a name to be recorded for the file.
c  form is either formatted or unformatted.
c  requst is the task to be done, it is one of
c    INIT: initialize the file (this will not open the file)
c    DUMP: print info about files
c    GET: get information about the files
c    DELETE: remove info about a file (this does not close the file)
c Error recovery: returns values from error.h
      integer function trunit(unit,name,form,requst)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer unit
      character*(*) name,form,requst

      integer istart, istop, i
      character*(10) names(maxunt)
      character*30 tmp
      integer forms(maxunt), iform
      save names
      save forms
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess
      data forms/maxunt*0/, names/maxunt*' '/


      if ((unit.gt.maxunt).or.(unit.lt.0)) then
        errcod = Earg
        trunit = errcod
        emess = 'Bad unit number (trunit)'
        return
        endif

      if (requst.eq.'INIT') then
        if (unit.eq.0) then
c         This is avoids conflicts with the old style unit assignment
c         during the transition
          if (form .eq. 'FORMATTED') then
            istart = 10
            istop  = 29
            iform  = 1
          else if (form .eq. 'UNFORMATTED') then
            istart = 30
            istop  = maxunt
            iform  = 2
            endif
          do 101 i = istart, istop
            if (forms(i) .eq. 0) then
              unit = i
              goto 901
              endif
  101       continue
  901     continue
        else if (forms(unit).ne.0) then
          errcod = Eopen
          trunit = errcod
          emess = 'Unit'
          call intchr(unit,tmp)
          call scatb(emess,tmp)
          call scatb(emess, 'is already in use (trunit)')
          return
          endif
        names(unit) = name
        forms(unit) = iform
      else if (requst.eq.'DELETE') then
        forms(unit) = 0
        names(unit) = ' '
      else if (requst.eq.'DUMP') then
        write(6,*) 'trunit: a DUMP has been requested'
        do 102 i = 1, maxunt
          if (forms(i).ne.0) then
            write(6,*) 'forms(',i,')=',forms(i)
            endif
 102      continue
      else if (requst.eq.'GET') then
        if (forms(unit).eq.1) then
          form = 'FORMATTED'
        else if (forms(unit).eq.2) then
          form = 'UNFORMATTED'
        else
          form = 'UNKNOWN'
          endif
        name = names(unit)
        endif

      errcod = EOK
      trunit = errcod
      emess = ' '

      return
      end
