c These routines deal with the current working keyword/working unit pairs.
c Since these routines are to be used by knowledgable programmers I am not
c doing much error checking.  However, I do avoid memory overwrites and other
c such potentially disastrous errors.
c---------------------------------------------------------------------------
c  The formatted i/o working keyword set defaults routine.
c This routine sets the current working keyword/unit to the defaults
c (k/u to progid/5 first and DEFAULT/5 second)
      subroutine fwkdef()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"

      integer lineno, linept
      integer cunit, icwk, ncwk, units(maxcwk)
      character*(lentok) keys(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /parsec/ keys

      icwk = 0
      ncwk = 2

c The default primary place for searching for data
      units(1) = 5
      call gprgid(keys(1))

c The default secondary place for searching for data
      units(2) = 5
      keys(2) = 'DEFAULT'

      return
      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword show working keywords
c This routine writes the working keywords to unit 6
      subroutine fwkswk()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"

      integer lineno, linept, i
      integer cunit, icwk, ncwk, units(maxcwk)
      character*(lentok) keys(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /parsec/ keys

      integer chrlen

      do 101 i = 1, ncwk
        write(6,'(I5,1X,I5,1X,A)') i,units(i),keys(i)(1:chrlen(keys(i)))
 101    continue
      return
      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword initialization routine.
c This routine is called when a search of all the inputs is to begin.
      subroutine fwkini()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"

      integer errcod
      integer lineno, linept
      integer cunit, icwk, ncwk, units(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units

c Reset the current working keyword.
      icwk = 0
      return
      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword increment routine.
c This is called to construct a new name and set the current unit to the
c next available input unit.
c  keywrd gives a string which is used with the current working keyword to
c construct the name.
c The return value is .true. when keywrd contains a new keyword and
c .false. when we are done parsing all of the inputs.
      logical function fwkinc(keywrd,name)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"
      character*(*) keywrd, name

      logical dupunt
      integer j
      integer lineno, linept
      integer cunit, icwk, ncwk, units(maxcwk)
      character*(lentok) keys(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /parsec/ keys
      integer chrlen, prcntr

      icwk = icwk + 1
      if (icwk .gt. ncwk) then
        fwkinc = .false.
        return
        endif

c Library routines might want to give an absolute working keyword
c This is denoted by a ':' as the first character of the keywrd
      if (keywrd(1:1) .eq. ':') then
        name = keywrd(2:len(keywrd))

c       Check to see if this unit has occured previously in the units list.
c       If it has, then it need not be considered again, since an absolute
c       keyword has been given.
        do 101 icwk = icwk, ncwk
c         Initialize duplicated unit to false
          dupunt = .false.
          do 201 j = 1, icwk - 1
            if (units(j) .eq. units(icwk)) dupunt = .true. 
 201        continue
          if (.not.dupunt) then
            fwkinc = .true.
            cunit = units(icwk)
            if (prcntr('IS_ON KEYWORD').eq.1) then
              write(6,1) 'Next keyword: ',name(1:chrlen(name)),
     &                            ' on unit ',cunit
              endif
            return
            endif
 101      continue
c       If we get here then we have exhausting the list of inputs.
        fwkinc = .false.
        return
        endif

c keywrd is not an absolute path, so name must be constructed
      name = keys(icwk)
      call scat(name,':')
      call scat(name,keywrd)
      cunit = units(icwk)

      fwkinc = .true.

c Print out debugging information
      if (prcntr('IS_ON KEYWORD').eq.1) then
        write(6,1) 'Next keyword: ',name(1:chrlen(name)),
     &                      ' on unit ',cunit
   1    format(1X,3A,I2)
        endif

      return
      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword add routine.
c This adds a unit and keyword to the end of the working directory path.
      integer function fwkadd(unit,keywrd)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"
#include "error.h"
      integer unit
      character*(*) keywrd

      integer lineno, linept
      integer cunit, icwk, ncwk, units(maxcwk)
      character*(lentok) keys(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /parsec/ keys
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      if (ncwk.ge.maxcwk) then
        errcod = Elimit
        fwkadd = errcod
        emess = 'Too many inputs (fwkadd)'
        return
        endif

      errcod = EOK
      fwkadd = errcod
      emess = ' '

      ncwk = ncwk + 1
      units(ncwk) = unit
      keys(ncwk) = keywrd

      return
      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword default user supplemented input routine.
c This lets an application easily allow users to specify alternate inputs.
c This searches the current inputs for the MORE_INPUT keyword.  If found
c the cwk are suplemented by the given unit/keyword pairs.  This
c routine works quietly, if an error occurs nothing happens.
c The MORE_INPUT keyword from PROGID is first used and then DEFAULT is
c used, even if MORE_INPUT for PROGID is found.
c Error recovery: returns a value from error.h
      integer function fwkdui()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"
#include "error.h"

      character*(lentok) progid,keywrd
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess
      character*11 morein

      integer fwkusi

      data morein/':MORE_INPUT'/

c First start with more_input input for the programid
      call gprgid(progid)

      keywrd = ':'
      call scat(keywrd,progid)
      call scat(keywrd,morein)

c This goes and gets the new inputs. EOD errors will be ignored
      fwkdui = fwkusi(keywrd)
      if (fwkdui.eq.EEOD) then
        errcod = EOK
        fwkdui = errcod
        emess = ' '
        endif
      if (fwkdui.ne.EOK) return

      keywrd = ':DEFAULT:MORE_INPUT'

c This goes and gets the new inputs. EOD errors will be ignored
      fwkdui = fwkusi(keywrd)
      if (fwkdui.eq.EEOD) then
        errcod = EOK
        fwkdui = errcod
        emess = ' '
        endif
      if (fwkdui.ne.EOK) return

      return
      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword user supplemented input routine.
c This take a keyword as an argument to let the application specifically
c give what input to use as extra inputs.
c  INPUT:
c keywrd = ( ("basis.data" lib basis) ("myown.basis" user basis) )
c first param is file
c second param is type
c   LIB means open up the input in the standard library directory.
c   USER means use the file name as is
c   SUSER means that a symbolic name is to be used.
c third param is keyword
c Error recovery: returns a value from error.h
      integer function fwkusi(keywrd)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"
#include "error.h"
      character*(*) keywrd

      integer unit, pos
      character*(lentok) token,fname,type
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdpos, frdlpn, frdchr, fioopn, fwkadd

      fwkusi = frdpos(keywrd)
      if (fwkusi.ne.EOK) return

      fwkusi = frdlpn()
      if (fwkusi.ne.EOK) return

 801  continue

        fwkusi = frdchr(token)
        if (fwkusi.eq.Ebvect) then
          pos = 0
c Default type
          type = 'USER'
c Default keyword
          call gprgid(keywrd)
          unit = 0
 802      continue
            fwkusi = frdchr(token)
            if (fwkusi.eq.Eevect) then
              if (pos.eq.0) goto 801
              fwkusi = fioopn(unit,fname,type,'READ')
              if (fwkusi.ne.EOK) return
              fwkusi = fwkadd(unit,keywrd)
              if (fwkusi.ne.EOK) return
              goto 801
              endif
            if (pos .eq. 0) then
              fname = token
            else if (pos .eq. 1) then
              type = token
            else if (pos .eq. 2) then
              keywrd = token
            else
              errcod = Esyntx
              fwkusi = errcod
              emess = 'Garbage after input specificiation: '''
              call scat(emess,token)
              call scat(emess,''' (fwkusi)')
              return
              endif
            pos = pos + 1
            goto 802
        else if (fwkusi.eq.Eevect) then
          errcod = EOK
          fwkusi = errcod
          emess = ' '
          return
          endif
      goto 801

      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword clear routine.
c This deletes all of the inputs.
      subroutine fwkclr()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"

      integer lineno, linept
      integer cunit, icwk, ncwk, units(maxcwk)
      character*(lentok) keys(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /parsec/ keys

      ncwk = 0

      return
      end
c---------------------------------------------------------------------------
c  The formatted i/o working keyword delete routine.
c This deletes one of the inputs.  If that input is unknown, then
c the list of current working keywords is unchanged.
      subroutine fwkdel(unit,keywrd)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"
      integer unit
      character*(*) keywrd

      integer i, delete
      integer lineno, linept
      integer cunit, icwk, ncwk, units(maxcwk)
      character*(lentok) keys(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /parsec/ keys

      delete = 0

      do 101 i = 1, ncwk
        if (delete .ne. 0) then
          units(i-delete) = units(i)
          keys(i-delete) = keys(i)
          endif
        if ((units(i).eq.unit).and.(keys(i).eq.keywrd)) then
          delete = delete + 1
          endif
  101   continue

      return
      end
