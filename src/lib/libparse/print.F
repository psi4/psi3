c Print option control routines.
c The names of print options are truncated to 8 characters.
c-----------------------------------------------------------------------------
c  Printing utility control routine
c This routine processes all print option functions
c  req is:
c    IS_ON followed by a print option names
c    ON followed one or more print option names
c    OFF followed by one or more print option names
c Error recovery: returns a value from error.h
c if the IS_ON request is given, a 0 or 1 is returned for true or false
      integer function prcntr(req)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) req

      integer i
      integer maxprt
      parameter (maxprt = 100)
      character*80 next
      character*8 word1, opt
      save names
      save nnames
      integer nnames,ndelet
      character*8 names(maxprt)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer chrlen

      data nnames/0/

c Default return action is EOK
      errcod = EOK
      prcntr = errcod
      emess = ' '

      word1 = req(1:chrlen(req))
      next  = req(chrlen(req)+1:len(req))
      call leftju(next)

      if ((len(req).gt.80).and.(word1.ne.'IS_ON')) then
        errcod = Elimit
        prcntr = errcod
        emess = 'Request string is too long (prcntr)'
        return
        endif

      if (word1.eq.'IS_ON') then
        do 101 i = 1, nnames
          if (names(i).eq.next) then
            prcntr = 1
            return
            endif
  101     continue
        prcntr = 0
        return
      else if ((word1.eq.'ON').or.(word1.eq.'OFF')) then

c       Loop thru all of the names in the next list
 801    continue
          opt = next(1:chrlen(next))
          next  = next(chrlen(next)+1:len(next))
          call leftju(next)
          if (opt.eq.' ') return
          if (word1.eq.'ON') then
c           If an option is given twice, ignore it.
            do 102 i = 1, nnames
              if (names(i).eq.opt) return
 102          continue
c           Bounds checking
            if (nnames.ge.maxprt) then
              errcod = Elimit
              prcntr = errcod
              emess = 'Too many print option names (prcntr)'
              return
              endif
            nnames = nnames + 1
            names(nnames) = opt
          else
c           This will delete all occurrances of the opt
            ndelet = 0
            do 103 i = 1, nnames
              names(i-ndelet) = names(i)
              if (names(i).eq.opt) ndelet = ndelet + 1
  103         continue
            nnames = nnames - ndelet
            endif
          goto 801
      else
        errcod = Earg
        prcntr = errcod
        emess = 'Bad request, '''
        call scat(emess,word1)
        call scat(emess,''' (prcntr)')
        return
        endif
c Execution cannot reach here
      end
c---------------------------------------------------------------------------
c  The Printing utility user print routine
c This finds the array of character specified by the keywrd argument.
c The form of the vector is as follows:
c  ((ON NAME_1 NAME_2 ...)
c   (OFF NAME_1 NAME_3 ... ))
c Printing is turned off for names following OFF.
c Printing is turned on for names following ON.
      integer function prup(keywrd)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "params.h"
#include "error.h"
      character*(*) keywrd

      character*(lentok) opt,val,req
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer frdpos, frdlpn, frdchr, prcntr

      prup = frdpos(keywrd)
      if (prup.ne.EOK) return

      prup = frdlpn()
      if (prup.ne.EOK) return

 801  continue
        prup = frdchr(val)
        if (prup.eq.Ebvect) then
          prup = frdchr(opt)
          if (prup.eq.Eevect) then
            goto 801
          else if (prup.ne.EOK) then
            return
            endif
 802      continue
            prup = frdchr(val)
            if (prup.ne.Eevect) then
              if (prup.ne.EOK) return
              req = opt
              call scatb(req,val)
              prup = prcntr(req)
              if (prup.ne.EOK) return
              goto 802
              endif
          goto 801
        else if (prup.eq.Eevect) then
          errcod = EOK
          prup = errcod
          emess = ' '
        else if (prup.eq.EOK) then
          errcod = Esyntx
          prup = errcod
          emess = 'Looking for a '')'' or ''('' but found a'''
          call scat(emess, val)
          call scat(emess, ''' (prup)')
          endif
      return
      end
c---------------------------------------------------------------------------
c  The Printing utility default user print routine
c This lets an application easily allow users to specify print options.
c This searches the current inputs for the PRINT_ON and PRINT_OFF keywords.
c   First the PRINT keyword from DEFAULT is examined.
c   Third the PRINT keyword from PROGID is examined.
      integer function prdup()
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

      integer prup

c Get the input from the defaults section
      keywrd = ':DEFAULT'
      call scat(keywrd,':PRINT')
      prdup = prup(keywrd)
      if (prdup.eq.EEOD) then
        errcod = EOK
        prdup = errcod
        emess = ' '
        endif
      if (prdup.ne.EOK) return

c Now get input for the programid
      call gprgid(progid)
      keywrd = ':'
      call scat(keywrd,progid)
      call scat(keywrd,':PRINT')
      prdup = prup(keywrd)
      if (prdup.eq.EEOD) then
        errcod = EOK
        prdup = errcod
        emess = ' '
        endif
      if (prdup.ne.EOK) return
      end

