c Position the unit to the keywrd
      integer function frdpos(keywrd)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) keywrd

      integer maxlvl
      parameter (maxlvl = 10)
      logical match,skip
      integer i,rdepth
      integer iseg, nseg, ntseg
      integer ilvl, lvltos(0:maxlvl)
      character*(lentok) name, seg(maxseg), tseg(maxseg), token
      character*(lentok) stoken
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      logical remtch
      integer askip, gtoken, segspl, segsp2
      logical fwkinc

c Initialize the working keyword and the current unit number. 
      call fwkini()

c Return to this position if the working keyword/working unit pair needs
c to be incremented.
 803  continue

c If there is no more data to be read, then return EEOD (End of Data)
c This is less severe than an EEOF which would indicated that an end of
c file occurred in an unexpected place, such as when an '=' is needed.
      if (.not. fwkinc(keywrd,name)) then
        errcod = EEOD
        frdpos = errcod
        emess = 'There is no more input to parse (frdpos)'
        return
        endif

c  This is where we restart if redirection for the desired for token is
c requesting in the input
      rdepth = 0
 802  continue

c ............... Break the name into segments (delimited by ':') ............
      frdpos = segspl(name,seg,nseg)
      if (frdpos .ne. EOK) return

c     write(6,*) 'frdpos: name= ', name
c     write(6,*) ' segments:'
c     write(6,*) (seg(i),i=1,nseg)

c ................ Begin parsing the file ..................
c Rewind the input file
      call frew()

c Set the current segment to the first segment
      iseg = 1
c The level of parentheses
      ilvl = 0
c When we go up a level of parenthesis, we must know what segment that
c we must then start comparing, lvltos (level to segment) gives this info.
      lvltos(ilvl) = iseg
c Match is false when we need a keyword segment (or ')') and not a ':'
      match = .false.

  801 continue

c Get the next token
        errcod = gtoken(token)
c An end of file error here is OK. It means that there is no more data
c to be read.  When this occurs we must check to see if there are any
c more working keyword/working unit pairs.
        if (errcod .eq. EEOF) then
          goto 803
        else if (errcod .ne. EOK) then
          frdpos = errcod
          return
          endif

c If token is of the form x:y:z then break it up into x, :, y, :, z and
c feed these pieces to the following code one at a time.
        frdpos = segsp2(token,tseg,ntseg)
        if (frdpos .ne. EOK) return
c       write(6,*) 'frdpos: token= ', token
c       write(6,*) ' segments:'
c       write(6,*) (tseg(i),i=1,ntseg)

        skip = .false.
        do 101 i = 1, ntseg
c If we skipped anything then we immediately terminitate this do loop
c and go back to gtoken
        if (skip) goto 801

        stoken = tseg(i)

        if (match) then
          if (stoken .eq. ':') then
            match = .false.
          else if (stoken .eq. '(') then
            errcod = Esyntx
            frdpos = errcod
            emess = 'Parser does not support function calls (frdpos)'
            return
          else
            errcod = Esyntx
            frdpos = errcod
            emess = 'Expected '':'', but found '''
            call scat(emess,stoken)
            call scat(emess,''' (frdpos)')
            endif
        else
          if (remtch(stoken,seg(iseg))) then
            match = .true.
            if (iseg .eq. nseg) then
c             If not all of the tokens in ntseg have processed then
c             this is a syntax error
              if (i .ne. ntseg) then
                errcod = Esyntx
                frdpos = errcod
                emess = 'Matched the initial part of a longer keyword'
                call scatb(emess,'for ''')
                call scat (emess,name)
                call scat (emess,''' (frdpos)')
                return
                endif
c             We have found the keyword, the next token must be '=' or '-'
              frdpos = gtoken(token)
              if (frdpos .ne. EOK) return
              if (token .eq. '-') then
c               The next token must be '>'
                frdpos = gtoken(token)
                if (frdpos .ne. EOK) return
                if (token .ne. '>') then
                  errcod = Esyntx
                  frdpos = errcod
                  emess = 'Missing ''>'' after the ''-'' for'
                  call scatb(emess,name)
                  call scatb(emess,'(frdpos)')
                  return
                  endif
c               This is a redirection. Use the following token as
c               the new name and start frdpos from the beginning.
                frdpos = gtoken(token)
                if (frdpos .ne. EOK) return
                name = token
                rdepth = rdepth + 1
c  The redirection depth is hardwired here because this is the only place
c it should be used. If we go below this depth, return with an error.
                if (rdepth .gt. 10) then
                  errcod = Erdrct
                  frdpos = errcod
                  emess = 'Redirection depth is too deep for '''
                  call scat(emess, name)
                  call scat(emess, ''' (frdpos)')
                  return
                  endif
                goto 802
              else if (token .ne. '=') then
                errcod = Esyntx
                frdpos = errcod
                emess = 'Missing ''='' for'
                call scatb(emess,name)
                call scatb(emess,'(frdpos)')
                return
              else
                errcod = EOK
                frdpos = errcod
                emess = ' '
                return
                endif
              endif
            iseg = iseg + 1
          else if (stoken .eq. '(') then
            ilvl = ilvl + 1
            if (ilvl .gt. maxlvl) then
              errcod = Eparen
              frdpos = errcod
              emess = 'Too many levels of parenthesis for'
              call scatb(emess,name)
              call scatb(emess,'(frdpos)')
              return
              endif
            lvltos(ilvl) = iseg
          else if (stoken .eq. ')') then
            ilvl = ilvl - 1
            if (ilvl.lt.0) then
              errcod = Esyntx
              frdpos = errcod
              emess = 'Too many '')'''
              return
              endif
c           If we pop up a level we must leave the current segment
            iseg = lvltos(ilvl)
          else 
            skip = .true.
            errcod = askip()
c Must reset the level after skipping
            iseg = lvltos(ilvl)
            if (errcod .ne. EOK) then
              frdpos = errcod
              return
              endif
            endif
          endif

101     continue

        goto 801
      end
c-------------------------------------------------------------------------
c The assignment skip routine.  This skips over keywords and their values.
      integer function askip()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"

      character*(lentok) token
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer gtoken, vskip, matchp

  801 continue

        askip = gtoken(token)
        if (askip .ne. EOK) return

        if (token .eq. ':') then
          goto 801
        else if (token .eq. '(') then
          askip = matchp()
          if (askip .ne. EOK) return
        else if (token .eq. '-') then
          askip = gtoken(token)
          if (askip .ne. EOK) return
          if (token .ne. '>') then
            errcod = Esyntx
            askip = errcod
            emess = 'A ''-'' was not followed by a ''>'' (askip)'
            return
            endif
          askip = vskip()
          if (askip .ne. EOK) return
        else if (token .eq. '=') then
          askip = vskip()
          if (askip .ne. EOK) return
          endif

      return
      end
c-------------------------------------------------------------------------
c The value skip routine.  This skips over values.
      integer function vskip()
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"

      character*(lentok) token
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer gtoken, matchp

      errcod = gtoken(token)
      if (errcod .ne. EOK) then
        vskip = errcod
        return
        endif

      if (token .eq. '(') then
        errcod = matchp()
        if (errcod .ne. EOK) then
          vskip = errcod
          return
          endif
        endif

      errcod = EOK
      emess = ' '
      vskip = errcod
      return
      end
c------------------------------------------------------------------------
c The array position routine.  After frdpos has been used to position
c the input to a value, this routine can be used to position the input
c to the desired array element (which may be an array).
      integer function arpos(nidx,idx)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer  nidx
      integer idx(nidx)

      integer i
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer arskip

c Position to the starting array location.
c If nidx == 0 then we have a scalar and the input stays where it is.
      do 102 i = 1, nidx
        errcod = arskip(idx(i))
        if (errcod .ne. EOK) then
          arpos = errcod
          return
          endif
  102   continue

      errcod = EOK
      arpos = errcod
      emess = ' '

      return
      end
c------------------------------------------------------------------------
c The array skip routine. This skips over nele - 1 elements of an array.
      integer function arskip(nele)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      integer  nele

      integer i
      character*(lentok) token
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      integer matchp, gtoken

      errcod = gtoken(token)
      if (errcod .ne. EOK) then
        arskip = errcod
        return
        endif

c The next token must be a '('

      if (token .ne. '(') then
        errcod = Esyntx
        arskip = errcod
        emess = 'Expected an array but did not find a ''('' (arskip)'
        return
        endif

      do 101 i = 1, nele - 1

        errcod = gtoken(token)
        if (errcod .ne. EOK) then
          arskip = errcod
          return
          endif

        if (token .eq. ')') then
c         The next token cannot be a ')'
          errcod = Elngth
          arskip = errcod
          emess = 'An array did not have enough elements (arskip)'
          return
        else if (token .eq. '(') then
c         If we are starting a new array then skip the whole thing.
          errcod = matchp()
          if (errcod .ne. EOK) then
            arskip = errcod
            return
            endif
          endif
  101   continue

      errcod = EOK
      arskip = errcod
      emess = ' '

      return
      end
c-----------------------------------------------------------------------
c The segment split routine.  This splits name into its segments and
c puts the number of segments into nseg and returns and error code as
c its value.  The elements of the segment array are the strings delimited
c by ':' in name.
      integer function segspl(name,seg,nseg)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) name, seg(maxseg)
      integer nseg

      integer i
      integer segpos(maxseg)
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

c First find the segment delimiters
      nseg = 1
      segpos(nseg) = 0
      do 101 i = 1, lentok
        if (name(i:i) .eq. ':') then
          nseg = nseg + 1
          if (nseg+1 .gt. maxseg) then
            errcod = Eseg
            segspl = errcod
            emess = 'Too many segments (segspl)'
            return
            endif
          segpos(nseg) = i
          endif
  101   continue
      segpos(nseg + 1) = lentok + 1

c Now get the segments
      do 102 i = 1, nseg
        seg(i) = name(segpos(i)+1:segpos(i+1)-1)
  102   continue
      do 103 i = nseg+1, maxseg
        seg(i) = ' '
  103   continue

      errcod = EOK
      segspl = errcod
      emess = ' '

      return
      end
c-----------------------------------------------------------------------
c The segment split routine version 2.  This is very similar to
c segspl, but this version retains the ':' and places them in the
c appropiate locations in the seg array.
      integer function segsp2(name,seg,nseg)
#ifdef sun
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) name, seg(maxseg)
      integer nseg

      integer i, start
      integer errcod
      character*80 emess
      common /lparer/ errcod, emess

      nseg = 0
      start = 1
      do 101 i = 1, len(name)
        if (name(i:i) .eq. ':') then
          if (start .lt. i) then
            nseg = nseg + 1
            if (nseg .gt. maxseg) goto 901
            seg(nseg) = name(start:i-1)
            endif
          nseg = nseg + 1
          if (nseg .gt. maxseg) goto 901
          seg(nseg) = ':'
          start = i + 1
          endif
  101   continue

      if (start .le. len(name)) then
        if (name(start:len(name)).ne.' ') then
          nseg = nseg + 1
          if (nseg .gt. maxseg) goto 901
          seg(nseg) = name(start:len(name))
          endif
        endif

      errcod = EOK
      segsp2 = errcod
      emess = ' '
      return

  901 continue
      errcod = Eseg
      segsp2 = errcod
      emess = 'Too many segments in a name (segsp2)'
      return
      end
