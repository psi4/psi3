c Get a token
      integer function gtoken(token)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include "error.h"
#include "params.h"
      character*(lentok) token

      logical quote
      integer errcod
      integer toknpt

      integer c1int,junk
      integer lineno,linept
      character*(lentok) line
      character*80 emess
      integer cunit, icwk, ncwk, units(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /curlin/ line
      common /lparer/ errcod, emess

      integer advlin
      logical tokchr, upprcs, lowrcs, lglch
      character*1 lctouc

      if (lineno .eq. 0) then
        errcod = advlin()
        if (errcod .ne. EOK) then
          gtoken = errcod
          return
          endif
        endif

c Default return code:
      gtoken = EOK
      errcod = EOK
      emess = 'The parser error message buffer is empty (gtoken)'

      quote = .false.
      toknpt = 1

      token = ' '

 801  continue

c     write(6,*) 'gtoken: top:',token

      if (linept .eq. 81) then
        if (quote) then
          errcod = Esyntx
          gtoken = errcod
          emess = 'Failure to provide closing ''"'' before end of 80'
          call scatb(emess,'character record')
          return
        else if (toknpt .ne. 1) then
          goto 1001
          endif
        errcod = advlin()
        if (errcod .ne. EOK) then
          gtoken = errcod
          return
          endif
        endif

      if (quote .and. (line(linept:linept) .eq. '"')) then
        linept = linept + 1
        goto 1001
      else if (.not. quote .and. (line(linept:linept) .eq. '"')) then
        linept = linept + 1
        quote = .true.
      else if (quote) then
        token(toknpt:toknpt) = line(linept:linept)
        linept = linept + 1
        toknpt = toknpt + 1
      else if (line(linept:linept) .eq. ' ') then
        if (toknpt .ne. 1) then
          linept = linept + 1
          goto 1001
          endif
        linept = linept + 1
c if a + or - is immediately followed by a digit or a '.', then
c it is treated as a part of the number following it.
c otherwise it will be treated as a token further below.
      else if (     (    (line(linept:linept) .eq. '+')
     &               .or.(line(linept:linept) .eq. '-'))
     &         .and.(linept.lt.80)) then
        if (    (c1int(line(linept+1:linept+1),junk).eq.0)
     &      .or.(line(linept+1:linept+1).eq.'.')) then
          token(toknpt:toknpt) = line(linept:linept)
          toknpt = toknpt + 1
          linept = linept + 1
          endif
      else if (tokchr(line(linept:linept))) then
        if (toknpt .eq. 1) then
          token(1:1) = line(linept:linept)
          linept = linept + 1
          goto 1001
        else
          goto 1001
          endif
clj   The comment character
      else if (line(linept:linept) .eq. '%') then
        errcod = advlin()
        if (errcod .ne. EOK) then
          gtoken = errcod
          return
          endif
      else if (upprcs(line(linept:linept))) then
        token(toknpt:toknpt) = line(linept:linept)
        toknpt = toknpt + 1
        linept = linept + 1
      else if (lowrcs(line(linept:linept))) then
        token(toknpt:toknpt) = lctouc(line(linept:linept))
        toknpt = toknpt + 1
        linept = linept + 1
      else if (lglch(line(linept:linept))) then
        token(toknpt:toknpt) = line(linept:linept)
        toknpt = toknpt + 1
        linept = linept + 1
      else
        gtoken = Esyntx
        errcod = Esyntx
        emess = 'Illegal character in token, current token:'
        call scatb(emess,token)
        call scatb(emess,', current char:')
        call scatb(emess,line(linept:linept))
        call scatb(emess,'(gtoken)')
        return
        endif

      goto 801

c This is where we branch for a successful find of a token
1001  continue
c     write(6,*) 'GTOKEN: ',token
      return
      end
c--------------------------------------------------------------
c Returns true if the character is a standalone token
      logical function tokchr(char)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      character*1 char
      integer i
      character*8 toks
      data toks/'=,()"->+'/
      tokchr = .false.
      do 101 i = 1, 8
        if (char .eq. toks(i:i)) then
          tokchr = .true.
          return
          endif
 101    continue
      return
      end
c--------------------------------------------------------------
c Returns true if the character is a nonletter legal character
      logical function lglch(char)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      character*1 char
      integer i
      character*26 lglchs
      data lglchs/'0123456789._*:$           '/
      lglch = .false.
      do 101 i = 1, 26
        if (char .eq. lglchs(i:i)) then
          lglch = .true.
          return
          endif
        if (lglchs(i:i).eq.' ') return
 101    continue
      return
      end
c--------------------------------------------------------------
c Returns true if the character is an uppercase letter
      logical function upprcs(char)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      character*(1) char
      integer i
      character*26 uplet
      data uplet/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/
      do 101 i = 1, 26
        if (char .eq. uplet(i:i)) then
          upprcs = .true.
          return
          endif
 101    continue
      upprcs = .false.
      return
      end
c--------------------------------------------------------------
c Returns true if the character is an uppercase letter
      logical function lowrcs(char)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      character*(1) char
      integer i
      character*26 lolet
      data lolet/'abcdefghijklmnopqrstuvwxyz'/
      do 101 i = 1, 26
        if (char .eq. lolet(i:i)) then
          lowrcs = .true.
          return
          endif
 101    continue
      lowrcs = .false.
      return
      end
c--------------------------------------------------------------
c Read in the next line of input (advance line)
      integer function advlin()
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include "error.h"
#include "params.h"

      integer lineno,linept
      integer errcod
      character*(lentok) line
      character*80 emess
      integer cunit, icwk, ncwk, units(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /curlin/ line
      common /lparer/ errcod, emess

      lineno = lineno + 1
      linept = 1
      read(cunit,'(A80)',end=901) line
clj     If we get a '#' then we are probably seeing an old style
clj     input file -- produce an end of file
clj     NOTE: the comment character has been changed from # to %
clj      to make the merger of the new style and old style inputs
clj      more transparent.
      if (line(1:1) .eq. '#') goto 901
      advlin = EOK
      errcod = EOK
      emess = 'The parser error message buffer is empty (advlin)'
      return

 901  continue
      advlin = EEOF
      errcod = EEOF
      emess = 'Reached the end of the input file (advlin)'
      line = 'End of File'
      return
      end
c----------------------------------------------------------------
c Push or pop info about where we are in the input into a stack.
c op can be:
c  PUSH
c  POP
c  DUPLICATE
c  PURGE
c where the meaning is as you would expect
      integer function fstck(op)
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include "error.h"
#include "params.h"
      character*(*) op

      integer lineno,linept
      integer errcod
      character*(lentok) line
      character*80 emess
      integer cunit, icwk, ncwk, units(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /curlin/ line
      common /lparer/ errcod, emess

      integer frstr

c     The maximum depth
      integer mxdpth
      parameter (mxdpth=10)
c     The stacks
      integer sicwk(mxdpth)
      integer slno(mxdpth)
      integer slpt(mxdpth)
crawdad
	save sicwk, slno, slpt
c     The depth of the stack
      integer dpth
      save dpth
      data dpth/0/

      if (op.eq.'PUSH') then
        dpth = dpth+1
        if (dpth.gt.mxdpth) goto 901
        sicwk(dpth) = icwk
        slno(dpth) = lineno
        slpt(dpth) = linept
      else if (op .eq. 'DUPLICATE') then
        if (dpth.eq.0) goto 902
        dpth = dpth+1
        if (dpth.gt.mxdpth) goto 901
        sicwk(dpth) = sicwk(dpth-1)
        slno(dpth) = slno(dpth-1)
        slpt(dpth) = slno(dpth-1)
      else if (op .eq. 'POP') then
        if (dpth.eq.0) goto 902
        icwk = sicwk(dpth)
        lineno = slno(dpth)
        linept = slpt(dpth)
        dpth = dpth-1
        fstck = frstr()
        if (fstck.ne.EOK) return
      else if (op .eq. 'PURGE') then
        dpth = 0
      else
        errcod = Earg
        fstck = errcod
        emess = 'Bad stack operation (fstck)'
        endif

c    NORMAL EXIT
      errcod = EOK
      fstck = errcod
      return

c    OVERFLOW ERROR
 901  continue
      errcod = Elimit
      fstck = errcod
      emess = 'Stack overflow (fstck)'
      return

c    UNDERFLOW ERROR
 902  continue
      errcod = Elimit
      fstck = errcod
      emess = 'Stack underflow (fstck)'
      return
      end
c----------------------------------------------------------------
c This restores the i/o state to be consistent with the values in
c the icwk, lineno, and linept common variables.
      integer function frstr()
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include "error.h"
#include "params.h"

      integer lineno,linept
      integer errcod
      character*(lentok) line
      character*80 emess
      integer cunit, icwk, ncwk, units(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /curlin/ line
      common /lparer/ errcod, emess

      integer i, savept, saveno

      integer advlin

      savept = linept
      saveno = lineno
      cunit = units(icwk)
      call frew()

      do 101 i = 1, saveno
        frstr = advlin()
        if (frstr.ne.EOK) then
          write(6,*) 'WARNING: UNEXPECTED ERROR IN FRSTR WHILE'
          write(6,*) 'RESTORING FORMATTED I/O STATE'
          return
          endif
 101    continue

      linept = savept
      return
      end
c----------------------------------------------------------------
c Rewind a file which will be parsed
      subroutine frew()
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include "params.h"

      integer lineno, linept
      integer cunit, icwk, ncwk, units(maxcwk)
      common /parse/ lineno, linept, cunit, icwk, ncwk, units

      rewind cunit
      lineno = 0

      return
      end
c----------------------------------------------------------------
c After a '(' is read, calling this will make the next token the token
c after the matching ')'. An error code is returned.
      integer function matchp()
#ifdef UNDEFINED
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
#include "error.h"
#include "params.h"

      integer lineno,linept
      character*(lentok) line
      character*80 emess
      integer cunit, icwk, ncwk, units(maxcwk)
      integer errcod
      common /parse/ lineno, linept, cunit, icwk, ncwk, units
      common /curlin/ line
      common /lparer/ errcod, emess

      integer i,j,level,start,advlin

      level = 1
 801  continue
        start = linept
        do 101 i = start, 80
c         Specify the action for the single most common character first
          if (line(i:i).eq.' ') then
            goto 101
          else if (line(i:i).eq.'%') then
            goto 901
          else if (line(i:i).eq.'#') then
            goto 902
          else if (line(i:i).eq.'(') then
            level = level + 1
          else if (line(i:i).eq.')') then
            level = level - 1
            if (level .eq. 0) then
              linept = i + 1
              errcod = EOK
              matchp = errcod
              emess = 'The error message buffer is empty (matchp)'
              return
              endif
          else if (line(i:i).eq.'"') then
            do 201 j = i+1, 80
              if (line(j:j).eq.'"') then
                linept = j+1
                goto 801
                endif
 201          continue
            linept = 81
            errcod = Esyntx
            emess = 'Missing a closing quote before line 81 (matchp)'
            matchp = errcod
            return
            endif

 101      continue

c Advance the line
 901  continue
      errcod = advlin()
      if (errcod.eq.EEOF) then
        goto 902
      else if (errcod.ne.EOK) then
        matchp = errcod
        return
        endif
      linept = 1
      goto 801

c Got to the end of the input
 902  continue
      errcod = EMRPAR
      matchp = errcod
      emess = 'Missing a right parenthesis (matchp)'
      return
      end
