      subroutine lexec(exec,nexec,mxexec,name)
#if defined(UNDEFINED)
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      integer nexec,mxexec,nrept
      character*(*) exec(mxexec),name
      integer rdepth
      integer errcod, frdpos, chrlen, var
      character*80 string,cnrep
      integer level
      integer i
#include <error.h>

      errcod = frdpos(name)
      if (errcod .ne. EOK) then
        write(6,*) 'ERROR: couldn''t locate ',name(1:chrlen(name)),
     &             ' in your input.'
        call pabort
        call qabort
        endif

c The number of things to execute
      nexec = 0
c The depth of parenthetical levels
      level = 0
c The depth of the REPEAT's
      rdepth = 0

 801  continue

      errcod = var('GET',string)
      if (errcod.eq.Ebvect) then
c       write(6,*) 'Ebvect'
        level = level + 1
        goto 801
      else if (errcod.eq.Eevect) then
c       write(6,*) 'Eevect'
        errcod = var('RESTORE',string)
        errcod = var('RESTORE',string)
        level = level - 1
c See if we are done
        if (level .eq. 0) return
        goto 801
      else if (errcod.ne.EOK) then
        write(6,*) 'ERROR: problem reading ',name(1:chrlen(name))
        call pabort
        call qabort
        endif

      if (string.eq.'REPEAT') then
        errcod = var('GET',cnrep)
        if (errcod .ne. EOK) then
          write(6,*) 'ERROR: couldn''t read a ''REPEAT'' argument'
          call pabort
          call qabort
          endif
        rdepth = rdepth + 1
        call scat(string,'(')
        call scat(string,cnrep)
        call scat(string,')')
        errcod = var('RESTORE',cnrep)
      else if (string.eq.'END') then
        rdepth = rdepth - 1
        if (rdepth.lt.0) then
          write(6,*) 'ERROR: too many ''END''s'
          call pabort
          call qabort
          endif
      else if (string.eq.'ERROR') then
        write(6,*) 'ERROR: found an ''ERROR'' program to execute'
        write(6,*) ' The most likely cause is an incompatible'
        write(6,*) ' combination of ''WFN'', ''DERTYPE'', and'
        write(6,*) '''REFERENCE''.'
        call qabort
        endif

      nexec = nexec + 1
      if (nexec .ge. mxexec) then
        write(6,*) 'ERROR: too many programs to execute'
        call qabort
        endif
      exec(nexec) = string
c     write(6,*) 'rdepth = ', rdepth
c     write(6,*) 'exec(',nexec,') = ', (' ',i=1,rdepth),
c    &           exec(nexec)(1:chrlen(exec(nexec)))

      errcod = var('RESTORE',string)

      goto 801
      end
c ---------------------------------
c This reads the next string, which may be of the form $xxx which
c is a variable.
      integer function var(task,string)
#if defined(UNDEFINED)
      implicit undefined(a-z)
#else
      implicit integer(a-z)
#endif
      character*(*) task,string
      integer mdepth
      parameter (mdepth=10)
      integer i,depth
      integer errcod, frdpos,frdchr,chrlen,fstck,npush(mdepth)
      save npush,depth
#include <error.h>
      data depth/0/
      if (task.eq.'GET') then
        depth = depth + 1
        npush(depth) = 0
c       write(6,*) 'GET: depth = ',depth
        if (depth.gt.mdepth) then
          write(6,*) 'ERROR: variable reference nested too deeply'
          call pabort
          call qabort
          endif
 802    continue
        errcod = frdchr(string)
        if (errcod.ne.EOK) then
          var = errcod
          return
          endif
        if (string(1:1) .eq. '$') then
          npush(depth) = npush(depth) + 1
          errcod = fstck('PUSH')
          if (errcod.ne.EOK) then
            write(6,*) 'ERROR: stack push failed'
            call pabort
            call qabort
            endif
          errcod = frdpos(string(2:chrlen(string)))
          if (errcod.ne.EOK) then
            write(6,*)'ERROR: problem finding ',string(2:chrlen(string))
            call pabort
            call qabort
            endif
          goto 802
          endif
      else if (task.eq.'RESTORE') then
c       write(6,*) 'RESTORE: depth, npush() = ',depth, npush(depth)
        if (depth.eq.0) then
          write(6,*) 'INTERNAL ERROR: did a restore with a depth of 0'
          call pabort
          call qabort
          endif
        do 101 i = 1, npush(depth)
          errcod = fstck('POP')
          if (errcod.ne.EOK) then
            write(6,*) 'ERROR: stack pop failed'
            call pabort
            call qabort
            endif
 101      continue
        npush(depth) = 0
        depth = depth - 1
      else
        write(6,*) 'INTERNAL ERROR: bad task to var'
        call qabort
        endif

      var = EOK
      return
      end
