      subroutine execut(exec,nexec)
#if defined(UNDEFINED)
      implicit undefined(a-z)
#else
      implicit integer (a-z)
#endif
      integer nexec
      character*(*) exec(nexec)
      integer mxexec
      parameter (mxexec = 100)
      integer i,j,pos,depth
      integer irep(mxexec),mirep(mxexec),prep(mxexec),pend(mxexec)
      character*40 cmrep
      integer errcod,chrlen,chrint,chrend

      if (nexec.gt.mxexec) then
        write(6,*) 'ERROR: too many commands to execute'
        call qabort
        endif

      do 102 i = 1, nexec
        irep(i) = 0
        mirep(i) = 0
        if (exec(i)(1:6).eq.'REPEAT') then
c Initialize all of the repeat indices to 0
          irep(i) = 0
c Compute the maximum index for each repeat we find.
          cmrep = exec(i)(8:chrlen(exec(i))-1)
          if (chrint(cmrep,mirep(i)).ne.0) then
            write(6,*) 'ERROR: invalid integer data as an argument'
     &                // ' to ''REPEAT'''
            write(6,*) '  Data:',cmrep
            call qabort
            endif
c Find the position of the END for this REPEAT (pend)
          depth = 0
          do 202 j = i, nexec
            if (exec(j)(1:6).eq.'REPEAT') then
              depth = depth + 1
            else if (exec(j).eq.'END') then
              depth = depth - 1
              if (depth.eq.0) then
c The position of the END for a REPEAT
                pend(i) = j
c The position of the REPEAT for an END
                prep(j) = i
                goto 902
                endif
              endif
  202       continue
  902     continue
          endif
  102   continue

c     write(6,*) 'pend = ',pend
c     write(6,*) 'prep = ',prep

c The position in the exec list
      pos = 1

 801  continue
        if (pos.gt.nexec) then
          write(6,*)
          write(6,*) 'Execution has completed.'
          return
          endif
        if (exec(pos)(1:6).eq.'REPEAT') then
          irep(pos) = irep(pos) + 1
          if (irep(pos).gt.mirep(pos)) then
c           go past the END statement
            pos = pend(pos) + 1
            goto 801
            endif
          write(6,*) exec(pos)(1:chrlen(exec(pos))),
     &              ' on iteration ',irep(pos)
c         Initialize the indices following this index to 0
          do 101 i = pos+1, nexec
            if (exec(i)(1:6).eq.'REPEAT') irep(i) = 0
  101       continue
        else if (exec(pos) .eq. 'END') then
          write(6,*) exec(pos)(1:chrlen(exec(pos)))
c         Goto the repeat statement
          pos = prep(pos)
          goto 801
        else
          write(6,*) 'EXECUTING: ',exec(pos)(1:chrend(exec(pos)))
c         write(6,*) 'nexec, pos = ',nexec,pos
          call runcmd(errcod,exec(pos),len(exec(pos)))
          if (errcod.ne.0) then
            write(6,*)
            write(6,*) 'COMMAND HAS RETURNED A NONZERO STATUS: ',errcod
            write(6,*) 'Execution is being terminated.'
            call qabort  
            endif
          endif

c Increment the position and go to the next statement
        pos = pos+1
        goto 801

      end
