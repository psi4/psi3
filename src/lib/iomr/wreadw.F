C Deck WReadW
      SUBROUTINE WREADW(file,array,nwords,fword1,nxtwrd)
c
c  NB. The first word of the file is considered to be "1" by R
c
      implicit integer(a-z)
c
      parameter (RECLEN=1024,numbuf=10)
c
      common /io/     inp,iout
      common /pointr/ wptr(128),tptr(128)
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &  untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
c
      dimension array(nwords)
c
c      write(6,1000) file,nwords,fword
 1000 format(' wread:    file',i3,', ',i12,' words,    first word',i8)
c
      fword=fword1
c
c     ---- determine which records need to be read ----
c
      frec=(fword-1)/reclen + 1
      lword=fword+nwords-1
      lrec=(lword-1)/reclen + 1
c
c      ---- loop over the records ----
c
      do 10 irec=frec,lrec
c
c       If this file does not have a buffer then grab one.
        if (untbuf(file) .eq. 0) then
c         Deallocate buffer from previous file.
          if (bufunt(nxtbuf) .ne. 0) then
c           if buffer has been written to then flush it
            call bioflu(bufunt(nxtbuf),'WREADW')
cc          if (wrtbuf(nxtbuf) .ne. 0) then
cc            write(bufunt(nxtbuf),rec=bufrec(nxtbuf),iostat=ierr)
cc   &           buffer(1,nxtbuf)
cc            if (ierr .ne. 0) then
cc              write(6,*) 'wreadw: file ',file,' rec ',bufrec(nxtbuf),
cc   &            ' ierr ', ierr
cc              call ioerr(ierr)
cc              call mabort
cc              endif
cc            wrtbuf(nxtbuf) = 0
cc            endif
            untbuf(bufunt(nxtbuf)) = 0
            endif
          untbuf(file) = nxtbuf
          bufunt(nxtbuf) = file
          bufrec(nxtbuf) = 0
          nxtbuf = nxtbuf + 1
          if (nxtbuf .gt. numbuf) nxtbuf = 1
          endif
 
c        If we don't have irec yet, then read it in.
         if (bufrec(untbuf(file)) .ne. irec) then
           call bioflu(file,'WREADW')
           read(file,rec=irec,iostat=ierr)
     &       (buffer(i,untbuf(file)),i=1,reclen)
c          write(6,*) 'wreadw: reading'
c          do 108 i = 1, 20
c            write(6,*) '   ',buffer(i,untbuf(file))
c108         continue
           if(ierr.ne.0) then
             write(6,1000) file,nwords,fword
             write(6,*) 'frec,lrec,irec=',frec,lrec,irec
             call ioerr(ierr)
             call mabort
             endif
           bufrec(untbuf(file)) = irec
           endif
c
c
c           -- transfer the appropriate portion --
c           --      of buffer to array          --
c
            if(irec.eq.frec) then
c
               offset=fword-(irec-1)*reclen - 1
               count=reclen-offset
               if(count.gt.nwords) count=nwords
               do 20 i=1,count
                  array(i)=buffer(i+offset,untbuf(file))
   20          continue
c
            else if(irec.eq.lrec) then
c
               count=lword-(irec-1)*reclen
               offset=nwords-count
               do 30 i=1,count
                  array(i+offset)=buffer(i,untbuf(file))
   30          continue
c
            else
c
              offset=reclen*(irec-1)+1-fword
              do 40 i=1,reclen
                 array(i+offset)=buffer(i,untbuf(file))
   40         continue
c
            end if
c
c
   10 continue
c
c      ---- update pointer ----
c
c
c     wptr is the next word after the last word transferred
      wptr(file)=lword+1
      nxtwrd=wptr(file)
c
      return
      END
