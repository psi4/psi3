*Deck WWritW
      SUBROUTINE WWRITW(file,array,nwords,fword1,nxtwrd)
c
c   NB.  The file is considered to start at word "1" in R
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
c     write(6,1000) file,nwords,fword
 1000 format(' wwrit:    file',i3,i12,' words   first word:',i8)
c
      fword=fword1
c
c     ---- determine which records need to be written ----
c
      frec=(fword-1)/reclen + 1
      lword=fword+nwords-1
      lrec=(lword-1)/reclen + 1
c
c      ---- loop over the records ----
c
      do 10 irec=frec,lrec
c
c       If file does not have a buffer, then allocate one.
        if (untbuf(file) .eq. 0) then
c         Deallocate buffer from previous file
          if (bufunt(nxtbuf) .ne. 0) then
c           Flush buffer if written to
            call bioflu(bufunt(nxtbuf),'WWRITW(PURGE)')
cc          if (wrtbuf(nxtbuf) .ne. 0) then
cc            write(bufunt(nxtbuf),rec=bufrec(nxtbuf),iostat=ierr)
cc   &              buffer(1,nxtbuf)
cc            if (ierr .ne. 0) then
cc              write(6,*) 'wwritw: file ',file,' rec ',bufrec(nxtbuf),
cc   &                     ' error ',ierr
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
c
c If we don't have irec yet, then read it in.
c       write(6,*) 'wwritw: file, irec ', file, irec
c       write(6,*) '        untbuf(file)', untbuf(file)
c       write(6,*) '        bufrec(untbuf(file))', bufrec(untbuf(file))
        if (bufrec(untbuf(file)) .ne. irec) then
          call bioflu(file,'WWRITW(NEW RECORD)')
          read(file,rec=irec,err=300,iostat=ierr)
     &        (buffer(i,untbuf(file)),i=1,reclen)
          if (ierr .ne. 0) then
            write(iout,8893) ierr
 8893       format(//' error code for read of buffer',i10)
            write(iout,*) 'file =',file,'irec =',irec
            call ioerr(ierr)
            call mabort
            endif
          bufrec(untbuf(file)) = irec
          endif
          go to 310
  300     continue
          bufrec(untbuf(file)) = irec
c         write(6,1300) irec,file,ierr
 1300     format('wwritw: first access of record',i10,' file ',i2,
     &            ' err ',i4)
  310     continue
c
c           -- update buffer --
c
            if(irec.eq.frec) then
c
               offset=fword-(irec-1)*reclen - 1
               count=reclen-offset
               if(count.gt.nwords) count=nwords
               do 20 i=1,count
                  buffer(i+offset,untbuf(file))=array(i)
   20          continue
c
            else if(irec.eq.lrec) then
c
               count=lword-(irec-1)*reclen
               offset=nwords-count
               do 30 i=1,count
                  buffer(i,untbuf(file))=array(i+offset)
   30          continue
c
            else
c
              offset=reclen*(irec-1)+1-fword
              do 40 i=1,reclen
                 buffer(i,untbuf(file))=array(i+offset)
   40         continue
c
            end if
c
c      ---- mark the buffer as written to ----
c
c889        continue
            wrtbuf(untbuf(file)) = 1
c
   10 continue
c
c      ---- update pointers ----
c
c     tptr is the highest record written
      if(lrec.gt.tptr(file)) tptr(file)=lrec
c
c     wptr is the next word after the last word transferred
      wptr(file)=lword+1
      nxtwrd=wptr(file)
clj022489
c     write(6,8893) lword,wptr(file),nxtwrd
c8893 format('wwritw: unit',i3,' lword',i10,' wptr',i10,' nxtwrd',i10)
c
c
      return
      end
