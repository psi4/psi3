*Deck WWritW
      SUBROUTINE iofw(file,array,fbyte,nbyte)
c
c   NB.  The file is considered to start at word "1" in R
c
      implicit integer(a-z)
c
      parameter (RECLEN=1024,numbuf=10)
c
      common /io/     inp,iout
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &  untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
c
      dimension array(1)
      fword1 = fbyte / 4 + 1
      nwords = nbyte / 4
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
            call iofflu(bufunt(nxtbuf),'WWRITW(PURGE)')
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
          call iofflu(file,'WWRITW(NEW RECORD)')
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
c
clj022489
c     write(6,8893) lword,wptr(file),nxtwrd
c8893 format('wwritw: unit',i3,' lword',i10,' wptr',i10,' nxtwrd',i10)
c
c
      return
      end
C Deck WReadW
      SUBROUTINE iofr(file,array,fbyte,nbyte)
c
c  NB. The first word of the file is considered to be "1" by R
c
      implicit integer(a-z)
c
      parameter (RECLEN=1024,numbuf=10)
c
      common /io/     inp,iout
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &  untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
c
      dimension array(1)
      fword1 = fbyte / 4 + 1
      nwords = nbyte / 4
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
            call iofflu(bufunt(nxtbuf),'WREADW')
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
           call iofflu(file,'WREADW')
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
c
      return
      END
C
      SUBROUTINE iofo(ITAPE)
C
C USE THIS TO OPEN ALL DIRECT ACCESS FILES
C
      IMPLICIT INTEGER (A-Z)
C
clj040389      parameter (reclen = 5120, numbuf = 10)
      parameter (reclen = 1024, numbuf = 10)
      COMMON /SECT/ SECTOR
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &   untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
      character*1 digit
      character*71 file
      character*65 prefix
      character*71 scrtch
      data first/1/
C
      ITAP3 = 3
      ITAP6 = 6
 
c Initialize on the first call to rfile.
      if (first .eq. 1) then
        first = 0
        do 101 i = 1, numbuf
          bufunt(i) = 0
          bufrec(i) = 0
          wrtbuf(i) = 0
  101     continue
        do 102 I = 1, 128
          untbuf(i) = 0
  102     continue
        nxtbuf = 1
        endif
c
clj022789 I think some programs require SECTOR = 1024
      SECTOR = 1024
c
 
clj040389#ifdef HITACHI
clj040389clj Now maxrec is read in from file28.  The allocation procedures
clj040389c automatically write the number of blocks available to file28.
clj040389c Each line, in free format, contains the unit and the number of
clj040389c 4096 byte records that the primary extent of that unit occupies.
clj040389c Only the primary extent will be filled in the present scheme.
clj040389      rewind 28
clj040389 922  read(28,*,end=923) itapet,mrect
clj040389      if (itapet .eq. itape) goto 925
clj040389      goto 922
clj040389 923  continue
clj040389      write(6,*) 'rfile: storage for unit ', itape,
clj040389     &           'has not been allocated'
clj040389      call mabort
clj040389 925  continue
clj040389      mrec = mrect
clj040389      if (mrec .ne. 0) then
clj040389        open (unit=itape,access='DIRECT',form='UNFORMATTED',
clj040389     #        iostat=ierr,status='UNKNOWN',maxrec=MREC,recl=reclen*4)
clj040389      else
clj040389        open (unit=itape,access='DIRECT',form='UNFORMATTED',
clj040389     #        iostat=ierr,status='OLD',recl=reclen*4)
clj040389        endif
clj040389#elif IBMMVS
clj040389      open (unit=itape,access='DIRECT',recl=reclen*4,form='UNFORMATTED',
clj040389     &      iostat=ierr,status='UNKNOWN')
clj040389#elif SUN
clj040389      open (unit=itape,access='DIRECT',recl=reclen*4,form='UNFORMATTED',
clj040389     &      iostat=ierr,status='UNKNOWN',file=file)
clj040389#elif APOLLO
      call locate(5,'# FILES ##',ierr)
clj   write(*,*) 'rfile: ierr=',ierr
clj   write(6,*) 'rfile: ierr=',ierr
      if (ierr .ne. 0) then
        iprfx = 0
      else
 804    continue
        prefix = ' '
        read(5,'(2I5,5X,A65)',err=904) fstu,lstu,prefix
        if (fstu .eq. 0) goto 904
        if (lstu .eq. 0) lstu = fstu
        if ((itape .lt. fstu).or.(itape .gt. lstu)) goto 804
 904    continue
        iprfx = 0
 805    continue
        iprfx = iprfx + 1
        if (iprfx .eq. 66) goto 905
        if (prefix(iprfx:iprfx).ne.' ') goto 805
 905    continue
        iprfx = iprfx - 1
        endif
      file='file'//digit(itape/10)//digit(itape-10*(itape/10))
      if (iprfx .ne. 0) then
        scrtch = file
        file=prefix(1:iprfx)//scrtch
        endif
clj   write(*,*) 'rfile:',file
clj   write(6,*) 'rfile:',file
      open (unit=itape,access='DIRECT',recl=reclen*4,form='UNFORMATTED',
     &      iostat=ierr,status='UNKNOWN',file=file)
clj040389#endif
clj Double check the record length
c     inquire(itape,recl=nbytes)
c     write(6,*) 'rfile: nbytes = ', nbytes
C
      IF (IERR .ne. 0) then
        WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        WRITE(6,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        write(*,*) ' rfile: filename=',file
        write(6,*) ' rfile: filename=',file
        call ioerr(ierr)
        call mabort
        endif
C
C PLACE THE POINTER OF ITAPE AT THE BEGINNING
C
      CALL SREW(ITAPE)
C
      RETURN
      END
      SUBROUTINE iofc(ITAPE,JCODE)
C
      IMPLICIT INTEGER (A-Z)
      parameter (numbuf=10,reclen=1024)
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &   untbuf(128), bufrec(numbuf), buffer(reclen, numbuf)
C
C JCODE = 4     CLOSE AND DELETE FILE
C JCODE = 3     CLOSE AND SAVE FILE
C
CTJL  WRITE(*,*) ' CLOSING FILE',ITAPE,'  JCODE = ',JCODE
      IF(JCODE.NE.3.AND.JCODE.NE.4) THEN
        WRITE(*,*) ' INVALID JCODE IN RCLOSE,  JCODE = ',JCODE
        WRITE(*,*) ' FILE ',ITAPE,'  CLOSED AND SAVED.'
        JCODE = 3
        END IF
      IF(ITAPE.EQ.6) then
        write(6,*) 'rclose: a close on file 6 was attempted - aborting'
        call mabort
        endif
c
c     If the the files buffer has been written to, it must be flushed.
      call iofflu(itape,'RCLOSE')
 
c
      if (jcode.eq.4) then
        close (unit=itape,status='DELETE',iostat=ierr)
      else
        close (unit=itape,status='KEEP',iostat=ierr)
        end if
c
      IF(IERR.NE.0) THEN
        WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE
        WRITE(*,*) ' IERR,JCODE = ',IERR,JCODE
CTJL    STOP
        END IF
c
      RETURN
      END
      SUBROUTINE iofflu(file,caller)
c
c  NB. The first word of the file is considered to be "1" by R
c
      implicit integer(a-z)
c
      parameter (RECLEN=1024,numbuf=10)
      character*(*) caller
c
      common /io/     inp,iout
      common /pointr/ wptr(128),tptr(128)
      common /iobufs/ nxtbuf, wrtbuf(numbuf), bufunt(numbuf),
     &  untbuf(128), bufrec(numbuf), buffer(reclen,numbuf)
 
      flubuf = untbuf(file)
c If the file does not have a buffer then return
      if (flubuf .eq. 0) return
c
c
c           if buffer has been written to then flush it
            if (wrtbuf(flubuf) .ne. 0) then
              write(bufunt(flubuf),rec=bufrec(flubuf),iostat=ierr)
     &           (buffer(i,flubuf),i=1,reclen)
c             write(6,*) 'bioflu: writing:'
c             do 103 i = 1, 10
c               write(6,*) '      ', buffer(i,flubuf)
c103            continue
              if (ierr .ne. 0) then
                write(6,*) 'bioflu: error in write '
                write(6,*) 'bioflu: called by ', caller
                write(6,*) ' nxtbuf ',nxtbuf,' file ',file,
     &                                           ' flubuf ',flubuf
                write(6,*) ' untbuf(file) ', untbuf(file)
                write(6,*) ' record ', bufrec(flubuf)
                write(6,*) ' byte ', (bufrec(flubuf)-1)*reclen*4 + 1
                write(6,2) ' buffer ',' wrtbuf ',' bufunt ',' bufrec '
   2            format(4a10)
                do 101 i = 1, numbuf
                  write(6,1) i,wrtbuf(i),bufunt(i),bufrec(i)
   1              format(4i10)
 101              continue
                call ioerr(ierr)
                call mabort
                endif
              wrtbuf(flubuf) = 0
              endif
c
      return
      END
