      SUBROUTINE bioflu(file,caller)
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
