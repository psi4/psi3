C
      SUBROUTINE RFILE(ITAPE)
C
C USE THIS TO OPEN ALL DIRECT ACCESS FILES
C
      IMPLICIT INTEGER (A-Z)
C
clj040389      parameter (reclen = 5120, numbuf = 10)
      parameter (reclen = 1024, numbuf = 10)
      common /pointr/ wptr(128),tptr(128)
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
