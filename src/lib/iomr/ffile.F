C
      SUBROUTINE FFILE(ITAPE,name,code)
C
C USE THIS TO OPEN ALL FORMATTED SEQUENTAL ACCESS FILES
C
      IMPLICIT INTEGER (A-Z)
C
      integer code
      character*1 digit
      character*12 file
      character*(*) name
      character*1 junk

      if (name .eq. ' ') then
        file='file'//digit(itape/10)//digit(itape-10*(itape/10))
     &       //'.dat'
      else
        file=name//'.dat'
        endif

      open (unit=itape,access='SEQUENTIAL',
     &      form='FORMATTED',
     &      iostat=ierr,status='UNKNOWN',file=file)

      IF (IERR .ne. 0) then
        WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        call ioerr(ierr)
        call mabort
        endif

      if (iand(code,1) .eq. 1) then
c       Go to the end of the file
 800    read(itape,'(A1)',end=900) junk
        goto 800
 900    continue
        endif
C
      RETURN
      END
