C///////////////////////////////////////////////////////////////////////

      SUBROUTINE OPENFF (ITAPE, IRECL)

 	IMPLICIT INTEGER (A-Z)
	character*1 digit
      character*71 file
      character*80 prefix
      character*71 scrtch


C Open formatted file.  Set IRECL=0 for sequential access.

      call locate(5,'# FILES ##',ierr)
      if (ierr .ne. 0) then
        iprfx = 0
      else
 804    continue
        prefix = ' '
        read(5,'(A80)',err=904) prefix

 904    continue
        iprfx = 0
 805    continue
        iprfx = iprfx + 1
        if (iprfx .eq. 81) goto 905
        if (prefix(iprfx:iprfx).ne.' ') goto 805
 905    continue
        iprfx = iprfx - 1
        endif

      if (itape.ne.3) then
         file='file'//digit(itape/10)//digit(itape-10*(itape/10))
         if (iprfx .ne. 0) then
              scrtch = file
              file=prefix(1:iprfx)//'.'//scrtch
              endif
                      else
         file = prefix(1:iprfx)//'.CHECK'
         endif

      if (irecl.eq.0) then
CSA   write (6,*) 'open seq',itape,irecl
	open (unit=itape, file=file, form='FORMATTED',
     $            access='SEQUENTIAL', iostat=ierr)

			     else
      write (6,*) 'open direct',itape,irecl
		open (unit=itape, file=file, form='FORMATTED',
     $            access='DIRECT', recl=irecl, iostat=ierr)
		end if

	if (ierr.ne.0) then

        WRITE(*,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        WRITE(6,*) ' ERROR ENCOUNTERED OPENING FILE ',ITAPE,' IERR =',
     &            ierr
        write(*,*) ' rfile: filename=',file
        write(6,*) ' rfile: filename=',file
        call ioerr(ierr)
        call mabort
        endif

	return

	end
