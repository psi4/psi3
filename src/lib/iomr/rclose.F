      SUBROUTINE RCLOSE(ITAPE,JCODE)
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
      call bioflu(itape,'RCLOSE')
 
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
