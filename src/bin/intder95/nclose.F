C///////////////////////////////////////////////////////////////////////
      SUBROUTINE NCLOSE(ITAPE)
C
      IF(ITAPE.EQ.6) STOP ' YOU CANNOT CLOSE A FILE ON UNIT 6'
      JCODE=4
      REWIND(ITAPE)
      READ(ITAPE,*,END=20)
      JCODE=3
 20   CONTINUE
C
      IF (JCODE.EQ.4) THEN
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='DELETE')
      ELSE
             CLOSE (UNIT=ITAPE,IOSTAT=IERR,STATUS='KEEP')
      END IF
C
      IF (IERR.NE.0) THEN
             WRITE(*,*) ' ERROR ENCOUNTERED CLOSING FILE',ITAPE
             WRITE(*,*) ' IERR,JCODE = ',IERR,JCODE
      END IF
C
      RETURN
      END
