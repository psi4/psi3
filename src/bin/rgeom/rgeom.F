      program rgeom
c Rotate the geometries in the geom.dat file.
c mtext = number of entries in the geom.dat file.
      implicit integer (a-h,o-z)
      parameter (mtext=2000)
      character*80 text(mtext)
      character*8 first8
      integer prcntr,errcod
      logical header

c Initialize
      call psinit('APPEND')
c     call tstart(6)

c Set up print control
      errcod = prcntr('ON BRIEF')
      call prdup()

c Open the geometry file
      call fioopn(7,'geom.dat','USER','WRITE')

c Read the geom.dat file into text
      rewind 7
      itext = 1
      nhead = 0
      header = .true.
      igeom = 0
      ngeom1 = 0
 801  continue
        if (itext.gt.mtext) then
          write(6,*) 'in RGEOM, the geom.dat rotation program:'
          write(6,*) 'ERROR: the geom.dat file is too large'
          call qabort
          endif
        read(7,1,end=901) text(itext)
   1    format(a80)
        first8 = text(itext)(1:8)
        call touppr(first8)
        if (text(itext)(1:2).eq.'%%') then
          if (header) then
            nhead = nhead + 1
          else
            write(6,*) 'in RGEOM, the geom.dat rotation program:'
            write(6,*) 'ERROR: the geom.dat file has the wrong format:'
            write(6,*) ' all the lines beginning ''%%'' must be at the'
            write(6,*) ' beginning of the file--these stay at the'
            write(6,*) ' beginning after the rotation of the geometries'
            call qabort
            endif
        else if (first8.eq.'GEOMETRY') then
          if (header) then
            header = .false.
            endif
          igeom = igeom + 1
          if (igeom.eq.2 .and. ngeom1.eq.0) then
            ngeom1 = itext - 1 - nhead
            endif
        else if (igeom.eq.1 .and. text(itext)(1:1) .eq. '%') then
            ngeom1 = itext - 1 - nhead
          endif
        itext = itext + 1
        goto 801
 901  continue

      itext = itext - 1

c If there isn't more than one geometry, then we are done.
      if (igeom.le.1) then
        goto 902
        endif

c Write out the geom.dat file with the geometries rotated
      rewind 7
      do 101 i = 1, nhead
        write(7,1) text(i)
 101    continue

      do 102 i = nhead + ngeom1 + 1, itext
        write(7,1) text(i)
 102    continue

      do 103 i = nhead + 1, ngeom1 + nhead
        write(7,1) text(i)
 103    continue

c Finish up
 902  continue
c     call tstop(6)
      call psdone()
      stop
      end

