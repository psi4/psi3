C   *** (MEMBER:DMPUNT  ) CONVERTED IN POTOPS-00-00. CONV.DATE:89-03-01
      subroutine dmpunt(unit,len)
      implicit integer (a-z)
      real*8 rbuf(512)
      integer ibuf(1024)
      common /pointr/ wptr(128), tptr(128)
      equivalence (rbuf,ibuf)
      write(6,*) 'dump for unit ', unit
      svwptr = wptr(unit)
      svtptr = tptr(unit)
      start = 1
      do 101 i = 1, len, 1024
        if (len - i + 1 .lt. 1024) then
          ndata = len - i + 1
        else
          ndata = 1024
          endif
        start1 = start
        call wreadw(unit,ibuf,ndata,start,start)
        do 201 j = 1, ndata
          write(6,1) start1+j-1, ibuf(j)
   1      format(1X,I7,1X,I15)
 201      continue
 101    continue
      wptr(unit) = svwptr
      tptr(unit) = svtptr
      return
      end
