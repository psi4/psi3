C   *** (MEMBER:DMPRFL  ) CONVERTED IN POTOPS-00-00. CONV.DATE:89-03-01
      subroutine dmprfl(unit,offset,len)
      implicit integer (a-z)
      real*8 rbuf(512)
      integer ibuf(1024)
      integer offset,len
      common /pointr/ wptr(128), tptr(128)
      equivalence (rbuf,ibuf)
      write(6,*) 'dump for unit ', unit
      svwptr = wptr(unit)
      svtptr = tptr(unit)
      start = 1 + intowp(offset)
      start1 = 1
      do 101 i = 1, len, 512
        if (len - i + 1 .lt. 512 ) then
          ndata = len - i + 1
        else
          ndata = 512
          endif
        call wreadw(unit,ibuf,intowp(ndata),start,start)
        do 201 j = 1, ndata
          write(6,1) start1+j-1, rbuf(j)
   1      format(1X,I7,1X,F22.15)
 201      continue
        start1 = start1 + ndata
 101    continue
      wptr(unit) = svwptr
      tptr(unit) = svtptr
      return
      end
