C//////////////////////////////////////////////////////////////////////
      subroutine bottom(iunit)
c  Skip to end of file 'iunit'
      character*1 dummy
1000  format(a1)
      k=1
10    read(iunit,1000,end=20) dummy
      k=k+1
      go to 10
20    continue
      rewind(iunit)
      do 30 i=1,k-1
30    read(iunit,1000) dummy
      return
      end
