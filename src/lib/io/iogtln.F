      subroutine iogtln(input,line)
      integer input
      character*80 line

      read(input,fmt='(a80)') line
      return
      end
