      subroutine rfile(iunit)
      implicit integer(a-z)
      common /pointr/ wptr(128), tptr(128)
      common /sect/ sector

      sector = 1024
      wptr(iunit) = 1
      call ioopen(iunit)
      return
      end
