      subroutine rclose(unit,stat)
      integer unit,stat
      call ioclos(unit,stat)
      return
      end
