      subroutine tstart(lu)
      integer lu
      integer istart, istop
      real*8 utime, stime
      character*80 host
      character*80 date
      save istart, host

      call timnow(istart)
      call cfdate(istart,date)
      call gthost(host)

      call wtstars(lu)
      write(lu,10) host
      write(lu,12) date
      write(lu,*)

      return
c
      entry tstop(lu)

      call timnow(iend)
      call cfdate(iend,date)
      call ctimes(utime,stime)

      call wtstars(lu)
      write(lu,11) host
      write(lu,12) date
      write(lu,*)
      write(lu,13) utime
      write(lu,14) stime
      write(lu,15) iend-istart

  10  format('tstart called on ',a40)
  11  format('tstop called on ',a40)
  12  format(a24)
  13  format('user time   = ',f10.2)
  14  format('system time = ',f10.2)
  15  format('total time  = ',i10)

      return
      end

      subroutine wtstars(lu)
      integer lu
      character*20 stars
      stars='********************'
      write(lu,10) stars,stars,stars,stars
  10  format(4(a20))
      return
      end
