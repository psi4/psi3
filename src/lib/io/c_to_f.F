c This is needed to interface the F77 gprgid routine to C.  The
c problem is that progid is character*(*) is cprgid, requiring that
c F77 hidden argument list is required.
      subroutine cprgid(progid)
      character*40 progid
      call gprgid(progid)
      return
      end

      integer function ffrdc(key,val)
      character*80 key
      character*80 val
      integer frdc
      ffrdc = frdc(key,val)
      return
      end
      integer function ffrdbo(key,val)
      character*80 key
      integer val
      integer frdboo
      ffrdbo = frdboo(key,val)
      return
      end
      integer function ffrdi(key,val)
      character*80 key
      integer val
      integer frdi
      ffrdi = frdi(key,val)
      return
      end
