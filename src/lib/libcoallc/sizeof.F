      integer function isizof(data)
#if defined(sun)
      implicit undefined (a-z)
#else
      implicit integer (a-z)
#endif
      character*(*) data

      if (data .eq. 'REAL') then
        isizof = 2
      else if (data .eq. 'REAL_IN_BYTES') then
        isizof = 8
      else if (data .eq. 'INT') then
        isizof = 1
      else if (data .eq. 'INT_IN_BYTES') then
        isizof = 4
      else
        write(6,*) 'isizof: illegal type: ', data
        call mabort
        endif

      return
      end
