      subroutine wwritw(unit,buffer,nwords,fword,nxtwrd)
      implicit integer (a-z)
      common /pointr/ wptr(128),tptr(128)
      dimension buffer(nwords)

      call iowrr(unit,buffer,fword,nwords)
      nxtwrd = fword + nwords
      if(nxtwrd.lt.fword) then
         write(*,*) 'Error: fortran I/O error'
         call qabort
      endif
      wptr(unit) = nxtwrd
      return
      end
