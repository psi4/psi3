clj This is the old find core routine.  It has been replaced by a version
clj that uses the new style input.
      subroutine fndcro(maxcor)
      integer maxcor
cets022091 this is a hack to make matout work, the common integer linewd has
c          to be set
c          see fndcor.f
      call psinit('TERMINAL')
cets022091
cets the next open statement was not in the mips version?
crawdad This is the second open call for unit 5 ???
crawdad      open(unit=5,file='input.dat',status='old')
c     call drum
      call locate(5,'# CORE ###',ierr)
      if (ierr .ne. 0) then
         call locate(5,'# MAXCOR #',ierr)
         if (ierr .ne. 0) then
           maxcor = 0
clj see comments about psdone below
clj	   call psdone()
           return
           endif
         endif
      read(5,*) maxcor

clj problems with multiple psinit calls should be handled in
clj psinit by examining the "inited" variable
clj calling psdone now will, in some cases, reset the output to fort.6
cljcrawdad
clj	call psdone()
cljcrawdad          close(5)
      return
      end
