      subroutine fentry(az,iz,maxcor)
c             INTDER
c*********************************************************************
c
c  1. Original version by Dr Wesley D Allen  (1984 - 1987 ?)
c
c  2. Modified to perform multiple sets of transformations
c                                  - Brian F Yates   21 December  1987
c
c  3. Modified to write force constants out to BMAT
c                                  - Brian F Yates    3 February  1988
c
c  4. Core dimension reduced from 1200000 to 800000
c                                  - Brian F Yates   11 May       1988
c
c  5. Reading of gradients when GEOM=1 fixed up.  Search cbfy24
c                                  - Brian F Yates   24 September 1988
c
c  6. Use subroutine 'bottom' to find end of file. Search cbfy7
c                                  - Brian F Yates    7 October   1988
c
c  7. Include Wesley Allen's modifications of January 16 1988:
c       - Program now does  nder=4  with  neq=1
c       - Some additions for the transformation of dipole moment
c         derivatives also included
c     Include Wesley's corrections for dipole moment transformations
c     of February 15 1988.
c                                  - Brian F Yates    3 February  1989
c
c  8. Take out all routines that are in IOMR FORTRAN.
c                                  - Brian F Yates    3 February  1989
c
c  9. Moved to PSI distribution disk on 2/4/89 - clj.
c
c 10. Put in run time core allocation on 3/3/89 - clj.
c
c 11. Write force constants to file fconst.dat rather than to bmat file
c        12/26/90 - Edward Seidl (cets122690)
c
c*********************************************************************
c 12.
c  July 1991 Edward Seidl(cets072991)
c  get rid of intco.dat file, put everything in input
c  split SYMM array into SYMM and ASYMM components
c
c  reworked by Edward Seidl (cets050491), May 1991
c  now read from intco.dat in format readable by new input parsing
c    routines.  intco.dat also used by bmat
c
c 13. George Vacek (croag) 
c	Wed Oct 21 16:08:19 EDT 1992
c	Added ZPVE determination in gfmat and normco subroutines
c
c 14. T. Daniel Crawford (Crawdad.7.18.93)
c	Sunday, July 18, 1993 
c	Added masses through U(92).  Some masses are omitted due 
c	inaccurate or just plain unavailable data.
c
c 15. Croag again with the zpve (3.30.94)
c	zpve was screwed up for a couple cases by the reordering 
c	of frequencies, so i bummed things around a bit in normco subr
***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include <libparse/error.h>
      integer errcod, fwkclr, fwkadd, fcnt, frdi
      integer nst,nbnd,ntor,nlin,nout,nspf
      integer coallc,isizof
      integer iz(1)
      integer type,u,ia,iu,lbls
      character junk
      DIMENSION AZ(maxcor)
 1    FORMAT(A10)
 2    FORMAT(2I5,25X,I5)
 3    FORMAT(1X,'NO INPUT FOUND.')
clj030389      MAXCOR= 700000
cets041389  IIN1=5
      iin1=7
      IOUT=6

      call psinit('APPEND')
C      CALL TSTART(3)
      CALL TSTART(IOUT)

croag      call ffile(9,'intder1',0)
croag      call locate(9,'# INTDER #',errcod)
croag      if(errcod.eq.0) then
croag        write(6,*) ' '
croag        write(6,*)'I see you have an old-style intder1.dat file.'
croag        write(6,*)'Please note: I am going to ignore that file.'
croag        write(6,*)'This version of intder only uses the new-style'
croag        write(6,*)
croag     + ' input format (ie., you must have intco: in your input).'
croag        write(6,*) ' '
croag        write(6,*)'If you really want to use the old-style input'
croag        write(6,*)' then you must use the INTDERV1 command.'
croag        write(6,*) ' '
croag      endif
croag      close(9,status='keep')

      call ffile(iin1,'intco',0)
      call ffile(8,'geom',0)
      call coinit(iz,maxcor*isizof('REAL'))

      errcod = fwkclr()
      errcod = fwkadd(5,'DUMMY')
      errcod = fwkadd(iin1,'DUMMY')

      nst=0
      nbnd=0
      ntor=0
      nout=0
      nlin=0
      nspf=0
      errcod = fcnt(':INTCO:STRE',nst)
      errcod = fcnt(':INTCO:BEND',nbnd)
      errcod = fcnt(':INTCO:TORS',ntor)
      errcod = fcnt(':INTCO:OUT',nout)
      errcod = fcnt(':INTCO:LIN2',nlin)
      errcod = fcnt(':INTCO:SPF',nspf)

      ns=nst+nbnd+ntor+nout+nlin+nspf

      iused = 0
      type = coallc(iz,ns*isizof('REAL'))
      ia = coallc(iz,6*ns*isizof('INT'))
      is = coallc(iz,ns*isizof('REAL'))
      iused = iused+ns*isizof('INT')
      iused = iused+2*ns*isizof('REAL')

      call rdsimp(iz(type),iz(ia),iz(is),ns,nst,nbnd,
     1  ntor,nlin,nout,nspf)

      nsym=0
      errcod = fcnt(':INTCO:SYMM',nsym)
      nasym=0
      errcod = fcnt(':INTCO:ASYMM',nasym)

      nsym=nsym+nasym

      lbls = coallc(iz,20*nsym*isizof('INT'))
      iu = coallc(iz,101*ns*isizof('INT'))
      u = coallc(iz,100*ns*isizof('REAL'))
      iused=iused+20*nsym*isizof('INT')
      iused=iused+100*ns*isizof('REAL')
      iused=iused+101*ns*isizof('INT')

      if(nsym.gt.0) 
     + call rdsymm(iz(lbls),iz(iu),iz(u),nsym,nasym,ns)

      errcod = fwkclr()
      errcod = fwkadd(8,'DUMMY')

      na=0
      errcod = fcnt(':GEOMETRY',na)
      if(errcod.ne.EOK.or.na.eq.0) then
        call ffile(11,' ',0)
        rewind(11)
        read(11,100) junk
        read(11,'(i5)') na
 100    format(a)
        close(11)
        if(na.eq.0) then
          errcod = frdi('NATOMS',na)
          if(errcod.ne.EOK) then
            write(6,*) 'need NATOMS'
            call pabort
            call qabort
          endif
        endif
      endif

      errcod = fwkclr()
      errcod = fwkadd(5,'INTDER')
      errcod = fwkadd(5,'DEFAULT')

      ndum=0
      errcod = fcnt('DUMMY_ATOMS',ndum)

      NAD=NA+ndum
      NC=3*NA
      NCD=3*NAD
      
      ixa = coallc(iz,nad*3*isizof('REAL'))
      ixmass = coallc(iz,na*isizof('REAL'))
      iss = coallc(iz,ns*isizof('REAL'))
      ib = coallc(iz,ns*nc*isizof('REAL'))
      ibs = coallc(iz,nc*nc*isizof('REAL'))
      iarr = coallc(iz,nc*nc*isizof('REAL'))
      iugf = coallc(iz,nc*nc*isizof('REAL'))
      ixs = coallc(iz,nc*nc*isizof('REAL'))
      ixt = coallc(iz,nc*nc*isizof('REAL'))
      ixki = coallc(iz,3*nc*isizof('REAL'))

      iused=iused+5*nc*nc*isizof('REAL')
      iused=iused+3*nc*isizof('REAL')
      iused=iused+ns*nc*isizof('REAL')
      iused=iused+ns*isizof('REAL')
      iused=iused+na*isizof('REAL')
      iused=iused+nad*3*isizof('REAL')
 
c it would be best to leave some room for error
      ileft = maxcor-iused
      icore = coallc(iz,ileft*isizof('REAL'))
      ileft = ileft*isizof('INT')/isizof('REAL')

      call intder(nad,nc,ns,iz(ixa),iz(ixmass),iz(is),
     + iz(iss),iz(ib),iz(ibs),
     1 iz(iarr),iz(iugf),iz(ixs),iz(ixt),iz(ixki),
     2 iz(type),iz(u),iz(ia),iz(iu),iflag,
     3 nst,nbnd,ntor,nout,nlin,nspf,nsym,na,
     4 iz(lbls),iz(icore),ileft)

C      CALL TSTOP(3)
      CALL TSTOP(IOUT)
      STOP
 1000 WRITE(IOUT,3)
      STOP
      END
