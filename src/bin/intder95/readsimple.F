      subroutine rdsimp(type,ia,s,ns)
      implicit integer(a-z)
#include <error.h>
      character*5 type(ns)
      character*80 keywrd
      real*8 s(ns)
      dimension ia(ns,6)
      integer nst,nbnd,ntor,nlin
      integer dum1

      nst=0
      errcod = fcnt(':INTCO:STRE',nst)
      if(nst.ne.0) then
        keywrd = ':INTCO:STRE'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 10 i=1,nst
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' STRE'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod = frdrpn()
  10      continue
        endif
      endif
      
      nbnd = 0
      errcod = fcnt(':INTCO:BEND',nbnd)
      if(nbnd.ne.0) then
        keywrd = ':INTCO:BEND'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 11 i=1,nbnd
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' BEND'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod = frdrpn()
  11      continue
        endif
      endif

      ntor = 0
      errcod = fcnt(':INTCO:TORS',ntor)
      if(ntor.ne.0) then
        keywrd = ':INTCO:TORS'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 12 i=1,ntor
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' TORS'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod = frdrpn()
  12      continue
        endif
      endif

      nlin = 0
      errcod = fcnt(':INTCO:LIN1',nlin)
      if(nlin.ne.0) then
        keywrd = ':INTCO:LIN1'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 13 i=1,nlin
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' LIN1'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod = frdrpn()
  13      continue
        endif
      endif

      nout = 0
      errcod = fcnt(':INTCO:OUT',nout)
      if(nout.ne.0) then
        keywrd = ':INTCO:OUT'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 14 i=1,nout
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)='  OUT'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod = frdrpn()
  14      continue
        endif
      endif

      nspf = 0
      errcod = fcnt(':INTCO:SPF',nspf)
      if(nspf.ne.0) then
        keywrd = ':INTCO:SPF'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 15 i=1,nspf
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)='  SPF'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdflt(s(ij))
            if(errcod.eq.Eevect) goto 15
            errcod = frdrpn()
  15      continue
        endif
      endif

      nlinx = 0
      errcod = fcnt(':INTCO:LINX',nlinx)
      if(nlinx.ne.0) then
        keywrd = ':INTCO:LINX'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 16 i=1,nlinx
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' LINX'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod = frdrpn()
  16      continue
        endif
      endif

      nliny = 0
      errcod = fcnt(':INTCO:LINY',nliny)
      if(nliny.ne.0) then
        keywrd = ':INTCO:LINY'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 17 i=1,nliny
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' LINY'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod = frdrpn()
  17      continue
        endif
      endif

      nrcom = 0
      errcod = fcnt(':INTCO:RCOM',nrcom)
      if(nrcom.ne.0) then
        keywrd = ':INTCO:RCOM'
        errcod = frdpos(keywrd)
        if(errcod.eq.EOK) then
          errcod = frdlpn()
          if(errcod.ne.EOK) goto 901
          do 18 i=1,nrcom
            errcod = frdlpn()
            if(errcod.ne.EOK) goto 901
            errcod = frdint(ij)
            type(ij)=' RCOM'
            errcod=frdint(ia(ij,1))
            errcod=frdint(ia(ij,2))
            errcod=frdint(ia(ij,3))
            errcod=frdint(ia(ij,4))
            errcod = frdrpn()
  18      continue
        endif
      endif

      return

  901 continue
      write(6,*) 'trouble parsing simple internal coordinates'
      call pabort
      call qabort
      end
