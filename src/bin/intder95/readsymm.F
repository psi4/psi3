      subroutine rdsymm(lbls,iu,u,ns,nsym)
      implicit integer(a-z)
#include <error.h>
      character*80 lbls(nsym)
      real*8 u(ns,100),cff
      real*8 norm
      dimension iu(ns,0:100)
      integer nsym, nasym, nsymco

      nasym=0
      errcod = fcnt(':INTCO:ASYMM',nasym)

      nsym=0
      errcod = fcnt(':INTCO:SYMM',nsym)
      if(nsym.ne.0) then
        errcod = frdpos(':INTCO:SYMM')
        if(errcod.ne.EOK) goto 901
      endif
      nsymco = nsym+nasym

      do 115  k=1,nsymco
 115  iu(k,0)=0

      errcod = frdlpn()
      do 10 symco=1,nsymco-nasym
        errcod = frdlpn()
        errcod = frdchr(lbls(symco))
        if(errcod.ne.EOK) goto 901
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 901
        nsimpl=0
 201    errcod = frdint(simno)
        if(errcod.ne.Eevect) then
          nsimpl=nsimpl+1
          if(simno.lt.0) then
            iu(symco,nsimpl) = -simno
            u(symco,nsimpl) = -1.0d0
          else
            iu(symco,nsimpl) = simno
            u(symco,nsimpl) = 1.0d0
          endif
          goto 201
        endif
        errcod = frdlpn()
        if(errcod.eq.EOK) then
          do 100 sim=1,nsimpl
            errcod = frdflt(cff)
            u(symco,sim)=u(symco,sim)*cff
 100      continue
          errcod=frdrpn()
          errcod=frdrpn()
        endif
        norm = 0.0d0
        do 101 sim=1,nsimpl
          norm = norm + u(symco,sim)*u(symco,sim)
 101    continue
        norm = 1.0d0/dsqrt(norm)
        do 102 sim=1,nsimpl
           u(symco,sim) = u(symco,sim) * norm
 102    continue
        iu(symco,0)=nsimpl
  10  continue

      errcod = frdpos(':INTCO:ASYMM')

      errcod = frdlpn()
      do 20 symco=nsymco-nasym+1,nsymco
        errcod = frdlpn()
        errcod = frdchr(lbls(symco))
        if(errcod.ne.EOK) goto 902
        errcod = frdlpn()
        if(errcod.ne.EOK) goto 902
        nsimpl=0
 221    errcod = frdint(simno)
        if(errcod.ne.Eevect) then
          nsimpl=nsimpl+1
          if(simno.lt.0) then
            iu(symco,nsimpl) = -simno
            u(symco,nsimpl) = -1.0d0
          else
            iu(symco,nsimpl) = simno
            u(symco,nsimpl) = 1.0d0
          endif
          goto 221
        endif
        errcod = frdlpn()
        if(errcod.eq.EOK) then
          do 200 sim=1,nsimpl
            errcod = frdflt(cff)
            u(symco,sim)=u(symco,sim)*cff
 200      continue
          errcod=frdrpn()
          errcod=frdrpn()
        endif
        iu(symco,0)=nsimpl

        norm = 0.0d0
        do 234 sim=1,nsimpl
          norm = norm + u(symco,sim)*u(symco,sim)
 234    continue
        norm = 1.0d0/dsqrt(norm)
        do 235 sim=1,nsimpl
           u(symco,sim) = u(symco,sim) * norm
 235    continue

  20  continue
      return 

 901  write(6,*) 'error parsing SYMM vector'
      call pabort
      call qabort
 902  write(6,*) 'error parsing ASYMM vector'
      call pabort
      call qabort
      end
