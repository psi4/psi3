      subroutine mxmb(a,ia,ja, b,ib,jb, c,ic,jc, nrow, nlnk, ncol)
      implicit real*8 (a-h, o-z)
      integer CUTOFF
c The CUTOFF parameter sets a maximum N*N*2 for which the faster
c matrix multiplies are used.  These malloc storage to rewrite the
c passed arrays.
#if defined(AIX)
      parameter(CUTOFF = 120000)
#else
      parameter(CUTOFF =  80000)
#endif
      dimension a(1), b(1), c(1)

c  Convert mxmb calls to mxmnn mxmnt mxmtn and mxmtt calls.


#if 1
      isize=nlnk*(nrow+ncol)
      if ( ic .eq. 1 ) then
          if ( ia .eq. 1 ) then
              if ( ib .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmnn(a,ja,b,jb,c,jc,nrow,nlnk,ncol)
                  else
                     call mmult(a,ja,0,b,jb,0,c,jc,0,nrow,nlnk,ncol,1)
                     endif
              else if ( jb .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmnt(a,ja,b,ib,c,jc,nrow,nlnk,ncol)
                  else
                     call mmult(a,ja,0,b,ib,1,c,jc,0,nrow,nlnk,ncol,1)
                     endif
              else
                  call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
                  endif
          else if ( ja .eq. 1 ) then
              if ( ib .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmtn(a,ia,b,jb,c,jc,nrow,nlnk,ncol)
                  else
                     call mmult(a,ia,1,b,jb,0,c,jc,0,nrow,nlnk,ncol,1)
                     endif
              else if ( jb .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmtt(a,ia,b,ib,c,jc,nrow,nlnk,ncol)
                  else
                     call mmult(a,ia,1,b,ib,1,c,jc,0,nrow,nlnk,ncol,1)
                     endif
              else
                  call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
                  endif
          else
              call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
              endif
      else if ( jc .eq. 1 ) then
          if ( ia .eq. 1 ) then
              if ( ib .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmtt(b,jb,a,ja,c,ic,ncol,nlnk,nrow)
                  else
                     call mmult(a,ja,0,b,jb,0,c,ic,1,nrow,nlnk,ncol,1)
                     endif
              else if ( jb .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmnt(b,ib,a,ja,c,ic,ncol,nlnk,nrow)
                  else
                     call mmult(a,ja,0,b,ib,1,c,ic,1,nrow,nlnk,ncol,1)
                     endif
              else
                  call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
                  endif
          else if ( ja .eq. 1 ) then
              if ( ib .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmtn(b,jb,a,ia,c,ic,ncol,nlnk,nrow)
                  else
                     call mmult(a,ia,1,b,jb,0,c,ic,1,nrow,nlnk,ncol,1)
                     endif
              else if ( jb .eq. 1 ) then
                  if(isize.gt.CUTOFF) then
                     call mxmnn(b,ib,a,ia,c,ic,ncol,nlnk,nrow)
                  else
                     call mmult(a,ia,1,b,ib,1,c,ic,1,nrow,nlnk,ncol,1)
                     endif
              else
                  call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
                  endif
          else
              call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
              endif
      else
          call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
          endif
#else

      isize=nlnk*(nrow+ncol)
      if ( ic .eq. 1 ) then
          if ( ia .eq. 1 ) then
             if ( ib .eq. 1 ) then
               call dgemm('n','n',nrow,ncol,nlnk,1.0,a,ja,b,jb,1.0,c,jc)
             else if ( jb .eq. 1 ) then
               call dgemm('n','t',nrow,ncol,nlnk,1.0,a,ja,b,ib,1.0,c,jc)
             else
               call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
             endif
          else if ( ja .eq. 1 ) then
            if ( ib .eq. 1 ) then
               call dgemm('t','n',nrow,ncol,nlnk,1.0,a,ia,b,jb,1.0,c,jc)
            else if ( jb .eq. 1 ) then
               call dgemm('t','t',nrow,ncol,nlnk,1.0,a,ia,b,ib,1.0,c,jc)
            else
               call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
               endif
          else
              call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
              endif
      else if ( jc .eq. 1 ) then
          if ( ia .eq. 1 ) then
            if ( ib .eq. 1 ) then
               call dgemm('t','t',ncol,nrow,nlnk,1.0,b,jb,a,ja,1.0,c,ic)
            else if ( jb .eq. 1 ) then
               call dgemm('n','t',ncol,nrow,nlnk,1.0,b,ib,a,ja,1.0,c,ic)
            else
               call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
               endif
          else if ( ja .eq. 1 ) then
           if ( ib .eq. 1 ) then
              call dgemm('t','n',ncol,nrow,nlnk,1.0,b,jb,a,ia,1.0,c,ic)
           else if ( jb .eq. 1 ) then
              call dgemm('n','n',ncol,nrow,nlnk,1.0,b,ib,a,ia,1.0,c,ic)
           else
              call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
              endif
          else
              call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
              endif
      else
          call mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
          endif
#endif
      return

      end

      subroutine mxmbol(a,ia,ja,b,ib,jb,c,ic,jc,nrow,nlnk,ncol)
      implicit real*8 (a-h, o-z)
      dimension a(ia,1), b(ib,1), c(ic,1)
      write(0,*) 'mxmbol: doesn''t exist yet'
      write(6,*) 'mxmbol: doesn''t exist yet'
      call lnkerr('Somebody should write a proper mxmb.')
      stop
      end
