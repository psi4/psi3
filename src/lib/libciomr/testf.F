#define NN 5
#define TRI 15
#define RSP 0
#define MXM 1
      program testf
#if RSP
      real*8 arr(TRI), res(NN,NN), vecs(NN,NN),vals(NN)
      real*8 scr1(NN),scr2(NN)
      real*8 k
#else
      real*8 arr(NN,NN),vecs(NN,NN),res(NN,NN)
#endif

      open(unit=6,file='test.out')
      rewind(6)

      n=NN
      nn=TRI
      k=1.0/100000.0

#if RSP
      do 10 i = 1,nn
         arr(i) = dble(i)
#else
   10 continue
       do 20 i=1,n
         do 30 j = 1,i
           arr(i,j) = k
           arr(j,i) = k+1.0
#if 0
           k=(k*100000.0d0 + 1.0d0)/10000.0d0
#else
           k=k+1.0d0
#endif
           vecs(i,j) = k
           vecs(j,i) = k+1.0
   30    continue
   20  continue
#endif
      call tstart(6)

#if RSP
      call rsp(n,n,nn,arr,vals,1,vecs,scr1,scr2)
      call eivout(vecs,vals,n,n,n,n,6)
#else
      call mxmb(arr,n,1,vecs,n,1,res,1,n,n,n,n)
      call matout(arr,n,n,n,n,6)
      call tstop(6)
      call mmult(arr,n,1,vecs,n,1,res,n,0,n,n,n,0)
      call matout(arr,n,n,n,n,6)
      call tstop(6)
      call dgemm('t','t',n,n,n,1,arr,n,vecs,n,0,res,n)
      call matout(arr,n,n,n,n,6)
#endif
      call tstop(6)

      end
