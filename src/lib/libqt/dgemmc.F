	subroutine dgemmc(transa,transb,m,n,k,alpha,a,nra,nca,b,nrb,
     +                    ncb,beta,c)

	character*1 tra,trb
        integer transa,transb
	integer m,n,k,nra,nrb,nca,ncb
	double precision b(ncb,nrb),a(nca,nra),c(n,m)
	double precision alpha,beta

        if(transa.eq.0) tra = 'n'
        if(transb.eq.0) trb = 'n'
        if(transa.eq.1) tra = 't'
        if(transb.eq.1) trb = 't'

 	call dgemm(trb,tra,n,m,k,alpha,b,ncb,a,nca,beta,c,n)

	return
	end
	
