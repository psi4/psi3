C  Fortran function to diagonalize a square nonsymmetric matrix using
C  the DGEEV BLAS routine.  The matrix must be square because of the
C  automatic matrix transposes which are performed below.  The eigenvalues
C  are returned in random order.

        subroutine dgeev_call(jobvlI, jobvrI, N, A, lda, wr,
     +   wi, VL, ldvl, VR, ldvr, work, lwork, info)

        character jobvl,jobvr
        integer jobvlI,jobvrI,i,j
        integer info,lda,ldvl,ldvr,lwork,N
        double precision A(lda,N),VL(ldvl,N),VR(ldvr,N)
        double precision wi(N),work(lwork),wr(N),temp

        do 100 i=1,N
          do 101 j=1,(i-1)
             temp = A(i,j)
             A(i,j) = A(j,i)
             A(j,i) = temp
          101 continue
        100 continue

        do 102 i=1,N
          do 103 j=1,(i-1)
             temp = VR(i,j)
             VR(i,j) = VR(j,i)
             VR(j,i) = temp
          103 continue
        102 continue

        if(jobvlI.eq.0) jobvl = 'N'
        if(jobvrI.eq.0) jobvr = 'N'
        if(jobvlI.eq.1) jobvl = 'V'
        if(jobvrI.eq.1) jobvr = 'V'

        call dgeev(jobvl, jobvr, N, A, lda, wr,
     +       wi, VL, ldvl, VR, ldvr, work, lwork, info)

        do 105 i=1,N
          do 106 j=1,(i-1)
             temp = A(i,j)
             A(i,j) = A(j,i)
             A(j,i) = temp
          106 continue
        105 continue

        do 200 i=1,N
           do 201 j=1,(i-1)
              temp = VR(i,j);
              VR(i,j) = VR(j,i);
              VR(j,i) = temp;
           201 continue
        200 continue

        return
        end

