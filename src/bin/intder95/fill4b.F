C///////////////////////////////////////////////////////////////////
	subroutine fill4b (NC, A)
	IMPLICIT REAL*8 (A-H,O-Z)
	dimension A(NC,NC,NC,NC)

	do 1 i = 1, NC
	do 1 j = 1, i
	do 1 k = 1, j
	do 1 l = 1, k
	A(i,j,l,k) = A(i,j,k,l)
        A(i,k,j,l) = A(i,j,k,l)
        A(i,k,l,j) = A(i,j,k,l)
        A(i,l,j,k) = A(i,j,k,l)
        A(i,l,k,j) = A(i,j,k,l)
        A(j,i,k,l) = A(i,j,k,l)
        A(j,i,l,k) = A(i,j,k,l)
        A(j,k,i,l) = A(i,j,k,l)
        A(j,k,l,i) = A(i,j,k,l)
        A(j,l,i,k) = A(i,j,k,l)
        A(j,l,k,i) = A(i,j,k,l)
        A(k,i,j,l) = A(i,j,k,l)
        A(k,i,l,j) = A(i,j,k,l)
        A(k,j,i,l) = A(i,j,k,l)
        A(k,j,l,i) = A(i,j,k,l)
        A(k,l,i,j) = A(i,j,k,l)
        A(k,l,j,i) = A(i,j,k,l)
        A(l,i,j,k) = A(i,j,k,l)
        A(l,i,k,j) = A(i,j,k,l)
        A(l,j,i,k) = A(i,j,k,l)
        A(l,j,k,i) = A(i,j,k,l)
        A(l,k,i,j) = A(i,j,k,l)
 1      A(l,k,j,i) = A(i,j,k,l)
	return
	END
