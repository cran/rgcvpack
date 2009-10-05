      subroutine dcpmut (x,ldx,nobs,npar,jpvt,job)
      integer ldx,nobs,npar,jpvt(npar),job
      double precision x(ldx,npar)
c
c Purpose: permute the columns of the matrix x according to the index 
c	vector jpvt (either forward or backward permutation).
c
c On Entry:
c   x(ldx,npar)		matrix whose columns are to be permuted
c   ldx			leading dimension of x as declared
c			in the calling program
c   nobs		number of rows of x used
c   npar		number of columns of x
c   jpvt(npar)		indices of the permutation
c   job			indicator of forward or backward permutation
c			if job = 0 forward permutation  
c				x(jpvt(i)) moved to x(i)
c			else backward permutation
c				x(i) moved to x(jpvt(i))
c On Exit:
c   x(ldx,npar)		matrix with columns permuted
c
c Subprograms Called Directly
c     Blas	- dswap
c
c  Written:	Yin Ling	U. of Maryland, August,1978
c
c $Header: dcpmut.f,v 2.100.1.1 86/10/07 12:47:45 lindstrom Exp $
c
      integer i,j,k
c
      if (npar .le. 1) then
         return
      endif
      do 10 j = 1,npar
         jpvt(j) = -jpvt(j)
   10 continue
      if (job .eq. 0) then
c		forward permutation
         do 30 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 30
            endif
            j = i
            jpvt(j) = -jpvt(j)
            k = jpvt(j)
c           while
   20       if (jpvt(k) .lt. 0) then
               call dswap (nobs,x(1,j),1,x(1,k),1)
               jpvt(k) = -jpvt(k)
               j = k
               k = jpvt(k)
               goto 20
c           endwhile
            endif
   30    continue
      endif
      if (job .ne. 0) then
c		backward permutation
         do 50 i = 1,npar 
            if (jpvt(i) .gt. 0) then
               goto 50
            endif
            jpvt(i) = -jpvt(i)
            j = jpvt(i)
c           while
   40       if (j .ne. i) then
               call dswap (nobs,x(1,i),1,x(1,j),1)
               jpvt(j) = -jpvt(j)
               j = jpvt(j)
               goto 40
c           endwhile
            endif
   50    continue
      endif
      return
      end
