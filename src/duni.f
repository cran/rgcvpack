      subroutine duni(x,ldx,nobs,ncx,xrep,xu,ldxu)
      integer ldx,nobs,ncx,ldxu,xrep(nobs)
      double precision x(ldx,*),xu(ldxu,*)
c
c Purpose: compute xu.
c
c On Entry:
c   x(ldx,ncx)		a matrix to be reduced to unique rows
c   ldx			leading dimension of x as declared in the
c			calling	program 
c   nobs		number of observations
c   ncx			number of columns in x
c   xrep(nobs)		xrep(i) contains 1 if ith row of x is a 
c			replicate row, 0 if not
c   ldxu		leading dimension of xu as declared in the 
c			calling	program 
c On Exit:
c   xu(ldxu,ncx) 	unique rows of x
c			may be identified with x in the calling sequence
c
c $Header: duni.f,v 2.100.1.1 86/10/07 12:59:15 lindstrom Exp $
c
      integer i,j,k
c
      j = 0
      do  20 i = 1,nobs
         if (xrep(i) .eq. 0) then
            j = j + 1
            do 10 k = 1,ncx
               xu(j,k) = x(i,k)
 10         continue
         endif
 20   continue
      return 
      end
