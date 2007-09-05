      subroutine dpdcr(fg,ldfg,nobs,nnull,fgaux,svals,npsing,u,ldu, 
     * nlamht,w1,g,pred,work)
      integer ldfg,nobs,nnull,npsing,ldu
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),nlamht,w1(nnull),g(npsing),pred(nobs),
     * work(*)
c
c Purpose: determine the predicted responses for a given value of 
c	nobs*lamhat and vectors g and w1.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder
c			transformations
c   svals(npsing)	singular values 
c   npsing		number of positive singular values 
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   nlamht		nobs*lambda hat
c   w1(nnull)		leading part of rotated response vector
c   g(npsing)		(D**2 + nlamht*I)*-1 Dz
c
c On Exit:
c   pred(nobs)		predicted responses
c
c Work Arrays:
c   work(nobs+npsing)	double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl 
c	Blas    - dcopy dgemv
c
c $Header: dpdcr.f,v 2.100.1.1 86/10/07 12:51:16 lindstrom Exp $
c
      integer i,locinf,nmh,np1
      double precision dummy(1)
c
c
      nlamht = 0.0d+0
      np1 = nobs + 1
      nmh = nobs - nnull
c			form the response vector
      call dcopy (nnull,w1,1,pred,1)
      call dcopy (npsing,g,1,work(np1),1)
      do 10 i = 1,npsing
         work(nobs+i) = work(nobs+i)*svals(i)
   10 continue
      call dgemv('N',nmh,npsing,1.0d0,u,ldu,work(np1),1,0.0d0,
     *  pred(nnull+1),1)
      call dqrsl (fg,ldfg,nobs,nnull,fgaux,pred,pred,dummy,dummy,dummy,
     * dummy,10000,locinf)
      return
      end
