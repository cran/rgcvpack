      subroutine drsap(fg,ldfg,nobs,nnull,fgaux,u,ldu,nmh,npsing,z,
     * ssqw2,addend,work)
      integer ldfg,nobs,nnull,ldu,nmh,npsing
      double precision fg(ldfg,nnull),fgaux(nnull),u(ldu,npsing),
     * z(nobs),ssqw2,addend,work(npsing)
c
c Purpose: apply Householder transformations to a response vector and
c	collect its inner product with u and the addend which are used 
c	to define the generalized cross validation function with a
c	semi-norm.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations 
c			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling program
c   nobs		number of rows in fg
c   nnull		number of columns in fg
c   fgaux(nnull)	auxiliary information on the fg	Householder 
c			transformations
c   u(ldu,npsing)		left singular vectors 
c   ldu			leading dimension of u as declared in the
c			calling program
c   nmh	     		number of rows in u. nmh = nobs - nnull
c   npsing		number of columns in u (maximum of npar - nnull)
c   z(nobs)		response vector
c
c On Exit:
c   z(nobs)		the first nnull positions contain w1 and the 
c			next npsing positions contain u'w2
c   ssqw2		the squared length of w2
c   addend		the squared length of z minus the squared length
c			of u'w2
c
c Work Arrays:
c   work(npsing) 	double precision work vector
c
c Subprograms Called Directly:
c      Linpack - dqrsl 
c      Blas    - ddot dcopy dgemv
c
c $Header: drsap.f,v 2.100.1.1 86/10/07 12:55:34 lindstrom Exp $
c
      integer locinf,hp1
      double precision dummy(1)
      double precision ddot
c			apply Householder transformations
c			which define f
      call dqrsl (fg,ldfg,nobs,nnull,fgaux,z,dummy,z,dummy,dummy,dummy,
     * 01000,locinf)
c			w1 in first nnull positions of z,w2 in
c			last nmh
      hp1 = nnull + 1
      addend = ddot(nmh,z(hp1),1,z(hp1),1)
      ssqw2=addend
      call dgemv('T',nmh,npsing,1.0d0,u,ldu,z(hp1),1,0.0d0,work,1)
c			u'w2 in positions nnull+1 to
c			nnull+npsing of z
      call dcopy (npsing,work,1,z(hp1),1)
      addend = addend - ddot(npsing,z(hp1),1,z(hp1),1)
      return
      end
