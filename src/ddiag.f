      subroutine ddiag(fg,ldfg,nobs,nnull,fgaux,svals,npsing,u,ldu,
     * nlamht,adiag,work)
      integer ldfg,nobs,nnull,npsing,ldu
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),nlamht,adiag(nobs),work(*)
c
c Purpose: determine the diagonal of hat matrix for nobs*lamhat
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
c
c On Exit:
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c
c Work Arrays:
c   work(nobs+npsing)	double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl 
c	Blas    - ddot dgemv
c	Other   - dset
c
c $Header: ddiag.f,v 2.100.1.1 86/10/07 13:28:57 lindstrom Exp $
c
      integer i,j,hp1,locinf,nmh,np1
      double precision dummy(1)
      double precision ddot
c
      np1 = nobs + 1
      hp1 = nnull + 1
      nmh = nobs - nnull
c			form adiag
      do 20 i = 1,nobs 
         call dset(nobs,0.0d0,work,1)
         work(i)=1.0d0
         call dqrsl(fg,ldfg,nobs,nnull,fgaux,work,dummy,work,dummy, 
     *    dummy,dummy,01000,locinf)
         adiag(i)=ddot(nnull,work,1,work,1)
         call dgemv('T',nmh,npsing,1.0d0,u,ldu,work(hp1),1,0.0d0,
     *    work(np1),1)
         do 10 j=1,npsing
            work(nobs+j)=work(nobs+j)*svals(j)/dsqrt(svals(j)**2+nlamht)
   10    continue
         adiag(i)=adiag(i) + ddot(npsing,work(np1),1,work(np1),1)
   20 continue
      return
      end
