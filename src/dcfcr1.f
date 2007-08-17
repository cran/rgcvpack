      subroutine dcfcr1(fg,ldfg,ncts1,fgaux,u,ldu,f1kf2,ldfkf,nuobs,
     * svals,npsing,nlamht,w1,z,coef,penlty,work,info)
      integer ldfg,ldfkf,ncts1,ldu,nuobs,npsing,info
      double precision fg(ldfg,ncts1),fgaux(ncts1),u(ldu,*),
     * f1kf2(ldfkf,*),svals(npsing),nlamht,w1(ncts1),z(npsing),coef(*),
     * penlty,work(*)
c
c Purpose: determine the coefficients for a given value of nlamht
c	and vectors z and w1.
c
c On Entry:
c   fg(ldfg,ncts1)	information on the Householder transformations
c   			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   ncts1		number of columns in g
c   fgaux(ncts1)	auxiliary information on the fg Householder
c			transformations
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu	    		leading dimension of u as declared in the 
c			calling	program
c   f1kf2(ldfkf,nuobs-ncts1) f1 k f2
c   ldfkf		leading dimension of f1kf2 as declared
c			in the calling program
c   nuobs		number of rows in fg
c   svals(npsing)	singular values of f2'k f2
c   npsing		number of positive singular
c   nlamht		nobs*(lambda hat)
c   w1(ncts1)		leading part of rotated response vector
c   z(npsing)		u'w2
c
c On Exit:
c   z(npsing)		g = [ (D**2 +nlamht)**-1 ] D z
c   coef(nuobs+ncts1)	estimated coefficients
c   penlty		smoothness penalty which equals	gamma'gamma
c   info		error indicator
c			  0 : successful completion
c			  1 : error in dtrco, g is singular
c
c Work Arrays:
c   work(nuobs-ncts1)  	double precision work vector
c
c   Subprograms Used:
c      Linpack - dqrsl dtrsl 
c      Blas    - ddot dcopy dgemv
c
c $Header: dcfcr1.f,v 2.100.1.1 86/10/07 12:47:31 lindstrom Exp $
c
      integer i,j,nmnct,nctp1,locinf
      double precision dummy,machpr,one,rcond
      double precision ddot
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      nmnct = nuobs - ncts1
      nctp1 = ncts1 + 1
c			form g and penalty
      do 20 j = 1,npsing
         z(j) = z(j)*svals(j)/(svals(j)**2 + nlamht)
   20 continue
      penlty = ddot(npsing,z,1,z,1)
c			z now contains g
c			compute xi 
      do 30 i = 1,npsing
         coef(i) = z(i)/svals(i)
   30 continue
      call dgemv('N',nmnct,npsing,1.0d0,u,ldu,coef,1,0.0d0,
     *  work,1)
      do 40 j = 1,ncts1
         coef(ncts1+j) = 0.0d0
   40 continue
      call dcopy(nmnct,work,1,coef(2*ncts1+1),1)
      call dqrsl(fg,ldfg,nuobs,ncts1,fgaux,coef(nctp1),coef(nctp1),
     * dummy,dummy,dummy,dummy,10000,locinf)
c			compute beta
      call dcopy(ncts1,w1,1,coef,1)
      call dgemv('N',ncts1,nmnct,-1.0d0,f1kf2,ldfkf,work,1,1.0d0,
     *  coef,1)
c			check condition number of g
      call dtrco(fg,ldfg,ncts1,rcond,work,1)
      if (rcond .le. machpr*100) then
         info = 1
         return
      endif
      call dtrsl (fg,ldfg,ncts1,coef,01,info)
      return
      end
