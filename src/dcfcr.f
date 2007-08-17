      subroutine dcfcr (fg,ldfg,nnull,qr,ldqr,npar,pmh,qraux,
     * sgpvt,capz,ldcapz,svals,npsing,v,ldv,nlamht,w1,z,coef,penlty,
     * work,info)
      integer ldfg,nnull,ldqr,npar,pmh,sgpvt(npar),ldcapz,npsing,ldv,
     * info
      double precision fg(ldfg,nnull),qr(ldqr,pmh),
     * qraux(pmh),capz(ldcapz,*),svals(npsing),v(ldv,npsing),nlamht,
     * w1(nnull),z(npsing),coef(npar),penlty,work(nnull)
c
c Purpose: determine the coefficients for a given value of nlamht
c	and vectors z and w1.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations
c   			that define f and g
c   ldfg		leading dimension of fg as declared in the 
c			calling	program
c   nnull		number of columns in g
c   qr(ldqr,pmh)	information on the Householder transformations 
c			that define q and r
c   ldqr		leading dimension of qr as declared
c			in the calling program
c   npar		number of rows in q
c   pmh			number of columns in r
c   qraux(pmh)		auxiliary information on the
c			qr Householder transformations
c   sgpvt(npar)		permuted indices from the pivoted Cholesky 
c			decomposition of the matrix which defines the 
c			semi-norm 
c   capz(ldcapz,pmh)	first part of the rotated design matrix
c   ldcapz		leading dimension of capz as declared
c			in the calling program
c   svals(npsing)	singular values 
c   npsing		number of positive singular values
c   v(ldv,npsing)	right singular vectors corresponding to svals
c   ldv			leading dimension of v as declared
c			in the calling program
c   nlamht		nobs*lambda hat
c   w1(nnull)		leading part of rotated response vector
c   z(npsing)		u'w2
c
c On Exit:
c   z(npsing)		g = [(D**2 +nlamht)**-1]Dz
c   coef(npar)		estimated coefficients
c   penlty		smoothness penalty which equals	gamma'gamma
c   info		error indicator
c			  0 : successful completion
c			  1 : error in dtrco, g is singular
c			  2 : error in dtrsl, r is singular
c Work arrays
c   work(nnull)		double precision work vector
c
c   Subprograms Called Directly:
c	Linpack - dqrsl dtrsl dtrco
c	Blas    - ddot dgemv
c	Other   - dprmut 
c
c $Header: dcfcr.f,v 2.100.1.1 86/10/07 12:47:09 lindstrom Exp $
c
      integer i
      double precision dummy(1),machpr,one,rcond
      double precision ddot
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c			form gamma-hat and penlty
      do 20 i = 1,npsing
         z(i) = z(i)*svals(i)/(svals(i)**2 + nlamht)
   20 continue
      call dgemv('N',pmh,npsing,1.0d0,v,ldv,z,1,0.0d0,coef,1)
      penlty = ddot(pmh,coef,1,coef,1)
c			form delta
      call dcopy(nnull,w1,1,coef(pmh+1),1)
      call dgemv('N',nnull,pmh,-1.0d0,capz,ldcapz,coef,1,1.0d0,
     *   coef(pmh+1),1)
c			check condition number of g
      if (nnull .ne. 0) then
          call dtrco(fg,ldfg,nnull,rcond,work,1)
          if (rcond .le. machpr*100) then
             info = 1
             return
          endif
          call dtrsl (fg,ldfg,nnull,coef(pmh+1),01,info)
      endif
c			form theta from gamma and delta
      call dtrsl (qr,ldqr,pmh,coef,11,info)
      if (info .ne. 0) then
         info = 2
         return
      endif
      call dqrsl (qr,ldqr,npar,pmh,qraux,coef,coef,dummy,dummy,dummy,
     * dummy,10000,info)
      call dprmut (coef,npar,sgpvt,1)
      return
      end
