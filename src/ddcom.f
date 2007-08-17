      subroutine ddcom (x,ldx,sigma,ldsigm,nobs,npar,nnull,tau,
     * npsing,svals,sgpvt,dcaux,ldcaux,normk,work,lwa,iwork,liwa,
     * job,info)
      integer ldx,ldsigm,nobs,npar,nnull,npsing,sgpvt(npar),
     * ldcaux,lwa,liwa,iwork(liwa),job,info
      double precision x(ldx,npar),sigma(ldsigm,npar),tau,svals(*),
     * dcaux(ldcaux),normk,work(lwa)
c
c Purpose: decompose sigma and the design matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx			leading dimension of x as declared in the 
c			calling program. Must be at least max(npar,nobs)
c   sigma(ldsigm,npar)	symmetric matrix that defines the semi-norm
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   tau			multiplier controlling the amount of truncation
c			if truncation is requested
c   ldcaux		length of auxiliary vector as declared in the 
c			calling program. 
c			Must be at least (npar-nnull)**2+2*npar-nnull
c   job 		job nonzero use truncation 
c			job = 0 no truncation
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c			must be	passed to dgcv intact
c   sigma(ldsigm,npar)	overwritten with the qr decomposition of the 
c			Cholesky factor of sigma
c			must be passed intact to dgcv
c   nnull		if input nnull is too small it is replaced by
c			larger value such that sigma has rank npar-nnull
c   npsing       	number of singular values calculated if info = 0
c			otherwise npsing contains nonzero info from 
c			dsvdc.
c   svals(npar-nnull)	singular values of the matrix j2 if info = 0
c			if info from dsvdc is nonzero then svals 
c			is as it was returned from dsvdc.
c   sgpvt(npar)		permuted indices from the Cholesky decomposition
c			with pivots of sigma
c			must be passed intact to dgcv 
c   dcaux(ldcaux)	auxiliary vector, must be passed intact to dgcv 
c   normk		frobenius norm of the k by k lower right sub 
c			matrix of j2
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			   1 : dimension error or tau < 0
c			   2 : ldcaux (length of dcaux) is too small
c			   3 : lwa (length of work) is too small
c			   4 : liwa (length of iwork) is too small
c			   10 < info < 20 : 10 + nonzero info returned
c					    from dsgdc
c			   20 < info < 30 : 20 + nonzero info returned 
c					    from dcrtz
c			   30 < info < 40 : 30 + nonzero info returned 
c					    from dzdc
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program Must be at least 
c			(npar-nnull)*(nobs-nnull+1)+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling
c			program 
c			must be at least npar-nnull
c
c Subprograms Called Directly:
c	Gcvpack - dsgdc dcrtz dzdc
c	Blas	- dcopy
c
c Subprograms Called Indirectly:
c	Gcvpack - dtsvdc
c	Linpack - dtrco dchdc dqrdc dqrsl dtrsl dsvdc
c	Blas    - dcopy dswap ddot
c	Other   - dprmut dcpmut dset
c
c $Header: ddcom.f,v 2.100.1.1 86/10/07 12:48:26 lindstrom Exp $
c
      integer pmh,pp1,pmhp1,wsize,nmh,sinfo
c
      sinfo = 0
      info = 0
c			check dimensions
      nmh = nobs - nnull
      pmh = npar - nnull
      pmhp1 = pmh + 1
      pp1 = npar + 1
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or. (pmh .le. 0) .or. (ldx .lt. nobs).or. (ldx .lt. npar)
     * .or. (tau .lt. 0)) then
         info = 1
         return
      endif
      wsize = pmh**2 + 2*npar -nnull
      if (ldcaux .lt. wsize) then
         info = 2
         return
      endif
      wsize = pmh*nmh + pmh + nobs
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      if (liwa .lt. pmh) then
         info = 4
         return
      endif
c			decompose sigma
      call dsgdc (sigma,ldsigm,npar,nnull,dcaux,sgpvt,info)
      if (info.gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0) then
	  sinfo = info
	  pmh = npar - nnull
          pmhp1 = pmh + 1
      endif
c			create z
      call dcrtz(x,ldx,nobs,npar,pmh,sigma,ldsigm,dcaux,sgpvt,work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
c			create j and decompose j2
      call dzdc(x,ldx,nobs,npar,pmh,tau,dcaux(pmhp1),dcaux(pp1),
     * npsing,dcaux(npar+pmhp1),pmh,normk,work,lwa,iwork,liwa,job,info)
c			copy svals
      call dcopy(pmh,dcaux(pp1),1,svals,1)
      if (info .gt. 0) then
         info = info + 30
         return
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
