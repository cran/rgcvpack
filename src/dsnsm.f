      subroutine dsnsm (x,ldx,y,sigma,ldsigm,nobs,npar,nnull,adiag,
     * tau,lamlim,cost,ntbl,dout,iout,coef,svals,tbl,ldtbl,auxtbl,
     * iwork,liwa,work,lwa,job,info)
      integer ldx,ldsigm,nobs,npar,nnull,ntbl,iout(3),ldtbl,liwa,
     * iwork(liwa),lwa,job,info
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * adiag(nobs),tau,lamlim(2),dout(5),coef(npar),svals(*),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa),cost
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized 
c	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	symmetric matrix that defines the semi-norm
c   ldsigm		leading dimension of sigma as declared
c			in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   adiag(nobs)		"true" y values on entry if computation of 
c			predictive mse is requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   tau			multiplier controlling the amount of truncation
c			if truncation is requested (try tau = 1 
c			to start then try 10 and 100)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abcd
c			if a is nonzero then truncation is used
c			if b is nonzero then predictive mse is computed
c			   using adiag as true y
c			if c is nonzero then user input limits on search
c			   for lambda hat are used
c			if d is nonzero then the diagonal of the hat 
c			   matrix is calculated
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c   y(nobs)		predicted values
c   sigma(ldsigm,npar)	overwritten with the QR decomposition of the 
c			Cholesky factor of sigma
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(5)		contains:
c  			1 lamhat   generalized cross validation 
c				   estimate of the smoothing parameter
c			2 penlty   smoothing penalty
c			3 rss	   residual sum of squares
c			4 tr(I-A)  trace of I - A
c			5 truncation ratio = 1/(1+(normk/(nobs*lamhat)))
c				   where normk = norm(R - R sub k)**2
c   iout(3)		contains:
c			1  npsing   number of positive singular
c				    values
c				    if info indicates nonzero info in 
c				    dsvdc then iout(1) contains info as
c				    it was returned from dsvdc
c			2  npar	    number of parameters
c			3  nnull    size of the null space of sigma
c   coef(npar)		coefficient estimates
c   svals(npar-nnull)	first npsing entries contain singular values 
c			of the matrix j2 
c			if info indicates nonzero info in dsvdc then 
c			svals is as it was returned from dsvdc
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c			  3     R(lambda) if requested
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lamhat), V(lamhat) and  
c			    R(lamhat) if requested
c			    where lamhat is the gcv estimate of lambda
c			2nd row contains:
c			    0, V(0) and  R(0) if requested
c			3rd row contains:
c			    0, V(infinity) and R(infinity) if requested
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			  -2 : log10(nobs*lamhat) >= lamlim(2) 
c			       (not fatal)
c			  -1 : log10(nobs*lamhat) <= lamlim(1)
c			       (not fatal)
c			   1 : dimension error 	
c			   2 : lwa (length of work) is too small
c			   3 : liwa (length of iwork) is too small
c			   4 : error in ntbl or tau
c			  100< info <200 : 100 + nonzero info returned
c					   from ddcom
c			  200< info <300 : 200 + nonzero info returned
c					   from dgcv
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			must be at least 
c			(npar-nnull)*(npar-2*nnull+2+nobs)+npar+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling 
c			program
c			must be at least 2*npar - nnull 
c
c Subprograms Called Directly:
c	Gcvpack - ddcom dgcv
c
c Subprograms Called Indirectly:
c	Gcvpack - dcrtz dsgdc dcfcr drsap dvlop dtsvdc
c		  dpmse dvmin dvl dzdc dpdcr ddiag
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc dtrco
c	Blas    - dcopy ddot dgemv dswap
c	Other 	- dcpmut dprmut dset 
c
c $Header: dsnsm.f,v 2.100.1.2 86/11/19 09:24:39 lindstrom Exp $
c
      integer bigp,pmh,pp1,lwa2,wsize,ldcaux,job1,job2,npsing,sinfo
      sinfo = 0
      info = 0
      iout(2) = npar
c			check dimensions
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or.(pmh .le. 0)) then
         info = 1
         return
      endif
      wsize = pmh**2 +2*npar - nnull + pmh*(nobs - nnull) + pmh + nobs
      if (lwa .lt. wsize) then
         info = 2
         return
      endif
      wsize = npar + pmh
      if (liwa .lt. wsize) then
         info = 3
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. (tau .lt. 0)) then
         info = 4
         return
      endif
c 			first p positions of iwork contain sgpvt
c 			first bigp-1  positions of work contain dcaux 
c			iwork vector starts at pp1
      pp1 = npar + 1
c 	work vector starts at bigp
      bigp = pmh**2 + 2*npar - nnull + 1
      ldcaux = pmh**2 + 2*npar - nnull
      lwa2 = lwa - ldcaux
      job1=job/1000
      job2=mod(job,1000)
c			decompose sigma and design matrix
      call ddcom (x,ldx,sigma,ldsigm,nobs,npar,nnull,tau,npsing,
     * svals,iwork,work,ldcaux,dout(5),work(bigp),lwa2,iwork(pp1),pmh,
     * job1,info)
      iout(1) = npsing
      iout(3) = nnull
      if (info .gt. 0) then
         info = info + 100
         return
      endif
      if (info .lt. 0) sinfo = info
c			compute lambda hat and other parameters
      call dgcv(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,iwork,work,ldcaux,
     * adiag,lamlim,cost,ntbl,dout,coef,tbl,ldtbl,auxtbl,work(bigp),
     * lwa2,job2,info)
      if (info .gt. 0) then
         info = info + 200
         return
      endif
      dout(5) = 1.0d0/(1.0d0+dout(5)/nobs/dout(1))
      if (sinfo .lt. 0) info = sinfo
      return
      end
