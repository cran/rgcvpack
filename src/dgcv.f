      subroutine dgcv(x,ldx,y,sigma,ldsigm,nobs,npar,nnull,sgpvt,dcaux,
     * ldcaux,adiag,lamlim,cost,ntbl,dout,coef,tbl,ldtbl,auxtbl,work,
     * lwa,job,info)
      integer ldx,ldsigm,nobs,npar,nnull,sgpvt(npar),ldcaux,ntbl,ldtbl,
     * lwa,info,job
      double precision x(ldx,npar),y(nobs),sigma(ldsigm,npar),
     * dcaux(ldcaux),adiag(nobs),lamlim(2),dout(4),coef(npar),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa),cost
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a penalized
c 	least squares problem with a semi-norm smoothing matrix.
c
c On Entry:
c   x(ldx,npar)		design matrix as returned by ddcom
c   ldx			leading dimension of x as declared in the 
c			calling program, must be at least max(nobs,npar)
c   y(nobs)		response vector
c   sigma(ldsigm,npar)	Cholesky factor of sigma as returned by ddcom
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   nobs		number of observations
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c   sgpvt(npar)		permuted indices from the pivoted Cholesky
c			decomposition of the semi-norm matrix
c   dcaux(ldcaux)	auxiliary vector from ddcom
c   ldcaux		length of ldcaux. Must be at least
c			(npar-nnull)**2+2*npar-nnull
c   adiag(nobs)	 	"true" y values on entry if predictive mse if 
c			requested 
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then diagonal of the hat matrix
c			   is calculated
c
c On Exit:
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   coef(npar)		coefficient estimates
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
c			  0 : successful completion
c			 -1 : log10(nobs*lamhat) <= lamlim(1) 
c			      (not fatal)
c			 -2 : log10(nobs*lamhat) >= lamlim(2) 
c			      (not fatal)
c			  1 : dimension error
c			  2 : error in ntbl
c			  3 : ldcaux (length of dcaux) is too small
c			  4 : lwa (length of work) is too small
c			  5 : lamlim(1) > lamlim(2)
c			 10 < info < 20 : 10 + nonzero info returned
c					  from dvlop
c			 20 < info < 30 : 20 +nonzero info returned
c					  from dcfcr
c
c Work Arrays:
c   work(lwa) 		double precision work vector
c   lwa			length of work as declared in the calling 
c			program must be at least (npar-nnull)+nobs
c
c Subprograms Called Directly:
c	Gcvpack - drsap dvlop dpmse dcfcr dpdcr ddiag
c
c Subprograms Called Indirectly:
c	Gcvpack - dvmin dvl
c	Linpack - dqrsl dtrsl
c	Blas    - dcopy ddot dgemv
c	Other   - dprmut dset
c
c $Header: dgcv.f,v 2.100.1.1 86/10/07 12:49:16 lindstrom Exp $
c
      double precision addend,nlamht,ssqw2
      integer mp1,mnh,pmh,minnp,pp1,pmhp1,wsize,npsing,i,jpmse,jlaml,
     * jadiag,sinfo
c
      jpmse = job/100
      jlaml = mod(job,100)/10
      jadiag = mod(job,10)
c			check dimensions
      sinfo = 0
      info = 0
      pmh = npar - nnull
      if ((nobs .le. 0) .or. (npar .le. 0) .or. (nnull .lt. 0) 
     * .or. (pmh .le. 0)) then
         info = 1
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl)) then
         info = 2
         return
      endif
      wsize = pmh**2 + 2*npar - nnull
      if (ldcaux .lt. wsize) then
         info = 3
         return
      endif
      wsize = pmh + nobs
      if (lwa .lt. wsize) then
         info = 4
         return
      endif
      if (jlaml .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 5
         return
      endif
      pmhp1 = pmh + 1
      pp1 = npar + 1
c			calculate npsing
      minnp = min(nobs-nnull,pmh)
      do 30 i=1,minnp
         if (dcaux(npar+i)**2 .gt. 0) npsing = i
   30 continue
c			apply rotations to y
      mp1 = nnull + 1
      mnh = nobs - nnull
      call drsap (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),x(mp1,1),ldx,
     * mnh,npsing,y,ssqw2,addend,work)
c			minimize V(lambda)
      call dvlop (y(mp1),dcaux(pp1),nobs,nnull,npsing,cost,addend,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout(3),jlaml,info)
      dout(1)=nlamht/nobs
      if (info .gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0) sinfo = info
c			calculate predictive mse
      if (jpmse .ne. 0) then
         call dpmse(x(1,pmhp1),ldx,nobs,nobs,nnull,dcaux(pmhp1),
     *    dcaux(pp1),npsing,x(mp1,1),ldx,y,y(mp1),ntbl,adiag,tbl,ldtbl,
     *    auxtbl,work)
      endif
c			calculate coefficients
      call dcfcr (x(1,pmhp1),ldx,nnull,sigma,ldsigm,npar,
     * pmh,dcaux,sgpvt,x,ldx,dcaux(pp1),npsing,dcaux(npar+pmhp1),pmh, 
     * nlamht,y,y(mp1),coef,dout(2),work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
c			calculate predicted values 
      call dpdcr (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),dcaux(pp1),
     * npsing,x(mp1,1),ldx,nlamht,y,y(mp1),y,work)
      if (jadiag .ne. 0) then
          call ddiag (x(1,pmhp1),ldx,nobs,nnull,dcaux(pmhp1),
     *        dcaux(pp1),npsing,x(mp1,1),ldx,nlamht,adiag,work)
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
