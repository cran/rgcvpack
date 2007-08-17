      subroutine dgcv1(fkf,ldfkf,y,nuobs,nobs,fg,ldfg,ncts1,fgaux,svals,
     * adiag,lamlim,cost,ssqrep,ntbl,dout,coef,tbl,ldtbl,auxtbl,work,
     * lwa,job,info)
      integer ldfkf,nuobs,nobs,ldfg,ncts1,ntbl,ldtbl,lwa,job,info
      double precision fkf(ldfkf,nuobs),y(nuobs),fg(ldfg,ncts1),
     * fgaux(ncts1),svals(*),adiag(nuobs),lamlim(2),ssqrep,dout(4),
     * coef(*),tbl(ldtbl,3),auxtbl(3,3),work(lwa),cost
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a semi-norm
c	thin plate spline model.
c
c On Entry:
c   fkf(ldfkf,nuobs) 	intermediate results as created by dsgdc1
c   ldfkf		leading dimension of fkf as declared in the 
c			calling program
c   y(nuobs)		B1'y 
c   nuobs		number of rows in fg
c   nobs		number of observations
c   fg(ldfg,ncts1)	qr decomposition of [t:s1]	
c   ldfg		leading dimension of fg as
c			declared in the calling program
c   ncts1		number of columns in [t:s1]	
c   fgaux(ncts1)	auxiliary information on the decomposition of 
c			[t:s1]
c   svals(nuobs-ncts1) 	positive singular values of the Cholesky factor
c			of f2'k f2
c   adiag(nuobs)	B1'(true y) if predictive mse is requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ssqrep		sum of squares for replication
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
c   y(nuobs)		B1'(predicted values)
c   adiag(nuobs)	diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   coef(nuobs+ncts1) 	estimated coefficients
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
c			  3 : lwa (length of work) is too small
c			  4 : lamlim(1) > lamlim(2)
c			 10 < info < 20 : 10 + nonzero info returned 
c					  from dvlop
c			 20 < info < 30 : 20 + nonzero info returned
c					  from dcfcr1
c
c Working Storage:
c   work(lwa) 		double precision work vector
c   lwa			length of work as declared in the calling 
c			program
c			must be at least nuobs-ncts1+nobs
c
c Subprograms Called Directly:
c       Gcvlib - drsap dvlop dpmse dcfcr1 dpdcr ddiag
c
c Subprograms Called Indirectly:
c	Gcvlib  - dvl vmin
c	Linpack - dqrsl dtrsl
c	Blas    - ddot dcopy dgemv
c
c $Header: dgcv1.f,v 2.100.1.1 86/10/07 12:49:43 lindstrom Exp $
c
      double precision addend,nlamht,ssqw2
      integer npsing,i,jpmse,jlaml,jadiag,nmnct,nctp1,sinfo
c
c
      sinfo = 0
      info = 0
      nmnct = nuobs - ncts1
      nctp1=ncts1+1
      jpmse = job/100
      jlaml = mod(job,100)/10
      jadiag = mod(job,10)
c			check dimensions
      if ((nuobs.le.0).or.(ncts1 .le. 0).or.(nmnct .le. 0)) then
         info = 1
         return
      endif
      if ((ntbl .lt. 0) .or. (ntbl .gt. ldtbl)) then
         info = 2
         return
      endif
      if (lwa .lt. nobs+nuobs - nmnct) then
         info = 3
         return
      endif
      if (jlaml .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 4
         return
      endif
c			calculate npsing
      do 30 i=1,nmnct
         if (svals(i)**2 .gt. 0.0d0) npsing = i
   30 continue
c			apply rotations to y
      call drsap(fg,ldfg,nuobs,ncts1,fgaux,fkf(nctp1,nctp1),ldfkf,nmnct,
     * npsing,y,ssqw2,addend,work)
      addend = addend + ssqrep
      ssqw2 = ssqw2 + ssqrep
c			minimize V(lambda)
      call dvlop (y(nctp1),svals,nobs,ncts1,npsing,cost,addend,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout(3),jlaml,info)
      dout(1)=nlamht/nobs
      if (info .gt. 0) then
         info = info + 10
         return
      endif
      if (info .lt. 0)  sinfo = info
c			calculate predictive mse
      if (jpmse .ne. 0) then
         call dpmse(fg,ldfg,nuobs,nobs,ncts1,fgaux,svals,npsing,
     *    fkf(nctp1,nctp1),ldfkf,y,y(nctp1),ntbl,adiag,tbl,ldtbl,
     *    auxtbl,work)
      endif
c			calculate coefficients
      call dcfcr1(fg,ldfg,ncts1,fgaux,fkf(nctp1,nctp1),ldfkf,
     * fkf(1,nctp1),ldfkf,nuobs,svals,npsing,nlamht,y,y(nctp1),coef,
     * dout(2),work,info)
      if (info .gt. 0) then
         info = info + 20
         return
      endif
      call dpdcr (fg,ldfg,nuobs,ncts1,fgaux,svals,npsing,
     * fkf(nctp1,nctp1),ldfkf,nlamht,y,y(nctp1),y,work)
      if (jadiag .ne. 0) then
          call ddiag(fg,ldfg,nuobs,ncts1,fgaux,svals,npsing,
     *     fkf(nctp1,nctp1),ldfkf,nlamht,adiag,work)
      endif
      if (sinfo .lt. 0) info = sinfo
      return
      end
