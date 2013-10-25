      subroutine dtpss(des,lddes,nobs,dim,m,s,lds,ncov,y,ntbl,adiag,
     * lamlim,cost,dout,iout,coef,svals,tbl,ldtbl,auxtbl,work,lwa,
     * iwork,liwa,job,info)
      integer lddes,nobs,dim,m,lds,ncov,ntbl,iout(4),ldtbl,lwa,
     * liwa,iwork(liwa),job,info
      double precision des(lddes,dim),s(lds,*),y(nobs),
     * adiag(nobs),lamlim(2),dout(5),coef(*),svals(*),
     * tbl(ldtbl,3),auxtbl(3,3),work(lwa),cost
c
c Purpose: determine the generalized cross validation estimate of the 
c 	smoothing parameter and fit model parameters for a thin plate
c 	smoothing spline.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c   lddes		leading dimension of des as declared in calling
c   			program 
c   nobs		number of observations
c   dim			number of columns in des
c   m			order of the derivatives in the penalty
c   s(lds,ncov) 	design for the covariates. The covariates
c			must duplicate the replication structure of des.
c			See dptpss to handle covariates which do not.
c   lds			leading dimension of s as declared in calling
c   			program 
c   ncov		number of covariates 
c   y(nobs)		response vector
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat 
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   adiag(nobs)	 	"true" y values on entry if predictive mse is 
c			requested
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda)	scale) if user input limits are 
c			requested if lamlim(1) = lamlim(2) then lamhat
c			is set to (10**lamlim(1))/nobs
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abdc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then adiag will be calculated
c			if d is nonzero then there are replicates in the
c			   design
c
c On Exit:
c   des(lddes,dim)	sorted unique rows of des if job indicates that
c			there are replicates otherwise not changed
c   s(lds,ncov)		unique rows of s sorted to correspond to des
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(5)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss	    residual sum of squares
c			4  tr(I-A)  trace of I - A
c   			5  ssqrep   sum of squares for replication
c   iout(4)		contains:
c			1  npsing   number of positive singular
c				    values (npsing = nuobs - ncts).
c				    if info indicates nonzero info in 
c				    dsvdc then npsing contains info as 
c				    it was returned from dsvdc.
c			2  npar	    number of parameters
c				    (npar = nuobs + ncts)
c			3  ncts     dimension of the polynomial space 
c				    plus ncov	
c				    ((m+dim-1 choose dim) + ncov)
c			4  nuobs    number of unique rows in des
c   coef(npar)		coefficient estimates [beta':alpha':delta']'
c			coef must have a dimension of at least nuobs+
c			ncts
c   svals(npar-nnull)	singular values of the matrix j2 if info = 0
c			if info indicates nonzero info from dsvdc then 
c			svals is as it was returned from dsvdc.
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
c			  2 : error in dreps, covariates do not 
c			      duplicate the replication structure of des
c			  3 : lwa (length of work) is too small
c			  4 : liwa (length of iwork) is too small
c			  10 < info < 20 : 10 + nonzero info returned 
c					   from dsetup
c			  100< info <200 : 100 + nonzero info returned
c					   from dsgdc1
c			  200< info <300 : 200 + nonzero info returned
c					   from dgcv1
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			Must be at least nuobs(2+ncts+nuobs)+nobs
c   iwork(liwa)		integer work vector
c   liwa		length of iwork as declared in the calling 
c			program
c			Must be at least 2*nobs + nuobs - ncts
c
c Subprograms Called Directly:
c       Gcvpack - dreps duni dsuy dsetup dsgdc1 dgcv1 
c       Other   - dprmut mkpoly
c
c Subprograms Called Indirectly:
c	Gcvpack - dcfcr1 drsap dvlop dsvtc dpdcr dpmse 
c		  dvmin dvl dmaketg dmakekg ddiag
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc
c	Blas    - ddot dcopy dgemv
c	Other 	- dprmut dset dftkf fact mkpoly 
c
c $Header: dtpss.f,v 2.100.1.1 86/10/07 12:57:52 lindstrom Exp $
c
      integer ncts,jrep,jadiag,p1,p1p1,p2,p2p1,p3,p3p1,ip1,ip1p1,ip2,
     * ip2p1,nuobs,p4,p4p1,j,i,lwa2,wsize,q1,q1p1,npsing,
     * jobgcv
      double precision ssqrep
      integer mkpoly
c
c
      info = 0
      ncts = mkpoly(m,dim) + ncov
      iout(3) = ncts
      jrep = mod(job,10)
      jadiag = mod(job,100)/10
      jobgcv = job/10
c			check dimensions
      if((nobs .le. 0) .or. (m .le. 0) .or. (dim .le. 0) .or.
     * (ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. 
     * (2*m-dim .le. 0)) then
         info = 1
         return
      endif
c			set up pointers for iwork vector
c  			first nobs positions of iwork contain order
      ip1 = nobs
c			next nobs positions of iwork contain xrep    
      ip1p1 = ip1 + 1
      ip2 = ip1 + nobs
c			rest of iwork is a integer work vector
      ip2p1 = ip2 + 1
      if (jrep .ne. 0) then
         call dreps(des,lddes,nobs,dim,s,lds,ncov,0,work,iwork,nuobs,
     *    iwork(ip1p1),1,info)
         if (info .ne. 0) then
            info = 2
            return
         endif
c			put unique row of des into des
         call duni(des,lddes,nobs,dim,iwork(ip1p1),des,lddes)
c			put unique row of s into s
         call duni(s,lds,nobs,ncov,iwork(ip1p1),s,lds)
      endif
      if (jrep .eq. 0) then
         work(1)= 0
         nuobs = nobs
         ssqrep = 0.0d0
      endif
      iout(2) = nuobs + ncts
      iout(4) = nuobs
c			check size of work vectors
      wsize = nuobs*(2+ncts+nuobs)+nobs
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      wsize = 2*nobs + nuobs - ncts
      if (liwa .lt. wsize) then
         info = 4
         return
      endif
c			set up pointers for c1,[tu:su1],ku,fgaux in work
c
c		c1  	 runs from 1    to p1,   p1 = nuobs
c		[tu:su1] runs from p1p1 to p2,   p2 = p1 + nuobs*ncts
c		ku	 runs from p2p1 to p3,   p3 = p2 + nuobs**2
c		fgaux    runs from p3p1 to p4,   p4 = p3 + ncts
c		the rest of work is a work vector

c  			after the call to dsetup 
c		f2'k f2  runs from q1p1 to p3, q1 = p2 + nuobs*ncts+ncts
c
      p1 = nuobs
      p1p1 = p1 + 1
      p2 = p1 + nuobs*ncts
      p2p1 = p2 + 1
      q1 = p2 + nuobs*ncts+ncts
      q1p1 = q1 + 1
      p3 = p2 + nuobs**2
      p3p1 = p3 + 1
      p4 = p3 + ncts
      p4p1 = p4 + 1
      lwa2 = lwa - (nuobs*(1+ncts+nuobs) + ncts)
c			set up structures needed for dsgdc1 and dgcv1
      call dsetup(des,lddes,s,lds,dim,m,ncov,nuobs,work,work(p1p1),
     * nuobs,ncts,work(p3p1),work(p2p1),nuobs,work(p4p1),iwork(ip2p1),
     * info)
      if (info .ne. 0) then
         info = info + 10
         return
      endif

c			decompose f2' k f2 
      npsing = nuobs - ncts
      call dsgdc1(work(q1p1),nuobs,npsing,svals,iwork(ip2p1),
     * work(p4p1),lwa2,info)
      iout(1) = npsing
      if (info .gt. 0) then
         info = info + 100
         return
      endif
c			setup y
      if (jrep .ne. 0) then
         call dsuy(y,nobs,nuobs,adiag,work,iwork,iwork(ip1p1),
     *    ssqrep,jadiag)
      endif
      dout(5) = ssqrep
c			compute lambda hat and other parameters
      call dgcv1(work(p2p1),nuobs,y,nuobs,nobs,work(p1p1),nuobs,ncts,
     * work(p3p1),svals,adiag,lamlim,cost,ssqrep,ntbl,dout,coef,tbl,
     * ldtbl,auxtbl,work(p4p1),lwa2,jobgcv,info)
      if (info .gt. 0) then
         info = info + 200
         return
      endif
c			if there are replicates then rescale the coef. 
c			vector,	the predicted values, and diagonal of A
      if (nuobs .ne. nobs) then
         do 10 i = 1,nuobs
            coef(ncts+i) = coef(ncts+i) * work(i)
 10      continue
         j = nuobs
         do 20 i = nobs,1,-1 
            y(i)=y(j)/work(j)
            if (jadiag .ne. 0) then
               adiag(i) = adiag(j)/work(j)**2
            endif
            if (iwork(ip1 + i) .eq. 0) then
               j = j - 1
            endif
 20      continue
      endif
c			undo the permutation of the predicted values and
c			adiag if necessary
      if (jrep .ne. 0) then
         call dprmut (y,nobs,iwork,1)
         if (jadiag .ne. 0) call dprmut (adiag,nobs,iwork,1)
      endif
      return
      end
