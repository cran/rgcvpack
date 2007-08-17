      subroutine dptpss(des,lddes,nobs,dim,m,s,lds,ncov1,ncov2,y,ntbl,
     * adiag,lamlim,cost,dout,iout,coef,svals,tbl,ldtbl,auxtbl,work,
     * lwa,iwork,liwa,job,info)
      integer lddes,nobs,dim,m,lds,ncov1,ncov2,ntbl,iout(4),ldtbl,lwa,
     * liwa,iwork(liwa),job,info
      double precision des(lddes,dim),s(lds,*),y(nobs),adiag(nobs),
     * lamlim(2),dout(4),coef(*),svals(*),tbl(ldtbl,3),
     * auxtbl(3,3),work(lwa),cost
c
c Purpose: determine the generalized cross validation estimate of the 
c	smoothing parameter and fit model parameters for a partial thin
c	plate spline model.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c   lddes		leading dimension of des as declared in the
c   			calling program 
c   nobs		number of observations
c   dim			number of columns in des
c   m			order of the derivatives in the penalty
c   s(lds,ncov1+ncov2)	design for the covariates
c			first ncov1 columns contain covariates which
c			  duplicate the replication structure of des
c			next ncov2 columns contain covariates which
c			  do not duplicate the replication structure of
c			  des
c   lds			leading dimension of s as declared in the
c   			calling program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of des
c   ncov2		number of covariates which do not duplicate the 
c			replication structure of des
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
c			requested. if lamlim(1) = lamlim(2) then lamhat 
c			is set to (10**lamlim(1))/nobs
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job			integer with decimal expansion abc
c			if a is nonzero then predictive mse is computed
c			   using adiag as true y
c			if b is nonzero then user input limits on search
c			   for lambda hat are used
c			if c is nonzero then adiag will be calculated
c On Exit:
c   des(lddes,dim)	unique rows of des
c   y(nobs)		predicted values
c   adiag(nobs)		diagonal elements of the hat matrix if requested
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   dout(4)		contains:
c  			1  lamhat   generalized cross validation 
c				    estimate of the smoothing parameter
c			2  penlty   smoothing penalty
c			3  rss      residual sum of squares
c			4  tr(I-A)  trace of I - A
c   iout(4)		contains:
c			1  npsing   number of positive singular values
c				    if info indicates nonzero info 
c				    from dsvdc then npsing contains
c				    info as it was returned from dsvdc
c			2  npar	    number of parameters 
c				    (npar = nuobs + nnull)
c			3  nnull    size of the null space of sigma
c				    (m+dim-1 choose dim)+ncov1+ncov2
c			4  nuobs    number of unique rows in des
c   coef(npar)		coefficient estimates [beta':alpha':delta']'
c			coef must have a dimension of at least 
c			nuobs+nnull
c   svals(npsing)	singular values, svals must have a dimension,
c			of at least nuobs-nnull.
c			if info indicates nonzero info in dsvdc then 
c			svals is as returned from dsvdc.
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
c			  -1 : log10(nobs*lamhat) <= lamlim(1) 
c			       (not fatal)
c			  -2 : log10(nobs*lamhat) >= lamlim(2) 
c			       (not fatal)
c			   1 : dimension error 	
c			   2 : error in dreps, the first ncov1 columns
c			       of s do not duplicate the replication 
c			       structure of des
c			   3 : lwa (length of work) is too small
c			   4 : liwa (length of iwork) is too small
c			   5 : error in dmaketg
c			   6 : sigma is rank deficient
c			  1000< info : 1000 + nonzero info returned from
c				       dsnsm
c
c Working Storage:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling 
c			program 
c			must be at least lwa1 + lwa2 where 
c			  lwa1 = (nnull-ncov2)*(nobs+nuobs+1)
c				 +npar*(nobs+npar)
c			  lwa2 = (npar-nnull)*(npar-2*nnull+2+nobs)
c			 	 +npar+nobs
c			
c			
c   iwork(liwa)		integer work vector
c   liwa		length of the iwork as declared in the calling 
c			program
c			must be at least 3*nobs - (nnull - ncov2)
c
c Subprograms Called Directly:
c 	Gcvpack - dreps dmaketg duni dmakekg dctsx dsnsm
c	Linpack - dqrdc dqrsl 
c	Blas    - dcopy
c	Other 	- dprmut dset prmut mkpoly
c
c Subprograms Called Indirectly:
c	Gcvpack - dcrtz ddcom dgcv dsgdc dtsvdc drsap ddiag
c		  dvlop dvlop dpmse dcfcr dpdcr dvmin dvl dzdc
c	Linpack - dchdc dqrdc dqrsl dtrsl dsvdc dtrco
c	Blas    - dcopy ddot dgemv dswap
c	Other 	- dcpmut dprmut dset dftkf fact mkpoly
c
c $Header: dptpss.f,v 2.100.1.2 86/11/21 12:22:18 lindstrom Exp $
c
      double precision dummy,dout2(5)
      integer i,ncts1,jadiag,p1,p1p1,p2,p2p1,p3,p3p1,p4,p4p1,p5,p5p1,
     * ip1,ip1p1,ip2,ip2p1,nuobs,lwa2,liwa2,npoly,nnull,snnpar,q1,
     * wsize,lwa1,sinfo,ldx,iout2(4)
      integer mkpoly
c
c
      info = 0
      sinfo = 0
      ncts1 = mkpoly(m,dim) + ncov1
      nnull = ncts1 + ncov2
      iout(3) = nnull
      jadiag = mod(job,10)
c			check dimensions
      if((nobs .le. 0) .or. (m .le. 0) .or. (dim .le. 0) .or.
     * (ntbl .lt. 0) .or. (ntbl .gt. ldtbl) .or. 
     * (2*m - dim .le. 0)) then
         info = 1
         return
      endif
c			set up pointers for iwork
c  			first nobs positions of iwork contain order
      ip1 = nobs
c			next nobs positions of iwork contain xrep
      ip1p1 = ip1 + 1
      ip2 = ip1 + nobs
c			rest of iwork is an integer work vector
      ip2p1 = ip2 + 1
      liwa2 = liwa - 2*nobs
      call dreps(des,lddes,nobs,dim,s,lds,ncov1,ncov2,dummy,iwork,nuobs,
     * iwork(ip1p1),0,info)
      iout(2) = nuobs+ncts1+ncov2
      iout(4) = nuobs
      if (info .ne. 0) then
         info = 2
         return
      endif
c			undo permutation of des, s and xrep
      do 10 i = 1,dim
          call dprmut(des(1,i),nobs,iwork,1)
   10 continue
      do 20 i = 1,ncov1+ncov2
          call dprmut(s(1,i),nobs,iwork,1)
   20 continue
      call prmut(iwork(ip1p1),nobs,iwork,1)
c
c			check size of work vectors
      snnpar = nuobs + ncov2
      lwa1 = ncts1*(nobs+nuobs+1)+snnpar*(nobs+snnpar)
      lwa2 = (snnpar-nnull)*(snnpar-2*nnull+2+nobs)+snnpar+nobs
      wsize = lwa1+lwa2
      if (lwa .lt. wsize) then
         info = 3
         return
      endif
      wsize = 3*nobs - ncts1
      if (liwa .lt. wsize) then
         info = 4
         return
      endif
c			set up pointers for work
c
c	     name     runs from  to		
c	     ----     ---------  --
c	     [t:s1]      1	 p1	  p1 = nobs*ncts1	
c	     [tu:s1u]    p1p1    p2       p2 = p1 + nuobs*ncts1
c	      x	         p2p1    p3       p3 = p2 + max(nobs,snnpar)
c								*snnpar
c	      sigma      p3p1    p4	  p4 = p3 + snnpar**2
c	      fgaux      p4p1    p5	  p5 = p4 + ncts1
c	      working    p5p1    p5+(snnpar-nnull)*(snnpar-2+nnull+
c					2+nobs)+snnpar+nobs
c
c		after the call to duni ku runs from q1 to p4
c
      p1 = nobs*ncts1
      p1p1 = p1 + 1
      p2 = p1 + nuobs*ncts1
      p2p1 = p2 + 1
      ldx=max(nobs,snnpar)
      p3 = p2 + ldx*snnpar
      p3p1 = p3 + 1
      q1 = p3 + ncov2*snnpar+ncov2+1
      p4 = p3 + snnpar**2
      p4p1 = p4 + 1
      p5 = p4 + ncts1
      p5p1 = p5 + 1
      lwa2 = lwa - ( ncts1*(nobs+nuobs+1) + snnpar*(nobs+snnpar))
c			make [t:s1] and k
      if (m .ne. 1) then
          call dmaketg(m,nobs,dim,des,lddes,s,lds,ncov1,npoly,work,nobs,
     *     iwork(ip2p1),info)
          if (info .ne. 0) then
             info = 5
             return
          endif
c			put unique rows of des into des
          call duni(des,lddes,nobs,dim,iwork(ip1p1),des,lddes)
c			put k into x
          call dmakekg(m,nobs,dim,work(1+nobs),nobs,nuobs,des,lddes,
     *     work(p2p1),ldx)
      else
c			put unique rows of des in 1st nuobs 
c			positions of work
          call duni(des,lddes,nobs,dim,iwork(ip1p1),work,1)
c			put k into x
          call dmakekg(m,nobs,dim,work,nobs,nuobs,des,lddes,
     *     work(p2p1),ldx)
c			copy unique rows of des from work to des 
          call dcopy(nuobs,work,1,des,1)
c			set t = 1
	  call dset(nobs,1.0d0,work,1)
      endif
c			put unique rows of k into (ncov2+1,ncov2+1)
c			position of sigma
      call duni(work(p2p1),ldx,nobs,nuobs,iwork(ip1p1),work(q1),snnpar)
c			compute unique rows of [t:s1] 
      call duni(work(1),nobs,nobs,ncts1,iwork(ip1p1),work(p1p1),nuobs)
c			make x and sigma
      call dctsx(work(1),nobs,ncts1,s(1,ncov1+1),lds,ncov2,nobs,
     * work(p1p1),nuobs, nuobs,work(p3p1),snnpar,work(p2p1),ldx,snnpar,
     * work(p4p1),work(p5p1))
c
      call dsnsm(work(p2p1),ldx,y,work(p3p1),snnpar,nobs,snnpar,nnull,
     * adiag,1.0d0,lamlim,cost,ntbl,dout2,iout2,coef,svals,tbl,ldtbl,
     * auxtbl,iwork(ip2p1),liwa2,work(p5p1),lwa2,job,info)
      call dcopy(4,dout2,1,dout,1)
      iout(3) = nnull
      iout(1) = iout2(1)
      if (info .eq. -3) then
	  info = 6
	  return
      endif
      if (info .gt. 0) then
	 info = info + 1000
	 return
      endif
      if (info .lt. 0) sinfo = info
c			form psi = f2 zeta
      call dset(ncts1,0.0d0,work(p5p1),1)
      call dcopy(nuobs-ncts1,coef(ncts1+ncov2+1),1,work(p5+ncts1+1),1)
      call dqrsl(work(p1p1),nuobs,nuobs,ncts1,work(p4p1),work(p5p1),
     * coef(ncts1+ncov2+1),dummy,dummy,dummy,dummy,10000,info)
      if (sinfo .lt. 0) info = sinfo
      return
      end
