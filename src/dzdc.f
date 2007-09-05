      subroutine dzdc(x,ldx,nobs,npar,pmh,tau,fgaux,svals,npsing,
     * v,ldv,normk,work,lwa,iwork,liwa,job,info)
      integer ldx,nobs,npar,pmh,npsing,ldv,lwa,iwork(liwa),liwa,info,job
      double precision x(ldx,npar),tau,fgaux(*),svals(*),v(ldv,pmh),
     * normk,work(lwa)
c
c Purpose: create and decompose j and decompose j2
c
c On Entry:
c   x(ldx,npar)		z = design matrix after sigma rotations 
c   ldx		     	leading dimension of x as declared in the 
c			calling	program
c   nobs   		number of observations
c   npar      		number of parameters
c   pmh    		npar minus the dimension of the null space of 
c			the semi-norm
c   tau			multiplier controlling the amount of truncation 
c			if truncation is requested
c   ldv    		leading dimension of v as declared in the 
c			calling	program
c   job			controls use of truncated singular value 
c			decomposition
c			   if job is nonzero then dtsvdc is used
c			   if job = 0 then dsvdc is called directly
c
c On Exit:
c   x(ldx,npar)		overwritten with many intermediate results
c   fgaux(npar-pmh)	auxiliary information on the qr decomposition 
c			of the matrix z2
c   svals(npsing)      	singular values of the matrix j2 if info = 0.
c			if info = 4,5,6 or 7 then svals is as 
c			it was returned from dsvdc.
c   npsing       	if info = 0 then npsing contains the number 
c			of singular values calculated 
c			if info = 4,5,6 or 7 then npsing contains info
c			as it was returned from dsvdc.
c   v(ldv,pmh) 		right singular vectors of the matrix j2
c   normk		frobenius norm of the k by k lower right sub 
c			matrix of j2
c   info    		error indicator
c			  0 : successful completion
c			  1 : lwa (length of work) is too small
c			  2 : tau < 0
c			  3 : transpose of j2 is necessary 
c				(npar .gt. nobs) and npar .gt. ldx
c			  4 : error in dtsvdc (using j2')
c			  5 : error in dsvdc (using j2')
c			  6 : error in dtsvdc
c			  7 : error in dsvdc
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa     		length of work as declared in the calling 
c			program	Must be at least 
c			pmh*(nobs - (npar-pmh) ) + pmh + nobs
c   iwork(liwa)		integer work vector
c   liwa		length of integer work vector, must be at least
c			pmh
c
c Subprograms Called Directly:
c	Gcvpack - dtsvdc 
c	Linpack - dqrdc dqrsl dsvdc
c	Blas    - dcopy dswap
c
c Subprograms Called Indirectly:
c	Linpack - dqrdc dsvdc dqrsl
c	Blas    - ddot dcopy
c	Other   - dset dprmut dswap
c
c $Header: dzdc.f,v 2.100.1.1 86/10/07 13:00:03 lindstrom Exp $
c
      double precision dummy,machpr,one,minrat
      integer i,j,idummy,nnull,hp1,nmh,pmhp1,wsize,locinf
c
      info = 0
c
      one = 1.0d0
      machpr = 1.0d0
 10   machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
      minrat = machpr*tau
c			check dimensions
      wsize = pmh**2 + pmh + nobs
      if (lwa .lt. wsize) then
         info = 1
         return
      endif
      if (tau .lt. 0) then
         info = 2
         return
      endif
c			z1 is the first pmh columns of x
c			z2 is next h columns
      pmhp1 = pmh + 1
      nnull = npar - pmh
      hp1 = nnull + 1
      nmh = nobs - nnull
c			decompose z2 into f and g
      call dqrdc (x(1,pmhp1),ldx,nobs,nnull,fgaux,idummy,dummy,0)
c			create j as f' z1 
      do 30 j = 1,pmh
         call dqrsl (x(1,pmhp1),ldx,nobs,nnull,fgaux,x(1,j),dummy,
     *    x(1,j),dummy,dummy,dummy,01000,locinf)
 30   continue
c			decompose j2 (last npar - nnull rows of j)
c			as udv'
      if (pmh .gt. nmh) then
c			transpose the j2 matrix
         if (npar .gt. ldx) then
            info = 3
            return
         endif
         do 40 i = 0, nmh-1
            call dcopy(pmh,x(hp1+i,1),ldx,work(i*pmh),1)
 40      continue
         do 50 i = 0, nmh-1
            call dcopy (pmh, work(i*pmh),1,x(hp1,i+1),1)
 50      continue
         if (job .ne. 0) then
            call dtsvdc (x(hp1,1),ldx,pmh,nmh,minrat,npsing,svals,
     *       x(hp1,1),ldx,v,ldv,normk,work,lwa,iwork,liwa,21,info)
            if (info .ne. 0) then
               npsing = info
               info = 4
               return
            endif
         else
            call dsvdc (x(hp1,1),ldx,pmh,nmh,svals,work(pmhp1),x(hp1,1),
     *       ldx,v,ldv,work,21,info)
            if (info .ne. 0) then
               npsing = info
               info = 5
               return
            endif
            npsing=nmh
            normk=0.0d0
         endif
         do 60 i = 1, npsing
            call dswap(pmh,x(hp1,i),1,v(1,i),1)
 60      continue
      else
         if (job .ne. 0) then
            call dtsvdc (x(hp1,1),ldx,nmh,pmh,minrat,npsing,svals,
     *       x(hp1,1),ldx,v,ldv,normk,work,lwa,iwork,liwa,21,info)
            if (info .ne. 0) then
               npsing = info
               info = 6
               return
            endif
         else
            call dsvdc (x(hp1,1),ldx,nmh,pmh,svals,work(nmh+1),x(hp1,1),
     *       ldx,v,ldv,work,21,info)
            if (info .ne. 0) then
               npsing = info
               info = 7
               return
            endif
            npsing=pmh
            normk=0.0d0
         endif
      endif
      return
      end
