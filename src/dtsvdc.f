      subroutine dtsvdc(x,ldx,n,p,minrat,k,s,u,ldu,v,ldv,normk,
     * work,lwa,iwork,liwa,job,info)
      integer ldx,n,p,k,ldu,ldv,lwa,iwork(liwa),liwa,job,info
      double precision x(ldx,p),minrat,s(p),u(ldu,*),v(ldv,p),
     * normk,work(lwa)
c
c Purpose: form the singular value decomposition of a truncated matrix 
c	obtained by first taking a qr decomposition with pivoting.
c
c On Entry:
c   x(ldx,p)   		matrix to be decomposed, x is destroyed
c   ldx     		leading dimension of x in the calling program
c   n       		number of rows in x
c   p       		number of columns in x   (p <= n)
c   minrat  		minimum ratio to determine truncation
c   ldu     		leading dimension of u as declared in the 
c			calling program
c   ldv     		leading dimension of v as declared in the 
c			calling program
c   job     		controls the computation of the singular vectors
c               	it has the decimal expansion ab with the meaning
c                            a = 0 no left singular vectors
c                            a = 1 all n left singular vectors
c                            a > 1 first p left singular vectors
c                            b = 0 no right singular vectors
c                            b > 0 all p right singular vectors
c
c On Exit:
c   k       		number of positive singular values 
c   s(p)		an approximation to the first k singular values
c			in decreasing order
c   u(ldu,m)   		first k singular vectors if job > 1 or all n 
c			left singular vectors if joba = 1 otherwise it 
c			is not accessed it may be identified with x in 
c			the subroutine call if joba > 1
c   v(ldv,p)   		the first k right singular vectors if job b > 0
c               	otherwise it is not accessed
c   normk		norm of the k by k lower right sub matrix of r
c   info    		error indicator
c		   	   0 : successful completion
c		  	   info > 0 : info was returned from dsvdc
c
c Work Arrays:
c   iwork(liwa) 	integer work vector, holds the pivot vector from
c			dqrdc
c   liwa		must be at least p
c   work(lwa)		work vector
c   lwa			must be at least p**2 + p + n if joba > 0, 
c			otherwise it must be at least 2*p
c
c Subprograms Called Directly
c	Linpack - dsvdc dqrdc dqrsl
c   	Blas    - ddot dcopy
c	Other   - dset dprmut
c
c $Header: dtsvdc.f,v 2.100.1.1 86/10/07 12:58:56 lindstrom Exp $
c
      integer i,j,jj,pp1,ppnp1,nmk,sjob,locinf
      double precision trxpx,accum,dummy(1),mintr
      double precision ddot
c
      info = 0
      call dset(p,0.0d0,s,1)
      pp1 = p+1
      ppnp1 = p+n+1
c                       calculate trace of x' x
      trxpx = 0.0d0
      do 10 j = 1,p
         trxpx = trxpx+ddot(n,x(1,j),1,x(1,j),1)
 10   continue
      mintr = trxpx * minrat
c                       qr decomposition of x
      do 20 j = 1,p
         iwork(j) = 0
 20   continue
      call dqrdc(x,ldx,n,p,work,iwork,work(pp1),1)
c                       calculate ratios for the truncated matrix
      accum = 0.0d0
      k = p
      do 30 i = p,1,-1
         accum = accum+ddot(pp1-i,x(i,i),ldx,x(i,i),ldx)
         if (accum .lt. mintr) then
            k = i
            normk = accum
         endif
 30   continue
      if (job/10.le.0) then
c                       no left singular vectors
c                       copy rk' (transpose of the first k rows of r
c			from qr decomposition) into x
         do 40 j = 1,k 
            call dcopy(p,x(j,1),ldx,x(1,j),1)
            call dset(j-1,0.0d0,x(1,j),1)
 40      continue
c                       svd of rk'
         sjob = 0
         if (mod(job,10).ne.0) sjob = 20
         call dsvdc(x,ldx,p,k,s,work,v,ldv,dummy,0,work(pp1),sjob,info)
         if (info.ne.0) return
         if (mod(job,10).ne.0) then
            do 50 j = 1,k
               call dprmut(v(1,j),p,iwork,1)
 50         continue
         endif
      else
c                       u or u1 is to be created
         i = n+1
c                       copy rk' to work
         do 60 j = 1,k 
            i = i+p
            call dcopy(p,x(j,1),ldx,work(i),1)
            call dset(j-1,0.0d0,work(i),1)
 60      continue
c                       create u2 if requested
         if (job/10.eq.1) then
            nmk = n-k
            do 70 i = 1,nmk 
               j = n+1-i
               call dset(n,0.0d0,work(pp1),1)
               work(p+j) = 1.0d0
               call dqrsl(x,ldx,n,min(j,p),work,work(pp1),work(pp1),
     *          dummy,dummy,dummy,dummy,10000,locinf)
               call dcopy(n,work(pp1),1,u(1,j),1)
 70         continue
         endif
         do 80 j = 1,k 
            jj = k+1-j
            call dset(n,0.0d0,work(pp1),1)
            i = p+jj
            work(i) = 1.0d0
            call dqrsl(x,ldx,n,jj,work,work(pp1),work(pp1),dummy,dummy,
     *       dummy,dummy,10000,locinf)
            call dcopy(n,work(pp1),1,u(1,jj),1)
 80      continue
c			  svd of rk'
         sjob = 1
         if (mod(job,10).ne.0) then
            sjob = 21
         endif
         call dsvdc(work(ppnp1),p,p,k,s,work,v,ldv,work(ppnp1),p,
     *    work(pp1),sjob,info)
         if (info.ne.0) return
         do 100 i = 1,n 
            call dcopy(k,u(i,1),ldu,work,1)
            jj = n+1
            do 90 j = 1,k 
               jj = jj+p
               u(i,j) = ddot(k,work,1,work(jj),1)
 90         continue
 100     continue
         if (mod(job,10).ne.0) then
c		 undo pivots on right singular vectors
            do 110 j = 1,k
               call dprmut(v(1,j),p,iwork,1)
 110        continue
         endif
      endif
      end
