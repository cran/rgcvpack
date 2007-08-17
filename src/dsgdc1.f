      subroutine dsgdc1(f2kf2,ldfkf,p,svals,iwork,work,lwa,info)
      integer ldfkf,p,iwork(p),lwa,info
      double precision f2kf2(ldfkf,p),svals(*),work(lwa)
c 
c Purpose: form the singular value decomposition of the Cholesky factor 
c 	of f2'k f2.
c
c On Entry:
c   f2kf2(ldfkf,p)	f2'k f2
c   ldfkf		leading dimension of f2'k f2 as declared in the 
c			calling	program 
c   p			number of rows and columns in f2'k f2
c
c On Exit:
c   f2kf2(p,p)  	overwritten with singular value decomposition 
c 			of Cholesky factor of f2'k f2
c   svals(p) 		the singular values of the Cholesky factor of 
c			f2'k f2 if info = 0.
c			if info = 3 then svals is as it was returned 
c			from dsvdc.
c   info 	   	error indicator
c			  0 : successful completion
c			  1 : lwa too small
c			  2 : f2'k f2 is not of full rank 
c			  3 : error in dsvdc
c   p			if info = 3 p contains info as it was returned 
c			from dsvdc (otherwise unchanged)
c
c Work Arrays:
c   work(lwa)		double precision work vector
c   lwa			length of work as declared in the calling
c			program
c			must be at least 2*p
c   iwork(p)		integer work vector
c
c Subprograms Called:
c	Linpack - dchdc dsvdc
c	Blas    - dcopy
c	Other   - dset dprmut
c
c $Header: dsgdc1.f,v 2.100.1.1 86/10/07 12:56:12 lindstrom Exp $
c
      integer i,j,pp1,locinf,k
      double precision dummy,one,machpr
c
      info = 0
      if (lwa .lt. 2*p) then
	  info = 1
	  return
      endif
      pp1 = p + 1
      call dset(p,0.0d0,svals,1)
c
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c			Cholesky decomposition of f2'k f2
      do 20 j = 1,p
         iwork(j) = 0
   20 continue
      call dchdc (f2kf2,ldfkf,p,work,iwork,1,locinf)
      do 30 i=1,locinf 
         if ((f2kf2(i,i)/f2kf2(1,1))**2 .gt. machpr) k = i
   30 continue
      if (k .lt. p) then
	 info = 2
	 return
      endif
c    		copy f2kf2' into f2kf2 
c    		svd of f2 k f2' = udv' return only u
      do 40 j = 1,p 
         call dcopy(1+p-j,f2kf2(j,j),ldfkf,f2kf2(j,j),1)
         call dset(j-1,0.0d0,f2kf2(1,j),1)
   40 continue
      call dsvdc(f2kf2,ldfkf,p,p,svals,work,f2kf2,ldfkf,dummy,1,
     * work(pp1),20,info)
      if (info .ne. 0) then
	  p = info
	  info = 3
	  return
      endif
      do 50 j=1,p
         call dprmut(f2kf2(1,j),p,iwork,1)
   50 continue
      return
      end
