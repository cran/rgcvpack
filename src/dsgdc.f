      subroutine dsgdc (sigma,ldsigm,npar,nnull,qraux,sgpvt,info)
      integer ldsigm,npar,nnull,sgpvt(npar),info
      double precision sigma(ldsigm,npar),qraux(npar)
c
c Purpose: decompose the semi-norm smoothing matrix into a QR 
c	decomposition of the transpose of the Cholesky factor.
c
c On Entry:
c   sigma(ldsigm,npar)	symmetric matrix which defines the semi-norm
c   ldsigm		leading dimension of sigma as declared in the
c			calling program
c   npar		number of parameters
c   nnull		dimension of the null space of sigma
c
c On Exit:
c   sigma(ldsigm,npar)	overwritten with the QR decomposition
c			of the Cholesky factor of sigma
c   nnull		if input nnull is too small it is replaced by
c			larger value such that sigma has rank npar-nnull
c   qraux(npar)		auxiliary information for the QR decomposition
c   sgpvt(npar) 	permuted indices from the Cholesky decomposition
c			with pivots of sigma
c   info		error indicator
c			   0 : successful completion
c			  -3 : nnull is too small (not fatal)
c			   1 : nnull is too large
c
c Subprograms Called Directly:
c	Linpack - dchdc dqrdc  
c	Blas    - dcopy
c	Other   - dset
c
c $Header: dsgdc.f,v 2.100.1.1 86/10/07 12:56:05 lindstrom Exp $
c
      integer locinf,i,j,nsm,idummy
      double precision dummy,machpr,one
c
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      locinf = 0
c			Cholesky decomposition of sigma
      do 20 j = 1,npar
         sgpvt(j) = 0
   20 continue
      call dchdc (sigma,ldsigm,npar,qraux,sgpvt,1,locinf)
      do 30 i=1,locinf 
         if (((sigma(i,i)/sigma(1,1))**2) .gt. machpr) nsm = i
   30 continue
      if (nsm .gt. npar - nnull) then
         info = 1
         return
      endif
      if (nsm .lt. npar - nnull) then
         nnull = npar - nsm
         info = -3
      endif
c			copy transpose of Cholesky factor to sigma
      do 40 i = 1,nsm 
         call dcopy (npar,sigma(i,1),ldsigm,sigma(1,i),1)
         j = i - 1
         call dset (j,0.0d0,sigma(1,i),1)
   40 continue
c			QR decomposition of Cholesky transpose
      call dqrdc (sigma,ldsigm,npar,nsm,qraux,idummy,dummy,0)
      return
      end
