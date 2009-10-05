      subroutine dcrtz (x,ldx,nobs,npar,pmh,qr,ldqr,qraux,sgpvt,work,
     * info)
      integer ldx,nobs,npar,pmh,ldqr,sgpvt(npar),info
      double precision x(ldx,npar),qr(ldqr,pmh),qraux(pmh),work(npar)
c
c Purpose: create z from x using permutations and householder 
c	   transformation from sigma
c
c On Entry:
c   x(ldx,npar)		design matrix
c   ldx		     	leading dimension of x as declared in the 
c			calling	program
c   nobs   		number of observations
c   npar      		number of parameters
c   pmh    		npar minus the dimension of the null space of 
c			the semi-norm
c   qr(ldqr,pmh)	information on the qr decomposition of the 
c			Cholesky factor of the semi-norm
c   ldqr    		leading dimension of qr as declared in the 
c			calling	program
c   qraux(pmh)		auxiliary information on the qr decomposition 
c			of the factor of the semi-norm
c   sgpvt(npar) 	permuted indices from the pivoted Cholesky
c			decomposition of the semi-norm matrix
c
c On Exit:
c   x(ldx,npar)		z (x after householder transformation)
c   info    		error indicator
c			  0 : successful completion
c			  1 : error in dtrsl, r is singular
c
c Work Arrays:
c   work(npar)		double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl dtrsl 
c	Blas    - dcopy 
c	Other   - dcpmut 
c
c Subprograms Called Indirectly:
c	Blas    - dswap
c
c $Header: dcrtz.f,v 2.100.1.1 86/10/07 12:47:59 lindstrom Exp $
c
      double precision dummy
      integer i,locinf
c
      info = 0
c			permute columns of x according to sgpvt
      call dcpmut (x,ldx,nobs,npar,sgpvt,0)
c			apply Householder transformations to rows of
c			permuted x to generate z1 and z2
      do 20 i = 1,nobs 
         call dcopy (npar,x(i,1),ldx,work,1)
         call dqrsl (qr,ldqr,npar,pmh,qraux,work,dummy,work,dummy,dummy,
     *    dummy,01000,locinf)
         call dtrsl (qr,ldqr,pmh,work,01,info)
         if (info .ne. 0) then
            info = 1
            return
         endif
         call dcopy (npar,work,1,x(i,1),ldx)
   20 continue
      return
      end
