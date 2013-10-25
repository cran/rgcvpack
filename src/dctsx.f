      subroutine dctsx(ts1,ldts1,ncts1,s2,lds2,ncov2,nobs,tbsb1,ldtbsb,
     * nb,sigma,ldsigm,x,ldx,npar,fgaux,work)
      integer ldts1,ncts1,lds2,ncov2,nobs,ldtbsb,nb,ldsigm,ldx,npar
      double precision ts1(ldts1,ncts1),s2(lds2,*),tbsb1(ldtbsb,ncts1),
     * sigma(ldsigm,*),x(ldx,*),fgaux(ncts1),work(nb)
c
c Purpose: compute x and sigma from t,s and k.
c
c On Entry:
c   ts1(nobs,ncts1)	[t:s1]
c   ldts1		leading dimension of ts1 as declared in the
c			calling	program 
c   ncts1	 	number of columns in [t:s1]
c   s2(lds2,ncov2) 	columns contain the covariates which do not 
c			duplicate the replication pattern of des
c   lds2		leading dimension of s2 as declared in the
c			calling	program 
c   ncov2		number of columns in s2
c   nobs		number of observations
c   tbsb1(ldtbsb,ncts1)	unique rows of t and s1 (or [t:s1] created from
c			basis functions)
c   ldtbsb		leading dimension of tbsb1 as declared in 
c			calling program
c   nb			number of basis functions
c   sigma		ku (ku must start in the ncov2+1 row and
c			the ncov2+1 column)
c   ldsigm		leading dimension of sigma as
c			declared in the calling program
c   x(ldx,npar)		k (k must start in the 1st row and column)
c   ldx			leading dimension of x as
c			declared in the calling program
c
c On Exit:
c   tbsb1(ldtbsb,ncts1)	qr decomposition of tbsb1
c   sigma(ldsigm,npar)  symmetric matrix that defines the semi-norm
c				[0:   0  ]
c				[0:f2'ku f2]
c   x(ldx,npar)		design matrix [t:s1:s2:kf2]
c   npar		number of parameters (nb+ncov2)
c   fgaux(ncts1)	auxiliary vector for qr decomposition of tbsb1
c
c Work Arrays:
c   work(nb)		double precision work vector
c
c Subprograms Called Directly
c	Linpack - dqrdc dqrsl
c	Blas    - dcopy
c	Other   - dset dftkf
c
c Subprograms Called Indirectly
c	Linpack - dqrsl
c	Blas    - dcopy
c
c $Header: dctsx.f,v 2.100.1.1 86/10/07 12:48:06 lindstrom Exp $
c
      double precision dummy
      integer i,locinf
c
      npar = nb + ncov2
      call dqrdc(tbsb1,ldtbsb,nb,ncts1,fgaux,0,0.0d0,0)
c			calculate k f2  put in last nb -ncts1 columns 
c			of x
      do 10 i=1,nobs 
         call dcopy(nb,x(i,1),ldx,work,1)
         call dqrsl(tbsb1,ldtbsb,nb,ncts1,fgaux,work,dummy,work,dummy,
     *    dummy,dummy,01000,locinf)
         call dcopy(nb-ncts1,work(ncts1+1),1,x(i,ncov2+ncts1+1),ldx)
   10 continue
c			copy [t:s1] into first ncts1 columns of x
      do 20 i=1,ncts1
         call dcopy(nobs,ts1(1,i),1,x(1,i),1)
   20 continue
c			copy s2 into next ncov2 columns of x
      do 30 i=1,ncov2
         call dcopy(nobs,s2(1,i),1,x(1,ncts1+i),1)
   30 continue
c			calculate sigma
      call dftkf(tbsb1,ldtbsb,nb,ncts1,fgaux,sigma(ncov2+1,ncov2+1),
     * ldsigm,work)
      do 40 i = 1,ncts1+ncov2
         call dset(npar,0.0d0,sigma(1,i),1)
   40 continue
      do 50 i = ncts1+ncov2+1,npar
         call dset(ncts1+ncov2,0.0d0,sigma(1,i),1)
   50 continue
      return
      end
