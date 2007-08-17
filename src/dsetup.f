      subroutine dsetup(des,lddes,su1,ldsu1,dim,m,ncov1,nuobs,c1,
     * tusu1,ldtu,ncts1,fgaux,ku,ldku,work,iwork,info)
      integer lddes,ldsu1,dim,m,ncov1,nuobs,ldtu,ncts1,ldku,
     * iwork(dim),info
      double precision des(lddes,dim),su1(ldsu1,*),c1(nuobs),
     * tusu1(ldtu,ncts1),fgaux(ncts1),ku(ldku,nuobs),work(ncts1)
c
c Purpose: set up [tu:su1] as f g and ku as f'c1 ku c1'f.
c
c On Entry:
c   des(lddes,dim)  	variables to be splined (unique rows)
c   lddes		leading dimension of des as declared in the 
c			calling	program 
c   su1(ldsu1,ncov1)	covariates (unique rows)
c   ldsu1		leading dimension of su1 as declared in the
c			calling	program 
c   dim 		dimension of the variables to be splined
c   m			order of the derivatives in the penalty
c   ncov1		number of covariates
c   nuobs		number of unique rows in des
c   c1(nuobs)		c1(i) contains the square root of the number of
c			replicates of the ith sorted design point
c   ldtu		leading dimension of tusu1 as declared in the 
c			calling	program 
c   ldku		leading dimension of ku as declared in the 
c			calling program 
c
c On Exit:
c   tusu1(ldtu,ncts1)	the qr decomposition of [tu:su1]
c   ncts1		number of columns in [tu:su1] = npoly + ncov1
c   fgaux(ncts1)	the auxiliary info on the qr decomposition of 
c			[tu:su1]
c   ku(p,p)  		f'ku f
c   info		error indicator
c			   0 : successful completion
c			   1 : error in dmaketg
c
c Work Arrays:
c   work(ncts1)		double precision work vector
c   iwork(dim)		integer work vector
c
c Subroutines called directly
c	Gcvpack - dmaketg dmakekg
c	Linpack - dqrdc
c	Other   - dftkf mkpoly
c
c Subroutines called indirectly
c	Blas    - dcopy
c	Gcvpack - dqrsl
c	Other   - fact mkpoly
c
c $Header: dsetup.f,v 2.100.1.1 86/10/07 12:55:48 lindstrom Exp $
c
      integer npoly,i,j
      integer mkpoly
c
      info = 0
      npoly=mkpoly(m,dim)
      ncts1=npoly+ncov1
c			make [tu:su1] and ku
      call dmaketg(m,nuobs,dim,des,lddes,su1,ldsu1,ncov1,npoly,tusu1,
     * ldtu,iwork,info)
      if (info .ne. 0) then
	 return
      endif
      call dmakekg(m,nuobs,dim,des,lddes,nuobs,des,lddes,ku,ldku)
      if (c1(1) .ne. 0) then
         do 30 i = 1,nuobs 
            do 10 j = 1,npoly+ncov1
               tusu1(i,j) = tusu1(i,j) * c1(i)
   10       continue
            do 20 j = 1,nuobs
               ku(i,j) = ku(i,j) * c1(i) * c1(j)
   20       continue
   30    continue
      endif
c			decompose [tu:su1] into fg
      call dqrdc(tusu1,ldtu,nuobs,ncts1,fgaux,0,0.0d0,0)
c      			calculate f'ku f
      call dftkf(tusu1,ldtu,nuobs,ncts1,fgaux,ku,ldku,work)
      return
      end
