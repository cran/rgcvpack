      subroutine dpred(pdes,ldpdes,npred,dim,m,desb,lddesb,ndesb,ps,
     * ldps,ncov1,ncov2,coef,npar,pred,work,lwa,iwork,info) 
      integer ldpdes,npred,dim,m,lddesb,ndesb,ldps,ncov1,ncov2,npar,lwa,
     * iwork(dim),info
      double precision pdes(ldpdes,dim),desb(lddesb,dim),ps(ldps,*),
     * coef(npar),pred(npred),work(lwa)
c
c   Purpose: determine predicted values at the locations in pdes and ps
c
c  On Entry:
c   pdes(ldpdes,dim) 	prediction design for splined variables
c   ldpdes		leading dimension of pdes as declared in the 
c			calling	program 
c   npred		number of rows in pdes
c   desb(lddesb,dim) 	locations for the basis functions
c			(returned from dtpss and dptpss in the 
c			variable des)
c   lddesb		leading dimension of desb as declared in the
c			calling	program 
c   ndesb		number of rows in desb
c   dim			number of columns in desb
c   m			order of the derivatives in the penalty
c   ps(ldps,ncov1+ncov2) prediction covariates corresponding to pdes
c   ldps		leading dimension of ps as declared in the
c			calling	program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of pdes
c   ncov2		number of covariates which do not duplicate the 
c			replication structure of pdes
c   coef(npar)		coefficient estimates  [delta':xi']'
c   npar		ndesb + (m+dim-1 choose dim) + ncov1 + ncov2
c
c On Exit:
c   pred(npred)		predicted values
c   info		error indicator
c			  0 : successful completion
c			  1 : dimension error
c			  2 : error in npar,ncov1,ncov2,m or dim
c			  3 : lwa too small
c			  4 : error in dmaketg
c			
c
c Working Storage:
c   work(lwa)		double precision work vector
c   lwa			length of work vector
c			must be at least npred*(nct+ndesb)
c			where nct = (m+dim-1 choose dim)
c   iwork(dim)		integer work vector
c
c Subprograms Called Directly:
c    Gcvpack - dmaketg dmakekg
c    Blas    - dgemv
c
c Subprograms Called Indirectly:
c    Blas    - dcopy
c    Other   - fact mkpoly
c
c $Header: dpred.f,v 2.100.1.1 86/10/07 12:51:40 lindstrom Exp $
c
      double precision dummy
      integer nct,p1,p1p1,npoly
      integer mkpoly 
c
      nct = mkpoly(m,dim) 
      if ((ndesb .le. 0) .or. (nct .le. 0) .or. (m .le. 0) .or. 
     * (dim .le. 0) .or. 2*m - dim .le. 0) then
	 info = 1
	 return
      endif
      if (npar .ne. ndesb + nct + ncov1 + ncov2) then
         info = 2
         return
      endif
      if (lwa .lt. npred*(nct+ndesb)) then
	 info = 3
	 return
      endif
c			first npred*nct positions of work contain t
      p1 = npred*nct
c			next npred*ndesb positions of work contain k
      p1p1 = p1 + 1
c
      call dmaketg(m,npred,dim,pdes,ldpdes,dummy,1,0,npoly,work(1),
     * npred,iwork,info)
      if (info .ne. 0) then
         info = 4
         return
      endif
      call dmakekg(m,npred,dim,pdes,ldpdes,ndesb,desb,lddesb,work(p1p1),
     * npred)
c compute predicted values
      call dgemv('N',npred,nct,1.0d0,work,npred,coef,1,0.0d0,
     * pred,1)
      call dgemv('N',npred,ncov1+ncov2,1.0d0,ps,ldps,coef(nct+1),1,
     * 1.0d0,pred,1)
      call dgemv('N',npred,ndesb,1.0d0,work(p1+1),npred,
     * coef(nct+ncov1+ncov2+1),1,1.0d0,pred,1)
      return
      end
