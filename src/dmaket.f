      subroutine dmaketg(m,n,dim,des,lddes,s1,lds1,ncov1,npoly,t,ldt,
     * wptr,info)
      integer m,n,dim,lddes,lds1,ncov1,npoly,ldt,wptr(dim),info
      double precision des(lddes,dim),s1(lds1,*),t(ldt,*)
c
c Purpose: create t matrix and append s1 to it.
c
c On Entry:
c   m			order of the derivatives in the penalty
c   n			number of rows in des
c   dim			number of columns in des
c   des(lddes,dim)	variables to be splined
c   lddes		leading dimension of des as declared in the 
c			calling program
c   s1(lds1,ncov1)	covariates which duplicate the replication 
c			structure of des
c   lds1		leading dimension of s1 as declared in the 
c			calling program
c   ncov1		number of columns in s1
c   ldt			leading dimension of t as declared in the 
c			calling program
c
c On Exit:
c   npoly		dimension of polynomial part of spline
c   t(ldt,npoly+ncov1)	[t:s1]
c   info 		error indication
c   			   0 : successful completion
c		 	   1 : error in creation of t
c Work Arrays:
c   wptr(dim)		integer work vector
c
c Subprograms Called Directly:
c	Blas  - dcopy
c	Other - mkpoly
c
c $Header: dmaketg.f,v 2.100.1.3 86/11/21 11:32:19 lindstrom Exp $
c
      integer i,j,k,tt,nt,bptr,eptr
      integer mkpoly
c
      info = 0
      npoly = mkpoly(m,dim)
      call dset(n,1.0d0,t(1,1),1)
      nt = 1
      if (npoly .gt. 1) then
          do 10 j=1,dim 
             nt = j + 1
             wptr(j) = nt
             call dcopy(n,des(1,j),1,t(1,nt),1)
   10     continue
c
c     get cross products of x's in null space for m>2
c
c     WARNING: do NOT change next do loop unless you fully understand:
c              This first gets x1*x1, x1*x2, x1*x3, then
c              x2*x2, x2*x3, and finally x3*x3 for dim=3,n=3
c              wptr(1) is always at the beginning of the current
c	       level of cross products, hence the end of the
c	       previous level which is used for the next.
c	       wptr(j) is at the start of xj * (previous level)
c
          do 50 k=2,m-1 
             do 40 j=1,dim 
                bptr = wptr(j)
                wptr(j) = nt + 1
                eptr = wptr(1) - 1
                do 30 tt=bptr,eptr 
                   nt = nt + 1
                   do 20 i=1,n
                      t(i,nt) = des(i,j) * t(i,tt)
   20              continue
   30           continue
   40        continue
   50     continue
          if (nt .ne. npoly) then
             info = 1
             return
          endif
      endif
c			append s1 to t
      do 60 i = 1,ncov1
         call dcopy(n,s1(1,i),1,t(1,nt+i),1)
   60 continue
      end
