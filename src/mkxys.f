      subroutine mkxys(bign,nobs,p,m,t,tpar,x,ldx,truey,
     *  sigma,ldsig,ftrig,ktrig,info)
      double precision t(nobs),tpar(4,*),x(ldx,*),truey(nobs),
     *  sigma(ldsig,*),ftrig(nobs),ktrig(nobs)
      integer bign,nobs,p,m,ldx,ldsig,info
 
c  Purpose: to create the design matrix, the observation vector and the
c	sigma matrix for the integral equation example.
c
c  On Entry:
c   bign		dimension of Hilbert space
c   nobs		number of observations
c   p			the number of parameters
c   m			order of the derivatives in the penalty
c   t(nobs)		locations for observations
c   tpar		values of the "true" parameters
c   ldx			leading dimension of x
c   ldsig 		leading dimension of sigma
c
c  On Exit:
c   x(nobs,p)		design matrix for dsnsm
c   truey(nobs)		"true" values for y (no noise added)
c   sigma(ldsig,p)	sigma matrix
c   ftrig(nobs)         trigonometric interpolant to f0
c   ktrig(nobs)         trigonometric interpolant to k0
c   info		error indicator. 0 indicates successful 
c			completion. Other values are:
c				1 - invalid testno
c				2 - bign <= 1 or nobs <= 1
c				3 - bign not even 
c  $Header: mkxys.f,v 2.100.1.1 86/10/07 13:30:50 lindstrom Exp $ 
      double precision twopi,sn,cs,tmp
      integer i,j,kb,r,flag
      
      twopi = 8.0d0 * atan(1.0d0)
      info = 1
      if (bign .le. 1 .or. nobs .le. 1) return
      info = 2
      if (int(bign/2) .ne. dble(bign)/2.0d0) return
      info = 3
      flag = 0
      if (p .eq. bign) then
	  flag = 1
      else
          if (int(p/2) .eq. dble(p)/2.0d0) return
      endif
      info = 0
      call dset(nobs,tpar(3,1),x(1,1),1)
      if (flag .eq. 1) then
	 kb = bign/2 + 1
	 r = bign/2 - 1
      else
         kb = (p - 1)/2 + 1
	 r = (p - 1)/2 
      endif
      do 20 i = 1, nobs
	  cs = cos( twopi*bign*t(i)/2.0d0)
	  truey(i) = tpar(1,1)*tpar(3,1)+ tpar(1,bign/2+1)*
     *	    tpar(3,bign/2+1)*cs*2.0d0
	  ftrig(i)= tpar(1,1) + tpar(1,bign/2+1)*cs/2.0d0
	  ktrig(i)= tpar(3,1) + tpar(3,bign/2+1)*cs/2.0d0
	  do 10 j = 1,bign/2-1
	      sn = sin(twopi*j*t(i))
	      cs = cos(twopi*j*t(i))
	      truey(i) = truey(i) + 2*((tpar(1,j+1)*tpar(3,j+1)-
     *	       tpar(2,j+1)*tpar(4,j+1))*cs+
     *	       (tpar(1,j+1)*tpar(4,j+1)+tpar(2,j+1)*tpar(3,j+1))*sn)
	      ftrig(i) = ftrig(i) + 2.0d0*tpar(1,j+1)*cs+
     *		    2.0d0*tpar(2,j+1)*sn
	      ktrig(i) = ktrig(i) + 2.0d0*tpar(3,j+1)*cs+
     *		    2.0d0*tpar(4,j+1)*sn
	      if (j .le. r) then
	          x(i,1+j) = 2.0d0*(tpar(3,j+1)*cs+tpar(4,j+1)*sn)
	          x(i,kb+j) = 2.0d0*(tpar(3,j+1)*sn-tpar(4,j+1)*cs)
	      endif
   10     continue
      if (flag .eq. 1) then
	  x(i,kb) = tpar(3,bign/2+1)*cos(twopi*bign*t(i)/2.0d0)/2.0d0
      endif
   20 continue
c			calculate sigma
      do 30 i = 1,p
	  call dset(p,0.0d0,sigma(1,i),1)
   30 continue
      do 40 i = 1,r
	  tmp = (twopi*i)**(2*m)
	  sigma(1+i,1+i) = tmp
	  sigma(kb+i,kb+i) = tmp
   40 continue
      if (flag .eq. 1) then
	  sigma(kb,kb) = ((twopi*dble(bign/2.0d0))**(2*m))/2.0d0
      endif
      return
      end
