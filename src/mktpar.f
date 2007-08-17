      subroutine mktpar(bign,nobs,testno,bandwd,tpar,t,f0,
     *   k0,info)
      double precision bandwd,tpar(4,*),t(nobs),f0(nobs),k0(nobs)
      integer bign,nobs,testno,info
c
c  Purpose: To construct the "true" parameter vectors for an
c	integral equation test case.
c
c  On Entry:
c   bign		dimension of Hilbert space
c   nobs		number of observations
c   testno		test number (determines the spread of the peaks)
c   bandwd		band width for the filter function
c   t(nobs)		the location of the observations
c
c  On Exit:
c   tpar(4,bign/2+1)	values of the "true" parameters
c   f0(nobs)		the "true" values of f at locations t
c   k0(nobs)		the "true" values of k at locations t
c   info		error indicator. 0 indicates successful 
c			completion. Other values are:
c				1 - invalid testno
c				2 - bign <= 1 or nobs <= 1
c				3 - bign not even 
c  $Header: mktpar.f,v 2.100.1.2 86/10/13 14:09:36 lindstrom Exp $ 
      double precision rt2pi,twopi,sn,cs,fo,f00,ko,k00,fdiff,
     *	iovern
      integer j,i
c     
      twopi = 8.0d0 * atan(1.0d0)
      rt2pi = sqrt(twopi)
      info = 1
      if ((testno .le. 0) .or. (testno .gt. 4)) return
      info = 2
      if ((bign .le. 1) .or. (nobs .le. 1)) return
      info = 3
      if (int(bign/2) .ne. dble(bign)/2.0d0) return
      info = 0
      fdiff = f00(0.0d+0,testno,rt2pi) - f00(1.0d+0,testno,rt2pi)
      call dset(4*(bign/2+1),0.0d0,tpar,1)
      do 20 j = 0, bign/2 
	  do 10 i = 1,bign
	      iovern = dble(i)/dble(bign)
	      sn = sin(twopi*j*iovern)
	      cs = cos(twopi*j*iovern)
 	      fo = f00(iovern,testno,rt2pi)+( iovern-0.5)*fdiff
	      tpar(1,j+1)= tpar(1,j+1)+cs*fo
	      tpar(2,j+1)= tpar(2,j+1)+sn*fo
	      ko = k00(iovern,bandwd,rt2pi)
	      tpar(3,j+1)= tpar(3,j+1)+cs*ko
	      tpar(4,j+1)= tpar(4,j+1)+sn*ko
   10     continue
   20 continue
      call dscal(4*(bign/2+1),1.0d0/dble(bign),tpar,1)
      do 30 j = 1,nobs
          f0(j) = f00(t(j),testno,rt2pi)+(t(j)-0.5)*fdiff
          k0(j) = k00(t(j),bandwd,rt2pi)
   30 continue
      return
      end
c
      double precision function f00(x,testno,rt2pi)
      double precision x,rt2pi
      integer testno
c
      double precision s1, s2, mu(4)
      data mu(1), mu(2), mu(3), mu(4) / 0.2, 0.15, 0.1, 0.05 /,
     *   s1, s2 / 0.015, 0.045 /
c
      f00 = (exp(-((x-0.3)/s1)**2/2.0d0)/s1+
     *    2.0d0*exp(-((x-0.3-mu(testno))/s2)**2/2.0d0)/s2)/(3.0d0*rt2pi)
      return
      end
c
      double precision function k00(x,bandwd,rt2pi)
      double precision x,bandwd,rt2pi
c
      k00 = (exp(-(x/bandwd)**2/2.0d0)+
     *    exp(-((1-x)/bandwd)**2/2.0d0))/(rt2pi*bandwd)
      return
      end
