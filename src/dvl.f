      double precision function dvl(lgnlam,svals,z,npsing,cost)
      integer npsing
      double precision lgnlam,svals(npsing),z(npsing),cost
c
c Purpose: evaluate the cross-validation function with a semi-norm.
c
c On Entry:
c   lgnlam		log10(nobs*lambda) where lambda is the value of
c			lambda for which V is evaluated
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive svals 
c
c On Exit:
c   dvl			V(lambda)
c
c $Header: dvl.f,v 2.100.1.1 86/10/07 12:59:22 lindstrom Exp $
c
      integer j
      double precision nlam,numrtr,denom,factor
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision rss,tria,addend
c     			see dvlop for definition of common block 
c			variables
c
      nlam = 10**lgnlam
      numrtr = addend
      denom = dble(n - h - npsing)
      do 10 j = 1,npsing 
         factor = 1.0d0/(1.0d0 + (svals(j)**2)/nlam)
         numrtr = numrtr + (factor*z(j))**2
         denom = denom + 1.0d+0 - cost*(1.0d+0 - factor)
 10   continue
      rss=numrtr
      tria=denom
      if (denom > 0.0d+0) then
         dvl=dble(n)*numrtr/denom**2
      else
         dvl=1.0d+20
      endif
      return
      end
      
