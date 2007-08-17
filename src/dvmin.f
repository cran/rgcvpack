      double precision function dvmin(lower,upper,svals,z,npsing,cost,
     * ntbl,tbl,ldtbl,vlamht,info)
      integer npsing,ntbl,ldtbl,info
      double precision lower,upper,svals(npsing),z(npsing),cost,
     * tbl(ldtbl,3),vlamht
c
c Purpose: evaluate V(lambda) for a grid of ln(nobs*lambda) values 
c	between	lower and upper, store these in the array tbl, and find
c	minimizer of v.
c
c On Entry:
c   lower		lower bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   upper		upper bound of interval (in nobs*ln(lambda)
c			scale) over which V(lambda) is to be minimized
c   svals(npsing)	singular values 
c   z(npsing)		data vector in canonical coordinates
c   npsing		number of positive elements of svals 
c   ntbl		number of evenly spaced values for 
c			ln(nobs*lambda)	to be used in the initial grid 
c			search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of ln(nobs*lambda) 
c			  2  	V(lambda)
c   vlamht		V(lambda hat)
c   dvmin		ln(nobs*lambda hat)
c   info		error indicator
c			  0 : successful completion
c			 -1 : dvmin <= lower (not fatal)
c			 -2 : dvmin >= upper (not fatal)
c
c Subprograms Called Directly:
c	Gcvpack - dvl
c
c $Header: dvmin.f,v 2.100.1.1 86/10/07 12:59:51 lindstrom Exp $
c

      double precision a,b,c,d,vc,vd,del,k1,k2,x,v
      integer j,jmin,k
      double precision dvl
c				null interval
      if (lower .eq. upper) then
	 dvmin = lower
	 info = -1
	 vlamht = dvl(lower,svals,z,npsing,cost)
	 do 10 j = 1, ntbl
	    tbl(j,1) = lower
	    tbl(j,2) = vlamht
   10    continue
	 return
      end if
c				non-null interval
      info = 0
      a = lower
      b = upper
      if (ntbl .eq. 1) then
	 x = (a + b)/2
	 tbl(1,1) = x
	 tbl(1,2) = dvl(x,svals,z,npsing,cost)
      else if (ntbl .ge. 2) then
c			do grid search
	 v=dvl(lower,svals,z,npsing,cost)*2.0d0
	 del=(upper-lower)/(ntbl-1)
	 do 20 j = 1, ntbl
	    tbl(j,1) = lower + (j - 1) * del
	    tbl(j,2) = dvl(tbl(j,1),svals,z,npsing,cost)
	    if (tbl(j,2) .le. v) then
	       jmin = j
	       v = tbl(j,2)
	    endif
   20    continue	      
	 a=tbl(jmin,1)-del
	 b=tbl(jmin,1)+del
      end if
c			do golden ratio search			
      k1=(3.0d0-dsqrt(5.0d0))/2.0d0
      k2=(dsqrt(5.0d0)-1)/2.0d0
      c = a + k1*(b - a)
      d = a + k2*(b - a)
      vc = dvl(c,svals,z,npsing,cost)
      vd = dvl(d,svals,z,npsing,cost)
      do 30 k=1,50
	 if (vd .lt. vc) then
	    a = c
	    c = d
	    d = a + k2*(b - a)
	    vc = vd
	    vd = dvl(d,svals,z,npsing,cost)
	 else
	    b = d
	    d = c
	    c = a + k1*(b - a)
	    vd = vc
	    vc = dvl(c,svals,z,npsing,cost)
	 end if
   30 continue
      x=(a+b)/2
      if (x .le. lower) info = -1
      if (x .ge. upper) info = -2
      vlamht=dvl(x,svals,z,npsing,cost)
      dvmin = x
      return
      end
      
