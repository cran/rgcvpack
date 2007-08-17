      integer function mkpoly(m,dim)
      integer m,dim
c
c  Purpose: compute the binomial coefficient of m + dim - 1 choose dim.
c  	This is the dimension of the space of polynomials which are in
c  	the null space of the smoothing penalty. Uses Andy Jaworski's
c	binomial coefficient algorithm that only requires integer
c	arithmetic.
c
c  On Entry:
c   m			order of derivatives in the penalty
c   dim	 		dimension of the variables to be splined
c
c  On Exit:
c   mkploy		(m + dim - 1) choose dim
c
c $Header: mkpoly.f,v 2.100.1.1 86/10/07 13:00:35 lindstrom Exp $
c
      integer i,j,k,k1,kcoef,n
c 			compute binomial coefficient 
c			m + dim - 1 choose dim
      n = m + dim - 1
      k1 = dim
      if (k1 .gt. n .or. k1 .lt. 0) then
         mkpoly = 0
         return
      endif
      k = k1
      if ((n - k1) .lt. k) then
         k = n - k1
      endif
      kcoef = 1
      j = n - k
      do 10 i = 1, k
         j = j + 1
         kcoef = (kcoef * j) / i
   10 continue
      mkpoly = kcoef
      return
      end
