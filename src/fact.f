      integer function fact(i)
      integer i
c 
c Purpose: quick factorial function for the bspline routine
c	returns zero for negative i.
c
c On Entry:
c   i			a non-negative integer
c On Exit:
c   fact		i factorial
c
c $Header: fact.f,v 2.100.1.1 86/10/07 13:00:27 lindstrom Exp $
c
      integer j
      fact = 0
      if (i .ge. 0) fact = 1
      if (i .le. 1) return
      do 10 j = 2,i
	 fact = fact*j
   10	continue
      return
      end
