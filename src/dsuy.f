      subroutine dsuy(y,nobs,nuobs,ytrue,c1,order,xrep,ssqrep,
     * job)
      integer nobs,nuobs,order(nobs),xrep(nobs),job
      double precision y(nobs),ytrue(nobs),c1(nuobs),ssqrep
c
c Purpose: compute B1'y, B1'ytrue and ssq for replication.
c
c On Entry:
c   y(nobs)  		response vector
c   nobs		number of observations
c   nuobs		number of unique design points
c   ytrue(nobs)		"true" response, if job is nonzero 1
c			not referenced if job = 0
c   c1(nuobs)		c1(i) contains the square root of the number of
c			replicates of the ith sorted design point
c   order(nobs)		order of sorted des
c   xrep(nobs)		xrep(i) = 1 if the ith sorted design point is a
c			replicate, 0 if not
c   job 		job is nonzero if B1'ytrue should be calculated 
c			job = 0 otherwise
c
c On Exit:
c   y(nuobs)  		B1'y
c   ytrue(nuobs)	B1'ytrue if job is nonzero
c   ssqrep		sum of squares for replication
c
c $Header: dsuy.f,v 2.100.1.1 86/10/07 12:57:32 lindstrom Exp $
c
      integer first,i,j
      double precision accum
c
      accum = 0.0d0
      ssqrep = 0.0d0
      call dprmut (y,nobs,order,0)
      if (job .ne. 0) call dprmut (ytrue,nobs,order,0)
c			compute ssq for replication
      first = 0
      do 20 i = 1,nobs
	  if (xrep(i) .eq. 1) then
	      accum = accum + y(i)
	  else if (first .eq. i - 1) then
	      first = i
	      accum = y(i)
	  else
	      accum = accum/(i-first)
	      do 10 j = first,i-1
	          ssqrep = (y(j)-accum)**2 + ssqrep
   10         continue
	      first = i
	      accum = y(i)
          endif
   20 continue
      if (xrep(nobs) .eq. 1) then
          accum = accum/(nobs + 1 - first)
          do 30 j = first,nobs
              ssqrep = (y(j)-accum)**2 + ssqrep
   30     continue
      endif
c			compute B1'y and B1'ytrue
      j = 0
      do 40 i = 1,nobs
 	  if (xrep(i) .eq. 0) then
 	      if (j .ne. 0) then
 		  y(j) = y(j) / c1(j)
 		  if (job .ne. 0) ytrue(j) = ytrue(j) / c1(j)
              endif
 	      j = j + 1
	      y(j) = y(i)
	      if (job .ne. 0) ytrue(j) = ytrue(i)
	  else
	      y(j) = y(j) + y(i)
	      if (job .ne. 0) ytrue(j) = ytrue(j) + ytrue(i)
 	  endif
   40 continue
      y(j) = y(j) / c1(j)
      if (job .ne. 0) ytrue(j) = ytrue(j) / c1(j)
      return 
      end
