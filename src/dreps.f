      subroutine dreps(des,lddes,nobs,dim,s,lds,ncov1,ncov2,c1,order,
     * nuobs,xrep,job,info)
      integer lddes,nobs,dim,lds,ncov1,ncov2,order(nobs),nuobs,
     * xrep(nobs),info,job
      double precision des(lddes,dim),s(lds,*),c1(*)
c
c Purpose: sort des and s, compute degrees of freedom for replication
c	and c1.
c
c On Entry:
c   des(lddes,dim) 	design for the variables to be splined
c			should be entered in lexicographical order 
c			(smallest to largest) if possible for efficient
c			computing
c   lddes		leading dimension of des as declared in the
c   			calling program 
c   nobs		number of observations
c   dim			number of columns in des
c   s(lds,ncov1+ncov2) 	design for the covariates
c   lds			leading dimension of s as declared in the
c   			calling program 
c   ncov1		number of covariates which duplicate the 
c			replication structure of des
c   ncov2		number of covariates which do not duplicate
c			replication structure of des
c   job			if job is nonzero then c1 is computed
c   			if job = 0 then c1 is not referenced
c
c On Exit:
c   des(lddes,dim) 	des sorted lexicographically
c   s(lds,ncov1+ncov2) 	s, sorted to correspond to des
c   c1(nuobs)		if job is nonzero then c1(i) the square root of
c			the number of replicates of the ith sorted 
c			design point 
c   order(nobs)		order of the sorted des
c   nuobs		number of unique rows in des
c   xrep(nobs)		xrep(i) = 1 if the ith sorted design point is a
c			replicate, 0 if not
c   info		error indicator
c			   0 : successful completion
c			   1 : ncov1 is incorrect
c
c
c $Header: dreps.f,v 2.100.1.1 86/10/07 12:55:24 lindstrom Exp $
c
      integer sw,oldsw,itemp,i,j,k,cont,dfrep
      double precision temp,diff,one,machpr,denom,wmin,wmax
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      do 20 i = 1,nobs
          order(i) = i
 	  xrep(i) = 0
   20 continue
      if (job .ne. 0) call dset(nobs,0.0d0,c1,1)
c			sort des and s 
      sw = nobs - 1
   30 if (sw .le. 0) goto 90 
          oldsw = sw
          sw = 0
          do 80 i = 1,oldsw 
	      cont = 1
	      k = 1
   40         if (cont .eq. 0) goto 70
                  if (k .le. dim) then
	              diff = des(i,k) - des(i+1,k)
                  else
	              diff = s(i,k-dim) - s(i+1,k-dim)
                  endif
		  if (diff .lt. 0.0d0) then
		      if (k .gt. dim) info = 1
		      cont = 0
		  else if (diff .gt. 0.0d0) then
		      if (k .gt. dim) info = 1
c		      switch the order of i and i+1
		      itemp = order(i)
		      order(i)=order(i+1)
		      order(i+1) = itemp
		      itemp = xrep(i)
		      xrep(i)= xrep(i+1)
		      xrep(i+1)= itemp
		      do 50 j = 1,dim
			  temp = des(i,j)
			  des(i,j) = des(i+1,j)
			  des(i+1,j) = temp
   50 		      continue
  		      do 60 j = 1,ncov1+ncov2
			  temp = s(i,j)
			  s(i,j) = s(i+1,j)
			  s(i+1,j) = temp
   60 		      continue
		      sw = i
		      cont = 0
                  else if (k .eq. dim + ncov1) then
		      xrep(i + 1) = 1
		      cont = 0
                  else 
	    	      k = k + 1
                  endif
              goto 40
   70         continue
   80     continue
      goto 30
   90 continue
c			compute range of design
      denom=0.0d0
      do 120 j=1,dim
          wmin = des(1,j)
	  wmax = des(1,j)
          do 110 i=1,nobs
	      if (des(i,j) .lt. wmin) wmin = des(i,j)
	      if (des(i,j) .gt. wmax) wmax = des(i,j)
  110     continue
	  denom = denom + (wmax-wmin)**2
  120 continue
	    
c			check for design points too close together
      do 140 i=1,nobs-1
	  if (xrep(i+1) .eq. 0) then
	     diff = 0.0d0
	     do 130 j=1,dim
	        diff = diff + (des(i,j)-des(i+1,j))**2
  130	     continue
	     if (abs(diff)/denom .lt. 100*machpr) xrep(i+1)=1
	  endif 
  140 continue
c			compute dfrep and c1
      dfrep = 0
      j = 0
       do 150 i = 1,nobs
	   j = j + 1 - xrep(i)
	   if (job .ne. 0) c1(j) = xrep(i)*c1(j) + 1.0d0
	   dfrep = dfrep + xrep(i)
  150 continue
      nuobs = nobs - dfrep
      if (job .eq. 0 ) return
      do 160 i = 1,nuobs
	  c1(i) = sqrt(c1(i))
  160 continue
      return 
      end
