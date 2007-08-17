      subroutine dvlop(z,svals,nobs,nnull,npsing,cost,inadd,ssqw2,
     * lamlim,ntbl,nlamht,tbl,ldtbl,auxtbl,dout,job,info)
      integer nobs,nnull,npsing,ntbl,ldtbl,job,info
      double precision z(npsing),svals(npsing),inadd,ssqw2,lamlim(2),
     * nlamht,tbl(ldtbl,3),auxtbl(3,3),dout(2),cost
c
c Purpose: determine the optimal lambda for the generalized cross 
c	validation function given singular values and the data vector 
c	in canonical coordinates.
c
c On Entry:
c   z(npsing)		data vector in canonical coordinates
c   svals(npsing)	singular values 
c   nobs		number of observations
c   nnull		dimension of the null space of sigma
c   npsing		number of positive elements of svals 
c   inadd		constant term in expression for V
c   ssqw2		squared length of w2
c   lamlim(2)		limits on lambda hat search (in log10(nobs*
c			lambda) scale) if user input limits are 
c			requested. if lamlim(1) = lamlim(2) then nlamht
c			is set to 10**lamlim(1)
c   ntbl		number of evenly spaced values for 
c			log10(nobs*lambda) to be used in the initial 
c			grid search for lambda hat
c			if ntbl = 0 only a golden ratio search will be 
c			done and tbl is not referenced, if ntbl > 0
c			there will be ntbl rows returned in tbl
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   job      		if job is nonzero then user input limits on 
c			lambda hat search are used
c
c On Exit:
c   lamlim(2)		limits on lambda hat search 
c			(in log10(nobs*lambda) scale)
c   nlamht		nobs*(lambda hat) where lambda hat is the gcv 
c			estimate of lambda
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  2  	V(lambda)
c   auxtbl(3,3)		auxiliary table
c			1st row contains:
c			    log10(nobs*lambda hat), V(lambda hat)
c			2nd row contains:
c			    0, V(0) 
c			3rd row contains:
c			    0, V(infinity) 
c   dout(2)		contains:
c			1  rss
c			2  tr(I-A)
c   info		error indicator
c			  0 : successful completion
c			 -1 : log10(nlamht) <= lamlim(1) (not fatal)
c			 -2 : log10(nlamht) >= lamlim(2) (not fatal)
c			  1 : svals(1) = 0.0d0
c			  2 : npsing is incorrect
c			  3 : lamlim(1) > lamlim(2)
c
c Subprograms Called Directly:
c	Gcvpack - dvmin 
c
c Subprograms Called Indirectly:
c	Gcvpack - dvl
c
c $Header: dvlop.f,v 2.100.1.1 86/10/07 12:59:38 lindstrom Exp $
c
      integer i,k
      double precision vlamht,w
      double precision dvmin
c
      common / gcvcom / addend,rss,tria,n,h
      integer n,h
      double precision addend,rss,tria,machpr,one
c
      info = 0
      one = 1.0d0
      machpr = 1.0d0
   10 machpr = machpr/2.0d0
      if (one .lt. 1.0d0 + machpr) goto 10
      machpr = machpr*2.0d0
c
      n=nobs
      h=nnull
      addend = inadd
      if (svals(1) .eq. 0.0d0) then
         info = 1
         return
      endif
      k = 0
      do 20 i = 1,npsing
         if (svals(i) .gt. 0) then
            k = i
         endif
   20 continue
      if (k .ne. npsing) then
	 info = 2
    	 return
      endif
      if (job .ne. 0 .and. (lamlim(1) .gt. lamlim(2))) then
         info = 3
         return
      endif
      if (job .eq. 0) then
         lamlim(2) = 2.0d0*dlog10(svals(1))+2.0d0
         lamlim(1) = 2.0d0*dlog10(svals(npsing))-2.0d0
      endif
      nlamht = dvmin (lamlim(1),lamlim(2),svals,z,npsing,cost,ntbl,tbl,
     * ldtbl,vlamht,info)
      dout(1) = rss
      dout(2) = tria
c			compute auxtbl
      auxtbl(1,1)=nlamht
      auxtbl(1,2)=vlamht
c			lambda = 0
      auxtbl(2,1)=0.0d0
      auxtbl(2,2)=0.0d0
      if ((nobs-nnull) .ne. npsing) then
         auxtbl(2,2)=inadd*(nobs)/(nobs-nnull-npsing)**2
      endif
      if ((nobs-nnull) .eq. npsing) then
         w=0.0d0
         do 30 i=npsing,1,-1
c           w=w+(z(i)*svals(npsing)**2/(svals(i)**2))**2
            w=w+(z(i)*(svals(npsing)/svals(i))**2)**2
   30    continue
         auxtbl(2,2)=nobs*w
         w=0.0d0
         do 40 i=npsing,1,-1
            w=w+(svals(npsing)/svals(i))**2
   40    continue
         auxtbl(2,2)=auxtbl(2,2)/(w**2)
      endif
c			lambda = infinity
      auxtbl(3,1)=0.0d0
      auxtbl(3,2)=ssqw2/(nobs - nnull)
      nlamht = 10**nlamht
      return
      end
