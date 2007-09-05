      subroutine dtpkm(des,lddes,nobs,dim,m,s,lds,ncov,knots,ldkt,nknt,
     * y,ntbl,adiag,lamlim,cost,dout,iout,coef,svals,tbl,ldtbl,auxtbl,
     * work,lwa,iwork,liwa,tau,job,info)
      
      integer lddes,nobs,dim,m,lds,ncov,ldkt,nknt,ntbl,iout(3),ldtbl,
     $     lwa,liwa,iwork(liwa),job,info
      double precision des(lddes,dim),s(lds,*),knots(ldkt,dim),y(nobs),
     $     adiag(nobs),lamlim(2),dout(5),coef(*),svals(*),tbl(ldtbl,3),
     $     auxtbl(3,3),work(lwa),tau,cost
      
c     
c Purpose: fit a thin plate spline with specified knots configuration. 
c     Smoothing parameter either chosen by GCV or given by user.
c
c On Entry:
c   des(lddes,dim)      design for the variables to be splined
c   lddes               leading dimension of des as declared in calling
c                       program 
c   nobs                number of observations
c   dim                 number of columns in des
c   m                   order of the derivatives in the penalty
c   s(lds,ncov)         not used
c   lds                 leading dimension of s as declared in calling
c                       program 
c   ncov                number of covariates 
c   knots(ldkt,dim)     knots placement for the variable to be splined
c   ldkt                leading dimension of the knots array
c   y(nobs)             response vector
c   ntbl                number of evenly spaced values for 
c                       log10(nobs*lambda) to be used in the initial 
c                       grid search for lambda hat 
c                       if ntbl = 0 only a golden ratio search will be 
c                       done and tbl is not referenced, if ntbl > 0
c   adiag(nobs)         "true" y values on entry if predictive mse is 
c                       requested
c   lamlim(2)           limits on lambda hat search (in log10(nobs*
c                       lambda) scale) if user input limits are 
c                       requested if lamlim(1) = lamlim(2) then lamhat
c                       is set to (10**lamlim(1))/nobs
c   ldtbl               leading dimension of tbl as declared in the 
c                       calling program 
c   job                 integer with decimal expansion abcd
c                       if a is nonzero then truncation is used
c                       if b is nonzero then predictive mse is computed
c                          using adiag as true y
c                       if c is nonzero then user input limits on search
c                          for lambda hat are used
c                       if d is nonzero then the diagonal of the hat 
c                          matrix is calculated
c
c On Exit:
c   y(nobs)             predicted values
c   adiag(nobs)         diagonal elements of the hat matrix if requested
c   lamlim(2)           limits on lambda hat search 
c                       (in log10(nobs*lambda) scale)
c   dout(5)             contains:
c                       1 lamhat   generalized cross validation 
c                                  estimate of the smoothing parameter
c                       2 penlty   smoothing penalty
c                       3 rss      residual sum of squares
c                       4 tr(I-A)  trace of I - A
c                       5 truncation ratio = 1/(1+(normk/(nobs*lamhat)))
c                                  where normk = norm(R - R sub k)**2
c   iout(3)             contains:
c                       1  npsing   number of positive singular
c                                   values
c                                   if info indicates nonzero info in 
c                                   dsvdc then iout(1) contains info as
c                                   it was returned from dsvdc
c                       2  npar     number of parameters
c                       3  nnull    size of the null space of sigma
c   coef(nknt+nnull)    coefficient estimates [beta'; delta']'
c   svals(nknt-nnull)   first npsing entries contain singular values 
c                       of the matrix j2 
c                       if info indicates nonzero info in dsvdc then 
c                       svals is as it was returned from dsvdc
c   tbl(ldtbl,3)        column  contains
c                         1     grid of log10(nobs*lambda) 
c                         2     V(lambda)
c                         3     R(lambda) if requested
c   auxtbl(3,3)         auxiliary table
c                       1st row contains:
c                           log10(nobs*lamhat), V(lamhat) and  
c                           R(lamhat) if requested
c                           where lamhat is the gcv estimate of lambda
c                       2nd row contains:
c                           0, V(0) and  R(0) if requested
c                       3rd row contains:
c                           0, V(infinity) and R(infinity) if requested
c   info                error indicator
c                          0 : successful completion
c                         -3 : nnull is too small (not fatal)
c                         -2 : log10(nobs*lamhat) >= lamlim(2) 
c                              (not fatal)
c                         -1 : log10(nobs*lamhat) <= lamlim(1)
c                              (not fatal)
c                          1 : dimension error  
c                          2 : lwa (length of work) is too small
c                          3 : liwa (length of iwork) is too small
c                          4 : error in ntbl or tau
c                         100< info <200 : 100 + nonzero info returned
c                                          from ddcom
c                         200< info <300 : 200 + nonzero info returned
c                                          from dgcv
c
c Work Arrays:
c   work(lwa)           double precision work vector, the first part from
c                       1 to nnull+nknt*(nnull+nobs+nknt) is used to hold
c                       input/output arrays for dsnsm, the rest of work
c                       array used as the working array of dsnsm
c   lwa                 length of work as declared in the calling 
c                       program 
c                       must be at least
c                       nnull*(1+nnull+nobs+nknt)+nknt+nobs+
c                       (nknt-nnull)*(2*nknt-nnull+2+2*nobs)
c   iwork(liwa)         integer work vector
c   liwa                length of iwork as declared in the calling 
c                       program
c                       must be at least 2*nknt - nnull 
c
c Subprograms Called Directly:
c       Gcvpack - mkpoly, dmaket, dmakek, dftkf, dsnsm
c       Linpack - dqrdc, dqrsl
c       Blas    - dcopy, dscal
c
c $Header: dtpkm.f,v 2.100.1.2 04/05/22 09:24:39 Xianhong Exp $    
c
      
      integer nnull,ptb,pqr,pk,pt0,px,pkb,psigma,pwk,llwa,iinfo
      double precision dummy
      
c     
c set the pointers of the variables in the work array
c     
      nnull = mkpoly(m,dim)      

c allocate space for qr decomposition of TB      
      ptb = 1
      pqr = ptb + nknt*nnull
c make K, T0, and X share the same starting addresses
      pk  = pqr + nnull
      pt0 = pk
      px  = pk
c make KB and Sigma share the same starting addresses
      pkb    = px + nobs*nknt
      psigma = pkb
c the true starting address of the work array and the length
      pwk  = psigma + nknt*nknt      
      llwa = (nknt-nnull)*(nknt-2*nnull+2+nobs)+nknt+nobs
      
c      
c generate the matrix TB
c
      call dmaketg(m,nknt,dim,knots,ldkt,s,lds,ncov,nnull,work(ptb),
     $     nknt,work(pwk),info)
      
c QR decomposition on TB (= [F1 F2]*[R; 0])
      call dqrdc(work(ptb),nknt,nknt,nnull,work(pqr),0,0.0d+0,0)
      
c generate the matrix K
      call dmakekg(m,nobs,dim,des,lddes,nknt,knots,ldkt,work(pk),nobs)
      
c calculate [K*F1 K*F2], store in K
      do 10, i=1,nobs
         call dcopy(nknt,work(pk+i-1),nobs,work(pwk),1)
         call dqrsl(work(ptb),nknt,nknt,nnull,work(pqr),work(pwk),dummy,
     $        work(pwk),dummy,dummy,dummy,01000,info)
         call dcopy(nknt,work(pwk),1,work(pk+i-1),nobs)
 10   continue
      
c generate matrix T0, overwrite K*F1      
      call dmaketg(m,nobs,dim,des,lddes,s,lds,ncov,nnull,work(pt0),
     $     nobs,work(pwk),info)
      
c generate the matrix KB
      call dmakekg(m,nknt,dim,knots,ldkt,nknt,knots,ldkt,work(pkb),nknt)
      
c calculate F'*KB*F, store in KB
      call dftkf(work(ptb),nknt,nknt,nnull,work(pqr),work(pkb),nknt,
     $     work(pwk))
      
c use 0 to overwrite the upper left, upper right, and lower left of KB
      do 20, i=1,nnull
         call dscal(nknt,0.0d+0,work(pkb+i-1),nknt)
 20   continue
      
      do 30, j=1,nnull
         call dscal(nknt,0.0d+0,work(pkb+(j-1)*nknt),1)
 30   continue
     
c call the generalized ridge regression subroutine to fit a tps
      call dsnsm(work(px),nobs,y,work(psigma),nknt,nobs,nknt,nnull,
     $     adiag,tau,lamlim,cost,ntbl,dout,iout,coef,svals,tbl,ldtbl,
     $     auxtbl,iwork,liwa,work(pwk),llwa,job,info)
      
c transform the coef \xi to the original space, i.e. calculate F2*\xi
      call dcopy(nknt,coef,1,work(pwk),1)
      call dscal(nnull,0.0d+0,work(pwk),1)
      call dqrsl(work(ptb),nknt,nknt,nnull,work(pqr),work(pwk),
     $     work(pwk),dummy,dummy,dummy,dummy,10000,iinfo)
      call dcopy(nknt,work(pwk),1,coef(nnull+1),1)
      
      return
      
c end of dtpkm
      end
