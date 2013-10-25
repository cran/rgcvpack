      subroutine dpmse(fg,ldfg,nuobs,nobs,nnull,fgaux,svals,npsing,u,
     * ldu,w1,z,ntbl,adiag,tbl,ldtbl,auxtbl,work)
      integer ldfg,nuobs,nobs,nnull,npsing,ldu,ntbl,ldtbl
      double precision fg(ldfg,nnull),fgaux(nnull),svals(npsing),
     * u(ldu,npsing),w1(nnull),z(npsing),adiag(nuobs),tbl(ldtbl,3),
     * auxtbl(3,3),work(npsing)
c
c Purpose: determine the predictive mean squared error for each lambda 
c	value in tbl.
c
c On Entry:
c   fg(ldfg,nnull)	information on the Householder transformations 
c			that define f and g
c   ldfg		leading dimension of fg as declared
c			in the calling program
c   nuobs		number of rows in f
c   nnull		number of columns in g
c   fgaux(nnull)	auxiliary information on the fg Householder 
c			transformations 
c   svals(npsing)	singular values 
c   npsing		number of singular values
c   u(ldu,npsing)	left singular vectors corresponding to svals
c   ldu			leading dimension of u as declared in the 
c			calling program
c   w1(nnull)		leading part of rotated response vector
c   z(npsing)		u'w2
c   ntbl		number of rows in tbl
c   adiag(nuobs)	"true" y values 
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c   ldtbl		leading dimension of tbl as declared in the 
c			calling program	
c   auxtbl(3,3)		auxiliary table
c			auxtbl(1,1) contains log10(nobs*lamhat) where
c			lamhat is the gcv estimate of lambda
c
c On Exit:
c   tbl(ldtbl,3)	column	contains
c			  1 	grid of log10(nobs*lambda) 
c			  3     R(lambda) 
c   auxtbl(3,3)		auxiliary table
c			3rd column contains:
c			    [R(lamhat) , R(0), R(infinity)]'
c
c Work Arrays:
c   work(npsing)	double precision work vector
c
c Subprograms Called Directly:
c      Linpack - dqrsl 
c      Blas    - ddot dgemv
c
c $Header: dpmse.f,v 2.100.1.1 86/10/07 12:51:29 lindstrom Exp $
c
      integer i,nmh,k,locinf
      double precision dummy,nlam,wrk1,addtru,wrk
      double precision ddot
c
      nmh = nuobs - nnull
      addtru = 0.0d0
      call dqrsl(fg,ldfg,nuobs,nnull,fgaux,adiag,dummy,adiag,dummy, 
     * dummy,dummy,01000,locinf)
c			the first nnull positions of adiag now contain 
c			w1 true the last nuobs-nnull positions contain 
c			w2 true
      do 10 i = 1,nnull
         addtru=addtru + (w1(i)-adiag(i))**2
   10 continue
      addtru = addtru + ddot(nmh,adiag(nnull+1),1,adiag(nnull+1),1)
      call dgemv('T',nmh,npsing,1.0d0,u,ldu,adiag(nnull+1),1,0.0d0,
     *  work,1)
      addtru = addtru - ddot(npsing,work,1,work,1)
c			addtru contains ||w1 - (w1 true)||**2 +
c			||w2 true||**2 - ||z true||**2
c			work contains z true
c
c			compute predictive mse for each lambda in tbl
      do 30 k = 1,ntbl 
         nlam = 10**tbl(k,1)
         wrk=0.0d0
         do 20 i=1,npsing 
            wrk1=(svals(i)**2)/(svals(i)**2+nlam)
            wrk = wrk + (work(i)-z(i)*wrk1)**2
   20    continue
         tbl(k,3)=(addtru+wrk)/nobs
   30 continue
c			add pred. mse for lambda hat to auxtbl
      wrk=0.0d0
      nlam=10**auxtbl(1,1)
      do 40 i=1,npsing 
         wrk1=(svals(i)**2)/(svals(i)**2+nlam)
         wrk = wrk + (work(i)-z(i)*wrk1)**2
   40 continue
      auxtbl(1,3)=(addtru+wrk)/nobs
c			add pmse for lambda = 0
      wrk=0.0d0
      do 50 i=1,npsing 
         wrk = wrk + (work(i)-z(i))**2
   50 continue
      auxtbl(2,3)=(addtru+wrk)/nobs
c			add pmse for lambda = infinity
      auxtbl(3,3)=(addtru+ddot(npsing,work,1,work,1))/nobs
      return
      end
