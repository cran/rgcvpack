      subroutine dftkf(fg,ldfg,nrf,ncg,fgaux,kk,ldkk,work)
      integer ldfg,nrf,ncg,ldkk
      double precision fg(ldfg,ncg),fgaux(ncg),kk(ldkk,nrf),work(nrf)
c
c Purpose: create f'k f.
c
c On Entry:
c   fg(ldfg,ncg)	qr decomposition of [t:s1]
c   ldfg		leading dimension of fg as declared in the
c   			calling program 
c   nrf 		number of rows in f
c   ncg			number of columns in g
c   fgaux(ncg)		auxiliary information on the qr decomposition
c			of [t:s1]
c   kk(ldkk,nrf) 	k
c   ldkk		leading dimension of kk as declared in the
c   			calling program 
c
c On Exit:
c   kk(ldkk,nrf)  	f'k f
c
c Work Array:
c   work(nrf)		double precision work vector
c
c Subprograms Called Directly:
c	Linpack - dqrsl
c	Blas    - dcopy
c
c $Header: dftkf.f,v 2.100.1.1 86/10/07 12:48:38 lindstrom Exp $
c
      double precision dummy
      integer i,locinf
c	  		calculate k f, store in kk
      do 10 i=1,nrf 
         call dcopy(nrf,kk(i,1),ldkk,work,1)
         call dqrsl(fg,ldfg,nrf,ncg,fgaux,work,dummy,work,dummy,dummy,
     *    dummy,01000,locinf)
         call dcopy(nrf,work,1,kk(i,1),ldkk)
   10 continue
c	  		calculate f'k f
      do 20 i=1,nrf 
         call dqrsl(fg,ldfg,nrf,ncg,fgaux,kk(1,i),dummy,kk(1,i),dummy,
     *    dummy,dummy,01000,locinf)
   20 continue
      return
      end
