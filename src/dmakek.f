      subroutine dmakekg(m,n,dim,des,lddes,nb,desb,lddesb,kk,ldkk)
      integer m,n,dim,lddes,nb,lddesb,ldkk
      double precision des(lddes,dim),desb(lddesb,dim),kk(ldkk,nb)
c
c Purpose: create the k matrix.
c   
c On Entry:
c   m			order of the derivatives in the penalty
c   n			number of rows in des
c   dim			dimension of the space to be splined 
c   des(lddes,dim)	variables to be splined		
c   lddes		leading dimension of des as declared in the 
c			calling	program
c   nb			number of rows in desb 
c   desb(lddesb,dim)	positions of unique design points or basis 
c			functions
c   lddesb		leading dimension of desb as declared in the 
c			calling	program
c   ldkk		leading dimension of kk as declared in the
c			calling	program
c On Exit:
c   kk(ldkk,nb)		k matrix
c
c Subprograms Called:
c	Other   - fact
c
c $Header: dmakekg.f,v 2.100.1.1 86/10/07 12:50:38 lindstrom Exp $
c
      integer i,j,k,fact
      double precision tauij,expo,theta,t,pi
c
c			t to be used in computation of theta
      pi = 4.0d0*atan(1.0d0)
      t = 2 ** (2*m) * pi**(dim/2.0d0) * fact(m-1)
c
c 			exponent for tauij
      expo = m - (dim / 2.0d0)
      if (dim .eq. 2*(dim/2)) then
c			1+dim odd
         theta = 1.0 / (0.5 * t * fact (m-dim/2))
         if ((2*m+dim) .eq. 4*((2*m+dim)/4)) theta = -theta
         do 30 i=1,n 
            do 20 j=1,nb 
               tauij = 0
               do 10 k=1,dim
                  tauij = tauij + (des(i,k)-desb(j,k))**2
   10          continue
               if (tauij .eq. 0.0d0) then
                  kk(i,j) = 0.0d0
               else
                  kk(i,j) = theta*tauij**expo * 0.5 * log(tauij)
               endif
   20       continue
   30    continue
      else
c			1+dim even
c			compute theta
c			compute gamma(dim/2 - m)
         j = (1 - (dim-2*m)) / 2
         theta = sqrt(pi)
         do 40 i=1,j
            theta = -theta / (i - 0.5d0)
   40    continue
         theta = theta / t

         do 70 i=1,n 
            do 60 j=1,nb 
               tauij = 0
               do 50 k=1,dim
                  tauij = tauij + (des(i,k)-desb(j,k))**2
   50          continue
               if (tauij .eq. 0.0d0) then
                  kk(i,j) = 0.0d0
               else
                  kk(i,j) = theta*tauij**expo
               endif
   60       continue
   70    continue
      endif
      end
