c.......................................................................
c
c     Initialize the initializable .............................
c
c.......................................................................
      subroutine peqina(xu, zw, au, cw, am, bm, cm, an, bn, cn,nx,ny,nz)
c
      implicit none
c
c     ...Parameters...
c
      integer nx,ny,nz
c     
c     ...Scalar Arguments...
c
      integer  mmx, mmz
c
c     ...Array Arguments...
c
      real             xu(nx), zw(nz)
      real             au(nx), cw(nz)
      double precision am(nx-2), bm(nx-2), cm(nx-2)
      double precision an(nz-2), bn(nz-2), cn(nz-2)
c
c     ...Local Scalars...
c
      integer  i, k
c.......................................................................
c     Executable Statements
c.......................................................................
c
c     create matrix for BLKTRI ......................................
c
      do i = 2, nx-1
         am(i-1) = DBLE(au(i-1))/DBLE(xu(i)-xu(i-1))
         cm(i-1) = DBLE(au(i))  /DBLE(xu(i)-xu(i-1))
         bm(i-1) = - am(i-1) - cm(i-1)
      end do
c
      do k = 2, nz-1
         an(k-1) = DBLE(cw(k-1))/DBLE(zw(k)-zw(k-1))
         cn(k-1) = DBLE(cw(k))  /DBLE(zw(k)-zw(k-1))
         bn(k-1) = - an(k-1) - cn(k-1)
      end do
c      
      return
      end
c
c.......................................................................
c
      subroutine peqinb(divv, dp, am, bm, cm, an, bn, cn,
     *                  ak, dy, bcpn, bcps, bcpw, bcpe,nx,ny,nz)
c
      implicit none
c
c     ...Parameters...
c
      integer nx,ny,nz
c     
c     ...Scalar Arguments...
c
      integer          bcpn, bcps, bcpw, bcpe
      double precision dy
c
c     ...Array Arguments...
c
      real             divv(nx,ny,nz), dp(nx,ny,nz)
      double precision ak(ny)
      double precision am(nx-2), bm(nx-2), cm(nx-2)
      double precision an(nz-2), bn(nz-2), cn(nz-2)
c
c     ...Local Scalars...
c
      integer          i, j, k, l, n2mh
      double precision dy2q, pi
c.......................................................................
c     Executable Statements
c.......................................................................
      pi = 4.0*atan(1.0)
c
c     Neuman boundary conditions for direct solver .................
c
c     ... south ....................................................
      if (bcps.eq.1) then
         k = 1
         an(k) = 0.
         bn(k) = -cn(k)
      else
         stop 'Dirichlet bound for p not implemented yet'
      end if
c     ... north ....................................................
      if (bcpn .eq. 1) then
         k = nz-2
         cn(k) = 0.
         bn(k) = -an(k)
      else
         k     = nz-2
         cn(k) = 0.
         stop 'Dirichlet bound for p not implemented yet'
      end if
c     ... west .....................................................
      if (bcpw.eq.1) then
         i = 1
         am(i) = 0.
         bm(i) = -cm(i)
      else
         stop 'Dirichlet bound for p not implemented yet'
      end if
c     ... east .....................................................
      if (bcpe.eq.1) then
         i = nx-2
         cm(i) = 0.
         bm(i) = -am(i)
      else
         stop 'Dirichlet bound for p not implemented yet'
      end if
c.......................................................................
c     Modified wave numbers
c.......................................................................
c     Contribute to the diagonal real coefficient from j direction
      dy2q = 1.0/(dy*dy)
      n2mh = (ny-2)/2
      do l=1, n2mh
        ak(l) = cos(2.*pi*l/real(ny-2))
        ak(l) = 2.*(1.-ak(l))*dy2q
c       print*, 'l',l,'    ak',ak(l)
      end do
      do l=n2mh+1, (ny-2)-1
        ak(l) = ak(ny-2-l)
c       print*, 'l',l,'    ak',ak(l)
      end do
c
      return
      end
c
c.......................................................................
