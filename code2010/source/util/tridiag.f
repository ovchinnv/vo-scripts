      subroutine tridiag(a,b,c,r,u,n,m)
c     V. Ovchinnikov 8/9/03. No warranty whatsoever.
      implicit none
      integer m,n,j
      real a(m),b(m),c(m),r(m),u(m)
      real dum,gam(n)
      gam(1)=b(1)
      u(1)=r(1)
      do j=2,n
       dum=a(j)/gam(j-1)
       gam(j)=b(j)-dum*c(j-1)
       u(j)=r(j)-dum*u(j-1)
      enddo
       u(n)=u(n)/gam(n)
      do j=n-1,1,-1
       u(j)=(u(j)-c(j)*u(j+1))/gam(j)
      enddo
      return
      end  
        