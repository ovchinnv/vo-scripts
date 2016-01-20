      subroutine pentdiag(a,b,c,d,e,r,n)
      implicit none
c     DO NOT MODIFY! V. Ovchinnikov, 2003.  No Warranty whatsoever.
      integer n
      complex a(n),b(n),c(n),d(n),e(n),r(n)
      integer i
      complex dummy
      do i=1,n-2
         dummy=b(i+1)/c(i)
         c(i+1)=c(i+1)-d(i)*dummy
         d(i+1)=d(i+1)-e(i)*dummy
         r(i+1)=r(i+1)-r(i)*dummy
         dummy=a(i+2)/c(i)
         b(i+2)=b(i+2)-d(i)*dummy
         c(i+2)=c(i+2)-e(i)*dummy
         r(i+2)=r(i+2)-r(i)*dummy
      enddo
      dummy=b(n)/c(n-1)
      c(n)=c(n)-d(n-1)*dummy
      r(n)=(r(n)-r(n-1)*dummy)/c(n)
      do i=n,3,-1
         dummy=r(i)
         r(i-1)=(r(i-1)-dummy*d(i-1))/c(i-1)
         r(i-2)=r(i-2)-dummy*e(i-2)
      enddo
      r(1)=(r(1)-r(2)*d(1))/c(1)
      return
      end
