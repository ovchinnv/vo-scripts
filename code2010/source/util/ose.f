      subroutine ose_(Re,omega,ky,U,d2U,phi,phip,y,d,dc,ny)
c     DO NOT MODIFY!     
c     purpose: to obtain an eigen-function from the continuous spectrum of the Orr-Sommerfeld eq.
c     Author: Victor Ovchinnikov 3/14/03. No warrany whatsoever.
      implicit none
      integer ny
      real pi,beta
      real omega,Re,ky,a_r,a_i
      complex c,alpha
      complex phi(ny),d2U(ny),U(ny),phip(ny)
      real y(ny),d(ny),dc(ny)
ccc   complex matrix coefficients
      complex a(ny),b(ny),cc(ny),dd(ny),e(ny)
      integer i
      complex coeff
      real Ap,Bp,Cp,Dp,BBp,AAp
c     
      pi=4.0*atan(1.)
cccc  initialize (do not remove unless you know what you are doing)
      phi(1)=cmplx(pi/3.,1)
      a(1)=cmplx(pi/4.,1)
cccc  initialize
c     
ccccccccccccccccccccccccccccccccc
c     Orr-Sommerfeld parameters
c      write(6,*) 'Re=',Re
c      do i=1,ny
c       write(11,*) y(i),real(U(i))
c      enddo 
c     ky=0.5*pi
      beta=0.5*(1+(2.*ky/Re)**2)
      a_r=0.5*Re*sqrt(sqrt(beta**2+4*(omega/Re)**2)-beta)
      a_i=0.5*Re*(omega/a_r-1)
      alpha=cmplx(a_r,a_i)
      c=cmplx(omega,0)/alpha
c     
c     now we can build the matrix
c     
      do i=3, ny-2
c     4th derivative:
         a(i)=d(i-2)*dc(i-1)*d(i-1)*dc(i)
         b(i)=-(dc(i)*d(i)+(d(i-1)*dc(i)+(d(i-1)+d(i-2))*dc(i-1)))*d(i-1)*dc(i)
         cc(i)=((d(i)*dc(i+1)+(d(i)+d(i-1))*dc(i))*d(i)+
     &        ((d(i)+d(i-1))*dc(i)+d(i-1)*dc(i-1))*d(i-1))*dc(i)
         dd(i)=-((d(i+1)+d(i))*dc(i+1)+d(i)*dc(i)+d(i-1)*dc(i))*d(i)*dc(i)
         e(i)=d(i)*d(i+1)*dc(i+1)*dc(i)
c     2nd derivative:
         coeff=-(2*alpha**2+(0,1)*alpha*Re*(U(i)-c))
         b(i)=b(i)+coeff*d(i-1)*dc(i)
         cc(i)=cc(i)-coeff*(d(i)+d(i-1))*dc(i)
         dd(i)=dd(i)+coeff*d(i)*dc(i)
c     0th derivative:
         coeff=alpha*((0,1)*alpha**2*Re*(U(i)-c)+(0,1)*Re*d2U(i)+alpha**3)
         cc(i)=cc(i)+coeff
         phi(i)=0.
      enddo
c     apply boundary conditions
c     assume phi'(1)=C*phi(1)+A*phi(2)+B*phi(3)
c     
      BBp=-d(2)**2/(d(1)+d(2))
      AAp=(1.-BBp*(1./d(1)+1./d(2)))*d(1)
c     C=-A-B
c     phi(2)=-B/A*phi(3)
      cc(3)=cc(3)-BBp/AAp*b(3)
      b(4)=b(4)-BBp/AAp*a(4)
c     free-stream (Jacobs & Durbin)  assume y(inf)=1 (normalization)
      coeff=exp(-alpha/d(ny-2))
      Dp=d(ny-1)*dc(ny-1)
      Cp=(ky**2-(d(ny-1)+d(ny-2))*dc(ny-1)-coeff*d(ny-2)*dc(ny-2))
      Bp=(d(ny-2)*dc(ny-1)+coeff*((d(ny-2)+d(ny-3))*dc(ny-2)-ky**2))
      Ap=-coeff*d(ny-3)*dc(ny-2)
c     
      phi(ny-2)=phi(ny-2)-e(ny-2)*1.
      cc(ny-2)=cc(ny-2)-Bp/Cp*dd(ny-2)
      b(ny-2)=b(ny-2)-Ap/Cp*dd(ny-2)
      phi(ny-2)=phi(ny-2)+Dp/Cp*dd(ny-2)
c     
      dd(ny-3)=dd(ny-3)-Bp/Cp*e(ny-3)
      cc(ny-3)=cc(ny-3)-Ap/Cp*e(ny-3)          
      phi(ny-3)=phi(ny-3)+Dp/Cp*e(ny-3)
c     
cccc  call pentadiagonal inversion routine
      call pentdiag(a(3),b(3),cc(3),dd(3),e(3),phi(3),ny-4)
      phi(1)=(0,0)
      phi(2)=-BBp/AAp*phi(3)
      phi(ny)=1.
      phi(ny-1)=(-Ap*phi(ny-3)-Bp*phi(ny-2)-Dp)/Cp
c     compute derivative (do not forget to multiply by 1./L to get d_phi/d_y)
      do i=2,ny
         phip(i)=(phi(i)-phi(i-1))*d(i-1)
      enddo
      phip(1)=0.                ! neumann at the wall by assumption
c     debug      
c      write(321,'(2 G25.15)') (y(i),real(phi(i)),i=1,ny)
c      close(321)
      return
      end
