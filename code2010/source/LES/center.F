c--------------------------------------------------------------------------
       subroutine center(u,v,w,uc,vc,wc,im,jm,km)
c
c      Tranfers variable u_i at cell center by 
c      linear interpolation
c      boundary values not defined
c
       include 'headers/common.h'
       real u(im,jm,km),uc(im,jm,km)
       real v(im,jm,km),vc(im,jm,km)
       real w(im,jm,km),wc(im,jm,km)
c
          do k=2,kz2
           do j=2,jy2
            do i=2,ix2
             uc(i,j,k)=0.5*( u(i-1,j,k)+u(i,j,k) )
             vc(i,j,k)=0.5*( v(i,j-1,k)+v(i,j,k) )
             wc(i,j,k)=0.5*( w(i,j,k-1)+w(i,j,k) )
            enddo
           enddo
          enddo
c
       return
       end

c--------------------------------------------------------------------------
       subroutine centerbc(u,v,w,uc,vc,wc,im,jm,km)
c
c      Tranfers variable u_i at cell center by 
c      linear interpolation
c      boundary values are set equal to u_i
c
       include 'headers/common.h'
       real u(im,jm,km),uc(im,jm,km)
       real v(im,jm,km),vc(im,jm,km)
       real w(im,jm,km),wc(im,jm,km)

c....streamwise
          do i=2,ix2
           do j=1,jy2+1
            do k=1,kz2+1
             uc(i,j,k)=0.5*( u(i-1,j,k)+u(i,j,k) )
            enddo
           enddo
          enddo
c....bc
           do j=1,jy2+1
            do k=1,kz2+1
             uc(1    ,j,k)= u(1    ,j,k)
             uc(ix2+1,j,k)= u(ix2+1,j,k)
            enddo
           enddo
c....spanwise          
          do i=1,ix2+1
           do j=2,jy2
            do k=1,kz2+1
             vc(i,j,k)=0.5*( v(i,j-1,k)+v(i,j,k) )
            enddo
           enddo
          enddo
c....bc
           do i=1,ix2+1
            do k=1,kz2+1
             vc(i,1    ,k)=v(i,1    ,k)
             vc(i,jy2+1,k)=v(i,jy2+1,k)
            enddo
           enddo
c....wall normal
          do i=1,ix2+1
           do j=1,jy2+1
            do k=2,kz2
             wc(i,j,k)=0.5*( w(i,j,k-1)+w(i,j,k) )
            enddo
           enddo
          enddo
c....bc
           do j=1,jy2+1
            do i=1,ix2+1
             wc(i,j,1    )=w(i,j,1    )
             wc(i,j,kz2+1)=w(i,j,kz2+1)
            enddo
           enddo
c
       return
       end


c--------------------------------------------------------------------------
       subroutine centerpol(u,v,w,uc,vc,wc,xf,yf,zf,im,jm,km)
c
c      Tranfers variable u at cell center by
c      linear interpolation
c      boundary values are calculated by 3rd order 
c      polynomial extrapolation
c
       include 'headers/common.h'
       real u(im,jm,km),uc(im,jm,km)
       real v(im,jm,km),vc(im,jm,km)
       real w(im,jm,km),wc(im,jm,km)
       real xf(im),yf(jm),zf(km)
c
       real wk1(10),wk2(10)

c.....inner values
          do i=2,ix2
           do j=2,jy2
            do k=2,kz2
             uc(i,j,k)=0.5*( u(i-1,j,k)+u(i,j,k) )
             vc(i,j,k)=0.5*( v(i,j-1,k)+v(i,j,k) )
             wc(i,j,k)=0.5*( w(i,j,k-1)+w(i,j,k) )
            enddo
           enddo
          enddo
c.....boundary values by extrapolation
c.......k=kz1
           do kk=1,3
            wk1(kk)=0.5*(zf(kk)+zf(kk+1))
           enddo
           do i=1,ix2+1
            do j=1,jy2+1
             do kk=1,3
              wk2(kk)=uc(i,j,kk+1)
             enddo
             call polint(wk1,wk2,3,zf(1),yy,dyy)
             uc(i,j,1)=yy
            enddo
           enddo
c......k=kz2+1
           do kk=1,3
            wk1(kk)=0.5*(zf(kz2+1-kk)+zf(kz2-kk))
           enddo
           do i=1,ix2+1
            do j=1,jy2+1
             do kk=1,3
              wk2(kk)=uc(i,j,kz2+1-kk)
             enddo
             call polint(wk1,wk2,3,zf(kz2),yy,dyy)
             uc(i,j,kz2+1)=yy
            enddo
           enddo
c.......i=ix1
           do ii=1,3
            wk1(ii)=0.5*(xf(ii)+xf(ii+1))
           enddo
           do k=1,kz2+1
            do j=1,jy2+1
             do ii=1,3
              wk2(ii)=uc(ii+1,j,k)
             enddo
             call polint(wk1,wk2,3,xf(1),yy,dyy)
             uc(1,j,k)=yy
            enddo
           enddo
c......i=ix2+1
           do ii=1,3
            wk1(ii)=0.5*(xf(ix2+1-ii)+xf(ix2-ii))
           enddo
           do k=1,kz2+1
            do j=1,jy2+1
             do ii=1,3
              wk2(ii)=uc(ix2+1-ii,j,k)
             enddo
             call polint(wk1,wk2,3,xf(ix2),yy,dyy)
             uc(ix2+1,j,k)=yy
            enddo
           enddo
c......j=1 and j=jy2+1
         do i=1,ix2+1
          do k=1,kz2+1
           uc(i,1,k)=uc(i,jy2,k)
           uc(i,jy2+1,k)=uc(i,2,k)
          enddo
         enddo

       return
       end


