c--------------------------------------------------------------
      subroutine filter2d(f,fb,wk1,wk2,im,jm,km,ifil)
c
cc box filter routine (trapezodial rule for UNIFORM grid)
cc no filtering in the wall normal direction (Z)
c
      include 'headers/common.h'
c
      real f(im,jm,km),fb(im,jm,km)
c
      real fx1d(1000),fy1d(1000)
      real wk1(im,jm),wk2(im,jm)
      real wh(-10:10)
      integer nfck,nf,imax,jmax,kmaxp
c
      if(ifil.eq.1) then
        nf=1
        vtot=64.
        wh(-1)=1.
        wh(0)=6.
        wh(1)=1.       
      elseif(ifil.eq.2) then
        nf=1
        vtot=16.
        wh(-1)=1.
        wh(0)=2.
        wh(1)=1. 
      else
       write(6,*) 'filter type is 1 or 2'
      endif
c
      do 1000 k=2,kz2
c intialize arrays
      do 1 i=1,ix2+1
      do 1 j=1,jy2+1
       wk1(i,j)=0.0
       wk2(i,j)=0.0
 1    continue
c--- filter in x ---
       do 100 j=2,jy2
c fill scratch 1d array
        do 10 i=1,ix2+1
         fx1d(i)=f(i,j,k)
 10   continue
c apply filter
        do 50 ii=-nf,nf
         do 50 i=2,ix2
          wk1(i,j)=wk1(i,j)+wh(ii)*fx1d(i+ii)
 50    continue
c
 100  continue
c
c--- filter in y ---
       do 200 i=2,ix2
c fill scratch array
        do 20 j=1,jy2+1
         fy1d(j)=wk1(i,j)
 20   continue
c apply filter
        do 60 jj=-nf,nf
         do 60 j=2,jy2
          wk2(i,j)=wk2(i,j)+wh(jj)*fy1d(j+jj)
 60    continue
c
 200  continue
c
       do i=2,ix2
       do j=2,jy2
       fb(i,j,k)=wk2(i,j)/vtot
       enddo
       enddo

 1000 continue
c
      return
      end
c-------------------------------------------------------------



c--------------------------------------------------------------
      subroutine filter3d(f,fb,wk1,wk2,im,jm,km,ifil)
c
cc box filter routine (trapezodial rule for NON-UNIFORM grid)
cc On the boundaries backward and forward formulas are used
cc fb is not defined on the boundaries on output
cc 
c
      include 'headers/common.h'
c
      real f(im,jm,km),fb(im,jm,km)
      real wk1(im,jm,km),wk2(im,jm,km)
      integer ifil,ii,jj,kk

c--- initialize work arrays
      do k=kz1,kz2+1
       do j=jy1,jy2+1
        do i=ix1,ix2+1
         wk1(i,j,k)=0.0
         wk2(i,j,k)=0.0
        enddo
       enddo
      enddo

c--- filter in x ---
      do k=kz1,kz2+1
       do j=jy1+1,jy2
        do i=3,ix2-1
           wk1(i,j,k)= whx(i,-1,ifil)*f(i-1,j,k)
     &                +whx(i, 0,ifil)*f(i  ,j,k)
     &                +whx(i, 1,ifil)*f(i+1,j,k)
        enddo
       enddo
      enddo 
c..boundary points
      do k=kz1,kz2+1
       do j=jy1+1,jy2
           wk1(2  ,j,k)=whx(2,-1,ifil)*f(2,j,k)
     &                 +whx(2, 0,ifil)*f(3,j,k)
     &                 +whx(2, 1,ifil)*f(4,j,k)
           wk1(ix2,j,k)=whx(ix2,-1,ifil)*f(ix2,j,k)
     &                 +whx(ix2, 0,ifil)*f(ix2-1,j,k)
     &                 +whx(ix2, 1,ifil)*f(ix2-2,j,k)
       enddo
      enddo    


c--- filter in y ---
      do k=kz1,kz2+1
       do j=3,jy2-1      
        do i=ix1+1,ix2
           wk2(i,j,k)= why(j,-1,ifil)*wk1(i,j-1,k)
     &                +why(j, 0,ifil)*wk1(i,j  ,k)
     &                +why(j, 1,ifil)*wk1(i,j+1,k)
        enddo
       enddo
      enddo
c..boundary points
      do k=kz1,kz2+1
        do i=ix1+1,ix2
           wk2(i  ,2,k)=why(2,-1,ifil)*wk1(i,2,k)
     &                 +why(2, 0,ifil)*wk1(i,3,k)
     &                 +why(2, 1,ifil)*wk1(i,4,k)
           wk2(i,jy2,k)=why(jy2,-1,ifil)*wk1(i,jy2,k)
     &                 +why(jy2, 0,ifil)*wk1(i,jy2-1,k)
     &                 +why(jy2, 1,ifil)*wk1(i,jy2-2,k)
c
       enddo
      enddo


c--- filter in z ---
      do k=3,kz2-1      
       do j=jy1+1,jy2
        do i=ix1+1,ix2
           fb(i,j,k)= whz(k,-1,ifil)*wk2(i,j,k-1)
     &               +whz(k, 0,ifil)*wk2(i,j,k  )
     &               +whz(k, 1,ifil)*wk2(i,j,k+1)
        enddo
       enddo
      enddo
c..boundary points
      IF(itype(5).GT.0) THEN    !BOTTOM BOUNDARY
         do j=jy1+1,jy2
            do i=ix1+1,ix2 
               fb(i,j,  2)=whz(2,-1,ifil)*wk2(i,j,2)
     &              +whz(2, 0,ifil)*wk2(i,j,3)
     &              +whz(2, 1,ifil)*wk2(i,j,4)
            enddo
         enddo
      ELSE   !GHOST LAYER
         do j=jy1+1,jy2
            do i=ix1+1,ix2 
               fb(i,j,2)=whz(2,-1,ifil)*wk2(i,j,1)
     &              +whz(2,0,ifil)*wk2(i,j,2)
     &              +whz(2,1,ifil)*wk2(i,j,3)
            enddo
         enddo
      ENDIF
      IF(itype(6).GT.0) THEN    !TOP BOUNDARY
         do j=jy1+1,jy2
            do i=ix1+1,ix2 
               fb(i,j,kz2)=whz(kz2,-1,ifil)*wk2(i,j,kz2)
     &              +whz(kz2, 0,ifil)*wk2(i,j,kz2-1)
     &              +whz(kz2, 1,ifil)*wk2(i,j,kz2-2)
            enddo
         enddo
      ELSE                      !GHOST LAYER
         do j=jy1+1,jy2
            do i=ix1+1,ix2 
               fb(i,j,kz2)=whz(kz2,-1,ifil)*wk2(i,j,kz2-1)
     &              +whz(kz2, 0,ifil)*wk2(i,j,kz2)
     &              +whz(kz2, 1,ifil)*wk2(i,j,kz2+1)
            enddo
         enddo
      ENDIF
      
c

c
      return
      end
c-------------------------------------------------------------


c--------------------------------------------------------------
      subroutine weights(xc,yc,zc,im,jm,km)
c
cc calculates the weights on each direction 
cc for a box filter on non-uniform grids
cc trapezodial rule is used
c
      include 'headers/common.h'
c
      real xc(im),yc(jm),zc(km)
      real dxl,dxr,dyl,dyr,dzl,dzr

c--------------------------------------------------------------------
c                   Weights for bar filter (DELTA)
c--------------------------------------------------------------------
c        
c--------------------------  x-direction  ---------------------------
       do i=3,ix2-1
        dxl=0.5*(xc(i)-xc(i-1))
        dxr=0.5*(xc(i+1)-xc(i))
        dxt=4.*(dxl+dxr)
        whx(i,-1, 1)=dxr/dxt
        whx(i, 0, 1)=3./4.
        whx(i, 1, 1)=dxl/dxt
       enddo
c...backward and forward formulas for bc points
        whx(2,-1, 1)= 1.0
        whx(2, 0, 1)= 0.0
        whx(2, 1, 1)= 0.0
c
        whx(ix2,-1, 1)=1.0 
        whx(ix2, 0, 1)=0.0
        whx(ix2, 1, 1)=0.0
c--------------------------  y-direction  ---------------------------
       do j=3,jy2-1
        dyl=0.5*(yc(j)-yc(j-1))
        dyr=0.5*(yc(j+1)-yc(j))
        dyt=4.*(dyl+dyr)
        why(j,-1, 1)=dyr/dyt
        why(j, 0, 1)=3./4.
        why(j, 1, 1)=dyl/dyt
       enddo
c...backward and forward formulas for bc points
        why(2,-1, 1)= 1.0
        why(2, 0, 1)= 0.0
        why(2, 1, 1)= 0.0
c
        why(jy2,-1, 1)= 1.0
        why(jy2, 0, 1)= 0.0
        why(jy2, 1, 1)= 0.0
c--------------------------  z-direction  ---------------------------
        
       do k=3,kz2-1
        dzl=0.5*(zc(k)-zc(k-1))
        dzr=0.5*(zc(k+1)-zc(k))
        dzt=4.*(dzl+dzr)
        whz(k,-1, 1)=dzr/dzt
        whz(k, 0, 1)=3./4.
        whz(k, 1, 1)=dzl/dzt
       enddo
c...backward and forward formulas for bc points
       IF(itype(5).GT.0) THEN !Bottom boundary
          whz(2,-1, 1)= 1.0
          whz(2, 0, 1)= 0.0
          whz(2, 1, 1)= 0.0
       ELSE
          k=2
          dzl=0.5*(zc(k)-zc(k-1))
          dzr=0.5*(zc(k+1)-zc(k))
          dzt=4.*(dzl+dzr)
          whz(k,-1, 1)=dzr/dzt
          whz(k, 0, 1)=3./4.
          whz(k, 1, 1)=dzl/dzt
       ENDIF
       
c     
       IF(itype(6).GT.0) THEN !TOP BOUNDARY
          whz(kz2,-1, 1)= 1.0
          whz(kz2, 0, 1)= 0.0
          whz(kz2, 1, 1)= 0.0
       ELSE
          k=kz2
          dzl=0.5*(zc(k)-zc(k-1))
          dzr=0.5*(zc(k+1)-zc(k))
          dzt=4.*(dzl+dzr)
          whz(k,-1, 1)=dzr/dzt
          whz(k, 0, 1)=3./4.
          whz(k, 1, 1)=dzl/dzt
       ENDIF


c--------------------------------------------------------------------
c                   Weights for hat filter (2*DELTA)
c--------------------------------------------------------------------
c        
c--------------------------  x-direction  ---------------------------
       do i=3,ix2-1
        dxl=xc(i)-xc(i-1)
        dxr=xc(i+1)-xc(i)
        dxt=2.*(dxl+dxr)
        whx(i,-1, 2)=dxr/dxt
        whx(i, 0, 2)=0.5
        whx(i, 1, 2)=dxl/dxt
       enddo
c...backward and forward formulas for bc points
        dxa =xc(3)-xc(2)
        dxb =xc(4)-xc(2)
        whx(2,-1, 2)= (dxa+dxb)/(2.*dxb)
        whx(2, 0, 2)= 0.5
        whx(2, 1, 2)=-dxa/(2.*dxb)
c
        dxa =xc(ix2)-xc(ix2-1)
        dxb =xc(ix2)-xc(ix2-2)
        whx(ix2,-1, 2)= (dxa+dxb)/(2.*dxb)
        whx(ix2, 0, 2)= 0.5
        whx(ix2, 1, 2)=-dxa/(2.*dxb)
c--------------------------  y-direction  ---------------------------
       do j=3,jy2-1
        dyl=yc(j)-yc(j-1)
        dyr=yc(j+1)-yc(j)
        dyt=2.*(dyl+dyr)
        why(j,-1, 2)=dyr/dyt
        why(j, 0, 2)=0.5
        why(j, 1, 2)=dyl/dyt
       enddo
c...backward and forward formulas for bc points
        dya =yc(3)-yc(2)
        dyb =yc(4)-yc(2)
        why(2,-1, 2)= (dya+dyb)/(2.*dyb)
        why(2, 0, 2)= 0.5
        why(2, 1, 2)=-dya/(2.*dyb)
c
        dya =yc(jy2)-yc(jy2-1)
        dyb =yc(jy2)-yc(jy2-2)
        why(jy2,-1, 2)= (dya+dyb)/(2.*dyb)
        why(jy2, 0, 2)= 0.5
        why(jy2, 1, 2)=-dya/(2.*dyb)
c--------------------------  z-direction  ---------------------------
       do k=3,kz2-1
        dzl=zc(k)-zc(k-1)
        dzr=zc(k+1)-zc(k)
        dzt=2.*(dzl+dzr)
        whz(k,-1, 2)=dzr/dzt
        whz(k, 0, 2)=0.5
        whz(k, 1, 2)=dzl/dzt
       enddo
c...backward and forward formulas for bc points
       IF(itype(5).GT.0) THEN !BOTTOM BOUNDARY
          dza =zc(3)-zc(2)
          dzb =zc(4)-zc(2)
          whz(2,-1, 2)= (dza+dzb)/(2.*dzb)
          whz(2, 0, 2)= 0.5
          whz(2, 1, 2)=-dza/(2.*dzb)
       ELSE  !GHOST LAYER
          k=2
          dzl=zc(k)-zc(k-1)
          dzr=zc(k+1)-zc(k)
          dzt=2.*(dzl+dzr)
          whz(k,-1, 2)=dzr/dzt
          whz(k, 0, 2)=0.5
          whz(k, 1, 2)=dzl/dzt
       ENDIF
c
       IF(itype(6).GT.0) THEN  !TOP BOUNDARY
          dza =zc(kz2)-zc(kz2-1)
          dzb =zc(kz2)-zc(kz2-2)
          whz(kz2,-1, 2)= (dza+dzb)/(2.*dzb)
          whz(kz2, 0, 2)= 0.5
          whz(kz2, 1, 2)=-dza/(2.*dzb)
       ELSE
          k=kz2
          dzl=zc(k)-zc(k-1)
          dzr=zc(k+1)-zc(k)
          dzt=2.*(dzl+dzr)
          whz(k,-1, 2)=dzr/dzt
          whz(k, 0, 2)=0.5
          whz(k, 1, 2)=dzl/dzt
       ENDIF
c
      

        return
        end


