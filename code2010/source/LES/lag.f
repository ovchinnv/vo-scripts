c---------------------------------------------------------------------
      subroutine lagrange(vin1,vin2,vout1,vout2,
     &                    uc,vc,wc,xc,yc,zc,dt,im,jm,km)
c
c     given intergal at i,j,k returns value at (n-1) along 
c     fluid pathline
c
      include 'headers/common.h'
c
      REAL vin1(im,jm,km),vout1(im,jm,km)
      REAL vin2(im,jm,km),vout2(im,jm,km)
      REAL uc(im,jm,km),vc(im,jm,km),wc(im,jm,km)
      REAL xc(im),yc(jm),zc(km)
      REAL dt
      
      integer  idir,jdir,kdir
c
      REAL    xpar,ypar,zpar
      REAL    csi,xni,eta
      REAL    fn1,fn2,fn3,fn4,fn5,fn6,fn7,fn8
      INTEGER ixpos,jypos,kzpos
c
      write(0,*) myrank,' entered lagrange'
      DO 10 i=2,ix2
       DO 20 j=2,jy2
        DO 30 k=2,kz2
c...particle possition at t-dt
      xpar=xc(i)-uc(i,j,k)*dt
      ypar=yc(j)-vc(i,j,k)*dt
      zpar=zc(k)-wc(i,j,k)*dt
c...find lower corner of the cell containing the particle
c....x-dir

c   commented for speed
c      if(xpar.le.xc(1).or.xpar.ge.xc(ix2+1)) then
c         write(6,*) 'x-uDt out of bounds',xpar
c      endif
c      if(ypar.le.yc(1).or.ypar.ge.yc(jy2+1)) then
c         write(6,*) 'y-uDt out of bounds',ypar
c       endif
c      if(zpar.le.zc(1).or.zpar.ge.zc(kz2+1)) then
c         write(6,*) 'z-uDt out of bounds',zpar
c      endif

c      goto 999
c    commented for speed

c       particle cannot travel more than one cell!
c        idir=int(uc(i,j,k)/abs(uc(i,j,k)))
	ixpos=i+1
c        jdir=int(vc(i,j,k)/abs(vc(i,j,k)))
	jypos=j+1
c        kdir=int(wc(i,j,k)/abs(wc(i,j,k)))
	kzpos=k+1
c
	do while (.not.(xc(ixpos).le.xpar.and.xc(ixpos+1).ge.xpar))
	 ixpos=ixpos-1
        enddo
	do while (.not.(yc(jypos).le.ypar.and.yc(jypos+1).ge.ypar))
	 jypos=jypos-1
        enddo
	do while (.not.(zc(kzpos).le.zpar.and.zc(kzpos+1).ge.zpar))
	 kzpos=kzpos-1
        enddo
c        if (abs(ixpos-i).gt.1) write(0,*) 'warning: lagrange: ', ixpos,i,idir
c        if (abs(jypos-j).gt.1) write(0,*) 'warning: lagrange: ', jypos,j,jdir
c        if (abs(kzpos-k).gt.1) write(0,*) 'warning: lagrange:,', kzpos,k,kdir
 999  continue
c        call locate(xc,im,xpar,ixpos)
c        call locate(yc,jm,ypar,jypos)
c        call locate(zc,km,zpar,kzpos)

c
c....particle is located between ixpos and ixpos+1
c                                jypos and jypos+1
c                                kzpos and jypos+1
c...interlopate between above values to find varout
c....math scace coordinates
         csi=(xpar-xc(ixpos))/(xc(ixpos+1)-xc(ixpos))
         xni=(ypar-yc(jypos))/(yc(jypos+1)-yc(jypos))
         eta=(zpar-zc(kzpos))/(zc(kzpos+1)-zc(kzpos))
c....weights
c     point 1: i,j,k
         fn1=(1.-csi)*(1.-xni)*(1.-eta)
c     point 2: i,j+1,k
         fn2=(1.-csi)*    xni *(1.-eta)
c     point 3: i,j+1,k+1  
         fn3=(1.-csi)*    xni*     eta
c     point 4: i,j,k+1
         fn4=(1.-csi)*(1.-xni)*    eta
c     point 5: i+1,j,k
         fn5=    csi* (1.-xni)*(1.-eta)
c     point 6: i+1,j+1,k
         fn6=    csi*     xni *(1.-eta)
c     point 7: i+1,j+1,k+1
         fn7=    csi*     xni *    eta
c     point 8: i+1,j,k+1
         fn8=    csi* (1.-xni)*    eta
c     values
         vout1(i,j,k)=vin1(ixpos  ,jypos  ,kzpos  )*fn1+
     &                vin1(ixpos  ,jypos+1,kzpos  )*fn2+
     &                vin1(ixpos  ,jypos+1,kzpos+1)*fn3+
     &                vin1(ixpos  ,jypos  ,kzpos+1)*fn4+
     &                vin1(ixpos+1,jypos  ,kzpos  )*fn5+
     &                vin1(ixpos+1,jypos+1,kzpos  )*fn6+
     &                vin1(ixpos+1,jypos+1,kzpos+1)*fn7+
     &                vin1(ixpos+1,jypos  ,kzpos+1)*fn8
c
         vout2(i,j,k)=vin2(ixpos  ,jypos  ,kzpos  )*fn1+
     &                vin2(ixpos  ,jypos+1,kzpos  )*fn2+
     &                vin2(ixpos  ,jypos+1,kzpos+1)*fn3+
     &                vin2(ixpos  ,jypos  ,kzpos+1)*fn4+
     &                vin2(ixpos+1,jypos  ,kzpos  )*fn5+
     &                vin2(ixpos+1,jypos+1,kzpos  )*fn6+
     &                vin2(ixpos+1,jypos+1,kzpos+1)*fn7+
     &                vin2(ixpos+1,jypos  ,kzpos+1)*fn8

 30   CONTINUE
 20   CONTINUE
 10   CONTINUE
      write(0,*) myrank,' quit lagrange'
         RETURN
         END


