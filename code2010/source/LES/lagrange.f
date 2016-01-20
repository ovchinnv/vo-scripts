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
c      write(0,*) myrank,' entered lagrange, dt=',dt
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
c         write(0,*) 'x-uDt out of bounds',xpar,uc(i,j,k),dt
c      endif
c      if(ypar.le.yc(1).or.ypar.ge.yc(jy2+1)) then
c         write(0,*) 'y-uDt out of bounds',ypar,vc(i,j,k),dt
c       endif
c      if(zpar.le.zc(1).or.zpar.ge.zc(kz2+1)) then
c         write(0,*) 'z-uDt out of bounds',zpar,wc(i,j,k),dt
c      endif

c    commented for speed

c       particle cannot travel more than one cell!
c        idir=int(uc(i,j,k)/abs(uc(i,j,k)))
	ixpos=i
c        jdir=int(vc(i,j,k)/abs(vc(i,j,k)))
	jypos=j
c        kdir=int(wc(i,j,k)/abs(wc(i,j,k)))
	kzpos=k
c
c       goto 999
	do while (.not.(xc(ixpos).le.xpar.and.xc(ixpos+1).ge.xpar))
	 ixpos=ixpos-1
c
c	 if ((ixpos.lt.1).or.ixpos.gt.ix2) write(0,'(5I4,2G20.10,I2)') myrank,ixpos,i,j,k,uc(i,j,k),dt,1
c
        enddo
	do while (.not.(yc(jypos).le.ypar.and.yc(jypos+1).ge.ypar))
	 jypos=jypos-1
	 
c	 if ((jypos.lt.1).or.jypos.gt.jy2) write(0,'(5I4,2G20.10,I2)') myrank,jypos,i,j,k,vc(i,j,k),dt,2
        enddo
	do while (.not.(zc(kzpos).le.zpar.and.zc(kzpos+1).ge.zpar))
	 kzpos=kzpos-1
c	 if ((kzpos.lt.1).or.kzpos.gt.kz2) write(0,'(5I4,2G20.10,I2)') myrank,kzpos,i,j,k,wc(i,j,k),dt,3
	 
        enddo
c        if (ixpos.lt.1) write(0,*) 'warning: lagrange: ixpos', ixpos,dt
c        if (jypos.lt.1) write(0,*) 'warning: lagrange: jypos', jypos,dt
c        if (kzpos.lt.1) write(0,*) 'warning: lagrange: kzpos', kzpos,dt
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
c	 if (kzpos.gt.km-1.or.kzpos.eq.0) then 
c	 write(0,*) 'kzpos:',kzpos,zpar-zc(k),k,dt,wc(i,j,k),1./cp(k),'***************************************'
c          do kk=2,kz2
c	   write(400+myrank,'(3I3,3G20.10)') i,j,kk,zc(kk),wc(i,j,kk),dt
c	  enddo
c	  close(400+myrank) 
c        stop
c         endif
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
c      write(0,*) myrank,' quit lagrange'
         RETURN
         END


