C
C 
C-----SUBROUTINE-ObjectGrdWann---------------------E. Balaras..04/14/02-----
C
      subroutine ObjectGrdWann(sb,xb,yb,coefx,coefy,nb,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
C
      include 'headers/immbd.h'
      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
c
      real    r,xc,zc,pi
      real    step,theta
      real    cofx(3),cofy(3)
      integer ns,nb      
c
c----------------------------------------------------------------------
      r = 1.   !  radious of the circle 
      xc = 3.5  !  x-coordinate of center
      zc = 3.5  !  z-coordinate of center
      ns = nimb   !  total number of points

      nb =ns-1  !  actual number of points (1 and ns coinside) 
      pi=3.14159265358979 
c----------------------------------------------------------------------
c....compute s,x(s),y(s)
      step = 2.*pi/float(ns-1)
      theta=-1.*step
      do i=1,ns
       theta= theta+step
       sb(i) = r*theta
       xb(i) = xc+r*sin(theta)
       yb(i) = zc+r*cos(theta)
c       write(123,*) i,xb(i),yb(i)
      enddo  
c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,nb
       call polcoe(sb(i-1),xb(i-1),3,cofx)
       call polcoe(sb(i-1),yb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c
      return
      end

C
C 
C-----SUBROUTINE-ObjectGrdCyl-------------------------E. Balaras..03/04/01-----
C
      subroutine ObjectGrdCyl(sb,xb,yb,coefx,coefy,nb,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
C
      include 'headers/immbd.h'
      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
c
      real    r,xc,zc,pi
      real    step,theta
      real    cofx(3),cofy(3)
      integer ns,nb      
c
c----------------------------------------------------------------------
c....input for the circle
      r = 0.5   !  radius of the circle 
      xc = 3.5  !  x-coordinate of center
      zc = 10.  !  z-coordinate of center
      ns = nimb !  total number of points

c      r = 1.   !  radius of the circle 
c      xc = 3.5  !  x-coordinate of center
c      zc = 3.5  !  z-coordinate of center
c      ns = 120.   !  total number of points

      nb =ns-1  !  actual number of points (1 and ns coinside) 
      pi=3.14159265358979 
c----------------------------------------------------------------------
c....compute s,x(s),y(s)
      step = 2.*pi/float(ns-1)
      theta=-1.*step
      do i=1,ns
       theta= theta+step
       sb(i) = r*theta
       xb(i) = xc+r*sin(theta)
       yb(i) = zc+r*cos(theta)
c       write(123,*) i,xb(i),yb(i)
      enddo  
c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,nb
       call polcoe(sb(i-1),xb(i-1),3,cofx)
       call polcoe(sb(i-1),yb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c
      return
      end
C
C 
C-----SUBROUTINE-ObjectOscCyl-------------------------E. Balaras..11/12/02-----
C
      subroutine ObjectOscCyl(sb,xb,yb,coefx,coefy,nb,nmax,
     &                        ubd,vbd,wbd,tlevel           )
C
C     finds location of the object at each time step
C     builds local quadratic functions
C     for an oscilating cylinder with prescribed motion
C
      include 'headers/immbd.h'

      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
c
      real    xref,zref
      real    r,xc,zc,pi
      real    xcr,zcr,amp,fr
      real    step,theta
      real    cofx(3),cofy(3)
      integer ns,nb      
c
c----------------------------------------------------------------------
c....input for the circle
      r = 0.5   !  radius of the circle 
      xref = 10.  !  x-coordinate of center at t=0
      zref = 15.  !  z-coordinate of center at t=0
      ns = 400   !  total number of points

c      r = 1.   !  radius of the circle 
c      xc = 3.5  !  x-coordinate of center
c      zc = 3.5  !  z-coordinate of center
c      ns = 120.   !  total number of points

      nb =ns-1  !  actual number of points (1 and ns coinside) 
      pi=3.14159265358979 


c....motion of the cylinder
      amp=2.*r/4.
      fr=0.1995
      xcr=0.0
c      zcr=amp*sin(2.*pi*fr*tlevel)
      zcr=amp*cos(2.*pi*fr*tlevel)
      xc=xcr+xref
      zc=zcr+zref

c----------------------------------------------------------------------
c....compute s,x(s),y(s)
      step = 2.*pi/float(ns-1)
      theta=-1.*step
      do i=1,ns
       theta= theta+step
       sb(i) = r*theta
       xb(i) = xc+r*sin(theta)
       yb(i) = zc+r*cos(theta)
c       write(123,*) i,xb(i),yb(i)
      enddo  
c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,nb
       call polcoe(sb(i-1),xb(i-1),3,cofx)
       call polcoe(sb(i-1),yb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c
c...body velocity
      ubd=0.0
      vbd=0.0
      wbd=2.*pi*fr*amp*sin(2.*pi*fr*tlevel)


      return
      end
C
C 
C-----SUBROUTINE-ObjectGrdHalfCyl-----------------E. Balaras..01/26/03-----
C
      subroutine ObjectGrdHalfCyl(sb,xb,yb,coefx,coefy,nb,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
C
      include 'headers/immbd.h'

      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
c
      real    r,xc,zc,pi
      real    step,theta
      real    cofx(3),cofy(3)
      integer ns,nb      
c
c----------------------------------------------------------------------
c....input for the circle
      r = 0.5   !  radius of the circle 
      xc = 10.  !  x-coordinate of center
      zc = 0.   !  z-coordinate of center
      ns = nimb   !  total number of points

      nb =ns  !  actual number of points 
      pi=3.14159265358979 
c----------------------------------------------------------------------
c....compute s,x(s),y(s)
      step = pi/float(ns-1)
      theta=-1.*step
      do i=1,ns
       theta= theta+step
       sb(i) = r*theta
       xb(i) = xc-r*cos(theta) 
       yb(i) = zc+r*sin(theta) 
      enddo  
c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,nb-1
       call polcoe(sb(i-1),xb(i-1),3,cofx)
       call polcoe(sb(i-1),yb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c
      return
      end
C
C
C 
C-----SUBROUTINE-ObjectGrdStn-------------------------E. Balaras..09/29/01-----
C
      subroutine ObjectGrdStn(sb,xb,yb,coefx,coefy,ntot,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
c
      include 'headers/common.h'
      include 'headers/immbd.h'

      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
      integer ntot
c..local arrays
      real   sbb(nmax),xbb(nmax),ybb(nmax)
      real   sbt(nmax),xbt(nmax),ybt(nmax)
      
      real    r,lx,Xref,Yref,dxref,dxx
      real    xc,zc,pi
      real    step,theta
      real    cofx(3),cofy(3)
      integer ns,nb      
c
c----------------------------------------------------------------------
c....geometry of stenosis (1)
c
c....................bottom wall......................................
      r = 0.5   !  pipe radius
      lx = 4.   !  lenght of stenosis
      Xref=13.0  !  starting point of stenosis from inflow (case4 & 5)
c      Xref=18.0  !  starting point of stenosis from inflow  (case3)
      Yref=0.0  !  starting point of stenosis from bottom wall
      ns = nimb ! VO 7/1/03   !  total number of points
      nb =ns-1  !  actual number of points (1 and ns coinside)
      pi=3.14159265358979

c...pass total number of points out
      ntot = 2*ns 

c....put marker points on curve
      dxref = lx/float(ns-1)
      dxx = (lx-2.*dxref)/float(ns-1)
      xbb(1) = dxx
      do i=2,ns
       xbb(i)=xbb(i-1)+dxx
      enddo
      do i=1,ns
       ybb(i)=0.5-0.25*(1.+cos(0.5*pi*xbb(i)))
      enddo
c..shift marker points to downstream location
      do i=1,ns
       xbb(i)=xbb(i)+Xref
       ybb(i)=ybb(i)+Yref
      enddo    
c..normalized arc-length
      sbb(1)=0.0
      do i=2,ns
       sbb(i)=sbb(i-1)+(1./(ns-1))
      enddo  
c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,ns-1
       call polcoe(sbb(i-1),xbb(i-1),3,cofx)
       call polcoe(sbb(i-1),ybb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c boundary points
       coefx(1,1) = coefx(2,1)
       coefx(1,2) = coefx(2,2)
       coefx(1,3) = coefx(2,3)
       coefy(1,1) = coefy(2,1)
       coefy(1,2) = coefy(2,2)
       coefy(1,3) = coefy(2,3)
c
       coefx(ns,1) = coefx(ns-1,1)
       coefx(ns,2) = coefx(ns-1,2)
       coefx(ns,3) = coefx(ns-1,3)
       coefy(ns,1) = coefy(ns-1,1)
       coefy(ns,2) = coefy(ns-1,2)
       coefy(ns,3) = coefy(ns-1,3)

c...store sb,xb,yb
      do i=1,ns 
       sb(i)=sbb(i)
       xb(i)=xbb(i)
       yb(i)=ybb(i)
      enddo

c....................top wall......................................
c...mirror the curve
      do i=1,ns
       ybb(i)=2.-ybb(i)
      enddo
c..mirror indexes
      ii=0
      do i=ns,1,-1
       ii=ii+1
       xbt(ii)=xbb(i)
       ybt(ii)=ybb(i)
      enddo
c..normalized arc-length
      sbt(1)=0.0
      do i=2,ns
       sbt(i)=sbt(i-1)+(1./(ns-1))
      enddo
c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,ns-1
       ii= ns+i
       call polcoe(sbt(i-1),xbt(i-1),3,cofx)
       call polcoe(sbt(i-1),ybt(i-1),3,cofy)
       coefx(ii,1) = cofx(3)
       coefx(ii,2) = cofx(2)
       coefx(ii,3) = cofx(1)
       coefy(ii,1) = cofy(3)
       coefy(ii,2) = cofy(2)
       coefy(ii,3) = cofy(1)
      enddo
c boundary points
       coefx(ns+1,1) = coefx(ns+2,1)
       coefx(ns+1,2) = coefx(ns+2,2)
       coefx(ns+1,3) = coefx(ns+2,3)
       coefy(ns+1,1) = coefy(ns+2,1)
       coefy(ns+1,2) = coefy(ns+2,2)
       coefy(ns+1,3) = coefy(ns+2,3)
c
       coefx(2*ns,1) = coefx(2*ns-1,1)
       coefx(2*ns,2) = coefx(2*ns-1,2)
       coefx(2*ns,3) = coefx(2*ns-1,3)
       coefy(2*ns,1) = coefy(2*ns-1,1)
       coefy(2*ns,2) = coefy(2*ns-1,2)
       coefy(2*ns,3) = coefy(2*ns-1,3)
c
c...store sb,xb,yb
      do i=ns+1,2*ns
       ii= i-ns
       sb(i)=sbt(ii)
       xb(i)=xbt(ii)
       yb(i)=ybt(ii)
      enddo
c
      return
      end

C
C 
C-----SUBROUTINE-ObjectGrdStnExp----------------------E. Balaras..05/29/03-----
C
      subroutine ObjectGrdStnExp(sb,xb,yb,coefx,coefy,ntot,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
c
      include 'headers/common.h'
      include 'headers/immbd.h'
c
      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
      integer ntot
c..local arrays
      real   sbb(nmax),xbb(nmax),ybb(nmax)
      real   sbt(nmax),xbt(nmax),ybt(nmax)
      
      real    r,lx,Xref,Yref,dxref,dxx
      real    xc,zc,pi
      real    step,theta
      real    cofx(3),cofy(3)
      integer ns,nb      

c...geometrical parameters
      r = 1.0     !  circle radius
      Xref=20.    !  reference point in the x-direction
      Yref=(r/2.) !  reference point in the y-direction
      ns = nimb/2 ! VO 7/1/03   !  number of points on each boundary
      ntot = 2*ns !  total number of points passed out
c
      pi=3.14159265358979 

c....................bottom wall......................................

c..normalized arc-length
      sbb(1)=0.0
      do i=2,ns
       sbb(i)=sbb(i-1)+(1./(ns-1))
      enddo

c....compute x(s),y(s)
      step = (2.*pi/3.)/float(ns-1)
      theta=-1.*step
      do i=1,ns
       theta= theta+step
       xbb(i) = xref-r*cos((pi/6.)+theta)
       ybb(i) = r*sin((pi/6.)+theta)-yref
      enddo

c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,ns-1
       call polcoe(sbb(i-1),xbb(i-1),3,cofx)
       call polcoe(sbb(i-1),ybb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c boundary points
       coefx(1,1) = coefx(2,1)
       coefx(1,2) = coefx(2,2)
       coefx(1,3) = coefx(2,3)
       coefy(1,1) = coefy(2,1)
       coefy(1,2) = coefy(2,2)
       coefy(1,3) = coefy(2,3)
c
       coefx(ns,1) = coefx(ns-1,1)
       coefx(ns,2) = coefx(ns-1,2)
       coefx(ns,3) = coefx(ns-1,3)
       coefy(ns,1) = coefy(ns-1,1)
       coefy(ns,2) = coefy(ns-1,2)
       coefy(ns,3) = coefy(ns-1,3)

c...store sb,xb,yb
      do i=1,ns 
       sb(i)=sbb(i)
       xb(i)=xbb(i)
       yb(i)=ybb(i)
      enddo

c....................top wall......................................

c..normalized arc-length
      sbt(1)=0.0
      do i=2,ns
       sbt(i)=sbt(i-1)+(1./(ns-1))
      enddo

c....compute x(s),y(s)
      do i=1,ns
       xbt(i) = xbb(ns-i+1)
       ybt(i) = -1.*ybb(ns-i+1)+2.*r
      enddo    

c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,ns-1
       ii= ns+i
       call polcoe(sbt(i-1),xbt(i-1),3,cofx)
       call polcoe(sbt(i-1),ybt(i-1),3,cofy)
       coefx(ii,1) = cofx(3)
       coefx(ii,2) = cofx(2)
       coefx(ii,3) = cofx(1)
       coefy(ii,1) = cofy(3)
       coefy(ii,2) = cofy(2)
       coefy(ii,3) = cofy(1)
      enddo
c boundary points
       coefx(ns+1,1) = coefx(ns+2,1)
       coefx(ns+1,2) = coefx(ns+2,2)
       coefx(ns+1,3) = coefx(ns+2,3)
       coefy(ns+1,1) = coefy(ns+2,1)
       coefy(ns+1,2) = coefy(ns+2,2)
       coefy(ns+1,3) = coefy(ns+2,3)
c
       coefx(2*ns,1) = coefx(2*ns-1,1)
       coefx(2*ns,2) = coefx(2*ns-1,2)
       coefx(2*ns,3) = coefx(2*ns-1,3)
       coefy(2*ns,1) = coefy(2*ns-1,1)
       coefy(2*ns,2) = coefy(2*ns-1,2)
       coefy(2*ns,3) = coefy(2*ns-1,3)
c
c...store sb,xb,yb
      do i=ns+1,2*ns
       ii= i-ns
       sb(i)=sbt(ii)
       xb(i)=xbt(ii)
       yb(i)=ybt(ii)
      enddo
c
      return
      end


C
C 
C-----SUBROUTINE-ObjectGrdWavy--------------------E. Balaras..02/15/02-----
C
      subroutine ObjectGrdWavy(sb,xb,yb,coefx,coefy,nb,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
c
      include 'headers/common.h'
      include 'headers/immbd.h'
c
      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)

c..local arrays
      real   sbb(nmax),xbb(nmax),ybb(nmax)
      real    r,l,a,lx,dxref,dxx
      real    xc,zc,pi
      real    step,theta
      real    cofx(3),cofy(3)
      integer ns,nb      
c
c----------------------------------------------------------------------
c....geometry a wavy wall located at the bottom of the channel
c
c      r = 0.1    !  2a/l
cc      r = 0.031   !  2a/l
c      l = 2.     !  wavelength
c      lx = 4.    !  domain size in x-direction
c      ns = 300   !  total number of points
c      nb =ns-1   !  actual number of points
c      pi=3.14159265358979
c      a=r*l/2.

c...Alberto's setup
      r = 0.1    !  2a/l
      l = 0.05     !  wavelength
      lx = 0.1    !  domain size in x-direction
      ns = nimb !VO 7/1/03   !  total number of points
      nb =ns-1   !  actual number of points
      pi=3.14159265358979
      a=r*l/2.

c....put marker points on curve
      dxref = lx/float(ns-1)
      dxx = (lx-2.*dxref)/float(ns-1)
      xbb(1) = dxx
      do i=2,ns
       xbb(i)=xbb(i-1)+dxx
      enddo
      do i=1,ns
c       ybb(i)=a*sin(2.*pi*xbb(i)/l - pi/2.) + a
       ybb(i)=a*(1.+cos(2*pi*xbb(i)/l))
      enddo  
c..normalized arc-length
      sbb(1)=0.0
      do i=2,ns
       sbb(i)=sbb(i-1)+(1./(ns-1))
      enddo  
c----------------------------------------------------------------------
c....generate functions: 
c       x(s)_i = coeffx(i,1)*s^2+coeffx(i,2)*s+coeffx(i,3)
c       y(s)_i = coeffy(i,1)*s^2+coeffy(i,2)*s+coeffy(i,3)
c....coefficients of the quadratic polinomials:
c       coefx(i,1)=ax, coefx(i,2)=bx,coefx(i,3)=cx
c       coefy(i,1)=ay, coefy(i,2)=by,coefy(i,3)=cy
      do i=2,ns-1
       call polcoe(sbb(i-1),xbb(i-1),3,cofx)
       call polcoe(sbb(i-1),ybb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c boundary points
       coefx(1,1) = coefx(2,1)
       coefx(1,2) = coefx(2,2)
       coefx(1,3) = coefx(2,3)
       coefy(1,1) = coefy(2,1)
       coefy(1,2) = coefy(2,2)
       coefy(1,3) = coefy(2,3)
c
       coefx(ns,1) = coefx(ns-1,1)
       coefx(ns,2) = coefx(ns-1,2)
       coefx(ns,3) = coefx(ns-1,3)
       coefy(ns,1) = coefy(ns-1,1)
       coefy(ns,2) = coefy(ns-1,2)
       coefy(ns,3) = coefy(ns-1,3)

c...store sb,xb,yb
      do i=1,ns 
       sb(i)=sbb(i)
       xb(i)=xbb(i)
       yb(i)=ybb(i)
      enddo
c
      return
      end
C
C 
C-----SUBROUTINE-bpoints-------------------------E. Balaras..03/04/01-----
C
      subroutine bpoints(xu,zu,sb,xb,yb,nb,coefx,coefy,
     &                   iim,jim,i2b,k2b,c4int,nim,
     &                   ibody,kbody,nbody,
     &                   xnim,ynim,snim,
     &                   nx,nz,nmax                    )
C
C     locates i,j,k points on the cartecian grid near
C     the immersed boundary to be modified     
C
      include 'headers/common.h'
      include 'mpif.h'
c
      real    xu(nx),zu(nz)
      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
      real    c4int(nmax,5)
      real    xnim(nmax),ynim(nmax),snim(nmax)
      integer iim(nmax),jim(nmax)
      integer i2b(nmax),k2b(nmax)
      integer ibody(nmax),kbody(nmax)
      integer nb,nim,nbody
c....local variables
      real    xx,zz,xn,yn,sa,sn
      real    a1,a2,a3,a4,tt,uu,d,h1,h2,h1h2
      integer i,k,ii,kk,n,inout
      real    wrk1(nx,nz),wrk2(nx,nz),wrk3(nx,nz)
      real    x2im(nmax),y2im(nmax)
      real    cofx(3),cofy(3),xpas(5),zpas(5)
      integer iwrk(nx,nz)
      integer ib(nmax),kb(nmax)
      real    deltax,deltaz,tlr,deltaxg,deltazg
      real    xcrd1,xcrd2,xcrd3,zcrd1,zcrd2,zcrd3
      real    xchk,zchk,h
      integer itest
      real c4sum
      integer ierr,type
c
#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
c

c-----------------------------------------------------------------------
c....locate the velocity cells that contain markers
c    indexes are stored in arrays ib and kb
c    marker is between ib(i) and ib(i)+1 and
c                      kb(i) and kb(i)+1
      do i=1,nb
       xx=xb(i)
       call locate(xu,nx,xx,ii)
       zz=yb(i)
       call locate(zu,nz,zz,kk)
       ib(i)=ii
       kb(i)=kk
      enddo

c-----------------------------------------------------------------------
c....compute average grid size near the immersed boundary.  Use this
c    info to calculate the tolerance that determines if a grid node
c    is located exactly on the IB
       deltax=0.0
       deltaz=0.0
       do i=1,nb
       if (ib(i).ne.0.and.ib(i).lt.nx-2.and.kb(i).ne.0.and.kb(i).lt.nz-1) then 
        deltax=deltax+(xu(ib(i)+1)-xu(ib(i)))
        deltaz=deltaz+(zu(kb(i)+1)-zu(kb(i)))
       endif	
       enddo
c      make tolerance the same in a  cpus & independent of number of cpus       
       call MPI_ALLREDUCE(deltax,deltaxg, 1, type, MPI_SUM, MPI_COMM_HYDRO,ierr)
       call MPI_ALLREDUCE(deltaz,deltazg, 1, type, MPI_SUM, MPI_COMM_HYDRO,ierr)
c
       deltaxg=deltaxg/float(nb)
       deltazg=deltazg/float(nb)
       tlr=sqrt(deltaxg**2.+deltazg**2.)/15.


c-----------------------------------------------------------------------
c...compute where normal intersects boundary for the cluster points
c   surrounding the marker (This version is for a 2x2 cluster)

c...initialize work array
      do i=1,ix2+1
       do k=1,kz2+1
        iwrk(i,k)=0
       enddo
      enddo

c...loop over marker points

      DO i=2,nb 

       if (ib(i).ne.0.and.ib(i).lt.nx-1.and.kb(i).ne.0.and.kb(i).lt.nz) then 
        do ii=1,2    
         do kk=1,2
          xx=xu(ib(i)+ii-1)
          zz=zu(kb(i)+kk-1)
          sa=sb(i) 
          do n=1,3
           cofx(n)=coefx(i,n)
           cofy(n)=coefy(i,n)
          enddo

          CALL Normal(cofx,cofy,xx,zz,sa,sn)
          xn=coefx(i,1)*sn**2.+coefx(i,2)*sn+coefx(i,3)
          yn=coefy(i,1)*sn**2.+coefy(i,2)*sn+coefy(i,3)

          CALL LocPnt(cofx,cofy,xx,zz,xn,yn,sn,tlr,inout)
          iwrk(ib(i)+ii-1,kb(i)+kk-1) = inout
          wrk1(ib(i)+ii-1,kb(i)+kk-1) = xn
          wrk2(ib(i)+ii-1,kb(i)+kk-1) = yn
          wrk3(ib(i)+ii-1,kb(i)+kk-1) = sn 
         enddo
        enddo
	else  
c	 write(0,*) 'i. b. ', ib(i),kb(i)
       endif
c
      ENDDO


c...mark all cells inside the body with +1 flag
c     commented out; VO 7/24/03
c      do k=2,kz2
c       do i=2,ix2
c        itest=iwrk(i,k)+iwrk(i+1,k)
c        if(itest.eq.1.or.itest.eq.2) then
c         iwrk(i+1,k) = 1
c        endif
c       enddo
c      enddo

c...find true boundary cells
      ii=0
      do k=2,kz2
       do i=2,ix2
c
         if(iwrk(i,k).ne.0) then !VO 7/24
c        if(iwrk(i,k).eq.-1) then
         a1=iwrk(i,k)+iwrk(i+1,k)
         a2=iwrk(i,k)+iwrk(i-1,k)
         a3=iwrk(i,k)+iwrk(i,k+1)
         a4=iwrk(i,k)+iwrk(i,k-1)
          if(a1.eq.0.or.a2.eq.0.or.a3.eq.0.or.a4.eq.0) then
           ii=ii+1
           iim (ii)=i
           jim (ii)=k
           xnim(ii)=wrk1(i,k)
           ynim(ii)=wrk2(i,k)
           snim(ii)=wrk3(i,k)
          endif         
        endif
c
       enddo
      enddo
      nim=ii

c...find interior cells and strore in 1D array
c  commented out; VO 7/24/03
      nbody=0
c      ii=0
c      do k=2,kz2
c       do i=2,ix2
c
c        if(iwrk(i,k).eq.1) then
c           ii=ii+1
c           ibody(ii)=i
c           kbody(ii)=k
c        endif
c
c       enddo
c      enddo
c      nbody=ii

c...arange all final flags into iwrk(i,k) array
      do i=2,ix2
       do k=2,kz2
         iwrk(i,k)=0
       enddo
      enddo
      do i=1,nim
         iwrk(iim(i),jim(i))=-1
      enddo
      do i=1,nbody
         iwrk(ibody(i),kbody(i))=1
      enddo


c....find 2nd point on the normal to the boundary used to linearize velocity

c...loop over true boundary points
      DO i=1,nim
      
c...find if the point is on the boundary.  In this case no interpolation
c   is required  (h is the distance of the point from the boundary)
      h=sqrt((xnim(i)-xu(iim(i)))**2.+(ynim(i)-zu(jim(i)))**2.)
      if(h.le.tlr) then
         i2b(i)=iim(i)
         k2b(i)=jim(i)
         x2im(i) = xu(iim(i))
         y2im(i) = zu(jim(i))
         c4int(i,1)=1.0
         c4int(i,2)=0.0
         c4int(i,3)=0.0
         c4int(i,4)=0.0
         c4int(i,5)=0.0
c         write(6,*) 'points on psi'

      else
c...find the point on the normal to the boundary at distance 2h, using the
c   midpoint formula
c                         h1       h2
c                    o--------o---------o
c                 1=bound  2=grid   3=virtual
c
       xcrd1=xnim(i)
       zcrd1=ynim(i)
       xcrd2=xu(iim(i))
       zcrd2=zu(jim(i))
       xcrd3=2.*xcrd2-xcrd1
       zcrd3=2.*zcrd2-zcrd1
c...find points around the virtual point 3 and check if other boundary 
c   points are involved in the interpolation stencil
 1000  continue
       call locate(xu,nx,xcrd3,ii)
       call locate(zu,nz,zcrd3,kk)
       itest=iwrk(ii,kk  )+iwrk(ii+1,kk  )+
     &       iwrk(ii,kk+1)+iwrk(ii+1,kk+1)
c....case corresponding to Figure 5a or 5c in the paper
c
       write(200+myrank,'(5 I5)') nb,ii,kk,nx,nz
c
         if(itest.eq.0) then
    
          i2b(i)=ii
          k2b(i)=kk
          x2im(i) = xcrd3
          y2im(i) = zcrd3

          h1=sqrt((xnim(i)-xu(iim(i)))**2.+(ynim(i)-zu(jim(i)))**2.)
          h2=sqrt((xu(iim(i))-x2im(i))**2.+(zu(jim(i))-y2im(i))**2.)
          h1h2=h1+h2

          tt = (x2im(i)-xu(i2b(i)))/(xu(i2b(i)+1)-xu(i2b(i)))
          uu = (y2im(i)-zu(k2b(i)))/(zu(k2b(i)+1)-zu(k2b(i)))
          a1=(1-tt)*(1-uu)
          a2=   tt *(1-uu)
          a3=   tt *   uu
          a4=(1-tt)*   uu

          c4int(i,1)=h2/h1h2
          c4int(i,2)=(h1/h1h2)*a1
          c4int(i,3)=(h1/h1h2)*a2
          c4int(i,4)=(h1/h1h2)*a3
          c4int(i,5)=(h1/h1h2)*a4



c....case corresponding to Figure 5b in the paper
         elseif(itest.eq.-1.and.iwrk(iim(i),jim(i)).eq.-1) then
          i2b(i)=ii
          k2b(i)=kk
          x2im(i) = xcrd3
          y2im(i) = zcrd3

          h1=sqrt((xnim(i)-xu(iim(i)))**2.+(ynim(i)-zu(jim(i)))**2.)
          h2=sqrt((xu(iim(i))-x2im(i))**2.+(zu(jim(i))-y2im(i))**2.)
          h1h2=h1+h2


          tt = (x2im(i)-xu(i2b(i)))/(xu(i2b(i)+1)-xu(i2b(i)))
          uu = (y2im(i)-zu(k2b(i)))/(zu(k2b(i)+1)-zu(k2b(i)))
          a1=(1-tt)*(1-uu)
          a2=   tt *(1-uu)
          a3=   tt *   uu
          a4=(1-tt)*   uu
          d=1.+(h1/h1h2)*( a1*float(iwrk(ii  ,kk  )) +
     &                     a2*float(iwrk(ii+1,kk  )) +
     &                     a3*float(iwrk(ii+1,kk+1)) +
     &                     a4*float(iwrk(ii  ,kk+1))   )

          c4int(i,1)=h2/h1h2/d
          c4int(i,2)=((h1/h1h2/d)*a1)*float(iwrk(ii  ,kk  )+1)
          c4int(i,3)=((h1/h1h2/d)*a2)*float(iwrk(ii+1,kk  )+1)
          c4int(i,4)=((h1/h1h2/d)*a3)*float(iwrk(ii+1,kk+1)+1)
          c4int(i,5)=((h1/h1h2/d)*a4)*float(iwrk(ii  ,kk+1)+1)
c	    VO 7.1.03
         else
c	 write(0,*) 'gimme a break, man'
          xcrd1=xcrd2
          zcrd1=zcrd2
          xcrd2=xcrd3
          zcrd2=zcrd3
          xcrd3=2.*xcrd2-xcrd1
          zcrd3=2.*zcrd2-zcrd1
c
	  if (xcrd3.gt.xu(ix2).or.zcrd3.gt.zu(kz2+1).or.xcrd3.lt.xu(ix1).or.zcrd3.lt.zu(kz1)) then 
c           cannot go farther out
             if (itest.eq.-3) then 
	         stop "bpoints: exception 1: ib point too close to logical boundary "
	     else if (itest.eq.-2) then 
	       write(0,*) "bpoints: warning : incomplete stencil: ib point too close to logical boundary "
c      proceed anyway -- interpolate based on distance using three (3) points

               c4sum=0. 
               i2b(i)=ii
               k2b(i)=kk
               x2im(i) = xcrd2
               y2im(i) = zcrd2
c  
c               tt = (x2im(i)-xu(i2b(i)))/(xu(i2b(i)+1)-xu(i2b(i)))
c               uu = (y2im(i)-zu(k2b(i)))/(zu(k2b(i)+1)-zu(k2b(i)))
c               a1=(1-tt)*(1-uu)
c               a2=   tt *(1-uu)
c               a3=   tt *   uu
c               a4=(1-tt)*   uu

               c4int(i,1)=1./(sqrt((xnim(i)-xu(iim(i)))**2+(ynim(i)-zu(jim(i)))**2))
               c4int(i,2)=sqrt((xu(i2b(i))-xu(iim(i)))**2+(zu(k2b(i))-zu(jim(i)))**2)
	        if (iwrk(ii,kk).eq.-1) c4int(i,2)=0.
               c4int(i,3)=sqrt((xu(i2b(i)+1)-xu(iim(i)))**2+(zu(k2b(i))-zu(jim(i)))**2)
	        if (iwrk(ii+1,kk).eq.-1) c4int(i,3)=0.
               c4int(i,4)=sqrt((xu(i2b(i)+1)-xu(iim(i)))**2+(zu(k2b(i)+1)-zu(jim(i)))**2)
	        if (iwrk(ii+1,kk+1).eq.-1) c4int(i,4)=0.
               c4int(i,5)=sqrt((xu(i2b(i))-xu(iim(i)))**2+(zu(k2b(i)+1)-zu(jim(i)))**2)
	        if (iwrk(ii,kk+1).eq.-1) c4int(i,5)=0.
               do nn=2,5
	         if (c4int(i,nn).ge.1.E-4) then 
	            c4int(i,nn)=1./c4int(i,nn)
	         else 
	            c4int(i,nn)=0.0
	         endif
	         c4sum=c4sum+c4int(i,nn)
	       enddo 	 
c
               do nn=1,5
	         c4int(i,nn)=c4int(i,nn)/c4sum
	       enddo 	 

             else 
	       write(0,*) 'bpoints: exception 2: itest = ',itest
	       stop
	     endif !(itest = -3)	  
	  else 
c           write(0,*) "bpoints: warning : incomplete stencil ",itest
c	   write(0,*) zu(kz1),zcrd3,zu(kz2+1)
           goto 1000
	  endif !xcrd            
         endif ! (itest = 0)
c
      endif !(h<tlr)
c       
 7654 continue
      ENDDO

c...sum of interpolation coef.
c       do i=1,nim
c        dd= c4int(i,1)+c4int(i,2)+c4int(i,3)+c4int(i,4)+c4int(i,5)
c        write(6,*) i,dd
c       enddo

c...grid
c      do i=2,ix2
c       do k=2,kz2
c        write(110,*) xu(i),zu(k)
c       enddo
c      enddo

c...all exterior control points
c      do i=1,nim
c       write(130,111) i,xu(iim(i)),zu(jim(i)),xnim(i),ynim(i),
c     &              x2im(i),y2im(i)
c      enddo

c...all interior control points
c      do i=1,nbody
c       write(140,112)i,xu(ibody(i)),zu(kbody(i))
c      enddo

c 111  format(i4,6(1x,e14.7))
c 112  format(i4,2(1x,e14.7))
      close(200+myrank)
      return
      end
C
C 
C-----SUBROUTINE-Normal-------------------------E. Balaras..03/04/01-----
C
      SUBROUTINE Normal(coefx,coefy,x,y,sa,sn)
c
c     computes the arclength coordinate at the interception
c     with the normal on (i,k)
c
      real coefx(3),coefy(3)
      real x,y,sn

c      double precision ax,bx,cx,ay,by,cy,f,fpr,sn1,snew,sntest,xx,yy
c
      ax=coefx(1)
      bx=coefx(2)
      cx=coefx(3)
      ay=coefy(1)
      by=coefy(2)
      cy=coefy(3)
c
      sn=sa
c
 50   continue
      f=(2.*ay**2.+2.*ax**2.)*sn**3. +
     &  (3.*ax*bx+3.*ay*by  )*sn**2. +
     &  (2.*ay*cy+2.*ax*cx+bx*bx+by*by-2.*ay*y-2.*ax*x)*sn+
     &  (cx*bx+by*cy-x*bx-y*by)
c
      fpr=3.*(2.*ay**2.+2.*ax**2.)*sn**2. +
     &    2.*(3.*ax*bx+3.*ay*by  )*sn     +
     &    2.*ay*cy+2.*ax*cx+bx*bx+by*by-2.*ay*y-2.*ax*x
c
      snew=sn-(f/fpr)
      sntest=abs(snew-sn)/snew
      sn=snew
c      write(123,*) sntest
      if(abs(sntest).lt.1.0e-06) then 
       goto 55
      else
       goto 50
      endif

 55   continue

      return
      end
C
C 
C-----SUBROUTINE-LocPnt-------------------------E. Balaras..03/04/01-----
C
      SUBROUTINE LocPnt(coefx,coefy,xi,yj,xn,yn,sn,tlr,inout)
c
c     computes the arclength coordinate at the interception
c     with the normal on (i,k)
c
      real    coefx(3),coefy(3)
      real    xi,yj,sn
      integer inout
      real    ax,bx,cx,ay,by,cy
      real    dxds,dyds
      real    wrk,vi,vj,ni,nj,delta
      real    tlr,h
c
      ax=coefx(1)
      bx=coefx(2)
      cx=coefx(3)
      ay=coefy(1)
      by=coefy(2)
      cy=coefy(3)
c
      wrk= sqrt((xn-xi)**2.+(yn-yj)**2.)
      vi = (xn-xi)/wrk
      vj = (yn-yj)/wrk
c
      dxds=2.*ax*sn+bx
      dyds=2.*ay*sn+by
      wrk= sqrt(dxds*dxds+dyds*dyds)
      ni = -1.*dyds/wrk
      nj = dxds/wrk
c
      delta = vi*ni+vj*nj
c      write(6,*) 'delta=',delta,vi,vj,ni,nj

c...mark point with a flag.  On the boundary (based on a specified tolerance)
c   or in the fluid inout=1.  In the solid inout=+1
c
      h=sqrt((xn-xi)**2.+(yn-yj)**2.)
      if(h.le.tlr) then
c        write(6,*) 'point on boundary'
        inout=-1
      elseif(delta.lt.0.0) then
        inout=-1
      else
        inout=1
      endif

      return
      end

C
C 
C-----SUBROUTINE-Forcing-------------------------E. Balaras..03/04/01-----
C
      SUBROUTINE Forcing(us,force,uim,c4int,i2u,k2u,iim,kim,
     &                   ibody,kbody,nbody,
     &                   dt,nim,im,jm,km,nmax           )
C
C          modify u*
C          because of the explicit nature of the algorithm, the 
C          application of forcing (as explained in the C&F paper) 
C          is equivalent in directly modifying (imposing bc) u*
C
      include 'headers/common.h'

      real    us(im,jm,km),force(im,jm,km)
      real    c4int(nmax,5)
      real    uim(nmax)
      integer iim(nmax),kim(nmax),i2u(nmax),k2u(nmax)
      integer ibody(nmax),kbody(nmax)
      real    dt
      real    u2,utarget
      integer nim,nbody



c...compute force on immerse boundary points
      DO j=jy1,jy2+1
       do i=1,nim
        us(iim(i),j,kim(i))= c4int(i,1)*uim(i) +
     &                       c4int(i,2)*us(i2u(i)  ,j,k2u(i)  ) +
     &                       c4int(i,3)*us(i2u(i)+1,j,k2u(i)  ) +
     &                       c4int(i,4)*us(i2u(i)+1,j,k2u(i)+1) +
     &                       c4int(i,5)*us(i2u(i)  ,j,k2u(i)+1)
       enddo
      ENDDO

c...set bc on interior points
c      DO j=jy1+1,jy2
c       do i=1,nbody
c         us(ibody(i),j,kbody(i)) = 0.0
c       enddo
c      ENDDO



c
      RETURN
      END

C
C 
C-----SUBROUTINE-RWim-------------------------E. Balaras..03/04/01-----
C
      SUBROUTINE RWim(nimu,nbodyu,iu,ku,i2u,k2u,ibodyu,kbodyu,
     &                c4uint,xnu,znu,snu,nmax,filestr,irw     )
C
C          read/write immersed boundary information
C

      integer nimu,nbodyu,irw

      integer iu(nmax),ku(nmax),i2u(nmax),k2u(nmax)
      integer ibodyu(nmax),kbodyu(nmax)

      real    c4uint(nmax,5)
      real    xnu(nmax),znu(nmax),snu(nmax)

      character*80 filestr

      OPEN(81,file=filestr,status='unknown',form='formatted')

       if(irw.eq.1) then
c..write to file
       write(81,1010) nimu,nbodyu
       write(81,1011) (iu (i),i=1,nimu)
       write(81,1011) (ku (i),i=1,nimu)
       write(81,1011) (i2u(i),i=1,nimu)
       write(81,1011) (k2u(i),i=1,nimu)
c
       write(81,1011) (ibodyu (i),i=1,nbodyu)
       write(81,1011) (kbodyu(i),i=1,nbodyu)
c
       write(81,1012) (c4uint(i,1),i=1,nimu)
       write(81,1012) (c4uint(i,2),i=1,nimu)
       write(81,1012) (c4uint(i,3),i=1,nimu)
       write(81,1012) (c4uint(i,4),i=1,nimu)
c
       write(81,1012) (xnu(i),i=1,nimu)
       write(81,1012) (znu(i),i=1,nimu)
       write(81,1012) (snu(i),i=1,nimu)
c
      else
c...read info from file
       read(81,1010) nimu,nbodyu
       read(81,1011) (iu (i),i=1,nimu)
       read(81,1011) (ku (i),i=1,nimu)
       read(81,1011) (i2u(i),i=1,nimu)
       read(81,1011) (k2u(i),i=1,nimu)
c
       read(81,1011) (ibodyu(i),i=1,nbodyu)
       read(81,1011) (kbodyu(i),i=1,nbodyu)
c
       read(81,1012) (c4uint(i,1),i=1,nimu)
       read(81,1012) (c4uint(i,2),i=1,nimu)
       read(81,1012) (c4uint(i,3),i=1,nimu)
       read(81,1012) (c4uint(i,4),i=1,nimu)
c
       read(81,1012) (xnu(i),i=1,nimu)
       read(81,1012) (znu(i),i=1,nimu)
       read(81,1012) (snu(i),i=1,nimu)

      endif
c
 1010 format(i4,2x,i4)
 1011 format(5(2x,i4))
 1012 format(5(1x,E14.7))

      close(81)

c
      RETURN
      END


C
C 
C-----SUBROUTINE-bpointspr-------------------------E. Balaras..04/03/03-----
C
      subroutine bpointspr(xu,zu,sb,xb,yb,nb,coefx,coefy,
     &                     iim,jim,i2b,k2b,c4int,nim,
     &                     ibody,kbody,nbody,
     &                     xnim,ynim,snim,
     &                     nx,nz,nmax                    )
C
C     locates i,j,k points on the cartecian grid near
C     the immersed boundary to be modified     
C
      include 'headers/common.h'
c
      real    xu(nx),zu(nz)
      real    sb(nmax),xb(nmax),yb(nmax)
      real    coefx(nmax,3),coefy(nmax,3)
      real    c4int(nmax,5)
      real    xnim(nmax),ynim(nmax),snim(nmax)
      integer iim(nmax),jim(nmax)
      integer i2b(nmax),k2b(nmax)
      integer ibody(nmax),kbody(nmax)
      integer nb,nim,nbody
c....local variables
      real    xx,zz,xn,yn,sa,sn
      real    a1,a2,a3,a4,tt,uu
      integer i,k,ii,kk,n,inout
      real    wrk1(nx,nz),wrk2(nx,nz),wrk3(nx,nz)
      real    x2im(nmax),y2im(nmax)
      real    cofx(3),cofy(3),xpas(5),zpas(5)
      integer iwrk(nx,nz)
      integer ib(nmax),kb(nmax)
      real    deltax,deltaz,tlr
      real    xcrd1,xcrd2,xcrd3,zcrd1,zcrd2,zcrd3
      real    xchk,zchk,h
      integer itest



c-----------------------------------------------------------------------
c....locate the velocity cells that contain markers
c    indexes are stored in arrays ib and kb
c    marker is between ib(i) and ib(i)+1 and
c                      kb(i) and kb(i)+1
      do i=1,nb
       xx=xb(i)
       call locate(xu,nx,xx,ii)
       zz=yb(i)
       call locate(zu,nz,zz,kk)
       ib(i)=ii
       kb(i)=kk
      enddo

c-----------------------------------------------------------------------
c....compute average grid size near the immersed boundary.  Use this
c    info to calculate the tolerance that determines if a grid node
c    is located exactly on the IB
       deltax=0.0
       deltaz=0.0
       do i=1,nb
        deltax=deltax+(xu(ib(i)+1)-xu(ib(i)))
        deltaz=deltaz+(zu(kb(i)+1)-zu(kb(i)))
       enddo
       deltax=deltax/float(nb)
       deltaz=deltaz/float(nb)
       tlr=sqrt(deltax**2.+deltaz**2.)/15.
c       tlr=sqrt(deltax**2.+deltaz**2.)/100.


c-----------------------------------------------------------------------
c...compute where normal intersects boundary for the cluster points
c   surrounding the marker (This version is for a 2x2 cluster)

c...initialize work array
      do i=1,ix2+1
       do k=1,kz2+1
        iwrk(i,k)=0
       enddo
      enddo

c...loop over marker points
      DO i=2,nb 

       do ii=1,2    
        do kk=1,2
         xx=xu(ib(i)+ii-1)
         zz=zu(kb(i)+kk-1)
         sa=sb(i) 
         do n=1,3
          cofx(n)=coefx(i,n)
          cofy(n)=coefy(i,n)
         enddo

         CALL Normal(cofx,cofy,xx,zz,sa,sn)
         xn=coefx(i,1)*sn**2.+coefx(i,2)*sn+coefx(i,3)
         yn=coefy(i,1)*sn**2.+coefy(i,2)*sn+coefy(i,3)

         CALL LocPnt(cofx,cofy,xx,zz,xn,yn,sn,tlr,inout)
         iwrk(ib(i)+ii-1,kb(i)+kk-1) = inout
         wrk1(ib(i)+ii-1,kb(i)+kk-1) = xn
         wrk2(ib(i)+ii-1,kb(i)+kk-1) = yn
         wrk3(ib(i)+ii-1,kb(i)+kk-1) = sn 
        enddo
       enddo
c
      ENDDO


c...mark all cells inside the body with +1 flag
      do k=2,kz2
       do i=2,ix2
        itest=iwrk(i,k)+iwrk(i+1,k)
        if(itest.eq.1.or.itest.eq.2) then
         iwrk(i+1,k) = 1
        endif
       enddo
      enddo

c...find true boundary cells
      ii=0
      do k=2,kz2
       do i=2,ix2
c
        if(iwrk(i,k).eq.-1) then
         a1=iwrk(i,k)+iwrk(i+1,k)
         a2=iwrk(i,k)+iwrk(i-1,k)
         a3=iwrk(i,k)+iwrk(i,k+1)
         a4=iwrk(i,k)+iwrk(i,k-1)
          if(a1.eq.0.or.a2.eq.0.or.a3.eq.0.or.a4.eq.0) then
           ii=ii+1
           iim (ii)=i
           jim (ii)=k
           xnim(ii)=wrk1(i,k)
           ynim(ii)=wrk2(i,k)
           snim(ii)=wrk3(i,k)
          endif         
        endif
c
       enddo
      enddo
      nim=ii

c...find interior cells and strore in 1D array
      ii=0
      do k=2,kz2
       do i=2,ix2
c
        if(iwrk(i,k).eq.1) then
           ii=ii+1
           ibody(ii)=i
           kbody(ii)=k
        endif
c
       enddo
      enddo
      nbody=ii

c...arange all final flags into iwrk(i,k) array
      do i=2,ix2
       do k=2,kz2
         iwrk(i,k)=0
       enddo
      enddo
      do i=1,nim
         iwrk(iim(i),jim(i))=-1
      enddo
      do i=1,nbody
         iwrk(ibody(i),kbody(i))=1
      enddo


c....find 2nd point on the normal to the boundary used to linearize velocity

c...loop over true boundary points
      DO i=1,nim
c...find if the point is on the boundary.  In this case no interpolation
c   is required  (h is the distance of the point from the boundary)
c      h=sqrt((xnim(i)-xu(iim(i)))**2.+(ynim(i)-zu(jim(i)))**2.)
c      if(h.le.tlr) then
c         i2b(i)=iim(i)
c         k2b(i)=jim(i)
c         x2im(i) = xu(iim(i))
c         y2im(i) = zu(jim(i))
c         c4int(i,1)=1.0
c         c4int(i,2)=0.0
c         c4int(i,3)=0.0
c         c4int(i,4)=0.0
c         c4int(i,5)=0.0
c         write(6,*) 'points on psi'

c      else
c...find the point on the normal to the boundary at distance 2h, using the
c   midpoint formula
c                         h1       h2
c                    o--------o---------o
c                 1=bound  2=grid   3=virtual
c
       xcrd1=xnim(i)
       zcrd1=ynim(i)
       xcrd2=xu(iim(i))
       zcrd2=zu(jim(i))
       xcrd3=2.*xcrd2-xcrd1
       zcrd3=2.*zcrd2-zcrd1
c...find points around the virtual point 3 and check if other boundary 
c   points are involved in the interpolation stencil
 1000  continue
       call locate(xu,nx,xcrd3,ii)
       call locate(zu,nz,zcrd3,kk)
       itest=iwrk(ii,kk  )+iwrk(ii+1,kk  )+
     &       iwrk(ii,kk+1)+iwrk(ii+1,kk+1)
c....case corresponding to Figure 5a or 5c in the paper
         if(itest.eq.0) then

          i2b(i)=ii
          k2b(i)=kk
          x2im(i) = xcrd3
          y2im(i) = zcrd3

          h1=sqrt((xnim(i)-xu(iim(i)))**2.+(ynim(i)-zu(jim(i)))**2.)
          h2=sqrt((xu(iim(i))-x2im(i))**2.+(zu(jim(i))-y2im(i))**2.)
          h1h2=h1+h2

          tt = (x2im(i)-xu(i2b(i)))/(xu(i2b(i)+1)-xu(i2b(i)))
          uu = (y2im(i)-zu(k2b(i)))/(zu(k2b(i)+1)-zu(k2b(i)))
          a1=(1-tt)*(1-uu)
          a2=   tt *(1-uu)
          a3=   tt *   uu
          a4=(1-tt)*   uu

          c4int(i,1)=0.0
          c4int(i,2)=a1
          c4int(i,3)=a2
          c4int(i,4)=a3
          c4int(i,5)=a4

c          write(6,*) 'points on with 5 point stencil',h1,h2

c....case corresponding to Figure 5b in the paper
c         elseif(itest.eq.-1.and.iwrk(iim(i),jim(i)).eq.-1) then

c          i2b(i)=ii
c          k2b(i)=kk
c          x2im(i) = xcrd3
c          y2im(i) = zcrd3

c          h1=sqrt((xnim(i)-xu(iim(i)))**2.+(ynim(i)-zu(jim(i)))**2.)
c          h2=sqrt((xu(iim(i))-x2im(i))**2.+(zu(jim(i))-y2im(i))**2.)
c          h1h2=h1+h2

c          tt = (x2im(i)-xu(i2b(i)))/(xu(i2b(i)+1)-xu(i2b(i)))
c          uu = (y2im(i)-zu(k2b(i)))/(zu(k2b(i)+1)-zu(k2b(i)))
c          a1=(1-tt)*(1-uu)
c          a2=   tt *(1-uu)
c          a3=   tt *   uu
c          a4=(1-tt)*   uu
c          d=1.+(h1/h1h2)*( a1*float(iwrk(ii  ,kk  )) +
c     &                     a2*float(iwrk(ii+1,kk  )) +
c     &                     a3*float(iwrk(ii+1,kk+1)) +
c     &                     a4*float(iwrk(ii  ,kk+1))   )

c          c4int(i,1)=h2/h1h2/d
c          c4int(i,2)=((h1/h1h2/d)*a1)*float(iwrk(ii  ,kk  )+1)
c          c4int(i,3)=((h1/h1h2/d)*a2)*float(iwrk(ii+1,kk  )+1)
c          c4int(i,4)=((h1/h1h2/d)*a3)*float(iwrk(ii+1,kk+1)+1)
c          c4int(i,5)=((h1/h1h2/d)*a4)*float(iwrk(ii  ,kk+1)+1)
c         write(6,*) 'points on 4 point stencil',h1,h2

         else
          xcrd1=xcrd2
          zcrd1=zcrd2
          xcrd2=xcrd3
          zcrd2=zcrd3
          xcrd3=2.*xcrd2-xcrd1
          zcrd3=2.*zcrd2-zcrd1 
          goto 1000         
         endif
c
c      endif
c       
      ENDDO

c...sum of interpolation coef.
c       do i=1,nim
c        dd= c4int(i,1)+c4int(i,2)+c4int(i,3)+c4int(i,4)+c4int(i,5)
c        write(6,*) i,dd
c       enddo

c...grid
c      do i=2,ix2
c       do k=2,kz2
c        write(110,*) xu(i),zu(k)
c       enddo
c      enddo

c...all exterior control points
c      do i=1,nim
c       write(130,111) i,xu(iim(i)),zu(jim(i)),xnim(i),ynim(i),
c     &              x2im(i),y2im(i)
c      enddo

c...all interior control points
c      do i=1,nbody
c       write(140,112)i,xu(ibody(i)),zu(kbody(i))
c      enddo

 111  format(i4,6(1x,e14.7))
 112  format(i4,2(1x,e14.7))

      return
      end
C
