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
      r = 1.   !  radius of the circle 
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
      real    sb(0:nmax),xb(0:nmax),yb(0:nmax)
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
      xc = 19.5  !  x-coordinate of center
      zc = 10.5  !  z-coordinate of center
      ns = nimb !  total number of points

      nb =ns-1  !  actual number of points (1 and ns coinside) 
      pi=3.14159265358979 
c----------------------------------------------------------------------
c....compute s,x(s),y(s)
      step = 2.*pi/float(ns-1)
      theta=-1.*step
      do i=0,ns
       theta= theta+step
       sb(i) = r*theta
       xb(i) = xc+r*cos(theta)
       yb(i) = zc+r*sin(theta)
c       write(123,*) i,xb(i),yb(i)
      enddo  
c----------------------------------------------------------------------
      do i=1,nb
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
C-----SUBROUTINE-ObjectBluntNose-------------------------V. Ovchinnikov---06/14/04-----
C
      subroutine ObjectBluntNose(sb,xb,yb,coefx,coefy,nb,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
C
      include 'headers/immbd.h'
      real    sb(0:nnmax),xb(0:nnmax),yb(0:nnmax)
      real    coefx(nnmax,3),coefy(nnmax,3)
c
      real    r,xc,zc,pi,m,ar
      real    theta,length
      real    cofx(3),cofy(3)
      integer ns,nb
c
c----------------------------------------------------------------------
ccc%%%%%%%%%%%%% from matlab

cclose all;
cAR=6;
cxc=13;
cr=0.5;
cm=4;% regular ellipse

c%%%%%%%  polar coords %%%%%%%

ctheta=linspace(0,1,50).^(m/2) * pi/2;

cxx=-AR*r*sin(theta).^(2./m);
cyy=r^(m/2.)*cos(theta);
cplot(xx,yy,'r*');

ccccc%%%%%%%
c....input
      m=4.      
      ar=6.0
      length=41.
      r = 0.75   !  radius 
      xc = 50.  !  x-coordinate of center
      zc = 54.5  !  z-coordinate of center
      ns = nimb/6. !  total number of points
      nb=nimb-1
c----------------------------------------------------------------------
      pi=3.14159265358979 

c....compute s,x(s),y(s)

c      do i=0,ns
c       ii=ns-i
c       theta= pi/2.*(real(ii)/real(ns))**(m/2.)
c       sb(i) = i
c       yb(i) = zc-r*ar*sin(theta)**(2./m)
c       xb(i) = xc+r*cos(theta)
c       write(123,'(I4, 3 G20.10)') i,theta,yb(i),xb(i)
c      enddo  


c.... 1 -- blunt nose
      
c     
      do i=0,ns
       ii=ns-i
       theta= pi/2.*(real(ii)/real(ns))**(m/2.)
       sb(i) = i
       yb(i) = zc-r*ar*sin(theta)**(2./m)
       xb(i) = xc+r*cos(theta)
      enddo  
      do i=ns+1,nimb/2
       sb(i)=i
       yb(i)=zc+(i-ns)*length/(nimb/2-ns)
       xb(i)=xc+r
c       write(501,'(I,2G20.10)') i,yb(i),xb(i)
c
      enddo 
      
      sb(nimb+2)=1
      yb(nimb+2)=yb(1)
      xb(nimb+2)=xb(1)
      
      do i=0,nimb/2
       ii=nimb+1-i
       sb(ii)=-i
       yb(ii)=yb(i)
       xb(ii)=2.*xc-xb(i)
      enddo 
c----------------------------------------------------------------------
      do i=1,nimb/2-1
       call polcoe(sb(i-1),xb(i-1),3,cofx)
       call polcoe(sb(i-1),yb(i-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
      do i=nimb/2, nimb-1
       ii=i+2
       call polcoe(sb(ii-1),xb(ii-1),3,cofx)
       call polcoe(sb(ii-1),yb(ii-1),3,cofy)
       coefx(i,1) = cofx(3)
       coefx(i,2) = cofx(2)
       coefx(i,3) = cofx(1)
       coefy(i,1) = cofy(3)
       coefy(i,2) = cofy(2)
       coefy(i,3) = cofy(1)
      enddo
c      return
c     shift to remove unphysical points
      do i=nimb/2,nimb-1
        sb(i)=sb(i+2)      
        xb(i)=xb(i+2)      
        yb(i)=yb(i+2)
      enddo	      
c
c   debug -- write stencil
c------ debug -----
c      do i=1,nimb-1
c        write(600,'(9G20.10)') sb(i),yb(i),xb(i), (coefx(i,k),k=1,3),(coefy(i,k),k=1,3)
c      enddo
c      close(600)

      return
      end
C
C-----SUBROUTINE-ObjectBluntNoseHalf-------------------------V. Ovchinnikov.06/14/04-----
C 
      subroutine ObjectBluntNoseHalf(sb,xb,yb,coefx,coefy,nb,nmax)
C
C     finds location of the object at each time step
C     builds local quadratic functions
C
      include 'headers/immbd.h'
      real    sb(0:nnmax),xb(0:nnmax),yb(0:nnmax)
      real    coefx(nnmax,3),coefy(nnmax,3)
c
      real    r,xc,zc,pi,m,ar
      real    theta,length
      real    cofx(3),cofy(3)
      integer ns,nb
c
c----------------------------------------------------------------------
ccc%%%%%%%%%%%%% from matlab

cclose all;
cAR=6;
cxc=13;
cr=0.5;
cm=4;% regular ellipse

c%%%%%%%  polar coords %%%%%%%

ctheta=linspace(0,1,50).^(m/2) * pi/2;

cxx=-AR*r*sin(theta).^(2./m);
cyy=r^(m/2.)*cos(theta);
cplot(xx,yy,'r*');

ccccc%%%%%%%
c....input
      m=4.      
      ar=6.0
      length=45.
      r = 0.75   !  radius  
      xc = 0.  !  x-coordinate of center
      zc = 54.5  !  z-coordinate of center
      ns = nimb/10. !  total number of points
      nb=nimb-2
c----------------------------------------------------------------------
      pi=3.14159265358979 

c....compute s,x(s),y(s)
c.... 1 -- blunt nose
      
     
      do i=0,ns
       ii=ns-i
       theta= pi/2.*(real(ii)/real(ns))**(m/2.)
       sb(i) = real(i)
       yb(i) = zc-r*ar*sin(theta)**(2./m)
       xb(i) = xc+r*cos(theta)
c       write(123,'(I4, 3 G20.10)') i,theta,yb(i),xb(i)
      enddo  
      do i=ns+1,nimb
       sb(i)=real(i)
c       yb(i)=zc+(i-ns)*length/(nimb-ns)
c       xb(i)=xc+r

       yb(i)=yb(i-1)+(yb(ns)-yb(ns-1)) ! 10/07 need to preserve uniform spacing because; quadratic interp very sensitive !!!
       xb(i)=xb(ns)
      enddo 
c----------------------------------------------------------------------
      do i=1,nimb-1
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
c   debug -- write stencil
c      do i=0,nimb
c        write(500,'(I5, 9G20.10)') i, sb(i), yb(i),xb(i), (coefx(i,j), j=1,3), (coefy(i,j), j=1,3)
c      enddo
c      close(500)	
c      stop
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
c
C-----SUBROUTINE-bpoints-------------------------E. Balaras..03/04/01-----
C-----Parallelization by V. Ovchinnikov-----------------------------3/04
      subroutine bpoints(x,zl,zg,sb,xb,yb,nb,coefx,coefy,
c
     &                   niml,iiml,jiml,ntop,tk,stk1,stk2,sti,nbot,bk,tib_ind,tib_c1,tib_c2,bib_ind,bib_c1,bib_c2,c4intl,
     &                   nbod, ibod, jbod, fname     )
C
C     locates i,j,k points on the cartecian grid near
C     the immersed boundary to be modified     
C
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/immbd.h'
      include 'mpif.h'
c

c   dbg
         character*5 fname
	 integer nbod, ibod((nx-2)*(nz-2)), jbod((nx-2)*(nz-2)) 
c   dbg

      real    x(nx),zg(nzg),zl(nz)
      real    sb(0:nnmax),xb(0:nnmax),yb(0:nnmax)
      real    coefx(nnmax,3),coefy(nnmax,3)
      real    c4int(nnmax,5)
      real    xnim(nnmax),ynim(nnmax),snim(nnmax)
      integer iim(nnmax),jim(nnmax)
      integer i2b(nnmax),k2b(nnmax)
      integer ibody((nzg-2)*(nx-2)),kbody((nzg-2)*(nx-2))
      integer nb,nim,nbody
c    for parallelization
      integer cpu,cpu2,cpu3,niml(0:ncpu-1),iiml(0:ncpu-1,nnmax),jiml(0:ncpu-1,nnmax),ntop(0:ncpu-1),tk(0:ncpu-1,nimb),
     &        stk1(0:ncpu-1,nnmax),stk2(0:ncpu-1,nnmax),nbot(0:ncpu-1),bk(0:ncpu-1,nimb),sti(0:ncpu-1,nnmax)
      real    tib_ind(0:ncpu-1,nimb),tib_c1(0:ncpu-1,nimb),tib_c2(0:ncpu-1,nimb),bib_ind(0:ncpu-1,nimb),
     &        bib_c1(0:ncpu-1,nimb),bib_c2(0:ncpu-1,nimb),
     &        c4intl(0:ncpu-1,nnmax,5)
      integer iii,kk1,kk2

c....local variables
      real    xx,zz,xn,yn,sa,sn
      real    a1,a2,a3,a4,tt,uu,d,h1,h2,h1h2
      integer i,k,ii,kk,n,inout
      real    wrk1(nx,nzg),wrk2(nx,nzg),wrk3(nx,nzg)
      real    x2im(nnmax),y2im(nnmax)
      real    cofx(3),cofy(3),xpas(5),zpas(5)
      integer iwrk(nx,nzg)
      integer ib(nnmax),kb(nnmax)
      real    deltax,deltaz,tlr,deltaxg,deltazg
      real    xcrd1,xcrd2,xcrd3,zcrd1,zcrd2,zcrd3
      real    xchk,zchk,h
      integer itest
      real c4sum
      integer ierr,type, step, o, p
      logical border,inside

#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
c
c  debug
c         write(0,*) myrank,' nb:',nb
c        do i=1,nb
c	 write(myrank+10,'(I,6G20.10)') (coefx(i,k),k=1,3),(coefy(i,k),k=1,3)
c        enddo
c   debug
c
c-----------------------------------------------------------------------
c....locate the velocity cells that contain markers
c    indexes are stored in arrays ib and kb
c    marker is between ib(i) and ib(i)+1 and
c                      kb(i) and kb(i)+1
      do i=1,nb
       xx=xb(i)
       call locate(x,nx,xx,ii)
       zz=yb(i)
       call locate(zg,nzg,zz,kk)
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
        if (ib(i).ne.0.and.ib(i).lt.nx-1.and.kb(i).ne.0.and.kb(i).lt.nzg-1) then 
         deltax=deltax+(x(ib(i)+1)-x(ib(i)))
         deltaz=deltaz+(zg(kb(i)+1)-zg(kb(i)))
        endif	
       enddo
       deltax=deltax/float(nb)
       deltaz=deltaz/float(nb)
       tlr=sqrt(deltax**2.+deltaz**2.)/80.
c        write(0,*) myrank,'delta ',deltax,deltaz, tlr

c-----------------------------------------------------------------------
c...compute where normal intersects boundary for the cluster points
c   surrounding the marker (This version is for a 2x2 cluster)

c...initialize work array
      do i=1,nx
       do k=1,nzg
        iwrk(i,k)=0
       enddo
      enddo

c...loop over marker points


c
c debug

      open(2,file=fname//'0',form='unformatted',status='unknown')

       write(2) nb*4
       iii=0
      DO i=1,nb 
       if (ib(i).ne.0.and.ib(i).lt.nx-1.and.kb(i).ne.0.and.kb(i).lt.nzg-1) then 
        do ii=1,2    
         do kk=1,2
          xx=x(ib(i)+ii-1)
          zz=zg(kb(i)+kk-1)
          sa=sb(i) 
          do n=1,3
           cofx(n)=coefx(i,n)
           cofy(n)=coefy(i,n)
          enddo
          iii=iii+1
          CALL Normal(iii,cofx,cofy,xx,zz,sa,sn)
          xn=coefx(i,1)*sn**2.+coefx(i,2)*sn+coefx(i,3)
          yn=coefy(i,1)*sn**2.+coefy(i,2)*sn+coefy(i,3)

      
          CALL LocPnt(cofx,cofy,xx,zz,xn,yn,sn,tlr,inout)
c         write(0,*) myrank,' locpnt'
          iwrk(ib(i)+ii-1,kb(i)+kk-1) = -inout 
          wrk1(ib(i)+ii-1,kb(i)+kk-1) = xn
          wrk2(ib(i)+ii-1,kb(i)+kk-1) = yn
          wrk3(ib(i)+ii-1,kb(i)+kk-1) = sn 
         write(2) xb(i),yb(i),xx,zz,coefx(i,1),coefx(i,2),coefx(i,3),coefy(i,1),coefy(i,2),coefy(i,3), real(-inout),sb(i), real(iii), xn, yn, sn
         enddo
        enddo
       endif
c
      ENDDO
      close(2)

cc  experimental 7/14/04
c...mark all cells inside the body with +1 flag

      step=-1 ! 1 -- start from bottom, -1 start from top
      do k=2,nzg-1
       inside=.false.
       border=.false.
C       do i=2,nx-1,step
       do i=(2+nx-1+(2-nx+1)*step)/2,(2+nx-1-(2-nx+1)*step)/2,step       
       
	if (inside) iwrk(i,k)=1

        o=iwrk(i,k)
	p=iwrk(i+step,k)

	if (border.and.(p.eq.0.or.p.eq.1)) then 
	 border=.false.
	 inside=.true.
	endif 

	if     (o.le.-1.and.p.eq.9) then
	 border=.true.
	elseif (o.eq.1.and.p.eq.9) then 
	 border=.true.
	 inside=.false.
	elseif (o.eq.9.and.p.le.-1) then
	 border=.false.
	 inside=.false.
        elseif (o.le.-1.and.p.eq.0) then 
	 inside=.false.
	 border=.false.
        elseif (o.le.-1.and.p.eq.1) then 
	 inside=.true.
	 border=.false.
        elseif (o.eq.1.and.p.le.-1) then 
	 inside=.false.
	 border=.false.
	endif 
       enddo
      enddo 	
cc
c     NOTE:  9 means inside the body, but effectively on the boundary
c           -9 means outside the body, but effectively on the boundary
c           -1 means inside the fluid
c            0 means not near an ib point
c            1 means inside the body      
c
cc remove +/-9
      do k=2,nzg-1
       do i=2,nx-1
        if(abs(iwrk(i,k)).eq.9) iwrk(i,k)=-1
       enddo
      enddo 	

cc    end of experimental	 	 
	  
ccccccccccccccccc obsolete -- will delete ccccccccccccccccccccc
c...mark all cells inside the body with +1 flag
c     commented out; VO 7/24/03

c      step=1
c      do k=2,nzg-1
c       inside=-1! flag :-1 -- outside, 1 -- inside
c       do i=2,nx-1,step
c        itest=iwrk(i,k)+iwrk(i+step,k)
c	if ((iwrk(i,k).eq.-iwrk(i+step,k)).or.(iwrk(i,k).eq.-1.and.iwrk(i+step,k).eq.0.)) then 
c	 inside=-inside
c	endif 
c        if (itest.eq.2..or.itest.eq.1.or.(iwrk(i,k).eq.-1.and.iwrk(i+step,k).eq.0.and.inside.eq.1)) then
c         iwrk(i+step,k) = 1
c        endif
c       enddo
c      enddo
cccccccccccccccccccccc end obsolete cccccccccccccccccc

c...find true boundary cells
      ii=0
      do k=2,nzg-1
       do i=2,nx-1
c
c         if(iwrk(i,k).ne.0) then !VO 7/24
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
c         write(0,*) myrank,' nim:',nim

c...find interior cells and strore in 1D array
c  commented out; VO 7/24/03
      nbody=0

      ii=0
      do k=2,nzg-1
       do i=2,nx-1

        if(iwrk(i,k).eq.1) then
           ii=ii+1
           ibody(ii)=i
           kbody(ii)=k
        endif
c
       enddo
      enddo
      nbody=ii

c...arrange all final flags into iwrk(i,k) array
      do i=1,nx-1 
       do k=1,nzg-1
         iwrk(i,k)=0
       enddo
      enddo
      do i=1,nim
         iwrk(iim(i),jim(i))=-1
      enddo
      do i=1,nbody
         iwrk(ibody(i),kbody(i))=1
      enddo

c debug     
c       do i=1,nim
c        write(123,'(2 G20.10)') xnim(i),ynim(i)
c       enddo
c       close(123)	


c....find 2nd point on the normal to the boundary used to linearize velocity

c...loop over true boundary points
      DO i=1,nim
      
c...find if the point is on the boundary.  In this case no interpolation
c   is required  (h is the distance of the point from the boundary)
      h=sqrt((xnim(i)-x(iim(i)))**2.+(ynim(i)-zg(jim(i)))**2.)
      if(h.le.tlr) then
         i2b(i)=iim(i)
         k2b(i)=jim(i)
         x2im(i) = x(iim(i))
         y2im(i) = zg(jim(i))
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
       xcrd2=x(iim(i))
       zcrd2=zg(jim(i))
       xcrd3=2.*xcrd2-xcrd1
       zcrd3=2.*zcrd2-zcrd1
c...find points around the virtual point 3 and check if other boundary 
c   points are involved in the interpolation stencil
 1000  continue
       call locate(x,nx,xcrd3,ii)
       call locate(zg,nzg,zcrd3,kk)
       itest=iwrk(ii,kk  )+iwrk(ii+1,kk  )+
     &       iwrk(ii,kk+1)+iwrk(ii+1,kk+1)
c....case corresponding to Figure 5a or 5c in the paper
c
         if(itest.eq.0) then
    
          i2b(i)=ii
          k2b(i)=kk
          x2im(i) = xcrd3
          y2im(i) = zcrd3

          h1=sqrt((xnim(i)-x(iim(i)))**2.+(ynim(i)-zg(jim(i)))**2.)
          h2=sqrt((x(iim(i))-x2im(i))**2.+(zg(jim(i))-y2im(i))**2.)
          h1h2=h1+h2

          tt = (x2im(i)-x(i2b(i)))/(x(i2b(i)+1)-x(i2b(i)))
          uu = (y2im(i)-zg(k2b(i)))/(zg(k2b(i)+1)-zg(k2b(i)))
          a1=(1-tt)*(1-uu)
          a2=   tt *(1-uu)
          a3=   tt *   uu
          a4=(1-tt)*   uu

          c4int(i,1)=h2/h1h2
          c4int(i,2)=(h1/h1h2)*a1
          c4int(i,3)=(h1/h1h2)*a2
          c4int(i,4)=(h1/h1h2)*a3
          c4int(i,5)=(h1/h1h2)*a4
ccc debug
c           if (c4int(i,2).le.0.or.c4int(i,3).le.0.or.c4int(i,4).le.0.or.c4int(i,5).le.0) then
c	       write(0,*) (c4int(i,jj),jj=1,4)
c	       write(0,*) x(i2b(i)+1),x2im(i),x(i2b(i))
c	       write(0,*) zg(k2b(i)+1),y2im(i),zg(k2b(i))
c               write(0,*) i2b(i),k2b(i)     
c	    stop
c	    endif        

c....case corresponding to Figure 5b in the paper
         elseif(itest.eq.-1) then
          i2b(i)=ii
          k2b(i)=kk
          x2im(i) = xcrd3
          y2im(i) = zcrd3

          h1=sqrt((xnim(i)-x(iim(i)))**2.+(ynim(i)-zg(jim(i)))**2.)
          h2=sqrt((x(iim(i))-x2im(i))**2.+(zg(jim(i))-y2im(i))**2.)
          h1h2=h1+h2


          tt = (x2im(i)-x(i2b(i)))/(x(i2b(i)+1)-x(i2b(i)))
          uu = (y2im(i)-zg(k2b(i)))/(zg(k2b(i)+1)-zg(k2b(i)))
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
         else
          xcrd1=xcrd2
          zcrd1=zcrd2
          xcrd2=xcrd3
          zcrd2=zcrd3
          xcrd3=2.*xcrd2-xcrd1
          zcrd3=2.*zcrd2-zcrd1
c
	  if (xcrd3.gt.x(nx-1).or.zcrd3.gt.zg(nzg-1).or.xcrd3.lt.x(1).or.zcrd3.lt.zg(1)) then 
c           cannot go farther out
             if (abs(itest).ge.2) then 
	       write(0,*)  myrank,': bpoints: exception: ib point too close to logical boundary: itest = ',itest
               write(0,*) 'xcrd3=',xcrd3,' zcrd3=',zcrd3
               write(0,*) 'xnim=',xnim(i),' ynim=',ynim(i)
	       
             else 
	       write(0,*)  myrank,': bpoints: unknown exception: itest = ',itest
	       stop
	     endif !abs(itest)	  
	  else 
           goto 1000
	  endif !xcrd            
         endif ! (itest = 0)
c
      endif !(h<tlr)
c       
 7654 continue
      ENDDO
cccccccccccccccccccc-BEGIN LOCAL DISTRIBUTION OF POINTS-cccccccccccccccccccc
c 1 - find local indices
      do n=0,ncpu-1
	niml(n)=0
      enddo
c 2 - sort ib points
c     initialize arrays
      do n=0,ncpu-1
        ntop(n)=0
        nbot(n)=0
      enddo
c   debug
c      write(0,*) myrank,'here..',nim
c       do i=1,nim
c        write(myrank+100,'(I,6G20.10)') i,(c4int(i,k),k=1,5),c4int(i,1)+c4int(i,2)+c4int(i,3)+c4int(i,4)+c4int(i,5)
c       enddo
c       close(myrank+100)

c    debug
      do i=1,nim
        cpu=(jim(i)-2)/(nz-2)
	niml(cpu)=niml(cpu)+1
	ii=niml(cpu)
	
	iiml(cpu,ii)=iim(i)
	jiml(cpu,ii)=jim(i)-cpu*(nz-2)
	
	cpu2=(k2b(i)-2)/(nz-2) ! cpu of 1st stencil point
	kk1=k2b(i)-cpu2*(nz-2) ! local index of 1st stencil point
	
	cpu3=(k2b(i)+1-2)/(nz-2) ! cpu of 2nd stencil point
	kk2=k2b(i)+1-cpu3*(nz-2) ! local index of 2nd stencil point
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc first point
         if (cpu2.eq.cpu+1) then
	 
	  ntop(cpu)=ntop(cpu)+1
	  iii=ntop(cpu)
	  tib_ind(cpu,iii)=ii
	  tib_c1(cpu,iii)=c4int(i,2)
	  tib_c2(cpu,iii)=c4int(i,3)
          tk(cpu,iii)=kk1 ! local index of point to pass
c         need to receive u(cpu,iii,1) and u(cpu,iii,2)	   
          c4intl(cpu,ii,2)=0.
	  c4intl(cpu,ii,3)=0.
	  stk1(cpu,ii)=1 ! set arbitrarily

         elseif (cpu2.eq.cpu-1) then
	 
	  nbot(cpu)=nbot(cpu)+1
	  iii=nbot(cpu)
	  bib_ind(cpu,iii)=ii
	  bib_c1(cpu,iii)=c4int(i,2)
	  bib_c2(cpu,iii)=c4int(i,3)
          bk(cpu,iii)=kk1 ! local index of point to pass
c         need to receive u(cpu,iii,1) and u(cpu,iii,2)	   
          c4intl(cpu,ii,2)=0.
	  c4intl(cpu,ii,3)=0.
	  stk1(cpu,ii)=1 ! set arbitrarily
         
	else
	  stk1(cpu,ii)=kk1
	  c4intl(cpu,ii,2)=c4int(i,2)
	  c4intl(cpu,ii,3)=c4int(i,3)
        endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc second point
         if (cpu3.eq.cpu+1) then
	 
	  ntop(cpu)=ntop(cpu)+1
	  iii=ntop(cpu)
	  tib_ind(cpu,iii)=ii
	  tib_c1(cpu,iii)=c4int(i,5)
	  tib_c2(cpu,iii)=c4int(i,4)
          tk(cpu,iii)=kk2 ! local index of point to pass
c         need to receive u(cpu,iii,1) and u(cpu,iii,2)	   
          c4intl(cpu,ii,4)=0.
	  c4intl(cpu,ii,5)=0.
	  stk2(cpu,ii)=1 ! set arbitrarily

         elseif (cpu3.eq.cpu-1) then
	 
	  nbot(cpu)=nbot(cpu)+1
	  iii=nbot(cpu)
	  bib_ind(cpu,iii)=ii
	  bib_c1(cpu,iii)=c4int(i,5)
	  bib_c2(cpu,iii)=c4int(i,4)
          bk(cpu,iii)=kk2 ! local index of point to pass
c         need to receive u(cpu,iii,1) and u(cpu,iii,2)	   
          c4intl(cpu,ii,4)=0.
	  c4intl(cpu,ii,5)=0.
	  stk2(cpu,ii)=1 ! set arbitrarily
         
	else
	  stk2(cpu,ii)=kk2
	  c4intl(cpu,ii,4)=c4int(i,4)
	  c4intl(cpu,ii,5)=c4int(i,5)
        endif

	  sti(cpu,ii)=i2b(i)
          c4intl(cpu,ii,1)=c4int(i,1)
      enddo
c	
c       experimental -- distribute body points
 	 nbod=0
       do i=1,nbody
        cpu=(kbody(i)-2)/(nz-2)
        if (cpu.eq.myrank) then 
	 nbod=nbod+1
	 ii=nbod
	
	 ibod(ii)=ibody(i)
	 jbod(ii)=kbody(i)-cpu*(nz-2)
	endif 
       enddo
c       experimental -- distribute body points

cccccccccccc debug
c...sum of interpolation coef.
c       do i=1,nim
c        dd= c4int(i,1)+c4int(i,2)+c4int(i,3)+c4int(i,4)+c4int(i,5)
c        write(6,*) i,dd
c       enddo

c...grid
      if (myrank.eq.0) then
        open(110,file=fname,form='unformatted',status='unknown')
        write(110) nx,nzg 
        write(110) (x(i),i=1,nx)
	write(110) (zg(i),i=1,nzg)
c...all exterior control points
        write(110) nim
         write(110) (i,x(iim(i)),zg(jim(i)),xnim(i),ynim(i),
     &              x2im(i),y2im(i),i=1,nim)
        write(110) nb,(xb(i),yb(i),i=1,nb)
        write(110) nbody,(x(ibody(i)),zg(kbody(i)),i=1,nbody)
       close (110)
c      write(0,*) myrank,': aborting..'
c      stop
      endif
cccccccccccc debug
c...all interior control points
c      do i=1,nbody
c       write(140,112)i,xu(ibody(i)),zu(kbody(i))
c      enddo

c 111  format(i4,6(1x,e14.7))
c 112  format(i4,2(1x,e14.7))
      return
      end
C
C 
C-----SUBROUTINE-Normal-------------------------E. Balaras..03/04/01-----
C     this version by VO 10/07; there is no apparent benefit over Elias's version
      SUBROUTINE Normal(i,coefx,coefy,x,y,sa,sn)
c
c     computes the arclength coordinate at the interception
c     with the normal on (i,k)
c
      real coefx(3),coefy(3)
      real x,y,sn
      real a,b,c,d
      integer n
      
      real ax,bx,cx,ay,by,cy
      real*8 coeff(0:3), ans(3)
      
c
c      write(0,*) x,y,sa,'normal'

      ax=coefx(1)
      bx=coefx(2)
      cx=coefx(3)
      ay=coefy(1)
      by=coefy(2)
      cy=coefy(3)

      coeff(3)=dble(2.*ay**2.+2.*ax**2.)
      coeff(2)=dble(3.*ax*bx+3.*ay*by  )
      coeff(1)=dble(2.*ay*cy+2.*ax*cx+bx*bx+by*by-2.*ay*y-2.*ax*x)
      coeff(0)=dble(cx*bx+by*cy-x*bx-y*by)
c
      call cubic(coeff, ans, n)
      if (n.eq.0) then 
       stop 'Unexpected error in Normal: no cubic roots found'
      endif
      
c     pick the closest root 
      sn=ans(1)
      do j=2,n
       if (abs(sn-sa) .gt. abs(ans(j)-sa) ) then 
        sn=ans(j)
       endif
      enddo	
      

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
c   or in the fluid inout=-1.  In the solid inout=+1
c
ccccc IMPORTANT: (VO) Actually, for the flipped case (x & z directions switched)  inout=-1 is in the solid.
c
      h=sqrt((xn-xi)**2.+(yn-yj)**2.)
c
      if(delta.lt.0.0) then
       if (h.le.tlr) then 
        inout=-9 ! so this corresponds to a point inside the body within tol of the bdry.
       else	
        inout=-1 !inside the body
       endif	
      else
       if (h.le.tlr) then
        inout=9 ! new: this corresponds to a point outside the body within tol of the boundary
       else 
        inout=1  !outside the body
       endif
      endif

      return
      end

C
C 
C-----SUBROUTINE-Forcing-------------------------E. Balaras..03/04/01-----
c     parallelization by V. Ovchinnikov 3/04
C
c
      SUBROUTINE Forcing(us,uim,
c
     & niml,iiml,jiml,ntop,tk,stk1,stk2,nbot,bk,sti,tib_ind,tib_c1,tib_c2,bib_ind,bib_c1,bib_c2,c4intl) 
C
C          modify u*
C          because of the explicit nature of the algorithm, the 
C          application of forcing (as explained in the C&F paper) 
C          is equivalent in directly modifying (imposing bc) u*
C
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/immbd.h'
      include 'mpif.h'

      real    us(nx,ny,nz)
      real    uim(nnmax)
c     parallelization variables
c 
      integer me,you,ierr,type,i,j,k,iii,status(MPI_STATUS_SIZE)
      integer niml(0:ncpu-1),iiml(0:ncpu-1,nnmax),jiml(0:ncpu-1,nnmax),ntop(0:ncpu-1),tk(0:ncpu-1,nimb),
     &        stk1(0:ncpu-1,nnmax),stk2(0:ncpu-1,nnmax),nbot(0:ncpu-1),bk(0:ncpu-1,nimb),sti(0:ncpu-1,nnmax)
      real    tib_ind(0:ncpu-1,nimb),tib_c1(0:ncpu-1,nimb),tib_c2(0:ncpu-1,nimb),bib_ind(0:ncpu-1,nimb),
     &        bib_c1(0:ncpu-1,nimb),bib_c2(0:ncpu-1,nimb),
     &        c4intl(0:ncpu-1,nnmax,5)
      real    u(2,ny,nimb),v(2,ny,nimb)

#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

      me=myrank
      do j=jy1,jy2+1
       do i=1,niml(me)
       
        us(iiml(me,i),j,jiml(me,i))= c4intl(me,i,1)*uim(i) +
     &                       c4intl(me,i,2)*us(sti(me,i)   ,j,stk1(me,i) ) +
     &                       c4intl(me,i,3)*us(sti(me,i)+1 ,j,stk1(me,i) ) +
     &                       c4intl(me,i,4)*us(sti(me,i)+1 ,j,stk2(me,i) ) +
     &                       c4intl(me,i,5)*us(sti(me,i)   ,j,stk2(me,i) )
       enddo
      enddo
c     nonlocal part      
c     top
c     step 1 -- send/receive data
      if (me.gt.0) then
       you=me-1
       do iii=1,ntop(you)
         i=sti(you,tib_ind(you,iii))
	 k=tk(you,iii)
	 do j=jy1,jy2+1 
           v(1,j,iii)=us(i,j,k)
	   v(2,j,iii)=us(i+1,j,k)
	 enddo
       enddo
       call mpi_send(v,ntop(you)*ny*2,type,you,100+me,MPI_COMM_HYDRO,ierr)
      endif 
      if (me.lt.ncpu-1) then 
        call mpi_recv(u,ntop(me)*ny*2,type,me+1,100+me+1,MPI_COMM_HYDRO,status,ierr) 
      endif
c
      do iii=1,ntop(me)
       i=iiml(me,tib_ind(me,iii))
       k=jiml(me,tib_ind(me,iii))
       do j=jy1,jy2+1
         us(i,j,k)=us(i,j,k)+tib_c1(me,iii)*u(1,j,iii)+tib_c2(me,iii)*u(2,j,iii) !i and i+1 points, resp
       enddo
      enddo
c     bottom
c     step 1 -- send/receive data
      if (me.lt.ncpu-1) then
       you=me+1
       do iii=1,nbot(you)
         i=sti(you,bib_ind(you,iii))
	 k=bk(you,iii)
	 do j=jy1,jy2+1 
           v(1,j,iii)=us(i,j,k)
	   v(2,j,iii)=us(i+1,j,k)
	 enddo
       enddo
       call mpi_send(v,nbot(you)*ny*2,type,you,100+me,MPI_COMM_HYDRO,ierr)
c       write(0,*) me,' sent to ',you
      endif 
      if (me.gt.0) then 
        call mpi_recv(u,nbot(me)*ny*2,type,me-1,100+me-1,MPI_COMM_HYDRO,status,ierr) 
c        write(0,*) me,' receiving from ',me-1
      endif 

      do iii=1,nbot(me)
       i=iiml(me,bib_ind(me,iii))
       k=jiml(me,bib_ind(me,iii))
       do j=jy1,jy2+1
         us(i,j,k)=us(i,j,k)+bib_c1(me,iii)*u(1,j,iii)+bib_c2(me,iii)*u(2,j,iii) !i and i+1 points, resp
       enddo
      enddo
c       write(0,*) me,': ib out'             

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
