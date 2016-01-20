C
C-----------------------------------------------------------------------
C                  ***************************                         
C                  *          plot.f         *                        
C                  ***************************                       
C----------------------------------------------------------------------- 
C
C This file contains all necessary routines to create plot files
C for EDDY (normally based on tecplot)
C
C	- TecHist : plots comp. history
C	- Tec2d   : plots 2-d slices (x,z)
C	- Tec3d   : plots instantaneous velocity field
C	- TecVar3d: plots one selected variable
C	- TecShot : inst. 2-d cut of the flow
C
C----------------------------------------------------------------------- 
C
C 
C-----SUBROUTINE-TecHist------------------------P. Flohr--20/02/1994----
C
      SUBROUTINE TecHist(TimeHist,TkeHist,DivHist,TwallHist,DisplHist,
     &                   FileName, Miter, Niter,Nstep) 
C
C-----------------------------------------------------------------------
C
C 
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C 
      CHARACTER*72  CTEXT
C
      DIMENSION      TimeHist(Miter),
     &               TkeHist(Miter), DivHist(Miter),
     &               TwallHist(Nstep,Miter),DisplHist(Nstep,Miter)
      INTEGER        Miter, Niter,Nstep
      CHARACTER*(*)  FileName
C
      INTEGER        i, j, FileUnit, ii
      CHARACTER*(60) Title, Variables,ZoneTitle
      DIMENSION      Data(20)
      INTEGER        iocheck
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA
C
C
C
C-----------------------------------------------------------------------
C 
      FileUnit = 90
      print*,'Simulation History'
      Title = 'simulation history'
      Variables = 'time,tke,div,twall1,twall2,twall3,disp1,disp2,disp3'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....begin zone record
C
      print*, 'begin zone record'
      ZoneTitle = 'history'
      iocheck = TECzneA(ZoneTitle, Niter, 1, 1, 'POINT', FileUnit)
C
C
C.....begin data record
C
      print*, 'begin data record'
      DO 10 i = 1, Niter
C
C *** Write History Ascii files ****
C        write(50,*) TimeHist(i),TkeHist(i),DivHist(i)
C        write(51,*) TimeHist(i),(TwallHist(j,i),j=1,Nstep)
C        write(52,*) TimeHist(i),(DisplHist(j,i),j=1,Nstep)
C ***           ****
C

        Data(1) = TimeHist(i)
        Data(2) = TkeHist(i)
        Data(3) = DivHist(i)
        ii = 3
        DO 101 j = 1,Nstep
          ii = ii + 1
          Data(ii) = TwallHist(j,i)
 101    CONTINUE
        DO 102 j = 1,Nstep
          ii = ii + 1
          Data(ii) = DisplHist(j,i)
 102    CONTINUE
        iocheck = TECdatA(ii, Data, FileUnit)
 10   CONTINUE
C
C..... close file
C
      print*,'close file'
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C 
C-----SUBROUTINE-Tec2d--------------------------P. FLOHR--20/02/1994----
C
      SUBROUTINE Tec2d(x,z,uo,vo,wo,uavt,p,rot,twavt,FileName,nx,ny,nz,
     &                 mx,my,mz) 
C
C-----------------------------------------------------------------------
C
C 
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C 
      CHARACTER*72  CTEXT
C
      REAL          x(mx), z(mz), 
     &              uo(mx,my,mz),vo(mx,my,mz),wo(mx,my,mz),
     &              rot(mx,my,mz),
     &              uavt(mx,my,mz),p(mx,my,mz),twavt(mx,my,mz)
      INTEGER       nx,ny,nz,mx,my,mz
      CHARACTER*(*) FileName
C
      INTEGER        i, j, k, FileUnit
      CHARACTER*(45) Title, Variables,ZoneTitle
      DIMENSION      Data(9)
      INTEGER        iocheck
      REAL           rotx, roty, rotz
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA
C
C-----------------------------------------------------------------------
C 
      FileUnit = 90
      Title = '2-d sliced data'
      Variables = 'x,z,ufl,vfl,wfl,uavt,p,rot,twavt'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....begin zone record
C
      ZoneTitle = 'statistics'
      iocheck = TECzneA(ZoneTitle, nz-2, nx-2, 1,'POINT', FileUnit)
C
C
C.....begin data record
C
      j = (jy2-jy1)/2
C
      DO 10 i = 2, nx-1
          DO 10 k = 2, nz-1
            Data(1) = 0.5*(x(i)+x(i-1))
            Data(2) = 0.5*(z(k)+z(k-1))
            Data(3) = 0.5*(uo(i,j,k)+uo(i-1,j,k))-uavt(i,j,k)
            Data(4) = 0.5*(vo(i,j,k)+vo(i,j-1,k))
            Data(5) = 0.5*(wo(i,j,k)+wo(i,j,k-1))
            Data(6) = uavt(i,j,k)
            Data(7) = p(i,j,k)
            Data(8) = rot(i,j,k)
            Data(9) = twavt(i,j,k)
            iocheck = TECdatA(9, Data, FileUnit)
 10   CONTINUE
C
C..... close file
C
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C-----SUBROUTINE-Tec3d--------------------------P. FLOHR--08/01/1994----
C
      SUBROUTINE Tec3d(x,y,z,uo,vo,wo,uavt,p,rot,FileName,nx,ny,nz,
     &                                                    mx,my,mz)
C
C-----------------------------------------------------------------------
C
C 
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C
      DIMENSION      x(mx), y(my), z(mz),
     &               uo(mx,my,mz), vo(mx,my,mz), wo(mx,my,mz),
     &               uavt(mx,my,mz), p(mx,my,mz), rot(mx,my,mz)
      CHARACTER*(*)  FileName
      INTEGER        nx,ny,nz,mx,my,mz
C
      INTEGER        i, j, k, FileUnit
      CHARACTER*(60) Title, Variables,ZoneTitle
      DIMENSION      Data(10)
      INTEGER        iocheck
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA 
C
C-----------------------------------------------------------------------
C
      FileUnit = 90
      Title = 'Instantaneous Velocity Field'
      Variables = 'x,y,z,ufl,vfl,wfl,p,rot'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....begin zone record
C
      ZoneTitle = 'u,v,w'
      iocheck = TECzneA(ZoneTitle, nz-2, ny-2, nx-2,'POINT',FileUnit)
C
C
C.....begin data record
C

      DO 10 i = 2, nx-1
        DO 10 j = 2, ny-1
          DO 10 k = 2, nz-1
            Data(1) = 0.5*(x(i)+x(i-1))
            Data(2) = 0.5*(y(j)+y(j-1))
            Data(3) = 0.5*(z(k)+z(k-1))
            Data(4) = 0.5*(uo(i,j,k)+uo(i-1,j,k))-uavt(i,j,k)
            Data(5) = 0.5*(vo(i,j,k)+vo(i,j-1,k))
            Data(6) = 0.5*(wo(i,j,k)+wo(i,j,k-1))
            Data(7) = p(i,j,k)
            Data(8) = rot(i,j,k)
C
            iocheck = TECdatA(8, Data, FileUnit)
 10   CONTINUE
C
C..... close file
C
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C 
C-----SUBROUTINE-TecMeanProf--------------------P. FLOHR--20/02/1994----
C
      SUBROUTINE TecMeanProf(x,z,u,tw,nx,ny,nz,FileName) 
C
C-----------------------------------------------------------------------
C 
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C 
      CHARACTER*72  CTEXT
C
      DIMENSION      x(nx),z(nz),u(nx,ny,nz),
     &               tw(nx,ny,nz)
      INTEGER        nx,ny,nz
      CHARACTER*(*)  FileName
C
      INTEGER        i, j, k
      INTEGER        FileUnit
      CHARACTER*(60) Title, ZoneTitle, Variables
      CHARACTER*(5)  CharIter
      DIMENSION      Data(35)
      INTEGER        iocheck
      REAL           utau, twall, zz, zplus, uplus
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA
C
C-----------------------------------------------------------------------
C
      FileUnit  = 10
      Title     = 'Velocity Profiles'
      Variables = 'z,z+,u,u+'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
      j    = (jy2+jy1)/2
      DO 70 i = ix1+2, ix2, 3
        WRITE(CharIter,3000) i
 3000   FORMAT(I4)
        ZoneTitle = 'xpos ='//CharIter
        iocheck = TECzneA(ZoneTitle,kz2-kz1-5,1,1,'POINT',FileUnit)
        twall = tw(i,j,kz1+1)
        utau  = SQRT(twall)
        DO 701 k = kz1+1, kz2-5
          zz    = 0.5*(z(k)+z(k-1))
          zplus = utau * zz / RU1
          uplus = u(i,j,k)/utau
c          Data(1) = zz
          Data(1) = zplus
c          Data(3) = u(i,j,k)
          Data(2) = uplus
          iocheck = TECdatA(2, Data, FileUnit)
 701    CONTINUE
 70   CONTINUE
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C 
C 
C-----SUBROUTINE-TecVar3d-----------------------P. FLOHR--28/02/1994----
C
      SUBROUTINE TecVar3d(x,y,z,nx,ny,nz,Var,FileName) 
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C
      DIMENSION      x(nx), y(ny), z(nz),
     &               Var(nx,ny,nz)
      CHARACTER*(*)  FileName
      INTEGER        nx, ny, nz
C
      INTEGER        i, j, k, FileUnit
      CHARACTER*(45) Title, Variables,ZoneTitle
      DIMENSION      Data(4)
      INTEGER        iocheck
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA
C
C-----------------------------------------------------------------------
C
      FileUnit = 90
      Title = 'Instantaneous Field (one variable)'
      Variables = 'x, y, z, var'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....begin zone record
C
      ZoneTitle = 'zone-1'
      iocheck = TECzneA(ZoneTitle, nz-2, ny-2, nx-2,'POINT',FileUnit)
C
C
C.....begin data record
C

      DO 10 i = 2, nx-1
        DO 10 j = 2, ny-1
          DO 10 k = 2, nz-1
            Data(1) = x(i)
            Data(2) = y(j)
            Data(3) = z(k)
            Data(4) = Var(i,j,k)
            iocheck = TECdatA(4, Data, FileUnit)
 10   CONTINUE
C
C..... close file
C
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C 
C-----SUBROUTINE-TecTurvis----------------------P. FLOHR--20/05/1994----
C
      SUBROUTINE TecTurvis(x,y,z,tnitx,tnity,tnitz,nx,ny,nz,FileName) 
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C
      REAL           x(nx), y(ny), z(nz),
     &               tnitx(nx,ny,nz),tnity(nx,ny,nz),tnitz(nx,ny,nz)
      CHARACTER*(*)  FileName
      INTEGER        nx, ny, nz
C
      INTEGER        i, j, k, FileUnit
      CHARACTER*(45) Title, Variables,ZoneTitle
      REAL           Data(6)
      INTEGER        iocheck
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA
C
C-----------------------------------------------------------------------
C                            plot turb. viscosity versus molecular visc.
C-----------------------------------------------------------------------
C
C
      FileUnit = 90
      Title = 'turbulent viscosities'
      Variables = 'x,y,z,turbx,turby,turbz'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....begin zone record
C
      ZoneTitle = 'zone-1'
      iocheck = TECzneA(ZoneTitle, nz-2, ny-2, nx-2,'POINT',FileUnit)
C
C
C.....begin data record
C
      DO 10 i = 2, nx-1
        DO 10 j = 2, ny-1
          DO 10 k = 2, nz-1
            Data(1) = 0.5 * (x(i)+x(i-1))
            Data(2) = 0.5 * (y(j)+y(j-1))
            Data(3) = 0.5 * (z(k)+z(k-1))
            Data(4) = (tnitx(i,j,k) + tnitx(i-1,j,k))/2/ru1
            Data(5) = (tnity(i,j,k) + tnity(i,j-1,k))/2/ru1
            Data(6) = (tnitz(i,j,k) + tnitz(i,j,k-1))/2/ru1
            iocheck = TECdatA(6, Data, FileUnit)
 10   CONTINUE
C
C..... close file
C
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C 
C-----SUBROUTINE-TecVar2d-----------------------P. FLOHR--28/02/1994----
C
      SUBROUTINE TecVar2d(x,y,z,nx,ny,nz,Var,FileName) 
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C
      DIMENSION      x(nx), y(ny), z(nz),
     &               Var(nx,ny,nz)
      CHARACTER*(*)  FileName
      INTEGER        nx, ny, nz
C
      INTEGER        i, j, k, FileUnit
      CHARACTER*(45) Title, Variables,ZoneTitle
      DIMENSION      Data(4)
      INTEGER        iocheck
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA 
C
C-----------------------------------------------------------------------
C
      FileUnit = 90
      Title = 'Instantaneous Field (one variable)'
      Variables = 'x, z, var'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....begin zone record
C
      ZoneTitle = 'zone-1'
      iocheck = TECzneA(ZoneTitle, nz, nx, 1,'POINT',FileUnit)
C
C
C.....begin data record
C
      j = 0.5 * (jy1 + jy2)
c
      DO 10 i = 1, nx
          DO 10 k = 1, nz
            Data(1) = x(i)
            Data(2) = z(k)
            Data(3) = Var(i,j,k)
            iocheck = TECdatA(3, Data, FileUnit)
 10   CONTINUE
C
C..... close file
C
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C-----SUBROUTINE-TecShot------------------------P. FLOHR--20/02/1994----
C
      SUBROUTINE TecShot(x,z,uo,vo,wo,p,nx,ny,nz,FileName) 
C
C-----------------------------------------------------------------------
C
C 
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C 
      CHARACTER*72  CTEXT
C
      DIMENSION      x(nx), z(nz),
     &               uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz),
     &               p(nx,ny,nz)
      INTEGER        nx,ny,nz
      CHARACTER*(*)  FileName
C
      INTEGER        i, j, k, FileUnit
      CHARACTER*(45) Title, Variables,ZoneTitle
      DIMENSION      Data(6)
      INTEGER        iocheck
      REAL           rotx,roty,rotz
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA
C
C
C
C-----------------------------------------------------------------------
C 
      FileUnit = 90
      Title = 'Movie'
c      Variables = 'x,z,u,w,p,rot'
      Variables = 'x,z,rot'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....begin zone record
C
      ZoneTitle = 'u_i_n_s_t'
      iocheck = TECzneA(ZoneTitle,kz2-1,ix2-1,1,'POINT',FileUnit)
C
C
C.....begin data record
C
      j = (jy2-jy1)/2
C
      DO 10 i = ix1+1, ix2
          DO 10 k = kz1+1, kz2
            Data(1) = x(i)
            Data(2) = z(k)
c            Data(3) = uo(i,j,k)
c            Data(4) = wo(i,j,k)
c            Data(5) = p(i,j,k)
            rotx = .5 * ( ( (av(i)/dx*(vo(i+1,j,k)-vo(i,j,k)))
     &                     -(bu(j)/dy*(uo(i,j+1,k)-uo(i,j,k))) )
     &                  + ( (av(i)/dx*(vo(i+1,j-1,k-1)-vo(i,j-1,k-1)))
     &                     -(bu(j-1)/dy*(uo(i,j,k-1)-uo(i,j-1,k-1)))) ) 
            roty = .5 * ( ( (bw(j)/dy*(wo(i,j+1,k)-wo(i,j,k)))
     &                     -(cv(k)/dz*(vo(i,j,k+1)-vo(i,j,k))) )
     &                  + ( (bw(j)/dy*(wo(i-1,j+1,k-1)-wo(i-1,j,k-1)))
     &                     -(cv(k-1)/dz*(vo(i-1,j,k)-vo(i-1,j,k-1)))) )
            rotz = .5 * ( ( (cu(k)/dz*(uo(i,j,k+1)-uo(i,j,k)))
     &                     -(aw(i)/dx*(wo(i+1,j,k)-wo(i,j,k))) )
     &                  + ( (cu(k)/dz*(uo(i-1,j-1,k+1)-uo(i-1,j-1,k)))
     &                     -(aw(i-1)/dx*(wo(i,j-1,k)-wo(i-1,j-1,k)))) )
            Data(3) = SQRT(rotx*rotx + roty*roty + rotz*rotz)
            iocheck = TECdatA(3, Data, FileUnit)
 10   CONTINUE
C
C..... close file
C
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
C 
C-----SUBROUTINE-TecStreak----------------------P. FLOHR--18/05/1994----
C
      SUBROUTINE TecStreak(xu,yv,zw,StrkPos,StrkIndx,Miter,icycle,
     &                     ipart,nx,ny,nz,npart,FileName) 
C
C-----------------------------------------------------------------------
C
C 
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C 
      CHARACTER*72  CTEXT
C
      DIMENSION      xu(nx), yv(ny), zw(nz),
     &               StrkPos(Miter,npart,3)
      INTEGER        StrkIndx(Miter,npart,3)
      INTEGER        nx,ny,nz,npart,Miter,icycle
      CHARACTER*(*)  FileName
C
      INTEGER        i, j, k, FileUnit
      CHARACTER*(45) Title, Variables,ZoneTitle
      DIMENSION      Data(3)
      INTEGER        iocheck
      INTEGER TECiniA, TECzneA, TECdatA ,TECnodA, TECendA
C
C-----------------------------------------------------------------------
C 
      FileUnit = 90
      Title = 'Streaks'
      Variables = 'x,y,z'
C
      iocheck = TECiniA(Title, Variables, FileName, FileUnit)
C
C.....zone: grid box
C
      ZoneTitle = 'grid box'
      iocheck = TECzneA(ZoneTitle,5,1,1,'POINT',FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
C
      ZoneTitle = 'grid box'
      iocheck = TECzneA(ZoneTitle,5,1,1,'POINT',FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
C
      ZoneTitle = 'grid box'
      iocheck = TECzneA(ZoneTitle,4,1,1,'POINT',FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
C
      ZoneTitle = 'grid box'
      iocheck = TECzneA(ZoneTitle,4,1,1,'POINT',FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz2-1)
      iocheck = TECdatA(3, Data, FileUnit)
      Data(1) = xu(ix2-1)
      Data(2) = yv(jy2-1)
      Data(3) = zw(kz1)
      iocheck = TECdatA(3, Data, FileUnit)
C
C.....zones: streaks
C
      DO 20 i = 1, npart
        Zonetitle = 'streak'
        iocheck = TECzneA(ZoneTitle,ipart,1,1,'POINT',FileUnit)
        DO 20 j = 1, ipart
          Data(1) = StrkPos(j,i,1)
          Data(2) = StrkPos(j,i,2)
          Data(3) = StrkPos(j,i,3)
          iocheck = TECdatA(3, Data, FileUnit)
 20   CONTINUE
C
C..... close file
C
      iocheck = TECendA(FileUnit)
C 
      RETURN 
      END 
C 
