      module grid
      use size
      use parser, only: command
      use boundary ! boundary conditions
!     contains all parameters for the computational grid
      private
#ifdef PARALLEL
       include 'mpif.h'
#endif
!
!
      public grid_initialize
      public grid_done
!
!      
!.....xu, yv, wz:		cell faces
!.....xc, yc, zc:               cell centers
      float, allocatable, public, save :: xu,xv,xw,xp,xc,xcf,
     &                                    yu,yv,yw,yp,yc,ycf,
     &                                    zu,zv,zw,zp,zc,zcf,
     &                                    zug,zvg,zwg,zpg,zcg,
! metrics
     &                                    app,ap,au,av,aw,
     &                                    bpp,bp,bu,bv,bw,
     &                                    cpp,cp,cu,cv,cw,
     &                                    cppg,cpg,cug,cvg,cwg
!
!     
      logical, public, save :: grid_initialized=.false.
!
!
      contains
       subroutine grid_initialize()
       use errors
       implicit none

       character*(80) :: xgrid_file, ygrid_file, zgrid_file
       int :: iostatus=0
       int :: i, j, k
       int :: bug
       float :: x, y, z
       float :: dxmax, dxmin, dymax, dymin, dzmax, dxmin
!
       character*(14) :: whoami
       data whoami /'GRID_INITIALIZE'/
!
       if (grid_initialized) call grid_done() ! clean up memory if needed
!
       if (.not.size_initialized) then 
        call error(me, whoami, ' SIZE NOT INITIALIZED. NOTHING DONE.')
        return
       endif
!
       if (.not.boundary_initialized) then 
        call error(me, whoami, ' BC NOT INITIALIZED. NOTHING DONE.')
        return
       endif
!
!     allocate arrays
       allocate(xu(nx),xv(nx),xw(nx),xp(nx),xc(nx),xcf(nx))
       allocate(yu(ny),yv(ny),yw(ny),yp(nx),yc(ny),ycf(ny))
       allocate(zu(nz),zv(nz),zw(nz),zp(nz),zc(nz),zcf(nz))
       allocate(zug(nzg),zvg(nzg),zwg(nzg),zpg(nzg),zcg(nzg))
       allocate(app(nx),ap(nx),au(nx),av(nx),aw(nx),
     &          bpp(ny),bp(ny),bu(ny),bv(ny),bw(ny),
     &          cpp(nz),cp(nz),cu(nz),cv(nz),cw(nz),
     &          cppg(nzg),cpg(nzg),cug(nzg),cvg(nzg),cwg(nzg))
!
! x-grid
       xgrid_file=command('XGRID') ! query parser module for grid filename 
       if (.not.iesmpty(xgrid_file) then 
        if (me.le.0) then 
         iostatus=0
         open(unit=100, file=xgrid_file, form='FORMATTED', status='OLD', IOSTAT=iostatus)
         i=0
         read(100, *) x
         do while (iostatus.ge.0)
          i=i+1
          if (i.le.nx) xu(i)=x
         enddo
         close(100)
        endif
       else
        call error(0, whoami, 'X-GRID FILE NOT SPECIFIED. ABORTING.')
        return
       endif
!
       ygrid_file=command('YGRID')
       if (.not.iesmpty(ygrid_file) then 
        if (me.le.0) then 
         iostatus=0
         open(unit=100, file=ygrid_file, form='FORMATTED', status='OLD', IOSTAT=iostatus)
         j=0
         read(100, *) y
         do while (iostatus.ge.0)
          j=j+1
          if (j.le.ny) yv(j)=y
         enddo
         close(100)
        endif
       else
        call error(0, whoami, 'Y-GRID FILE NOT SPECIFIED. ABORTING.')
        return
       endif
!
       zgrid_file=command('ZGRID')
       if (.not.iesmpty(zgrid_file) then 
        if (me.le.0) then 
         iostatus=0
         open(unit=100, file=zgrid_file, form='FORMATTED', status='OLD', IOSTAT=iostatus)
         k=0
         read(100, *) z
         do while (iostatus.ge.0)
          k=k+1
          if (k.le.nz) zw(k)=z
         enddo
         close(100)
        endif
       else
        call error(0, whoami, 'Z-GRID FILE NOT SPECIFIED. ABORTING.')
        return
       endif
!
#if PARALLEL
       if (ncpu.gt.1) then 
        call MPI_BCAST(i,1,mpi_int,0,communicator,bug)
        call MPI_BCAST(j,1,mpi_int,0,communicator,bug)
        call MPI_BCAST(k,1,mpi_int,0,communicator,bug)
!
        call MPI_BCAST(xu,nx,mpi_float,0,communicator,bug)
        call MPI_BCAST(yv,ny,mpi_float,0,communicator,bug)
        call MPI_BCAST(zw,nz,mpi_float,0,communicator,bug)
       endif
#endif	 
       if     (i.gt.nx) then 
        call warning(0, whoami, 'X-GRID FILE HAS MORE LINES THAN EXPECTED.')
       elseif (i.lt.nx) then 
        call error(0, whoami, 'X-GRID FILE DOES NOT CONTAIN ENOUGH ENTRIES.')
        return
       elseif (j.gt.ny) then 
        call warning(0, whoami, 'Y-GRID FILE HAS MORE LINES THAN EXPECTED.')
       elseif (j.lt.ny) then 
        call error(0, whoami, 'Y-GRID FILE DOES NOT CONTAIN ENOUGH ENTRIES.')
        return
       elseif (k.gt.nzg) then 
        call warning(0, whoami, 'Z-GRID FILE HAS MORE LINES THAN EXPECTED.')
       elseif (k.lt.nzg) then 
        call error(0, whoami, 'Z-GRID FILE DOES NOT CONTAIN ENOUGH ENTRIES.')
        return
       endif
!
!    offset for staggered components
!
       do i=ix1+1,ix2+1
         xv(i)=(xu(i)+xu(i-1))*0.5
       enddo
       xv(ix1)=xu(ix1)-0.5*(xu(ix1+1)-xu(ix1))
       xw=xv; xp=xv
!
       do j=jy1+1,jy2+1
         yu(j)=(yv(j)+yv(j-1))*0.5
       enddo
       yu(jy1)=yv(jy1)-(yv(jy1+1)-yv(jy1))*0.5
       yw=yu; yp=yu
!
       do k=2,nzg
         zug(k)=(zwg(k)+zwg(k-1))*0.5
       enddo
       zug(1)=zwg(1)-0.5*(zwg(2)-zwg(1))
       zvg=zug; zpg=zug;
!
!      (formerly in centercord)
!      coordinates of cell centers; 
!      end-points are on the actual boundaries
!      not sure this is really needed
!
! global z
       zcg=zug
!
       xcf=xv; xc=xv; xc(1)=xu(1); xc(ix2+1)=xu(ix2);
       ycf=yu; yc=yu; yc(1)=yv(1); yc(jy2+1)=yv(jy2);
! local z later
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    metrics
!
       do I=IX1+1,IX2
        AU(I)=1./(XU(I+1)-XU(I-1))*2
        AV(I)=1./(XV(I+1)-XV(I-1))*2
        AW(I)=1./(XW(I+1)-XW(I-1))*2
        APP(I)=1./(XP(I+1)-XP(I-1))*2
        AP(I)=1./(XU(I)-XU(I-1))
       enddo
       AU(IX1)=1./(-3.*XU(IX1)+4.*XU(IX1+1)-XU(IX1+2))*2
       AV(IX1)=1./(-3.*XV(IX1)+4.*XV(IX1+1)-XV(IX1+2))*2
       AW(IX1)=1./(-3.*XW(IX1)+4.*XW(IX1+1)-XW(IX1+2))*2
       APP(IX1)=1./(-3.*XP(IX1)+4.*XP(IX1+1)-XP(IX1+2))*2
       AP(IX1)=1./(XU(IX1+1)-XU(IX1))
       AU(IX2+1)=AU(IX2)
       AP(IX2+1)=AP(IX2)
       AV(IX2+1)=1./(3.*XV(IX2+1)-4.*XV(IX2)+XV(IX2-1))*2
       AW(IX2+1)=1./(3.*XW(IX2+1)-4.*XW(IX2)+XW(IX2-1))*2
       APP(IX2+1)=1./(3.*XP(IX2+1)-4.*XP(IX2)+XP(IX2-1))*2
!
       do J=JY1+1,JY2
        BU(J)=1./(YU(J+1)-YU(J-1))*2
        BV(J)=1./(YV(J+1)-YV(J-1))*2
        BW(J)=1./(YW(J+1)-YW(J-1))*2
        BPP(J)=1./(YP(J+1)-YP(J-1))*2
        BP(J)=1./(YV(J)-YV(J-1))
       enddo
       BU(JY1)=1./(-3.*YU(JY1)+4.*YU(JY1+1)-YU(JY1+2))*2
       BV(JY1)=1./(-3.*YV(JY1)+4.*YV(JY1+1)-YV(JY1+2))*2
       BW(JY1)=1./(-3.*YW(JY1)+4.*YW(JY1+1)-YW(JY1+2))*2
       BPP(JY1)=1./(-3.*YP(JY1)+4.*YP(JY1+1)-YP(JY1+2))*2
       BP(JY1)=1./(YV(JY1+1)-YV(JY1))
       BU(JY2+1)=1./(3.*YU(JY2+1)-4.*YU(JY2)+YU(JY2-1))*2
       BV(JY2+1)=1./(3.*YV(JY2+1)-4.*YV(JY2)+YV(JY2-1))*2
       BW(JY2+1)=1./(3.*YW(JY2+1)-4.*YW(JY2)+YW(JY2-1))*2
       BPP(JY2+1)=1./(3.*YP(JY2+1)-4.*YP(JY2)+YP(JY2-1))*2
!
       do K=2,NZG-1
        CUG(K) = 1./(ZUG(K+1)-ZUG(K-1))*2
        CVG(K) = 1./(ZVG(K+1)-ZVG(K-1))*2
        CWG(K) = 1./(ZWG(K+1)-ZWG(K-1))*2
        CPPG(K)= 1./(ZPG(K+1)-ZPG(K-1))*2
        CPG(K) = 1./(ZWG(K)-ZWG(K-1))
       enddo
       CUG(1)   = 1./(-3.*ZUG(1)+4.*ZUG(2)-ZUG(3))*2
       CVG(1)   = 1./(-3.*ZVG(1)+4.*ZVG(2)-ZVG(3))*2
       CWG(1)   = 1./(-3.*ZWG(1)+4.*ZWG(2)-ZWG(3))*2
!       CWG(1)   = CWG(2)
       CPPG(1)  = 1./(-3.*ZP(1)+4.*ZPG(2)-ZP(3))*2
       CPG(1)   = 1./(ZWG(2)-ZWG(1))
       CUG(NZG) = 1./(3.*ZU(NZG)-4.*ZU(NZG-1)+ZU(NZG-2))*2
       CVG(NZG) = 1./(3.*ZV(NZG)-4.*ZV(NZG-1)+ZV(NZG-2))*2
       CPPG(NZG)= 1./(3.*ZP(NZG)-4.*ZP(NZG-1)+ZP(NZG-2))*2
!
!     SET THE LOCAL Z-METRIC
!
       DO K=1,KZ2+1
         KK=(ME)*(NZG-2)/NCPU+K
         ZW(K)=ZWG(KK)
         CPP(K)=CPPG(KK)
         CP(K)=CPG(KK)
         CU(K)=CUG(KK)
         CV(K)=CVG(KK)
         CW(K)=CWG(KK)
!
         zc(k)=zcg(kk)
         zcf(k)=zcg(kk)
       ENDDO
!
       if (bound(top).gt.0) then
          zc(kz2+1)=zw(kz2)
       else
          zc(kz2+1)=(zw(kz2+1)+zw(kz2))*0.5
       endif
!
!   compute min and max cell size
!
       dxmax=maxval(xu(2:nx)-xu(1:nx-1))
       dxmin=minval(xu(2:nx)-xu(1:nx-1))
!
       dymax=maxval(yv(2:ny)-yv(1:ny-1))
       dymin=minval(yv(2:ny)-yv(1:ny-1))
!
       dzmax=maxval(zwg(2:nzg)-zwg(1:nzg-1))
       dzmin=minval(zwg(2:nzg)-zwg(1:nzg-1))
      
!     
       if (me.eq.0) then 
        write(output,'(2A)')      whoami, '============== GRID INITIALIZED =============='
        write(output,'(2A, 3I5)') whoami, ' Local domain size: ',(ix2+1),(jy2+1),(kz2+1)
        write(output,'(2A, 3I5)') whoami, ' Global domain size: ',nx,ny,nzg
!
        write(output,'(2A,G11.5)') whoami, ' Dx max    =',dxmax
        write(output,'(2A,G11.5)') whoami, ' Dx min    =',dxmin
        write(output,'(2A,G11.5)') whoami, ' Dy max    =',dymax
        write(output,'(2A,G11.5)') whoami, ' Dy min    =',dymin
        write(output,'(2A,G11.5)') whoami, ' Dz max    =',dzmax
        write(output,'(2A,G11.5)') whoami, ' Dz min    =',dzmin 
        write(output,'(2A,G11.5)') whoami, ' dx/dy max =',dxmax/dymin
        write(output,'(2A,G11.5)') whoami, ' dx/dy min =',dxmin/dymax
        write(output,'(2A,G11.5)') whoami, ' dx/dz max =',dxmax/dzmin
        write(output,'(2A,G11.5)') whoami, ' dx/dz min =',dxmin/dzmax
        write(output,'(2A,G11.5)') whoami, ' dy/dz max =',dymax/dzmin
        write(output,'(2A,G11.5)') whoami, ' dy/dz min =',dymin/dzmax 
        write(output,'(2A)')      whoami,  '============================================='
       endif
!
!
       grid_initialized=.true.
!
       end subroutine grid_initialize 
!
       subroutine grid_done
       implicit none
       deallocate(xu,xv,xw,xp,xc,xcf,
     &            yu,yv,yw,yp,yc,ycf,
     &            zu,zv,zw,zp,zc,zcf,
     &            zug,zvg,zwg,zpg,zcg,
! metrics
     &            app,ap,au,av,aw,
     &            bpp,bp,bu,bv,bw,
     &            cpp,cp,cu,cv,cw,
     &            cppg,cpg,cug,cvg,cwg)
!
       grid_initialized=.false.
!
       end subroutine grid_done
!       
