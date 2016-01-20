      module ib_pressure
!     compute and output pressure at immersed boundary points
       integer presfreq, nsurf
       character*80 presname
       parameter (presfreq=20)
       parameter (nsurf=20)
       data presname/'surf_pres.dat'/












Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      SUBROUTINE CP
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      subroutine pressure_setup(nb,sb,xb,yb,coefx,coefy,x,zl,zg,npresl,ipresl,kpresl,spxl,spzl, me)
      subroutine pressure_setup(nb,sb,xb,yb,coefx,coefy,x,zl,zg,npres,npresl,ipresl,kpresl,npres_array,npres_disp,me,comm)
      implicit none
c      
      include 'headers/immbd.h'
      include 'headers/dimension.h'
      include 'mpif.h'
      include 'headers/spress.h'
c
      real    sb(0:nnmax),xb(0:nnmax),yb(0:nnmax) ! body point coordinates
      real    coefx(nnmax,3),coefy(nnmax,3) ! polynomial coefficients
      real    x(nx), zg(nzg), zl(nz) 
      real    tlr
      real    cofx(3),cofy(3)
      real    sa, sn, xx, zz, xn, yn, dist
      integer nb, rad, indx, indz, i, k, n, m, inout, cpu, me
      integer*1 pres_out(nx,nzg)
      real sp_x(nx,nzg), sp_z(nx,nzg)
      integer npres, ipres(nimb), kpres(nimb), npresl, ipresl(nimb), kpresl(nimb),
     &        npres_array(ncpu), npres_disp(ncpu)
      real spx(nimb), spz(nimb), spxl(nimb), spzl(nimb)
      real     zpres(nimb), xpres(nimb)
      real     zpresl(nimb), xpresl(nimb)
c
      integer cell1(2,nimb), cell2(2,nimb), cell3(2,nimb), cell4(2,nimb), cell5(2,nimb)! cell1 -- closest, cell5 -- farthest
      real dist1(nimb), dist2(nimb), dist3(nimb), dist4(nimb), dist5(nimb)
      real sp1(2,nimb), sp2(2,nimb), sp3(2,nimb), sp4(2,nimb), sp5(2,nimb)
      real deltax, deltaz, vx1, vx2, vy1, vy2, mag1, mag2, dot
      integer ib(nnmax), kb(nnmax)
      
      integer type, ierr, comm ! mpi variables
      integer nlen
      
c
#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
c
      nlen=index(presname,' ')-1
      if(nlen.le.0) nlen=80
c
      do n=1,nb
       dist1(n)=999999 ! put infinite distance initially
       dist2(n)=999999
       dist3(n)=999999
       dist4(n)=999999
       dist5(n)=999999
      enddo

      rad=5 ! search distance in cells
c
      deltax=0.0
      deltaz=0.0
c-----------------------------------------------------------------------
      do n=1,nb
       call locate(x,nx,xb(n),indx) !xb wall normal
       call locate(zg,nzg,yb(n),indz) !yb streamwise
c  
       ib(n)=indx ! save indices
       kb(n)=indz
c       
        if (ib(n).ne.0.and.ib(n).lt.nx-1.and.kb(n).ne.0.and.kb(n).lt.nzg-1) then 
         deltax=deltax+(x(ib(n)+1)-x(ib(n)))
         deltaz=deltaz+(zg(kb(n)+1)-zg(kb(n)))
        endif	
      enddo
c
       deltax=deltax/float(nb)
       deltaz=deltaz/float(nb)
       tlr=sqrt(deltax**2.+deltaz**2.)/2.
c        write(0,*) 'deltas ',deltax,deltaz, tlr
c
      do n=1,nb ! loop over ib points
c  
       indx=ib(n)
       indz=kb(n)
c       
       do i=-rad,rad,1
        do k=-rad,rad,1
         if ((indx+i).gt.0.and.(indz+k).gt.0.and.(indx+i).le.nx.and.(indz+k).le.nzg) then 
c    check if the point is inside or outside
          xx=x(indx+i)
          zz=zg(indz+k)
          sa=sb(n) ! the (arclength) coordinate
c
c     make local copies of polynomial coeffs	   
          do m=1,3
           cofx(m)=coefx(n,m)
           cofy(m)=coefy(n,m)
          enddo
c
          call Normal(cofx,cofy,xx,zz,sa,sn)
          xn=cofx(1)*sn**2.+cofx(2)*sn+cofx(3)
          yn=cofy(1)*sn**2.+cofy(2)*sn+cofy(3)
c
          call LocPnt(cofx,cofy,xx,zz,xn,yn,sn,tlr,inout)
c     NOTE: inout =  9 means inside the body, (although effectively on the boundary)
c                    1 means inside the fluid
c                   -1 means inside the body      
          if (inout.eq.1) then ! inside the fluid
c     dot product weighting function for distance
           vx1=xx-xn ! normal vector at xx
	   vy1=zz-yn	  
	   mag1=sqrt(vx1**2+vy1**2) ! its magnitude
	   vx1=vx1/mag1 !normalize to unit length
	   vy1=vy1/mag1 
c     now vector to the marker point
           vx2=xx-xb(n)
	   vy2=zz-yb(n)	  
	   mag2=sqrt(vx2**2+vy2**2) ! magnitude
	   vx2=vx2/mag2 !normalize to unit length
	   vy2=vy2/mag2 
c     dot product
           dot=abs(vx1*vx2+vy1*vy2)	   
	  
	   dist=sqrt( (xb(n)-xx)**2 + (yb(n)-zz)**2 ) ! distance to the boundary marker
           dist=dist/dot ! weight points that lie close to the normal more
c     find 5 nearest-neighbor cells
           if (dist1(n).gt.dist) then
            cell5(1,n)=cell4(1,n) ! i-indx
            cell5(2,n)=cell4(2,n) ! k-indx
            sp5(1,n)=sp4(1,n) ! x-coord of surface point along normal 
            sp5(2,n)=sp4(2,n) ! z-coord of surface point along normal
              dist5(n)=dist4(n)   ! distance to the marker point
            cell4(1,n)=cell3(1,n)
            cell4(2,n)=cell3(2,n)
            sp4(1,n)=sp3(1,n)  
            sp4(2,n)=sp3(2,n) 
              dist4(n)=dist3(n)
            cell3(1,n)=cell2(1,n)
            cell3(2,n)=cell2(2,n)
            sp3(1,n)=sp2(1,n)  
            sp3(2,n)=sp2(2,n) 
              dist3(n)=dist2(n)
            cell2(1,n)=cell1(1,n)
            cell2(2,n)=cell1(2,n)
            sp2(1,n)=sp1(1,n)  
            sp2(2,n)=sp1(2,n) 
              dist2(n)=dist1(n)
            cell1(1,n)=indx+i
            cell1(2,n)=indz+k
            sp1(1,n)=xn
            sp1(2,n)=yn
              dist1(n)=dist   
           elseif (dist2(n).gt.dist) then
            cell5(1,n)=cell4(1,n) ! i-indx
            cell5(2,n)=cell4(2,n) ! k-indx
            sp5(1,n)=sp4(1,n) ! x-coord of surface point along normal 
            sp5(2,n)=sp4(2,n) ! z-coord of surface point along normal
              dist5(n)=dist4(n)   ! distance to the marker point
            cell4(1,n)=cell3(1,n)
            cell4(2,n)=cell3(2,n)
            sp4(1,n)=sp3(1,n)  
            sp4(2,n)=sp3(2,n) 
              dist4(n)=dist3(n)
            cell3(1,n)=cell2(1,n)
            cell3(2,n)=cell2(2,n)
            sp3(1,n)=sp2(1,n)  
            sp3(2,n)=sp2(2,n) 
              dist3(n)=dist2(n)
            cell2(1,n)=indx+i
            cell2(2,n)=indz+k
            sp2(1,n)=xn
            sp2(2,n)=yn
              dist2(n)=dist
           elseif (dist3(n).gt.dist) then
            cell5(1,n)=cell4(1,n)
            cell5(2,n)=cell4(2,n)
            sp5(1,n)=sp4(1,n) ! x-coord of surface point along normal 
            sp5(2,n)=sp4(2,n) ! z-coord of surface point along normal
              dist5(n)=dist4(n)
            cell4(1,n)=cell3(1,n)
            cell4(2,n)=cell3(2,n)
            sp4(1,n)=sp3(1,n)  
            sp4(2,n)=sp3(2,n) 
              dist4(n)=dist3(n)
            cell3(1,n)=indx+i
            cell3(2,n)=indz+k
            sp3(1,n)=xn
            sp3(2,n)=yn
              dist3(n)=dist
           elseif (dist4(n).gt.dist) then
            cell5(1,n)=cell4(1,n)
            cell5(2,n)=cell4(2,n)
            sp5(1,n)=sp4(1,n) ! x-coord of surface point along normal 
            sp5(2,n)=sp4(2,n) ! z-coord of surface point along normal
              dist5(n)=dist4(n)
            cell4(1,n)=indx+i
            cell4(2,n)=indz+k
              dist4(n)=dist
           elseif (dist5(n).gt.dist) then
            cell5(1,n)=indx+i
            cell5(2,n)=indz+k
            sp5(1,n)=xn
            sp5(2,n)=yn
               dist5(n)=dist
           endif     
          endif ! inout
         endif !indx
	enddo ! rad
       enddo ! rad
      enddo !nb
ccccccccccccccccccccccc
c     AT THIS POINT WE HAVE 5 NEAREST NEIGHBORS FOR EACH MARKER POINT
C     NEXT, WE FILL A 2D ARRAY WITH 1 CORRESPONDING TO THE POINTS WE WANT TO SAVE (O ELSEWHERE)      
c     initialize to 0
      do i=1,nx
       do k=1,nzg
        pres_out(i,k)=0
        sp_x(i,k)=0
        sp_z(i,k)=0
       enddo
      enddo        	
      do n=1,nb
       pres_out(cell1(1,n), cell1(2,n))=1
       pres_out(cell2(1,n), cell2(2,n))=1
       pres_out(cell3(1,n), cell3(2,n))=1
       pres_out(cell4(1,n), cell4(2,n))=1
       pres_out(cell5(1,n), cell5(2,n))=1
c
       sp_x(cell1(1,n), cell1(2,n))=sp1(1,n)
       sp_x(cell2(1,n), cell2(2,n))=sp2(1,n)
       sp_x(cell3(1,n), cell3(2,n))=sp3(1,n)
       sp_x(cell4(1,n), cell4(2,n))=sp4(1,n)
       sp_x(cell5(1,n), cell5(2,n))=sp5(1,n)
c
       sp_z(cell1(1,n), cell1(2,n))=sp1(2,n)
       sp_z(cell2(1,n), cell2(2,n))=sp2(2,n)
       sp_z(cell3(1,n), cell3(2,n))=sp3(2,n)
       sp_z(cell4(1,n), cell4(2,n))=sp4(2,n)
       sp_z(cell5(1,n), cell5(2,n))=sp5(2,n)
      enddo
c     now fill the array      
      npres=0.
      do i=1,nx
       do k=1,nzg
        if (pres_out(i,k).eq.1) then 
	 npres=npres+1
         ipres(npres)=i
	 kpres(npres)=k
         spx(npres)=sp_x(i,k)
         spz(npres)=sp_z(i,k)
	endif  
       enddo
      enddo        	
cccccccccccccccccccc NOW DISTRUBITE POINTS BETWEEN CPUS ccccccccccccccc      
c     find local indices
      npresl=0
      do i=1,npres
        cpu=(kpres(i)-2)/(nz-2)
        if (cpu.eq.me) then 
	 npresl=npresl+1
 	 ipresl(npresl)=ipres(i)
	 kpresl(npresl)=kpres(i)-cpu*(nz-2)
         spxl(npresl)=spx(i)
         spzl(npresl)=spz(i)
        endif
      enddo

ccccc open surface pressure file & write grid
      call mpi_reduce(npresl, npres, 1, MPI_INTEGER, MPI_SUM, ncpu-1, comm, ierr)
      call mpi_gather(npresl, 1, MPI_INTEGER, npres_array, 1, MPI_INTEGER, ncpu-1, comm, ierr)
c
c     form local arrays
      do i=1,npresl
       xpresl(i)=x(ipresl(i))
       zpresl(i)=zl(kpresl(i))
c       presl(i)=p(ipresl(i), ny/2+1, kpresl(i)) ! pressure
      enddo
c     collect data              
      npres_disp(1)=0
      do i=2,ncpu
       npres_disp(i)=npres_array(i-1)
      enddo 
c      write(0,*) npresl, npres, npres_array, npres_disp
      call mpi_gatherv(zpresl, npresl, type, zpres, npres_array, npres_disp, type, ncpu-1, comm, ierr)
      call mpi_gatherv(xpresl, npresl, type, xpres, npres_array, npres_disp, type, ncpu-1, comm, ierr)
      call mpi_gatherv(spxl, npresl, type, spx, npres_array, npres_disp, type, ncpu-1, comm, ierr)
      call mpi_gatherv(spzl, npresl, type, spz, npres_array, npres_disp, type, ncpu-1, comm, ierr)
c      call mpi_gatherv(presl, npresl, type, pres, npres_array, npres_disp, type, ncpu-1, comm, ierr)
c
      if (me.eq.ncpu-1) then 
        open(700,file=presname(1:nlen),form='unformatted', status='unknown')
	write(700) ny-2, npres, (xpres(i), zpres(i), spx(i), spz(i), i=1,npres)
      endif		
c
      return
      end
C 
c
      subroutine surface_pressure(p, tlevel, npres, npresl, ipresl, kpresl, npres_array, npres_disp, me, comm)
      implicit none
      include 'headers/dimension.h'
      include 'mpif.h'
      include 'headers/immbd.h'
      
      integer comm, npresl, npres_array(ncpu), npres_disp(ncpu), npres, type, ierr, i, j, n
      integer npres_array2(ncpu), npres_disp2(ncpu)
      real pres(ny-2, nimb), presl(ny-2, nimb), tlevel, p(nx,ny,nz) 
      integer ipresl(nimb), kpresl(nimb), me
c
#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
c
c      create arrays for 2D (spanwise) data
      do i=1,ncpu
        npres_array2(i)=npres_array(i)*(ny-2)
        npres_disp2(i) =npres_disp(i)*(ny-2)
      enddo
c
c     form local pressure array
      do n=1,npresl
       do j=2,ny-1
        presl(j-1,n)=p(ipresl(n), j, kpresl(n)) ! pressure
       enddo
      enddo
c     collect data              
c
      call mpi_gatherv(presl, npresl*(ny-2), type, pres, npres_array2, npres_disp2, type, ncpu-1, comm, ierr)

c      write(0,*) npresl, me

c
      if (me.eq.ncpu-1) then 
	write(700) tlevel, ((pres(j,i), i=1,npres),j=1,ny-2)
      endif		
c
      return
      end
C 
      subroutine pressure_close(me)
      implicit none
      include 'headers/dimension.h'
      integer me
       if (me.eq.ncpu-1) close(700)
      return
      end
       