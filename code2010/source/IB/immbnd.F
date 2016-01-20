      module immbnd
      private
#ifdef PARALLEL
       include 'mpif.h'
#endif
!
	integer, public, save :: nimb=0, nnmax=0
!	parameter (nimb=4000)! number of points to describe the surface boundary points
!        parameter (nnmax=3*nimb)
	
	
	
	
!     IMMERSED BOUNDARY DATA
      real sb(0:nnmax),yb(0:nnmax),xb(0:nnmax)
      real coefx(nnmax,3),coefy(nnmax,3)
      real uim(nnmax),vim(nnmax),wim(nnmax),tvim(nnmax),Ilmim(nnmax)
      integer nb

!     surface pressure computation
      include 'headers/spress.h'      
      integer npres ! for surface pressure computation
      integer npresl, ! for surface pressure computation (local)
     &        npres_array(ncpu), npres_disp(ncpu),
     &        ipresl(nimb), kpresl(nimb)
!      real     zpres(nimb), xpres(nimb), pres(nimb)              
!      real     zpresl(nimb), xpresl(nimb), presl(nimb)
!      real spx(nimb), spxl(nimb), spz(nimb), spzl(nimb)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef #PARALLEL
!cccc VARIABLES FOR THE PARALLELIZATION OF IMMERSED BOUNDARY cccccccccccccccc
! U
      int, allocatable ::     
     &        uniml(0:ncpu-1),uiiml(0:ncpu-1,nnmax),ujiml(0:ncpu-1,nnmax),untop(0:ncpu-1),utk(0:ncpu-1,nimb),
     &        ustk1(0:ncpu-1,nnmax),ustk2(0:ncpu-1,nnmax),unbot(0:ncpu-1),ubk(0:ncpu-1,nimb),usti(0:ncpu-1,nnmax)
      float, allocatable ::   
     &        utib_ind(0:ncpu-1,nimb),utib_c1(0:ncpu-1,nimb),utib_c2(0:ncpu-1,nimb),ubib_ind(0:ncpu-1,nimb),
     &        ubib_c1(0:ncpu-1,nimb),ubib_c2(0:ncpu-1,nimb),uc4intl(0:ncpu-1,nnmax,5)
      int, allocatable ::
              uibod((nz-2)*(nx-2)),ujbod((nz-2)*(nx-2))
      int     unbod!(0:ncpu-1)
! V
      int     vniml(0:ncpu-1),viiml(0:ncpu-1,nnmax),vjiml(0:ncpu-1,nnmax),vntop(0:ncpu-1),vtk(0:ncpu-1,nimb),
     &        vstk1(0:ncpu-1,nnmax),vstk2(0:ncpu-1,nnmax),vnbot(0:ncpu-1),vbk(0:ncpu-1,nimb),vsti(0:ncpu-1,nnmax)
      float   vtib_ind(0:ncpu-1,nimb),vtib_c1(0:ncpu-1,nimb),vtib_c2(0:ncpu-1,nimb),vbib_ind(0:ncpu-1,nimb),
     &        vbib_c1(0:ncpu-1,nimb),vbib_c2(0:ncpu-1,nimb),vc4intl(0:ncpu-1,nnmax,5)
      int     vibod((nz-2)*(nx-2)),vjbod((nz-2)*(nx-2))
      int     vnbod!(0:ncpu-1)
! W
      int     wniml(0:ncpu-1),wiiml(0:ncpu-1,nnmax),wjiml(0:ncpu-1,nnmax),wntop(0:ncpu-1),wtk(0:ncpu-1,nimb),
     &        wstk1(0:ncpu-1,nnmax),wstk2(0:ncpu-1,nnmax),wnbot(0:ncpu-1),wbk(0:ncpu-1,nimb),wsti(0:ncpu-1,nnmax)
      float   wtib_ind(0:ncpu-1,nimb),wtib_c1(0:ncpu-1,nimb),wtib_c2(0:ncpu-1,nimb),wbib_ind(0:ncpu-1,nimb),
     &        wbib_c1(0:ncpu-1,nimb),wbib_c2(0:ncpu-1,nimb),wc4intl(0:ncpu-1,nnmax,5)
      int     wibod((nz-2)*(nx-2)),wjbod((nz-2)*(nx-2))
      int     wnbod!(0:ncpu-1)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccc  VARIABLES FOR PRESSURE RECOVERY
      int     ptrashn, ptrashi(0:ncpu-1,nnmax), ptrashj(0:ncpu-1,nnmax)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif












#ifdef immbnd 
      if (ibswitch) then
!.... obstacle location on the grid
               CALL ObjectBluntNoseHalf(sb,xb,yb,coefx,coefy,nb,nnmax)
!.... locate boundary points and interpolation coeff. for u
               CALL bpoints(xu,zc,zcg,sb,xb,yb,nb,coefx,coefy,
     &                   uniml,uiiml,ujiml,untop,utk,ustk1,ustk2,usti,unbot,ubk,utib_ind,utib_c1,utib_c2,ubib_ind,
     &                   ubib_c1,ubib_c2,
     &                   uc4intl,unbod,uibod,ujbod,'u.dat')
!.... locate boundary points and interpolation coeff. for v
               CALL bpoints(xc,zc,zcg,sb,xb,yb,nb,coefx,coefy,
     &                   vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vsti,vnbot,vbk,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                   vbib_c1,vbib_c2,
     &                   vc4intl,vnbod,vibod,vjbod,'p.dat')
!.... locate boundary points and interpolation coeff. for w
               CALL bpoints(xc,zw,zwg,sb,xb,yb,nb,coefx,coefy,
     &                   wniml,wiiml,wjiml,wntop,wtk,wstk1,wstk2,wsti,wnbot,wbk,wtib_ind,wtib_c1,wtib_c2,wbib_ind,
     &                   wbib_c1,wbib_c2,
     &                   wc4intl,wnbod,wibod,wjbod,'w.dat')
!ccccccccccc experimental cccccccccccc surface pressure computation cccccccccccccccccc

      if (pres_out) then
        call pressure_setup(nb,sb,xb,yb,coefx,coefy,xc,zc,zcg,npres,npresl,ipresl,kpresl,npres_array,npres_disp,myrank,
     &	                    mpi_comm_hydro)
      endif		

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       for pressure recovery -- experimental

!      me=myrank
!      do i=1,vnbod
!        write(300+myrank,'(3 I)') i,vibod(i),vjbod(i)
!      enddo
!      close (300+myrank)
!cccccccccccccccccccccccccccc EXPERIMENTAL AND NON_ESSENTIAL ccccccccccccccccccccccccccc
      me=myrank
      k=0               !unphysical pressure pts
      do i=1,uniml(me)
       j=1               !v counter
       found=.false.
       do while (j.le.vniml(me).and..not.found)
        if (uiiml(me,i).eq.viiml(me,j)) then 
	 if (ujiml(me,i).eq.vjiml(me,j)) then 
	  found=.true.
	 endif
	endif
	j=j+1
       enddo
       if (found) then 
        j=1               !w counter
        found=.false.
        do while (j.le.wniml(me).and..not.found)
         if (uiiml(me,i).eq.wiiml(me,j)) then 
	  if (ujiml(me,i).eq.wjiml(me,j)) then
	   found=.true.
          endif
	 endif
	 j=j+1
        enddo
       endif
       if (found) then
        k=k+1
	ptrashi(me,k)=uiiml(me,i)
	ptrashj(me,k)=ujiml(me,i)
       endif
      enddo 	
      ptrashn=k
!     
!      do i=1, ptrashn
!        write(me+250,'(3 I)') i,ptrashi(me,i),ptrashj(me,i)  
!      enddo
!        close (me+250)
!
!      do i=1, uniml(me)
!        write(me+300,'(7 I)') i,uiiml(me,i),ujiml(me,i),viiml(me,i),vjiml(me,i),wiiml(me,i),wjiml(me,i)  
!      enddo
!        close (me+300)
!ccccccccccccccccccccc      for pressure recovery experimental	
     
!
!      debug
!      write(myrank+200,'(I,6G20.10)') (i,(uc4intl(myrank,i,k),k=1,5),uc4intl(myrank,i,1)+uc4intl(myrank,i,2)+uc4intl(myrank
!     &     ,i,3)+uc4intl(myrank,i,4)+uc4intl(myrank,i,5)
!     &     ,i=1,uniml(myrank))
!      close(myrank+200)
!
!
!      write(myrank+100,'(2I3,7G20.10)') (utib_ind(myrank,i),utk(myrank,i),(uc4intl(myrank,utib_ind(myrank,i),k),k=1,5),
!     & utib_c1(myrank,i),utib_c2(myrank,i),i=1,untop(myrank))
!      write(myrank+100,*) '********************************************************'
!      write(myrank+100,'(2I3,7G20.10)') (ubib_ind(myrank,i),ubk(myrank,i),(uc4intl(myrank,ubib_ind(myrank,i),k),k=1,5),
!     &  ubib_c1(myrank,i),ubib_c2(myrank,i),i=1,unbot(myrank))
!      close(myrank+100)
!      return
!      debug
      
 9875  continue                  
!...  number of points to be forced
               write(0,'("ib points: (me)",I2,": u-grid=",2 I5)') myrank,uniml(myrank)
               write(0,'("ib points: (me)",I2,": v-grid=",2 I5)') myrank,vniml(myrank)
               write(0,'("ib points: (me)",I2,": w-grid=",2 I5)') myrank,wniml(myrank)
!...  set up immersed boundary velocities
               do i=1,uniml(myrank)
                uim(i)=0.0
               enddo
               do i=1,vniml(myrank)
                vim(i)=0.0
#ifdef LES
                tvim(i)=0.0
#endif		
               enddo
               do i=1,wniml(myrank)
                wim(i)=0.0
               enddo
      endif
!      if (ifield.ne.2) then 
!       me=myrank
!       do j=jy1,jy2+1
!        do i=1,vnbod
!         p (vibod(i),j,vjbod(i))=1.
!         vo(vibod(i),j,vjbod(i))=0.
!        enddo
!        do i=1,unbod
!         uo(uibod(i),j,ujbod(i))=0.
!        enddo
!        do i=1,wnbod
!         wo(wibod(i),j,wjbod(i))=0.
!        enddo
!       enddo
!      endif 
#endif      





