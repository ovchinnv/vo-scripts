      module precursor_output
      
       integer savenum
       real saveloc, savecut
       integer savefreq
       character*80 savename
       parameter (saveloc=65.0)
       parameter (savecut=0.74)
       parameter (savefreq=5)
       data savenum/0/
       data savename/'b2/t3bh.new'/


!
      if (save) call save_setup(cpusave,ksave,xu,yv,zwg,tlevel,t_prev,index)



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
      subroutine save_setup(cpusave,ksave,xu,yv,zwg,tlevel,t_prev,idx)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/save.h'
c
      real xu(nx),yv(ny),zwg(nzg),tlevel,t_prev
      integer ksave, cpusave, k,i,j,isave
      character*3 idx(0:999)
c
      integer nlen
      nlen=index(savename,' ')-1
      if(nlen.le.0) nlen=80
c     
      if (myrank.eq.0) then 
   	 if (lcfl) write(0,*) 'WARNING: saving slices with cfl ON' 
      endif
c
c     find part of domain that gets saved
      isave=1
       do while (xu(isave).lt.savecut)
       isave=isave+1
      enddo 
c
      ksave=1
      cpusave=0 
      k=1
        do while (zwg(k).lt.(saveloc))
        k=k+1
      enddo
      if (k.gt.1) then 
           cpusave=(k-2)/(nz-2)
	   ksave=mod(k-2,nz-2)+1
      endif
      if (cpusave.gt.mysize-1) then 
        cpusave=cpusave-1
         ksave=nz
      endif  
      
      if (myrank.eq.cpusave) then 
        open(59,file=savename(1:nlen)//'.'//idx(savenum),form='binary',status='new')
        write(59) tlevel,nx-isave+1,ny,(xu(i),i=isave,nx),(yv(j),j=1,ny)
        t_prev=tlevel
      endif
      return
      end 
cccccccccccccccccccccccccccccccccccccccccccccc
      subroutine save_slice(uo,vo,wo,cpusave,ksave,xu,yv,tlevel,t_prev,icycle,idx)
      implicit none
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/save.h'
c
      real xu(nx),yv(ny),tlevel,t_prev
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      integer ksave, cpusave, k, i, j, icycle, isave
      character*3 idx(0:999)
      save savenum
c
      integer nlen
      nlen=index(savename,' ')-1
      if(nlen.le.0) nlen=80
c
c          find part of domain that gets saved
      isave=1
       do while (xu(isave).lt.savecut)
       isave=isave+1
      enddo 
c
        if (mod(icycle,savefreq).eq.1) then  
         if (myrank.eq.cpusave) then 
      	    if (icycle-itini.gt.1.and.mod(icycle,ievery*savefreq).eq.1) then
	      savenum=savenum+1 
      	      close(59)
      	      open(59,file=savename(1:nlen)//'.'//idx(savenum),form='binary',status='new')
              write(59) t_prev,nx-isave+1,ny,(xu(i),i=isave,nx),(yv(j),j=1,ny)
	    endif !icycle
            write(59) tlevel-t_prev
	    t_prev=tlevel
            do i=isave+1,nx-1
               write(59) (uo(i,j,ksave), vo(i,j,ksave), wo(i,j,ksave),j=2,ny-1)
            enddo
         endif !myrank
        endif !mod	  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine save_close(cpusave)
      implicit none
      include 'headers/common.h'
      integer cpusave
                    
        if (myrank.eq.cpusave) close(59)
      
      return
      end
       