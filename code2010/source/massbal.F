!
!-----Formerly SUBROUTINE-MasBlas-----------------------E. Balaras--17/01999--
!-----V.Ovchinnikov 5/10--
!
      SUBROUTINE correct_mass(ust,vst,wst,xf,yf,zf,im,jm,km,where)
!
!       Balances the mass by adjusting the normal component
!       of the velocity at the boundary specified in 'where'
!
!-----------------------------------------------------------------------
!
       use size, only: me, communicator, ncpu
       use boundary, only: left, right, top, bottom, back, front
       use fout, only: output
!
       implicit none
!
#ifdef PARALLEL
       include 'mpif.h' 
#endif
!-----------------------------------------------------------------------
!
       float, dimension(im,jm,km) :: ust, vst, wst
       float :: xf(im),yf(jm),zf(km) ! grids
!
       float, dimension(6) :: mold=0.0, moldg 
     &          alpha, dA
       int :: ierr, where
! 
!    compute global mass balance
!    boundary itype(1/2)
       do k=2,kz2
        do j=2,jy2
         dA=(zf(k)-zf(k-1))*(yf(j)-yf(j-1))
         mold(left ) = mold(left ) - ust(1,j,k  )*dA
         mold(right) = mold(right) + ust(ix2,j,k)*dA
        enddo 
       enddo
!     boundary itype(3/4)
       do k=2,kz2
         do i=2,ix2
           dA=(zf(k)-zf(k-1))*(xf(i)-xf(i-1))
           mold(front) = mold(front) - vst(i,1,k  )*dA
           mold(back ) = mold(back ) + vst(i,jy2,k)*dA
         enddo 
       enddo
!..   boundary itype(5/6)
       do j=2,jy2
         do i=2,ix2
           dA=(yf(j)-yf(j-1))*(xf(i)-xf(i-1))
           mold(bottom) = mold(mbottom) - wst(i,j,1)*dA
           mold(top) = mold(top) + wst(i,j,kz2)*dA
         enddo 
       enddo
       if (me.gt.0)      mold(bottom)=0.0
       if (me.ne.ncpu-1) mold(top)=0.0
#IF PARALLEL       
       call MPI_ALLREDUCE(mold,moldg,6,MPI_float,MPI_SUM,communicator,ierr)
#ENDIF
!...correction factor to be applied at boundary 'where'
       alpha   = (-sum(moldg)+moldg(where))/moldg(where)
!
       select case(where)
        case(right);
         ust(ix2,2:jy2,2:kz2)=ust(ix2,2:jy2,2:kz2)*alpha
        case(left);
         ust(1,2:jy2,2:kz2) = ust(1,2:jy2,2:kz2)*alpha
        case(front);
         vst(2:ix2,1,2:kz2) = vst(2:ix2,1,2:kz2)*alpha
        case(back);
         vst(2:ix2,jy2,2:kz2) = vst(2:ix2,jy2,2:kz2)*alpha
        case(top);
         if (me.eq.ncpu-1) wst(2:ix2,2:jy2,kz2)= wst(2:ix2,2:jy2,kz2)*alpha
        case(bottom);
         if (me.le.0) wst(2:ix2,2:jy2,1)= wst(2:ix2,2:jy2,1)*alpha
        case(default)
!       nothing
       end select
!     


HERE

      if (mod(icycle,infnum).eq.0) then 
        if (myrank.eq.0) then
          write(6,1234)  moutletg ,minletg, mtopg, mbottomg,alpha
        endif
      endif	
 1234 format ('moutlet,minlet, mtop, mbottom, alpha ', 5 G23.17)     
!
      RETURN
      END
