#ifdef PARALLEL
      SUBROUTINE REFRESHBC(X,IM,KM,MYSIZE,MYRANK,COMMM,QPER)
!     Original author: A. Scotti
!     handles both the periodic/non-periodic versions via QPER flag
!     VO 5/10
       include 'mpif.h'
       int :: IM,KM,MYSIZE,MYRANK,STATUS(MPI_STATUS_SIZE)
       int :: COMM, IERR
       bool :: QPER
!      
       int :: MTYPE,DEST,SOURCE
!      
       float :: X(IM,*)
       MTYPE=MPI_float
!
       if (mysize.gt.1) then 
         DEST=modulo(MYRANK+1,MYSIZE)
         SOURCE=modulo(MYRANK-1,MYSIZE)
         
         if (.not.QPER.and.(myrank.ne.0.and.myrank.ne.mysize-1)) then   
            CALL MPI_SENDRECV(X(1,KM+1),IM,MTYPE,
     &           DEST,10,X(1,1),IM,MTYPE,SOURCE,10,COMM,
     &           STATUS,IERR)
         else 
            if (myrank.eq.mysize-1) then 
               call MPI_RECV(X(1,1),IM,MTYPE,SOURCE,10,COMM,STATUS,IERR) 	  
            else 
               call MPI_SEND(X(1,KM+1),IM,MTYPE,DEST,10,COMM,IERR)
            endif
         endif	  
!
         DEST=SOURCE
         SOURCE=MOD(MYRANK+1,MYSIZE)
!
         if (.not.QPER.and.(myrank.ne.0.and.myrank.ne.mysize-1)) then 
            CALL MPI_SENDRECV(X(1,2),IM,MTYPE,
     &           DEST,11,X(1,KM+2),IM,MTYPE,SOURCE,11,COMM,
     &           STATUS,IERR)
         else
            if (myrank.eq.0) then 
               call MPI_RECV(X(1,KM+2),IM,MTYPE,SOURCE,11,COMM,STATUS,IERR) 	  
            else 
               call MPI_SEND(X(1,2),IM,MTYPE,DEST,11,COMM,IERR)
            endif
         endif	  
       endif
       RETURN
      END SUBROUTINE REFRESHBC
#else
! nothing : this subroutine should not be called in serial mode
#endif
