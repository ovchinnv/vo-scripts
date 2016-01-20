C
C
C-----SUBROUTINE-BOUNDARY-----------------------E. Balaras--23/05/1998----
C
      SUBROUTINE BOUNDARY(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    u2d,v2d,w2d,t,t_old,ksub,nt,
     &                    uavneu,wavneu,
     &                    xf,yf,zf,imax,jmax,kmax,
     &                    Ltstep,uinf,duinf,dtm,xref)
C
C-----------------------------------------------------------------------
C
C     PURPOSE:    - calls all corresponding routines
C                   variables are defined at lower level routines
C                              
C-----------------------------------------------------------------------
C 
      include 'headers/common.h'
C 
      integer imax,jmax,kmax,nt,ksub
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),vec3(imax,jmax,kmax),
     & vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),vec33(imax,jmax,kmax),
     & u2d(nt,jmax,kmax),v2d(nt,jmax,kmax),w2d(nt,jmax,kmax),t,t_old,
     & uavneu(imax,kmax),wavneu(imax,kmax),xf(imax),yf(jmax),zf(kmax)
      real uinf(imax),duinf(imax),dtm,xref
      logical Ltstep

C-----------------------------------------------------------------------
C                                  select appropriate boundary condition
C-----------------------------------------------------------------------
C
      DO 10 ibnd=1,6
C
C-----------------------------------------------------------------------
C                                                          non-slip wall
C-----------------------------------------------------------------------
C
        IF( itype(ibnd).eq.50 ) THEN
           CALL BOUND050(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
C
C-----------------------------------------------------------------------
C                                           part slip, part no-slip wall
C-----------------------------------------------------------------------
C
        ELSEIF( itype(ibnd).eq.55 ) THEN
           CALL BOUND055(xf,vec1,vec2,vec3,imax,jmax,kmax,ibnd)
        ELSEIF( itype(ibnd).eq.550 ) THEN
           CALL BOUND0550(zf,vec1,vec2,vec3,imax,jmax,kmax,ibnd)
C
C-----------------------------------------------------------------------
C                                                            slip walls
C-----------------------------------------------------------------------
C
        ELSEIF( itype(ibnd).eq.60 ) THEN
           CALL BOUND060(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
C
        ELSEIF( itype(ibnd).eq.62 ) THEN
           CALL BOUND062(vec1,vec2,vec3,vec11,vec22,vec33,imax,jmax,kmax,ibnd,Ltstep)
C-----------------------------------------------------------------------
C                                              laminar INFLOW conditions 
C-----------------------------------------------------------------------
C
c.....Backward facing step
        ELSE IF( itype(ibnd).eq.800 ) THEN
c          CALL BOUND050(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
          CALL BOUND800(vec1,vec2,vec3,imax,jmax,kmax,zf,ibnd)
        ELSE IF( itype(ibnd).eq.8000 ) THEN
          CALL BOUND8000(vec1,vec2,vec3,imax,jmax,kmax,xf)
c.....Blasius bl
        ELSE IF( itype(ibnd).eq.801 ) THEN
          CALL BOUND801(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf)
        ELSE IF( itype(ibnd).eq.8010 ) THEN
          CALL BOUND8010(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf,xref)
        ELSE IF( itype(ibnd).eq.8030 ) THEN
          CALL BOUND8030(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf,xref)
c.....Developing channel flow
        ELSE IF( itype(ibnd).eq.802 ) THEN
          CALL BOUND802(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,zf)
        ELSE IF( itype(ibnd).eq.8020 ) THEN
          CALL BOUND8020(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax)
c.....Developed channel flow
        ELSE IF( itype(ibnd).eq.803 ) THEN
          CALL BOUND803(vec1,vec2,vec3,imax,jmax,kmax,zf,ibnd)
c.....Developed channel flow
        ELSE IF( itype(ibnd).eq.804 ) THEN
          CALL BOUND804(vec1,vec2,vec3,imax,jmax,kmax,zf,ibnd)

c.....Blasius bl + random noise
        ELSE IF( itype(ibnd).eq.805 ) THEN
          CALL BOUND805(vec1,vec2,vec3,vec11,vec22,vec33,
     &                  Ltstep,imax,jmax,kmax,xf,zf,ibnd )
        ELSE IF( itype(ibnd).eq.806 ) THEN
          CALL BOUND806(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf)
C
C-----------------------------------------------------------------------
C                                            turbulent INFLOW conditions 
C-----------------------------------------------------------------------
C
c....impose at inflow profiles from an extermal simulation
c    Carefull! ...set-up at input file is required...
        ELSE IF(ITYPE(IBND) .EQ. 650) THEN
            CALL BOUND650(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    u2d,v2d,w2d,t,t_old,ksub,
     &                    nt,imax,jmax,kmax,ibnd,Ltstep)
c....inflow bc for tbl. A rescaling procedure from a downstream
c     station is applied
        ELSE IF(ITYPE(IBND) .EQ. 660) THEN
            CALL BOUND660(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,uavneu,wavneu,
     &                    imax,jmax,kmax,ibnd,Ltstep,uinf,duinf)
C
C-----------------------------------------------------------------------
C                                                     OUTFLOW conditions 
C-----------------------------------------------------------------------
C
c.... Homogeneous Newmman (du/dx=0, dv/dx=0, dw/dx=0)
        ELSE IF( itype(ibnd).eq.710 ) THEN
          CALL BOUND710(vec1,vec2,vec3,vec11,vec22,vec33,
     &                      xf,yf,zf,imax,jmax,kmax,Ltstep)
c....extrapolation from the interior
        ELSE IF( itype(ibnd).eq.720 ) THEN
          CALL BOUND720(vec1,vec2,vec3,vec11,vec22,vec33,
     &                   xf,yf,zf,imax,jmax,kmax,Ltstep,dtm)
c....convective condition (with mass balance)
        ELSE IF( itype(ibnd).eq.730 ) THEN
          CALL BOUND730(vec1,vec2,vec3,vec11,vec22,vec33,
     &                   xf,yf,zf,t,imax,jmax,kmax,Ltstep,dtm)
c....convective condition (constant speed)
        ELSE IF( itype(ibnd).eq.740 ) THEN
          CALL BOUND740(vec1,vec2,vec3,vec11,vec22,vec33,
     &                   xf,yf,zf,t,imax,jmax,kmax,Ltstep,dtm)
c....convective condition (constant speed)
        ELSE IF( itype(ibnd).eq.7300 ) THEN
          CALL BOUND7300(vec1,vec2,vec3,vec11,vec22,vec33,
     &                   xf,yf,zf,t,imax,jmax,kmax,Ltstep,dtm)

C
C-----------------------------------------------------------------------
C                                               OPEN boundary conditions 
C-----------------------------------------------------------------------
c
c....Freestream boundary condition for laminar bl with dpdx=0
        ELSE IF(ITYPE(IBND) .EQ. 240) THEN
            CALL BOUND240(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
        ELSE IF(ITYPE(IBND) .EQ. 2400) THEN
            CALL BOUND2400(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep,dtm,xref    )
c     strictly for ELLIPTICAL LEADING EDGE FLOW -- BLOCKAGE CORRECTION
        ELSE IF(ITYPE(IBND) .EQ. 2401) THEN
            CALL BOUND2401(vec1,vec2,vec3,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd    )

c....Freestream boundary condition for tbl with dpdx=0
c    (used in the fringe computations) GRESHO and SANI approach
        ELSE IF(ITYPE(IBND) .EQ. 250) THEN
            CALL BOUND250(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )

c....Freestream boundary condition for turbulent bl with dp/dx=0
c    Normal component from cont.       
        ELSE IF(ITYPE(IBND) .EQ. 255) THEN
            CALL BOUND255(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )

c....Freestream boundary condition for tbl with dp/dx=0
c    Lund et. al. approach  (linear regrasion on running averages)
        ELSE IF(ITYPE(IBND) .EQ. 260) THEN
            CALL BOUND260(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,uavneu,
     &                    imax,jmax,kmax,ibnd,Ltstep        )

c....Freestream boundary condition for tbl with dp/dx=0
c    w computed from power laws (Dirichlet)
        ELSE IF(ITYPE(IBND) .EQ. 265) THEN
            CALL BOUND265(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )


c....Freestream boundary condition for turbulent bl with dpdx
c    Lund et. al. approach
        ELSE IF(ITYPE(IBND) .EQ. 270) THEN
            CALL BOUND270(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,uavneu,
     &                    imax,jmax,kmax,ibnd,Ltstep,uinf,duinf)


c....Freestream boundary condition for turbulent bl with dpdx
c    Lund et. al. approach - Pressure gradient is introduced
C    gradualy
        ELSE IF(ITYPE(IBND) .EQ. 275) THEN
            CALL BOUND275(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,uavneu,
     &                    imax,jmax,kmax,ibnd,Ltstep,uinf,duinf,t)


c....Freestream boundary condition for wake-boundary-layer interaction
c    GRESHO and SANI approach
c    Correction is applied by MasWake routine
        ELSE IF(ITYPE(IBND) .EQ. 280) THEN
            CALL BOUND280(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )

C
C-----------------------------------------------------------------------
C                                                Analytical solution bc 
C-----------------------------------------------------------------------
c
        ELSE IF(ITYPE(IBND) .EQ. 901) THEN
            CALL BOUND901(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,t,
     &                    imax,jmax,kmax,ibnd,Ltstep        )

        ELSE IF(ITYPE(IBND) .EQ. 902) THEN
            CALL BOUND902(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,t,
     &                    imax,jmax,kmax,ibnd,Ltstep        )

        ELSE IF(ITYPE(IBND) .EQ. 905) THEN
            CALL BOUND905(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,t,
     &                    imax,jmax,kmax,ibnd,Ltstep        )

        ELSE IF(ITYPE(IBND) .EQ. 906) THEN
            CALL BOUND906(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,t,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
c
        ENDIF
C
 10   CONTINUE
C
C-----------------------------------------------------------------------
C                                           periodic boundary conditions
C-----------------------------------------------------------------------
C
      DO 20 ibnd=1,6
        IF( itype(ibnd).eq.500 ) THEN
          CALL BOUND500(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
        elseif( itype(ibnd).eq.5000 ) THEN
          CALL BOUND5000(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
        ENDIF
 20   CONTINUE
C
C
      RETURN
      END


*//
*//
*//////// WALL BOUNDARY CONDITIONS  ///////////////////////////////////
*//
*//
C
C
C---Subroutine Bound050------------------------I. Balaras -19/2/99-----
C
      SUBROUTINE BOUND050(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
C
C     No-slip wall boundary conditions
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C
      REAL   vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &       vec3(imax,jmax,kmax)
C
C-----------------------------------------------------------------------
C
*
****  Starting and ending indexes for ibnd 
*
      ibx=ibegin(ibnd)
      iex=iend(ibnd)
      jby=jbegin(ibnd)
      jey=jend(ibnd)
      kbz=kbegin(ibnd)
      kez=kend(ibnd)
c     write(6,*) ibx, iex, jby, jey, kbz, kez
*     
****  Find orientation of boundary : 
*     i=const (isx=1,jsy=0,ksz=0)
*     j=const (isx=0,jsy=1,ksz=0)
*     k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez
ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  HARDWIRING NECESSARY IN CASE OF 2D FLOW:    VO 1/2006
      if (ibnd.eq.1.or.ibnd.eq.2.or.ibnd.eq.5.or.ibnd.eq.6.and.jey.eq.jby) jsy=0 !i.e. east/west/top/bottom bdrs need fixing
ccccccccccccccccccccccccccccccccccccccccccccccccc      
c     write(6,*) isx, jsy, ksz
*     
****  Find possition of boundary
*     ipx=1 : begining, ipx=-1 : end, etc.
*     position indexes for the other surfaces are 0
*     
      ipx=isign( 1,(ix2-ibx) ) * isx
      jpy=isign( 1,(jy2-jby) ) * jsy
      kpz=isign( 1,(kz2-kbz) ) * ksz
c     write(6,*) ipx, jpy, kpz
      
*     
****  Indexes pointing the first inner values
*     
      icx=( ((ix2-2)*ibx+ix2+2)/ix2 - ibx ) * isx
      jcy=( ((jy2-2)*jby+jy2+2)/jy2 - jby ) * jsy
      kcz=( ((kz2-2)*kbz+kz2+2)/kz2 - kbz ) * ksz
c     write(6,*) icx, jcy, kcz
*     
****  Main loop
*     
      DO 10 i=ibx,iex
         ii=i+icx
         DO 20 j=jby,jey
            jj=j+jcy
            DO 30 k=kbz,kez
               kk=k+kcz
               vec1(i,j,k) =  - vec1(ii    ,jj    ,kk    )*(jsy+ksz)
               vec2(i,j,k) =  - vec2(ii    ,jj    ,kk    )*(isx+ksz)
               vec3(i,j,k) =  - vec3(ii    ,jj    ,kk    )*(isx+jsy)
*     
****  Special treatment of bondaries with position -1 for
*     the impermeability condition. i.e. w(i,j,kz2)=0.
*     
*     East boundary
               IF(isx.eq.1 .and. ipx.eq.-1) THEN
                  vec1(ix2,j,k)=0.0
*     Back boundary
               ELSE IF(jsy.eq.1 .and. jpy.eq.-1) THEN
                  vec2(i,jy2,k)=0.0
*     Top boundary
               ELSE IF(ksz.eq.1 .and. kpz.eq.-1) THEN
                  vec3(i,j,kz2)=0.0
               ENDIF
               
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
      
      RETURN
      END
C
C
C---Subroutine Bound055------------------------ U Piomelli 10/5/02
C
      SUBROUTINE BOUND055(x,vec1,vec2,vec3,imax,jmax,kmax,ibnd)
C
C     Part sli, part no-slip wall boundary conditions
C     Works only for ibnd = 5 and 6
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C
      REAL   vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &       vec3(imax,jmax,kmax)
      REAL   x(imax)
C
C-----------------------------------------------------------------------
C
      if(ibnd .lt. 5) then
         write(6,*) ' Error on bc, ibnd .ne. 5'
         stop
      endif
c     
c***  Main loop
c     
c     xle is the leading edge
c     width is the point where the tanh is +/-0.9
c
      xle   = 10.
      width = 0.4
*
****  Starting and ending indexes for ibnd 
*
      ibx=ibegin(ibnd)
      iex=iend(ibnd)
      jby=jbegin(ibnd)
      jey=jend(ibnd)
      kbz=kbegin(ibnd)
      kez=kend(ibnd)
*     
****  Find orientation of boundary : 
*     i=const (isx=1,jsy=0,ksz=0)
*     j=const (isx=0,jsy=1,ksz=0)
*     k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez
*     
****  Find possition of boundary
*     ipx=1 : begining, ipx=-1 : end, etc.
*     position indexes for the other surfaces are 0
*     
      ipx=isign( 1,(ix2-ibx) ) * isx
      jpy=isign( 1,(jy2-jby) ) * jsy
      kpz=isign( 1,(kz2-kbz) ) * ksz
*     
****  Indexes pointing the first inner values
*     
      icx=( ((ix2-2)*ibx+ix2+2)/ix2 - ibx ) * isx
      jcy=( ((jy2-2)*jby+jy2+2)/jy2 - jby ) * jsy
      kcz=( ((kz2-2)*kbz+kz2+2)/kz2 - kbz ) * ksz
*     
****  Main loop
      DO 10 i=ibx,iex
         ii=i+icx
         smooth = TANH(1.1*width*(x(i)-xle))
         DO 20 j=jby,jey
            jj=j+jcy
            DO 30 k=kbz,kez
               kk=k+kcz
               vec1(i,j,k) =  - smooth *
     $              vec1(ii    ,jj    ,kk    )*(jsy+ksz)
               vec2(i,j,k) =  - smooth *
     $              vec2(ii    ,jj    ,kk    )*(isx+ksz)
               vec3(i,j,k) =  - vec3(ii    ,jj    ,kk    )*(isx+jsy)
*     
****  Special treatment of bondaries with position -1 for
*     the impermeability condition. i.e. w(i,j,kz2)=0.
*     
*     East boundary
               IF(isx.eq.1 .and. ipx.eq.-1) THEN
                  vec1(ix2,j,k)=0.0
*     Back boundary
               ELSE IF(jsy.eq.1 .and. jpy.eq.-1) THEN
                  vec2(i,jy2,k)=0.0
*     Top boundary
               ELSE IF(ksz.eq.1 .and. kpz.eq.-1) THEN
                  vec3(i,j,kz2)=0.0
               ENDIF
               
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
c      
      RETURN
      END
C
      SUBROUTINE BOUND0550(z,vec1,vec2,vec3,imax,jmax,kmax,ibnd) 
c     for flipped wake-bl calculation (bottom wall)
c     U. Piomelli 2002, V. Ovchinnikov (adaptation), 07/3/03
      implicit none
      include 'headers/common.h'
c      
      real xle,width,smooth
      integer imax,jmax,kmax,ibnd,i,j,k
      real z(kmax)
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),vec3(imax,jmax,kmax)

      xle=10.5
      width=0.8
c      write(2000+myrank,'(G20.10)') z
c      close(2000+myrank) 
      
      do k=1,kmax
        smooth = TANH(1.1*width*(z(k)-xle))
        do j=1,jmax
          vec3(1,j,k)=-smooth*vec3(2,j,k)
          vec2(1,j,k)=-smooth*vec2(2,j,k)
          vec1(1,j,k)=0.
	enddo
c	write(1000+myrank,*) k,smooth
      enddo
c      close(1000+myrank)
      return
      end	  
	  
      
C
C---Subroutine Bound060------------------------I. Balaras -21/07/2000---
C
      SUBROUTINE BOUND060(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
C
C     Slip-wall boundary conditions
C
C       - normal component is set to 0
C       - homogeneous Newmann for the other two
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C
      REAL   vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &       vec3(imax,jmax,kmax)
C
C-----------------------------------------------------------------------
C
*
****  Starting and ending indexes for ibnd 
*
      ibx=ibegin(ibnd)
      iex=iend(ibnd)
      jby=jbegin(ibnd)
      jey=jend(ibnd)
      kbz=kbegin(ibnd)
      kez=kend(ibnd)
*
****  Find orientation of boundary : 
*      i=const (isx=1,jsy=0,ksz=0)
*      j=const (isx=0,jsy=1,ksz=0)
*      k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez
ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  HARDWIRING NECESSARY IN CASE OF 2D FLOW:    VO 1/2006
      if (ibnd.eq.1.or.ibnd.eq.2.or.ibnd.eq.5.or.ibnd.eq.6.and.jey.eq.jby) jsy=0 !i.e. east/west/top/bottom bdrs need fixing
ccccccccccccccccccccccccccccccccccccccccccccccccc      
*
****  Find possition of boundary
*     ipx=1 : begining, ipx=-1 : end, etc.
*     position indexes for the other surfaces are 0
*
      ipx=isign( 1,(ix2-ibx) ) * isx
      jpy=isign( 1,(jy2-jby) ) * jsy
      kpz=isign( 1,(kz2-kbz) ) * ksz
     
*
****  Indexes pointing the first inner values
*
      icx=( ((ix2-2)*ibx+ix2+2)/ix2 - ibx ) * isx
      jcy=( ((jy2-2)*jby+jy2+2)/jy2 - jby ) * jsy
      kcz=( ((kz2-2)*kbz+kz2+2)/kz2 - kbz ) * ksz
*
**** Main loop
*
      DO 10 i=ibx,iex
         ii=i+icx
         DO 20 j=jby,jey
            jj=j+jcy
            DO 30 k=kbz,kez
               kk=k+kcz
               vec1(i,j,k) =  vec1(ii    ,jj    ,kk    )*(jsy+ksz)
               vec2(i,j,k) =  vec2(ii    ,jj    ,kk    )*(isx+ksz)
               vec3(i,j,k) =  vec3(ii    ,jj    ,kk    )*(isx+jsy)
*     
****  Special treatment of bondaries with position -1 for
*     the impermeability condition. i.e. w(i,j,kz2)=0.
*     
*     East boundary
               IF(isx.eq.1 .and. ipx.eq.-1) THEN
                  vec1(ix2,j,k)=0.0
*     Back boundary
               ELSE IF(jsy.eq.1 .and. jpy.eq.-1) THEN
                  vec2(i,jy2,k)=0.0
*     Top boundary
               ELSE IF(ksz.eq.1 .and. kpz.eq.-1) THEN
                  vec3(i,j,kz2)=0.0
               ENDIF
               
 30         CONTINUE
 20      CONTINUE
 10   CONTINUE
c      
      RETURN
      END
c
c------- Subroutine Bound062 ----- V. Ovchinnikov 8/02
c      
      SUBROUTINE BOUND062(vec1,vec2,vec3,vec11,vec22,vec33,imax,jmax,kmax,ibnd,Ltstep)
c     0 - neumann conditions for top boundary for all components      
      include 'headers/common.h'
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),vec3(imax,jmax,kmax),
     &     vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),vec33(imax,jmax,kmax)      
      integer i,j,ibnd!,imax,jmax,kmax
c      
      if (Ltstep) then 
       if (ibnd.eq.6) then
        do i=ix1,ix2+1
	  do j=jy1,jy2+1
	    vec1(i,j,kz2+1)=uinfty!vec1(i,j,kz2)
	    vec2(i,j,kz2+1)=vec2(i,j,kz2)
	    vec3(i,j,kz2)=vec33(i,j,kz2)
c	    write(4,*) vec3(i,j,kz2)
	  enddo
	enddo
       elseif (ibnd.eq.5) then
        do i=ix1,ix2+1
	  do j=jy1,jy2+1
	    vec1(i,j,kz1)=vec1(i,j,kz1+1)
	    vec2(i,j,kz1)=vec2(i,j,kz1+1)
	    vec3(i,j,kz1)=vec33(i,j,kz1)
	  enddo
	enddo
       endif	
c       
      else
       if (ibnd.eq.6) then	     
        do i=ix1,ix2+1
	  do j=jy1,jy2+1
c	    vec1(i,j,kz2+1)=vec1(i,j,kz2)
c	    vec2(i,j,kz2+1)=vec2(i,j,kz2)
	    vec3(i,j,kz2)=vec3(i,j,kz2-1)*1.
c            write(3,*) vec3(i,j,kz2)  
c	    vec3(i,j,kz2+1)=vec3(i,j,kz2)
	  enddo
	enddo
       elseif (ibnd.eq.5) then 
         do i=ix1,ix2+1
	   do j=jy1,jy2+1
	    vec3(i,j,kz1)=vec3(i,j,kz1+1) 	
	   enddo
	 enddo
       endif     	       
      endif
c      
      return
      end
c
*//
*//
*////////////////////  LAMINAR INFLOW CONDITIONS  ///////////////////////
*//
*//
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BOUND800(vec1,vec2,vec3,imax,jmax,kmax,zf,ibnd)
*                                                             *
*         Poiseille Inflow Backward Facing Step
*         valid only for i=1 !!
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      include 'headers/common.h'

      DIMENSION vec1(imax,jmax,kmax),zf(kmax)
      DIMENSION vec2(imax,jmax,kmax),vec3(imax,jmax,kmax)
C
C      DO j=jy1+1,jy2
C       DO k=kz1+1,kz2
CC        vec1(ix1,j,k) = -3*zf(k)*(zf(k)-1.)
CC        vec1(ix1,j,k) = 1.0
C        vec2(ix1,j,k) = -vec2(ix1+1,j,k)
C        vec3(ix1,j,k) = -vec3(ix1+1,j,k)
C
C       ENDDO
C      ENDDO
C      DO j=jy1,jy2+1
C      vec1(ix1,j,kz1) = 0.0
C      vec1(ix1,j,kz2+1) = 0.0
C      ENDDO
C
       DO j=jy1+1,jy2
        DO k= (kz2+1)/2+1,kz2
         vec1(ix1,j,k) = -16*zf(k)**2+24.*zf(k)-8.
C         vec1(ix1,j,k) = 1.0
        ENDDO
       ENDDO
C
      RETURN
      END
c
c     SUBROUTINE BOUND8000 Vitaly Ovchinnikov --- 5/29/03 
      SUBROUTINE BOUND8000(vec1, vec2, vec3, imax, jmax, kmax, xf)
      include 'headers/common.h'
      
      real u,x
      integer i,j,k
      
      DIMENSION vec1(imax,jmax,kmax),xf(imax)
      DIMENSION vec2(imax,jmax,kmax),vec3(imax,jmax,kmax)

      k=1
      do i=2,ix2
        x=0.5*(xf(i)+xf(i-1))/xf(ix2)
        u=4.*uinfty*x*(1-x)
	do j=2,jy2
	 vec1(i,j,k)=0.
	 vec2(i,j,k)=0.
	 vec3(i,j,k)=u
	enddo
      enddo
      return
      end	     

C
C
C---Subroutine Bound801------------------------E. Balaras-19/1/99------
      SUBROUTINE BOUND801(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf)
C
C
C       Imposes bc at inflow for a laminar blasius bl (i=1!!!)
C         - u=blasius
C         - v=0
C         - w=blasius
C
C         - Calls  subroutine Blasius
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &     vec3(imax,jmax,kmax),vec11(imax,jmax,kmax),
     &     vec22(imax,jmax,kmax),vec33(imax,jmax,kmax)
      real xf(imax),zf(kmax)
c
      real u1da(mmz),w1da(mmz),u1db(mmz),w1db(mmz)
c
      real   xloc,xref,dstaref
C
      logical compute
      data compute /.true./
    
      save u1da, w1db, compute
C
C----------------------------------------------------------------------
C
      if (compute) then
      compute=.false.
c..compute u blasius profile at inflow
      xloc   = 0.0
      xref   = -10.0
      dstaref= 0.
      CALL BLASIUS(zf,u1da,w1da,kmax,kmax,mmz,
     &             ru1,uinfty,xloc,xref,dstaref)

c$$$      do k=1,kmax
c$$$         write(6,*) zf(k),u1da(k),w1da(k)
c$$$      enddo
c$$$      stop

c..compute w blasius profile at inflow
      xloc   =-0.5*(xf(2)-xf(1))
      xref   = -10.0
      dstaref= 0.
      CALL BLASIUS(zf,u1db,w1db,kmax,kmax,mmz,
     &             ru1,uinfty,xloc,xref,dstaref)
      endif
c..impose bc for u_i*

      if (.not.ltstep) then 
       DO j=jy1+1,jy2
        DO k= kz1+1,kz2
         vec1(ix1,j,k) =  u1da(k)        
         vec2(ix1,j,k) = 0.!-vec2(ix1+1,j,k)
         vec3(ix1,j,k) =  w1db(k)
        ENDDO
       ENDDO
      else
       DO j=jy1+1,jy2
        DO k= kz1+1,kz2
         vec1(ix1,j,k) =  vec11(ix1,j,k)        
         vec2(ix1,j,k) =  vec22(ix1,j,k)
         vec3(ix1,j,k) =  vec33(ix1,j,k)
        ENDDO
       ENDDO
      endif        
      RETURN
      END
C----------------------------------------------------------------------
C
C---Subroutine Bound8010------------------------E. Balaras-19/1/99--V.Ovchinnikov 05/03----
      SUBROUTINE BOUND8010(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf,xref)
C
C
C       Imposes bc at inflow for a laminar blasius bl (i=1!!!)
C         - u=blasius
C         - v=0
C         - w=blasius
C
C         - Calls  subroutine Blasius
C
C----------------------------------------------------------------------
C
      implicit none
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      integer imax,jmax,kmax,i,j,k
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &     vec3(imax,jmax,kmax),vec11(imax,jmax,kmax),
     &     vec22(imax,jmax,kmax),vec33(imax,jmax,kmax)
      real xf(imax),zf(kmax)
c
      real u1da(mmx),w1da(mmx),u1db(mmx),w1db(mmx)
c
      real   xloc,xref,dstaref
C
      logical compute,ltstep
      data compute /.true./
    
      save u1da, w1db, compute
C
C----------------------------------------------------------------------
C
      if (compute) then
      compute=.false.
c..compute u blasius profile at inflow
      xloc   = 0.0
      dstaref= 0.
      CALL BLASIUS(xf,u1da,w1da,imax,imax,imax,
     &             ru1,uinfty,xloc,xref,dstaref)

c$$$      do k=1,kmax
c$$$         write(6,*) zf(k),u1da(k),w1da(k)
c$$$      enddo
c$$$      stop

c..compute w blasius profile at inflow
      xloc   =-0.5*(zf(2)-zf(1))
      dstaref= 0.
      CALL BLASIUS(xf,u1db,w1db,imax,imax,imax,
     &             ru1,uinfty,xloc,xref,dstaref)
      endif
c..impose bc for u_i*

      if (.not.ltstep) then 
       DO i=ix1+1,ix2
        DO j= jy1+1,jy2
         vec3(i,j,kz1) =  u1da(i)        
         vec2(i,j,kz1) = 0.!-vec2(i,j,kz1+1)
         vec1(i,j,kz1) =  w1db(i)
        ENDDO
       ENDDO
      else
       DO i=ix1+1,ix2
        DO j= jy1+1,jy2
         vec1(i,j,kz1) =  vec11(i,j,kz1)        
         vec2(i,j,kz1) =  vec22(i,j,kz1)
         vec3(i,j,kz1) =  vec33(i,j,kz1)
        ENDDO
       ENDDO
      endif        
      RETURN
      END
C----------------------------------------------------------------------
C---Subroutine Bound8030------------------------E. Balaras-19/1/99--V.Ovchinnikov 05/03----
      SUBROUTINE BOUND8030(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf,xref)
C
C
C       Imposes bc at OUTFLOW (test) for a laminar blasius bl (i=1!!!)
C         - u=blasius
C         - v=0
C         - w=blasius
C
C         - Calls  subroutine Blasius
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &     vec3(imax,jmax,kmax),vec11(imax,jmax,kmax),
     &     vec22(imax,jmax,kmax),vec33(imax,jmax,kmax)
      real xf(imax),zf(kmax)
c
      real u1da(mmx),w1da(mmx),u1db(mmx),w1db(mmx)
c
      real   xloc,xref,dstaref
C
      logical compute
      data compute /.true./
    
      save u1da, w1db, compute
C
C----------------------------------------------------------------------
C
      if (compute) then
      compute=.false.
c..compute u blasius profile at outflow
      xloc   = zf(kz2)
      dstaref= 0.
      CALL BLASIUS(xf,u1da,w1da,imax,imax,imax,
     &             ru1,uinfty,xloc,xref,dstaref)

c$$$      do k=1,kmax
c$$$         write(6,*) zf(k),u1da(k),w1da(k)
c$$$      enddo
c$$$      stop

c..compute w blasius profile at outflow
      xloc   =0.5*(3.*zf(kz2)-zf(kz2-1))
      dstaref= 0.
      CALL BLASIUS(xf,u1db,w1db,imax,imax,imax,
     &             ru1,uinfty,xloc,xref,dstaref)
      endif
c..impose bc for u_i*

      if (.not.ltstep) then 
       DO i=ix1+1,ix2
        DO j= jy1+1,jy2
         vec3(i,j,kz2) =  u1da(i)        
         vec2(i,j,kz2+1) = -vec2(i,j,kz2)
         vec1(i,j,kz2+1) =  w1db(i)
        ENDDO
       ENDDO
      else
       DO i=ix1+1,ix2
        DO j= jy1+1,jy2
         vec1(i,j,kz2+1) =  vec11(i,j,kz2+1)        
         vec2(i,j,kz2+1) =  vec22(i,j,kz2+1)
         vec3(i,j,kz2) =  vec33(i,j,kz2)
        ENDDO
       ENDDO
      endif        
      RETURN
      END
C----------------------------------------------------------------------
C
C---Subroutine Bound802------------------------E. Balaras-14/6/99------
C
      SUBROUTINE BOUND802(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,zf)

C     Uniform inflow for a Developing Channel (i=1 !!!)
C
C         - u   = 1.
C         - v   = 0.
C         - dwdx= 0
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      DIMENSION vec1(imax,jmax,kmax),zf(kmax)
      DIMENSION vec2(imax,jmax,kmax),vec3(imax,jmax,kmax),vec11(imax,jmax,kmax),
     &          vec22(imax,jmax,kmax),vec33(imax,jmax,kmax) 
      logical ltstep
C
      if (ltstep.ne..true.) then 
       DO j=jy1+1,jy2
         DO k= kz1+1,kz2
            vec1(ix1,j,k) =  uinfty
            vec2(ix1,j,k) = -vec2(ix1+1,j,k)
            vec3(ix1,j,k) =  vec3(ix1+1,j,k)
         ENDDO
       ENDDO
      else
       DO j=jy1+1,jy2
         DO k= kz1+1,kz2
            vec1(ix1,j,k) =  vec11(ix1,j,k)
            vec2(ix1,j,k) =  vec22(ix1,j,k)
            vec3(ix1,j,k) =  vec33(ix1,j,k)
         ENDDO
       ENDDO
      endif 
c     
      RETURN
      END
C

C---Subroutine Bound8020------------------------E. Balaras-14/6/99------
C
      SUBROUTINE BOUND8020(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax)

C     Uniform inflow for a Developing Channel (i=1 !!!)
C
C         - w   = 1.
C         - v   = 0.
C         - dudz= 0
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      DIMENSION vec1(imax,jmax,kmax)
      DIMENSION vec2(imax,jmax,kmax),vec3(imax,jmax,kmax),vec11(imax,jmax,kmax),
     &          vec22(imax,jmax,kmax),vec33(imax,jmax,kmax) 
      logical ltstep
C
      if (ltstep.ne..true.) then 
       DO i=ix1+1,ix2
         DO j= jy1+1,jy2
            vec3(i,j,kz1) = uinfty
            vec2(i,j,kz1) = 0!-vec2(i,j,kz1+1)
            vec1(i,j,kz1) = 0!vec1(i,j,kz1+1)
         ENDDO
       ENDDO
      else
       DO i=ix1+1,ix2
         DO j=jy1+1,jy2
            vec1(i,j,kz1) =  vec11(i,j,kz1)
            vec2(i,j,kz1) =  vec22(i,j,kz1)
            vec3(i,j,kz1) =  vec33(i,j,kz1)
         ENDDO
       ENDDO
      endif 
c     
      RETURN
      END
C
C---Subroutine Bound803------------------------E. Balaras-02/4/99------
C
      SUBROUTINE BOUND803(vec1,vec2,vec3,imax,jmax,kmax,zf,ibnd)

C        Poiseuille Inflow Channel
C
C         - u= 2.*zf(k)*(1.-zf(k))
C         - v=0
C         - w=0
C
C         - Calls  subroutine Blasius
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      DIMENSION vec1(imax,jmax,kmax),zf(kmax)
      DIMENSION vec2(imax,jmax,kmax),vec3(imax,jmax,kmax)
C
       DO j=jy1+1,jy2
        DO k= kz1+1,kz2
         vec1(ix1,j,k) = 2.*zf(k)*(1.-zf(k))
         vec2(ix1,j,k) = -vec2(ix1+1,j,k)
         vec3(ix1,j,k) = -vec3(ix1+1,j,k)
        ENDDO
       ENDDO
C
      RETURN
      END
C
C
C---Subroutine Bound804------------------------E. Balaras-14/6/99------
C
      SUBROUTINE BOUND804(vec1,vec2,vec3,imax,jmax,kmax,zf,ibnd)

C     Uniform inflow for a Developing Channel (i=1 !!!) at 45^o
C
C         - u= 1.
C         - v= 1.
C         - dwdx=0
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      DIMENSION vec1(imax,jmax,kmax),zf(kmax)
      DIMENSION vec2(imax,jmax,kmax),vec3(imax,jmax,kmax)
C
       DO j=jy1+1,jy2
        DO k= kz1+1,kz2
         vec1(ix1,j,k) =  1.0
c         vec2(ix1,j,k) =  1.0
         vec2(ix1,j,k) =  2.-vec2(ix1+1,j,k)
c         vec2(ix1,j,k) =  -vec2(ix1+1,j,k)
         vec3(ix1,j,k) =  vec3(ix1+1,j,k)
        ENDDO
       ENDDO
C
      RETURN
      END




C
C
C---Subroutine Bound805------------------------E. Balaras-1/8/2000------
C
      SUBROUTINE BOUND805(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    Ltstep,imax,jmax,kmax,xf,zf,ibnd        )
C
C       Imposes bc at inflow for a laminar blasius bl (i=1!!!)
C         - u=blasius
C         - v=0
C         - w=blasius
C
C         - super-imposes rondom noise  (set amp in the routine)
C
C         - Calls  subroutine Blasius
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &     vec3(imax,jmax,kmax)
      real vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &     vec33(imax,jmax,kmax)
      real xf(imax),zf(kmax)
c
      real u1da(mmz),w1da(mmz),u1db(mmz),w1db(mmz)
      real f1(mmz)
c
      real   xloc,xref,dstaref
      real   amp,a1,a2,z2
      integer ii
c
      logical Ltstep
C
C----------------------------------------------------------------------
C
c..local scalars
c      amp = 0.05/0.29
      amp = 0.15/0.29
      a1=1.
      a2=0.003
      z2=200.
      ii=10

c...create function
      do k=2,kz2
       f1(k) = 1.-exp(-a1*zf(k)**2)
       if(zf(k).gt.z2) then
         f1(k)=exp(-a2*(zf(k)-z2)**2)
       endif
       f1(k)=f1(k)*amp
c       write(123,*) k,zf(k),f1(k)
      enddo

c..compute u blasius profile at inflow
      xloc   =0.0
      xref   =0.0
      dstaref=1.0
      CALL BLASIUS(zf,u1da,w1da,kz2,kmax,mmz,
     &             ru1,uinfty,xloc,xref,dstaref)

c..compute w blasius profile at inflow
      xloc   =-0.5*(xf(2)-xf(1))
      xref   =0.0
      dstaref=1.0
      CALL BLASIUS(zf,u1db,w1db,kz2,kmax,mmz,
     &             ru1,uinfty,xloc,xref,dstaref)

c...bc for u^(n+1)
      IF(Ltstep) THEN  
       i=ix1
       DO j=jy1+1,jy2
        DO k=kz1+1,kz2
         vec1(i,j,k) = vec11(i,j,k)
         vec2(i,j,k) = f1(k)*(ran1(ii)-0.5)
         vec3(i,j,k) = w1db(k)+ f1(k)*(ran1(ii)-0.5)
        ENDDO
       ENDDO
c...bc for u*
      ELSE
       DO j=jy1+1,jy2
        DO k= kz1+1,kz2
         vec1(ix1,j,k) =  u1da(k) + f1(k)*(ran1(ii)-0.5)
        ENDDO
       ENDDO
c
      ENDIF
C
      RETURN
      END
c
C---Subroutine Bound806------------------------adopted from 801------
      SUBROUTINE BOUND806(vec1,vec2,vec3,vec11,vec22,vec33,ltstep,imax,jmax,kmax,xf,zf)
C
C
C       Imposes bc at inflow for a laminar blasius bl (i=1!!!)
C         - u=blasius
C         - v=0
C         - w=blasius
C
C         - Calls  subroutine Blasius
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &     vec3(imax,jmax,kmax),vec11(imax,jmax,kmax),
     &     vec22(imax,jmax,kmax),vec33(imax,jmax,kmax)
      real xf(imax),zf(kmax)
c
      real u1db(mmz),w1db(mmz)
c
      real   xloc,xref,dstaref
C
      logical compute
      data compute /.true./
    
      save w1db, compute
C
C----------------------------------------------------------------------
C
      if (compute) then
      compute=.false.
c..compute w blasius profile at inflow
      xloc   =-0.5*(xf(2)-xf(1))
      xref   = -10.0
      dstaref= 0.
      CALL BLASIUS(zf,u1db,w1db,kmax,kmax,mmz,
     &             ru1,uinfty,xloc,xref,dstaref)
      endif
c..impose bc for u_i*

      if (.not.ltstep) then 
       DO j=jy1+1,jy2
        DO k= kz1+1,kz2
         vec3(ix1,j,k) =  w1db(k)
        ENDDO
       ENDDO
      else
       DO j=jy1+1,jy2
        DO k= kz1+1,kz2
         vec1(ix1,j,k) =  vec11(ix1,j,k)
         vec2(ix1,j,k) =  vec22(ix1,j,k)
         vec3(ix1,j,k) =  vec33(ix1,j,k)
        ENDDO
       ENDDO
      endif        
      RETURN
      END
C----------------------------------------------------------------------

*//
*//
*////////////////////  TURBULENT INFLOW CONDITIONS  ///////////////////////
*//
*//
c
c
c---Subrotine bound650----------------------------------------------
c
      SUBROUTINE BOUND650(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    u2d,v2d,w2d,t,t_old,ksub,
     &                    nt,imax,jmax,kmax,ibnd,Ltstep)
c
c
c      Imposes prescribed values for u,v,w at the inflow plane
c      Data are coming from an auxiliary simulation
c
c       - reads a block of slice data and finds velocity
c         components at current time-step by linear 
c         interpolation
c
c       - appropriate set-up is required for this routine
c             a. parameter in the main
c             b. input file
c
c             Applied at i=ix1 !!!!!
c
c--------------------------------------------------------------------
c
      include 'headers/common.h'
c
c--------------------------------------------------------------------
c
      PARAMETER(iwrk=200)
      REAL scratchu(iwrk,iwrk),scratchv(iwrk,iwrk),
     &     scratchw(iwrk,iwrk)
      logical Ltstep
      REAL t,t_old,tt,tt1,tt2,coef1,coef2
      INTEGER idim,ipos,icomp,ksub

      DIMENSION vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &          vec3(imax,jmax,kmax),
     &          vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &          vec33(imax,jmax,kmax),
     &          u2d(nt,jmax,kmax),v2d(nt,jmax,kmax),
     &          w2d(nt,jmax,kmax)

c..predictor step  
      IF(.not.Ltstep) THEN   
        icomp=nslice/nsub
        ksub_old=ksub
        CALL AbsTime(t,t_old,ksub)
c...read a new block
         IF(ksub.ne.ksub_old) then
c....save last plane from previous block
          do j=jy1,jy2+1
           do k=kst,kz2
            scratchu(j,k)=u2d(nsub,j,k)
            scratchv(j,k)=v2d(nsub,j,k)
            scratchw(j,k)=w2d(nsub,j,k)
           enddo
          enddo
c....rewind file if required
           if(ksub.eq.1) then
            rewind(111)
           endif
          read(111) u2d,v2d,w2d
          write(6,*) 'read...ksub=',ksub
         ENDIF
c....find velocity at t with linear interpolation
         i=ix1
         ti=(isub-1)*dtslice
         coef1=(ti+dtslice-tsub)/dtslice
         coef2=(tsub-ti)/dtslice
c...interpolate between last ksub and first ksub+1
        IF(isub.eq.nsub) then
         DO j=jy1,jy2+1
          DO k=kst,kz2
           vec1(i,j,k) =  coef1*scratchu(j,k)
     &                  + coef2*u2d(1,j,k)
c
           vec2(i,j,k) =  coef1*scratchv(j,k)
     &                  + coef2*v2d(1,j,k)
c
           vec3(i,j,k) =  coef1*scratchw(j,k)
     &                  + coef2*w2d(1,j,k)
          ENDDO
         ENDDO
c
        ELSE
c...interpolate between to values of ksub
         DO j=jy1,jy2+1
          DO k=kst,kz2
           vec1(i,j,k) =  coef1*u2d(isub,j,k)
     &                  + coef2*u2d(isub+1,j,k)
c
           vec2(i,j,k) =  coef1*v2d(isub,j,k)
     &                  + coef2*v2d(isub+1,j,k)
c
           vec3(i,j,k) =  coef1*w2d(isub,j,k)
     &                  + coef2*w2d(isub+1,j,k)
c
          ENDDO
         ENDDO
        ENDIF
c
      ELSE
c...corrector step
       i=ix1
       DO j=jy1,jy2+1
        DO k=kst,kz2
         vec1(i,j,k) = vec11(i,j,k)
         vec2(i,j,k) = vec22(i,j,k) 
         vec3(i,j,k) = vec33(i,j,k)
        ENDDO
       ENDDO
c
      ENDIF
c
      RETURN
      END
c
c
c---Subrotine bound660----------------------------------------------
c
      SUBROUTINE BOUND660(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xu,yv,zw,uavneu,wavneu,
     &                    imax,jmax,kmax,ibnd,Ltstep,uinf,duinf)
c
c         Imposes u,v,w at the inflow plane by rescaling
c         data from a downstream station
c         For dp/dx=0 turbulent bl
c
c         - Calls linint
c         - The routine that calculates runing averages should be
c           called. Provides the arrays:  uavneu,wavneu
c
c         Applied at i=ix1 only !!!!!
c
c--------------------------------------------------------------------
c
      include 'headers/common.h'
c
c--------------------------------------------------------------------
c
      logical Ltstep

      REAL  vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &      vec3(imax,jmax,kmax),
     &      vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &      vec33(imax,jmax,kmax)
      REAL  xu(imax),yv(jmax),zw(kmax)
      real    uinf(imax), duinf(imax)
      REAL  uavneu(imax,kmax),wavneu(imax,kmax)
c     ...External Functions...
      external g
      real     g
c
      integer counter, idump, irc, jj, npl
      real alfa, beta, eta, fjp, gamma, lambda, weta
      real ddstai, dstin, delin, deltar, dstinn, dstarr
      real conpar, etabeg, etaend  
      real hshapi, hshapr
      real cappi, kappa, theinn, ueutau, und, wgth1, wgth2
      real smother(mmz)

      real urc(mmz),unl(mmz),uinner(mmz),uouter(mmz)
      real wrc(mmz),wnl(mmz),winner(mmz),wouter(mmz)
      real urcn(mmz),unln(mmz)
      real uold(mmz),wold(mmz)
      real u1dinner(mmz),u1douter(mmz)
      real v1dinner(mmz),v1douter(mmz)
      real w1dinner(mmz),w1douter(mmz)
      real ubinle(mmz), vbinle(mmz), wbinle(mmz)
      real uprecy(mmz), vprecy(mmz), wprecy(mmz)

      real ucord(mmz)
      real zplusrc(mmz), zplusnl(mmz)
      real wzplusrc(mmz), wzplusnl(mmz)
      real etarc(mmz), etanl(mmz), wetarc(mmz), wetanl(mmz)
      real uspace(mmz),wspace(mmz)
      real utest(mmz)
      real eps,step
      integer itermax

c------------------------------------------------------------------------
c                                                                set-up  
c------------------------------------------------------------------------
c
      eps  =0.001
      step =0.01 
      itermax=400

c...recycle station
c      irc =40
      irc =120

c...Power law index (npl = 7)
      npl = 6.0

c...Coefficients for the  Law of the wall and law of the wake
      kappa  = 0.41
      ueutau = 25.0
      cappi  = 0.55

c...Weighting function coefficients
      alfa = 4.0
      beta = 0.2

c...spanwise number of samples
      fjp   = 1.0/real(jy2-1)

c------------------------------------------------------------------------
c                                                        predictor step  
c------------------------------------------------------------------------
      IF(.not.Ltstep) THEN 
c
c...........................Recycle Station............................
c
        ucord(1) = 0.0
        urc(1)   = 0.0
        wrc(1)   = 0.0
c
        do k=2,kz2
          ucord(k)= 0.5*(zw(k)+zw(k-1))
c...<u>
          urc(k)  = uavneu(irc,k)
c...<v>
c     spanwise component is zero
c...<w>
          wrc(k)  = wavneu(irc,k)
        end do

c...Uinf at recycling station
        uinfrc=urc(kz2)

c...<u_t>
        h1 = ucord(2)
        h2 = ucord(3)
        coef1 = h2/(h1*(h2-h1))
        coef2 = h1/(h2*(h2-h1))
        dudz= abs(coef1*urc(2) - coef2*urc(3))
        utrc = sqrt(ru1*dudz)/uinfrc

c...Compute <d*>, <theta> and shape factor
        dstarrc  = 0.0
        thetarc  = 0.0
        do k = 2, kz2
          und     = urc(k)/uinfrc
          urcn(k) = urc(k)/uinfrc
            dstarrc = dstarrc +     (1.0 - und)*(zw(k) - zw(k-1))
            thetarc = thetarc + und*(1.0 - und)*(zw(k) - zw(k-1))
        end do

c...compute <d>
       uval=0.995
       call locatemine(urcn,mmz,kz2,uval,kpos)
       alpha1 = (ucord(kpos+1)-ucord(kpos))/(urcn(kpos+1)-urcn(kpos))
       beta1  = ucord(kpos)-alpha1*urcn(kpos)
       deltarc = alpha1*uval+beta1

c...normalized coordinates at recycle station
        do k=1, kz2
          etarc(k)    = ucord(k)/deltarc       !  n=z/delta
          wetarc(k)   = zw(k)/deltarc
          zplusrc(k)  = ucord(k)*utrc/ru1     !  z+   
          wzplusrc(k) = zw(k)*utrc/ru1
        end do
c
c.............................Inlet Station..............................
c
c..........Adjast theta at inflow through an iterative process...........

        ii=0
        thetaref = 1.0   !  adjust theta
cc        dstarref = 0.80     !  adjust d*

c...compute <d*>, <theta>
        uinfnl   = uavneu(1,kz2)
        thetanl  = 0.0
        dstarnl  = 0.0
        do k = 2, kz2
          und     = uavneu(1,k)/uinfnl
          unln(k) = uavneu(1,k)/uinfnl
            dstarnl = dstarnl +     (1.0 - und)*(zw(k) - zw(k-1))
            thetanl = thetanl + und*(1.0 - und)*(zw(k) - zw(k-1))
        end do

c...compute delta
       uval=0.995
       call locatemine(unln,mmz,kz2,uval,kpos)
       alpha1 = (ucord(kpos+1)-ucord(kpos))/(unln(kpos+1)-unln(kpos))
       beta1 = ucord(kpos)-alpha1*unln(kpos)
       deltanl = alpha1*uval+beta1

       utnl = utrc*(thetarc/thetanl)**(1.0/(npl + 1.0))
       gamma = utnl/utrc

c...iterative process to control d* 
 10     continue     
        ii=ii+1

c...normalized coordinates at inlet station
        do k=1, kz2
          etanl(k)    = ucord(k)/deltanl
          wetanl(k)   = zw(k)/deltanl
          zplusnl(k)  = ucord(k)*utnl/ru1
          wzplusnl(k) = zw(k)*utnl/ru1
        end do

c...Interpolate mean flow
c..<u>
        idump = 0
        call linint(zplusrc,zplusnl,urc,uinner,kz2,mmz,kz2,mmz,
     *              idump)

        call linint(etarc,  etanl,  urc,uouter,kz2,mmz,kz2,mmz,
     *              idump)

c..<v> (mean of the v is zero: no interpolation required)

c..<w>
        call linint(wzplusrc,wzplusnl,wrc,winner,kz2,mmz,kz2,mmz,
     *              idump)
c
        call linint(wetarc,  wetanl,  wrc,wouter,kz2,mmz,kz2,mmz,
     *              idump)

c
c...compute resulting <u>
c
        do k = 2, kz2
          eta  = etanl(k)
          wgth1 = 1.0
          if (eta.le.1.0) then
            wgth1 = 0.5*( 1.0 + ( tanh( (alfa*(eta - beta))/
     *             ((1.0 - 2.0*beta)*eta + beta) ) ) / tanh(alfa) )
          end if
          utest(k)= gamma*uinner(k)*(1.0-wgth1)+
     &             (gamma*uouter(k)+(1.0-gamma)*uinfnl)*wgth1
        end do

c...Compute resulting <theta>
        dstarnl    = 0.0
        thetanl    = 0.0
        uinfnl   = utest(kz2)
        do k = 2, kz2
          und = utest(k)/uinfnl
            dstarnl = dstarnl +     (1.0 - und)*(zw(k) - zw(k-1))
            thetanl = thetanl + und*(1.0 - und)*(zw(k) - zw(k-1))
        end do

c...change delta until desired thata is achived
        tstatus = abs(thetanl-thetaref)  !adjust theta
cc        tstatus = abs(dstarnl-dstarref)    !adjust d*
        if(tstatus.ge.eps.and.ii.lt.itermax) then
         signstep = abs(thetaref-thetanl)/(thetaref-thetanl) !adjust theta
cc         signstep = abs(dstarref-dstarnl)/(dstarref-dstarnl)   !adjust d*
         deltanl = deltanl + signstep*step
         utnl = utrc*(thetarc/thetanl)**(1.0/(npl + 1.0))  !adjust theta
cc         utnl = utrc*(dstarrc/dstarnl)**(1.0/(npl + 1.0))    !adjust d*
         gamma = utnl/utrc
c
c         write(6,*) ii,deltanl,dstarnl
         write(6,*) ii,deltanl,thetanl
         goto 10
        endif
c
c...scale final <u>, <w> & compute smother function
c
        etabeg = 1.0
        etaend = etanl(kz2)
        conpar = 5
c
        do k = 1, kz2
          uinner(k) = gamma*uinner(k)
          uouter(k) = gamma*uouter(k)+(1.0-gamma)*uinfnl
          winner(k) = winner(k)*uinfnl
          wouter(k) = wouter(k)*uinfnl
c
          eta         = etanl(k)
          smother(k)  = g(eta,etabeg,etaend,conpar) 
        end do
c
c...Interpolate u', v' & w' (do one k-line at the time)
c
        DO j = 2, jy2      
c
          uprecy(1) = 0.0    
          vprecy(1) = 0.0
          wprecy(1) = 0.0
          do k = 2, kz2
            uprecy(k) = vec11(irc,j,k) - urc(k)       
            vprecy(k) = vec22(irc,j,k)      
            wprecy(k) = vec33(irc,j,k) - wrc(k)
          end do
c
          idump = 1
c
c...u'
          call linint(zplusrc,zplusnl,uprecy,u1dinner,kz2,mmz,
     *                kz2,mmz,idump)
          call linint(etarc,  etanl,  uprecy,u1douter,kz2,mmz,
     *                kz2,mmz,idump)

c...v'
          call linint(zplusrc,zplusnl,vprecy,v1dinner,kz2,mmz,
     *                kz2,mmz,idump)
          call linint(etarc,  etanl,  vprecy,v1douter,kz2,mmz,
     *                kz2,mmz,idump)
c...w'
          call linint(wzplusrc,wzplusnl,wprecy,w1dinner,kz2,mmz,
     *                kz2,mmz,idump)
          call linint(wetarc,  wetanl,  wprecy,w1douter,kz2,mmz,
     *                kz2,mmz,idump)
c
c...Compose final velocity at the inflow plane
c
          do k = 2, kz2
c..Weighting function 
            eta  = etanl(k)
            weta = wetanl(k)
c
            wgth1 = 1.0
            if (eta.le.1.0) then
              wgth1 = 0.5*( 1.0 + ( tanh( (alfa*(eta - beta))/
     *                ((1.0 - 2.0*beta)*eta + beta) ) ) / tanh(alfa) )
            end if
c
            wgth2 = 1.0
            if (weta.le.1.0) then
              wgth2 = 0.5*( 1.0 + ( tanh( (alfa*(weta - beta))/
     *                ((1.0 - 2.0*beta)*weta + beta) ) ) / tanh(alfa) )
            end if
c
c...Smooth out the perturbation in the free-stream
c
            u1douter(k) = u1douter(k)*smother(k)
            v1douter(k) = v1douter(k)*smother(k)
            w1douter(k) = w1douter(k)*smother(k)
c
            vec1(ix1,j,k) = (uinner(k)+gamma*u1dinner(k))*(1.0-wgth1)+
     *                      (uouter(k)+gamma*u1douter(k))*     wgth1
            vec2(ix1,j,k) = gamma*v1dinner(k)*(1.0-wgth1)+
     *                      gamma*v1douter(k)*     wgth1
            vec3(ix1,j,k) = (winner(k)+gamma*w1dinner(k))*(1.0-wgth2)+
     *                      (wouter(k)+gamma*w1douter(k))*     wgth2
c
          end do
c
        END DO
c

      ELSE
c------------------------------------------------------------------------
c                                                        corrector step      
c------------------------------------------------------------------------
        i=ix1
        do k=2,kz2
          do j=2,jy2
            vec1(i,j,k) = vec11(i,j,k)
            vec2(i,j,k) = vec22(i,j,k) 
            vec3(i,j,k) = vec33(i,j,k)
          end do
        end do
c
      ENDIF
c
      RETURN
      END




*//
*//
*/////////// OUTFLOW BOUNDARY CONDITIONS  //////////////////////////
*//
*//
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BOUND710(vec1,vec2,vec3,vec11,vec22,vec33,
     &                       xf,yf,zf,imax,jmax,kmax,Ltstep)
*                                                             *
*            Outflow BC                                       *
*            Homogeneous Newmann bc for u,v,w                 *
*            (!!!Valid only for ibnd=2)                       *
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      include 'headers/common.h'

      DIMENSION vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &          vec3(imax,jmax,kmax)
      DIMENSION vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &          vec33(imax,jmax,kmax)
      DIMENSION xf(imax),yf(jmax),zf(kmax)
      LOGICAL Ltstep
C
C...BC for U^(n+1)
      if(Ltstep) then
c     
         DO j=jy1+1,jy2 
            DO k=kz1+1,kz2
c...  U^(n+1)=U^* 
               vec1(ix2,j,k)   = vec11(ix2,j,k)
c...  V^(n+1)_b=V^n_(b-1)
               vec2(ix2+1,j,k) = vec2(ix2,j,k)
c...  Same as for V
               vec3(ix2+1,j,k) = vec3(ix2,j,k)
            ENDDO
         ENDDO
C     
      else
C...  BC for U^*
C.... comp. normal to the boundary : U^*_b  = U^n_(b-1)
C     
         DO j=jy1+1,jy2
            DO k=kz1+1,kz2
               vec1(ix2,j,k) = vec11(ix2-1,j,k) 
c     c         vec1(ix2,j,k) = vec1(ix2-1,j,k) 
            ENDDO
       ENDDO

      endif
C
      RETURN
      END


* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BOUND720(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,imax,jmax,kmax,Ltstep,dtm)
*                                                             *
*       Outflow BC                                            *
*       Extrapolation from the interior for u,v,w             *
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *


      include 'headers/common.h'

      DIMENSION vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &          vec3(imax,jmax,kmax)
      DIMENSION vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &          vec33(imax,jmax,kmax)
      DIMENSION xf(imax),yf(jmax),zf(kmax)
      LOGICAL Ltstep
c
      REAL par1,par2,dum1,dum2

C
C...BC for U^(n+1)
      if(Ltstep) then
C
       dum1=(xf(ix2+1)-xf(ix2-1))/2.
       dum2=(xf(ix2)-xf(ix2-2))/2.
       par2=(dum1+dum2)/dum2
       DO j=jy1+1,jy2 
        DO k=kz1+1,kz2
         vec1(ix2,j,k)   = vec11(ix2,j,k)
         vec2(ix2+1,j,k) = par2     *vec2(ix2  ,j,k)+
     &                     (1.-par2)*vec2(ix2-1,j,k)
         vec3(ix2+1,j,k) = par 2    *vec3(ix2  ,j,k)+
     &                     (1.-par2)*vec3(ix2-1,j,k)
        ENDDO
       ENDDO
C
      else

C...BC for U^* (component normal to the boundary)
       par1=(xf(ix2)-xf(ix2-2))/(xf(ix2-1)-xf(ix2-2))
       DO j=jy1+1,jy2
        DO k=kz1+1,kz2
         vec1(ix2,j,k)= par1    *vec11(ix2-1,j,k)+
     &                 (1.-par1)*vec11(ix2-2,j,k)
        ENDDO
       ENDDO
C 
      endif
C
      RETURN
      END
C
C
C---Subroutine Bound730------------------------E. Balaras-19/1/99------
C
      SUBROUTINE BOUND730(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,t,imax,jmax,kmax,Ltstep,dtm)
C
C            Outflow BC           
C            Convective bc (du_i/dt=U_con du_i/dx_i) 
C            (!!!! Valid only for ibnd=2)
C            Convective speed is computed at outflow plane
C            MassBla routine should be called to ensure
C            mass conservation
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
      include 'mpif.h'
C
C----------------------------------------------------------------------
C
      REAL    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax)
      REAL    vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      REAL    xf(imax),yf(jmax),zf(kmax)
      REAL    uav1d(mmz)
      LOGICAL Ltstep
      REAL    deltax,coeff,sum,t
      real qstar1,area
      integer type
#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
C
C...BC for U^(n+1)
C
      if(Ltstep) then
C
         DO j=jy1+1,jy2 
            DO k=kz1+1,kz2
               vec1(ix2,j,k)   = vec11(ix2,j,k)
               vec2(ix2+1,j,k) = vec22(ix2+1,j,k)
               vec3(ix2+1,j,k) = vec33(ix2+1,j,k)
            ENDDO
         ENDDO
C     
      else
C     
C...  BC for U^*
C...  du_i/dt=U_con du_i/dx_i
c...  calculate uav
         samples=real(jy2-1)
         do k=2,kz2
            uav1d(k)=0.0
            do j=2,jy2
               uav1d(k)=uav1d(k)+vec11(ix2,j,k)
            enddo
            uav1d(k)=uav1d(k)/samples
         enddo
c     
         sum=0.0
         do k=2,kz2
            sum = sum +(zf(k)-zf(k-1))*uav1d(k)
         enddo
         call mpi_reduce(sum,qstar1,1,type,MPI_SUM,mysize-1,
     &        MPI_COMM_HYDRO,IERR)
c     
c     the last processor computes the cross sectional area
c     
         if (myrank.eq.mysize-1) then
            AREA= ( ZF(KZ2))
            qstar1=qstar1/AREA
         endif
         
         call MPI_BCAST(qstar1,1,type,mysize-1,MPI_COMM_HYDRO,IERR)
         uconv=qstar1
         deltax=xf(ix2)-xf(ix2-1)
         coeff  =uconv*dtm/deltax
         DO j=jy1+1,jy2
            DO k=kz1+1,kz2
               vec1(ix2,j,k)   = vec11(ix2,j,k)
     &              -coeff*(vec11(ix2,j,k)-vec11(ix2-1,j,k))
               vec2(ix2+1,j,k) = vec22(ix2+1,j,k)
     &              -coeff*(vec22(ix2+1,j,k)-vec22(ix2,j,k))
               vec3(ix2+1,j,k) = vec33(ix2+1,j,k)
     &              -coeff*(vec33(ix2+1,j,k)-vec33(ix2,j,k))
            ENDDO
         ENDDO
c     
      endif
C     
      RETURN
      END
C
C
C---Subroutine Bound740------------------------E. Balaras-27/07/2000------
C
      SUBROUTINE BOUND740(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,t,imax,jmax,kmax,Ltstep,dtm)
C
C            Outflow BC           
C            Convective bc (du_i/dt=U_con du_i/dx_i) 
C            (!!!! Valid only for ibnd=2)
C            Convective speed is computed at inflow plane 
C            and is kept constant
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
      include 'mpif.h'
C
C----------------------------------------------------------------------
C
      REAL    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax)
      REAL    vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      REAL    xf(imax),yf(jmax),zf(kmax)
      REAL    uav1d(mmz)
      LOGICAL Ltstep
      REAL    deltax,coeff,sum,t
      real qstar1,area
      integer type
#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

C
C...BC for U^(n+1)
      if(Ltstep) then
C
       DO j=jy1+1,jy2 
        DO k=kz1+1,kz2
         vec1(ix2,j,k)   = vec11(ix2,j,k)
         vec2(ix2+1,j,k) = vec22(ix2+1,j,k)
         vec3(ix2+1,j,k) = vec33(ix2+1,j,k)
        ENDDO
       ENDDO
C
      else
C
C...BC for U^*
C...du_i/dt=U_con du_i/dx_i
c...calculate uav
       samples=real(jy2-1)
       do k=2,kz2
        uav1d(k)=0.0
        do j=2,jy2
         uav1d(k)=uav1d(k)+vec11(1,j,k)
        enddo
        uav1d(k)=uav1d(k)/samples
       enddo
c
       sum=0.0
       do k=2,kz2
        sum = sum +(zf(k)-zf(k-1))*uav1d(k)
       enddo

         call mpi_reduce(sum,qstar1,1,type,MPI_SUM,mysize-1,
     &     MPI_COMM_HYDRO,IERR)



c     the last processor computes the cross sectional area

      if (myrank.eq.mysize-1) then
         AREA= ( ZF(KZ2))
         qstar1=qstar1/AREA
      endif

      call MPI_BCAST(qstar1,1,type,mysize-1,MPI_COMM_HYDRO,IERR)
      
      uconv=qstar1

       write(6,*) 'uconv=',uconv       

       deltax=xf(ix2)-xf(ix2-1)
       coeff  =uconv*dtm/deltax
       DO j=jy1+1,jy2
        DO k=kz1+1,kz2
         vec1(ix2,j,k)   = vec11(ix2,j,k)
     &                     -coeff*(vec11(ix2,j,k)-vec11(ix2-1,j,k))
         vec2(ix2+1,j,k) = vec22(ix2+1,j,k)
     &                     -coeff*(vec22(ix2+1,j,k)-vec22(ix2,j,k))
         vec3(ix2+1,j,k) = vec33(ix2+1,j,k)
     &                     -coeff*(vec33(ix2+1,j,k)-vec33(ix2,j,k))
        ENDDO
       ENDDO
c 
      endif
C
      RETURN
      END
c
c
C
C---Subroutine Bound7300------------------------E. Balaras-19/1/99--V.Ovchinikov 05/03
C
      SUBROUTINE BOUND7300(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,t,imax,jmax,kmax,Ltstep,dtm)
C
C            Outflow BC           
C            Convective bc (du_i/dt=U_con du_i/dx_i) 
C            (!!!! Valid only for ibnd=6)
C            Convective speed is computed at outflow plane
C            MassBlas routine should be called to ensure
C            mass conservation
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
      include 'mpif.h'
C
C----------------------------------------------------------------------
C
      REAL    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax)
      REAL    vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)!,dp(imax,jmax,kmax)
      REAL    xf(imax),yf(jmax),zf(kmax)
      REAL    uav1d(mmx)
      LOGICAL Ltstep
      REAL    deltaz,coeff,sum,t
      real qstar1,area
      integer type
#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
C
C...BC for U^(n+1)
C
      if(Ltstep) then
C
         DO i=ix1+1,ix2 
            DO j=jy1+1,jy2
               vec3(i,j,kz2)   = vec33(i,j,kz2)!+dp(i,j,kz2)/(zf(kz2)-zf(kz2-1))*dtm
               vec1(i,j,kz2+1) = vec11(i,j,kz2+1)
               vec2(i,j,kz2+1) = vec22(i,j,kz2+1)
            ENDDO
         ENDDO
C     
      else
C     
C...  BC for U^*
C...  du_i/dt=U_con du_i/dx_i
c...  calculate uav
         samples=real(jy2-1)
         do i=2,ix2
            uav1d(i)=0.0
            do j=2,jy2
               uav1d(i)=uav1d(i)+vec33(i,j,kz2)
            enddo
            uav1d(i)=uav1d(i)/samples
         enddo
c     
         sum=0.0
         do i=2,ix2
            sum = sum +(xf(i)-xf(i-1))*uav1d(i)
         enddo
c
            AREA= ( XF(IX2))
            qstar1=sum/AREA
c         
         uconv=qstar1
         deltaz=zf(kz2)-zf(kz2-1)
         coeff  =uconv*dtm/deltaz
         DO i=ix1+1,ix2
            DO j=jy1+1,jy2
               vec3(i,j,kz2)   = vec33(i,j,kz2)
     &              -coeff*(vec33(i,j,kz2)-vec33(i,j,kz2-1))
               vec2(i,j,kz2+1) = vec22(i,j,kz2+1)
     &              -coeff*(vec22(i,j,kz2+1)-vec22(i,j,kz2))
               vec1(i,j,kz2+1) = vec11(i,j,kz2+1)
     &              -coeff*(vec11(i,j,kz2+1)-vec11(i,j,kz2))
            ENDDO
         ENDDO
c     
      endif
C     
      RETURN
      END


*//
*//
*///////////////   OPEN BOUNDARY CONDITIONS  //////////////////////////
*//
*//
C
C
C---Subroutine Bound240------------------------E.Balaras 99------
C
      SUBROUTINE BOUND240(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
C
C       Freestream boundary condition for laminar bl
C       with dp/dx=0
C
C         - dudz=0
C         - dvdz=0
C         - w = prescribed
C
C                valid for k=kz2 !!
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)
      real    dsdx(mmx)
      real    xw,X0, xstart
      real    dstar
      real    dref
C
C----------------------------------------------------------------------
       xstart=0.0 ! this is the position at which the boundary layer is imposed(inflow)      
       dstar=0. ! this is the displacement thickness we impose at the position xstart
       dref=1.21982081876391! this is the integral of 1-f' over [0,inf) 
c      X0   = 0.0 - 0.5*1.0*((1.0/1.2168)**2)/ru1
c       X0   = xstart - 0.5*uinfty*((dstar/dref)**2)/ru1
       X0   =-10.0
c       write(6,*) x0 !(should be negative) -- this is where the plate starts
c      if (.not.Ltstep) then 
      do i=1,ix2
       xw = 0.5*(xf(i)+xf(i-1))
       dsdx(i) = 0.5*sqrt(2.*ru1/(xw-x0)/uinfty)*dref
      enddo
c     
c      wavg=0.0
      do i=ix1,ix2+1
         do j=jy1,jy2+1
            vec1(i,j,kz2+1)=uinfty!vec1(i,j,kz2)
            vec2(i,j,kz2+1)=vec2(i,j,kz2)
            vec3(i,j,kz2  )=uinfty*dsdx(i)*1.
c	    wavg=wavg+vec3(i,j,kz2)
	enddo
c	write(98,*) i,dsdx(i)*uinfty
      enddo
c       wavg=wavg/(ix2-ix1)/(jy2-jy1)
c       write(6,*) 'average w on top: ', wavg
c      else 
c	do i=ix1,ix2+1
c         do j=jy1,jy2+1
c            vec1(i,j,kz2+1)=vec1(i,j,kz2)
c            vec2(i,j,kz2+1)=vec2(i,j,kz2)
c            vec3(i,j,kz2  )=vec33(i,j,kz2)
c         enddo
c	enddo 
c      endif
      RETURN
      END
C
C
C---Subroutine Bound2400------------------------E.Balaras 99-V. Ovchinnikov-5/03----
C
      SUBROUTINE BOUND2400(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep,dtm,x0        )
C
C       Freestream boundary condition for laminar bl
C       with dp/dx=0
C
C         - dudz=0
C         - dvdz=0
C         - w = prescribed
C
C                valid for k=kz2 !!
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)!,dp(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)
      real    dsdx(kmax)
      real    xw,X0
      real    dstar
      real    dref
      real    dtm
C
C----------------------------------------------------------------------
c       x0 is the position at which the boundary layer is imposed(inflow)      
       dstar=0. ! this is the displacement thickness we impose at the position xstart
       dref=1.21982081876391! this is the integral of 1-f' over [0,inf) 
c      X0   = 0.0 - 0.5*1.0*((1.0/1.2168)**2)/ru1
c       X0   = xstart - 0.5*uinfty*((dstar/dref)**2)/ru1
c       write(6,*) x0 !(should be negative) -- this is where the plate starts
c      if (.not.Ltstep) then
 
       if (myrank.ne.0) then 
         xw = 0.5*(3.*zf(1)-zf(2))
         dsdx(1) = 0.5*sqrt(2.*ru1/(xw-x0)/uinfty)*dref
       endif
c
       do k=2,kz2+1
         xw = 0.5*(zf(k)+zf(k-1))
         dsdx(k) = 0.5*sqrt(2.*ru1/(xw-x0)/uinfty)*dref
       enddo
c     
      if (.not.Ltstep) then 
       do j=jy1,jy2+1
         do k=kz1,kz2+1
            vec3(ix2+1,j,k)=uinfty!vec3(ix2,j,k)
            vec2(ix2+1,j,k)=vec2(ix2,j,k)
            vec1(ix2,j,k  )=uinfty*dsdx(k)*1.
	enddo
       enddo
      else
       do j=jy1,jy2+1
         do k=kz1,kz2+1
            vec3(ix2+1,j,k)=uinfty!vec3(ix2,j,k)
            vec2(ix2+1,j,k)=vec2(ix2,j,k)
            vec1(ix2,j,k  )=uinfty*dsdx(k)*1.
	enddo
       enddo
      endif
      RETURN
      END
C
C---Subroutine Bound2401-------------------------V. Ovchinnikov-6/04----
C
      SUBROUTINE BOUND2401(vec1,vec2,vec3,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd)
C
C       Freestream boundary condition for blockage flow
c       only for elliptic/superelliptic LE
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
      include 'mpif.h'
C
C----------------------------------------------------------------------
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax),zc(kmax)
      real    xw
      integer ibnd,k,k1,k2,kk
c      parameters
      real ar, x0, r, m, h, dh(mmz), dhg(mmz), sum, sum2, sumg,sum2g, zg(mmz)
      integer ierr, type
      
      data init /.true./
      save dh,init
C
#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

      
      if (init) then 
      
c       init=.false. 
       
       ar=6.
       x0=10.
       r=1.
       m=4. !2.- ellipse, 4.- super-ellipse
      
       zc(1) = 0.5*(3.*zf(1)-zf(2))
       k1=1
       do while (zc(k1).le.x0.and.k1.lt.kmax)
         k1=k1+1
 	 zc(k1)=0.5*(zf(k1)+zf(k1-1))
         dh(k1)=0
       enddo	
       k2=k1
       do while (zc(k2).le.x0+ar*r.and.k2.lt.kmax)
         k2=k2+1
 	 zc(k2)=0.5*(zf(k2)+zf(k2-1))
         z=(zc(k2)-x0)/ar-1
         h=sqrt(r**m-(z)**m)
         dh(k2)=-0.5*m*(z)**(m-1.)/h/ar
       enddo	
      
       do k=k2+1,kz2+1
         zc(k)=0.5*(zf(k)+zf(k-1))
         dh(k)=0
       enddo	

        sum=0.
        do k=kz1+1,kz2
         sum=sum+dh(k)*(zf(k)-zf(k-1)) ! total actual displacement
        enddo
	call MPI_ALLREDUCE(sum, sumg, 1, type, MPI_SUM, MPI_COMM_HYDRO, ierr)
ccccccccccccccccccccccccccccccccccccccccc
c       gather before smoothing
        call MPI_GATHER(dh(2), kmax-2, type, dhg(2), kmax-2, type, mysize-1, MPI_COMM_HYDRO, ierr)
        call MPI_GATHER(zc(2), kmax-2, type, zg(2), kmax-2, type, mysize-1, MPI_COMM_HYDRO, ierr)
      

        if (myrank.eq.mysize-1) then 
c       filter to smooth
        kk=2
	do while (zg(kk).le.(x0+ar).and.kk.lt.mmz-1)
	 kk=kk+1
	enddo 
         do n=1,300
          do k=2,kk
           dhg(k)=0.25*(dhg(k-1)+2.*dhg(k)+dhg(k+1) )
          enddo
          do k=kk,2,-1
           dhg(k)=0.25*(dhg(k-1)+2.*dhg(k)+dhg(k+1) )
          enddo
         enddo  	
c
         do n=1,10   ! smooth corner at kk
          do k=2,mmz-1
           dhg(k)=0.25*(dhg(k-1)+2.*dhg(k)+dhg(k+1) )
          enddo
         enddo
        endif
c	
       call MPI_SCATTER(dhg(2), kmax-2, type, dh(2), kmax-2, type, mysize-1, MPI_COMM_HYDRO, ierr)
	 
cccccccccccccccccccccccccccccccccccccccccc
c        do k=kz1+1,kz2
c         dh2(k)=-tanh(0.08*(zc(k)-1.15*x0))**2+1.
c        enddo

        sum2=0.
        do k=kz1+1,kz2
         sum2=sum2+dh(k)*(zf(k)-zf(k-1)) ! total actual displacement
        enddo
c
	call MPI_ALLREDUCE(sum2, sum2g, 1, type, MPI_SUM, MPI_COMM_HYDRO, ierr)
        fact=sumg/sum2g
c
        do k=kz1+1,kz2
         dh(k)=dh(k)*fact
        enddo
c
c     debug -- check velocity field
       do k=kz1+1,kz2
         write(myrank+200,'(2G20.10)') zc(k),dh(k)
       enddo
       close(myrank+200)	 
c     debug
c       
      endif
 
       if (ibnd.eq.2) then      
        do k=kz1,kz2+1
          do j=jy1,jy2+1
             vec3(ix2+1,j,k)=vec3(ix2,j,k)
             vec2(ix2+1,j,k)=vec2(ix2,j,k)
             vec1(ix2,j,k  )=uinfty*dh(k)
	 enddo
        enddo
       elseif (ibnd.eq.1) then 
        do k=kz1,kz2+1
          do j=jy1,jy2+1
             vec3(ix1,j,k)=vec3(ix1+1,j,k)
             vec2(ix1,j,k)=vec2(ix1+1,j,k)
             vec1(ix1,j,k)=-uinfty*dh(k)
	 enddo
        enddo 
       endif	
c       

      RETURN
      END
C
C
C---Subroutine Bound250------------------------------------------------
C
      SUBROUTINE BOUND250(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
C
C                         
C       Freestream boundary condition for tbl with dpdx=0
C       (used in the fringe computations)
C 
C       GRESHO and SANI approach
C       Imposes du*/dn=0  + correction for mass balance.                   
C       Other componets are obtained by applying omega=0 
C       ----- the routine works only for k=kz2!! ----    
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)

c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------

      IF(.not.Ltstep) THEN
c..impose homogeneous Newmann to u*_n  
      do i=ix1,ix2+1
       do j=jy1,jy2+1
        vec3(i,j,kz2)=vec3(i,j,kz2-1)
       enddo
      enddo
c..corect with uniform velocity to conserve mass
      q1=0.0
      k=kz2
      do i=ix1+1,ix2
       do j=jy1+1,jy2
        q1=q1+vec3(i,j,k)*(xf(i)-xf(i-1))*(yf(j)-yf(j-1))
       enddo
      enddo
c
c      write(6,*) 'qtop=',q1
c
      surf=(xf(ix2)-xf(ix1))*(yf(jy2)-yf(jy1))
      wuni=-1.*(q1/surf)
c
c      write(6,*) 'uni w=',wuni
c
      do i=ix1+1,ix2
       do j=jy1+1,jy2
        vec3(i,j,k)=vec3(i,j,k)+wuni
c        vec3(i,j,k)=0.0
       enddo
      enddo
c-------------------------------------------------------------------
c                                                 corrector step      
c-------------------------------------------------------------------
      ELSE
c..u^(n+1)=u* for normal component
      do i=ix1,ix2+1
       do j=jy1,jy2+1
        vec3(i,j,kz2)=vec33(i,j,kz2)
       enddo
      enddo
c..derive bc for other two components from vorticity
      dzodx= DZ/DX
      dzody= DZ/DY
      do i = ix1+1, ix2
       do j = jy1+1, jy2
c         vec1(i,j,kz2+1) = vec1(i,j,kz2) +
c     &                     dzodx*(au(i)/cw(kz2))*
c     &                     (vec3(i+1,j,kz2)-vec3(i,j,kz2))

c         vec2(i,j,kz2+1) = vec2(i,j,kz2) +
c     &                     dzody*(bv(j)/cw(kz2))*
c     &                     (vec3(i,j+1,kz2)-vec3(i,j,kz2))
          vec1(i,j,kz2+1) = vec1(i,j,kz2)  
          vec2(i,j,kz2+1) = vec2(i,j,kz2) 
       enddo
      enddo 
c
      ENDIF
c
      RETURN
      END
C
C
C---Subroutine Bound255------------------------------------------------
C
      SUBROUTINE BOUND255(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
C            
C       Freestream boundary condition for tbl with dpdx=0
C 
C       Normal component from cont.      
C
C       ----- the routine works only for k=kz2!! ----    
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'

      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)
      real    cellx,celly,cellz
      real    fluxX,fluxY

c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------

      IF(.not.Ltstep) THEN
c
      DO 10 i=ix1+1,ix2
       DO 20 j=jy1+1,jy2
       cellx=(yf(j)-yf(j-1))*(zf(kz2)-zf(kz2-1))
       celly=(xf(i)-xf(i-1))*(zf(kz2)-zf(kz2-1))
       cellz=(xf(i)-xf(i-1))*(yf(j)-yf(j-1))
       fluxX=(vec11(i,j,kz2)-vec11(i-1,j,kz2))*cellx
       fluxY=(vec22(i,j,kz2)-vec22(i,j-1,kz2))*celly
       vec3(i,j,kz2)= vec33(i,j,kz2-1) + (fluxX+fluxY)/cellz
20     CONTINUE
10    CONTINUE


c-------------------------------------------------------------------
c                                                 corrector step      
c-------------------------------------------------------------------
      ELSE
c..u^(n+1)=u* for normal component
      do i=ix1,ix2+1
       do j=jy1,jy2+1
        vec1(i,j,kz2+1)=vec1(i,j,kz2)
        vec2(i,j,kz2+1)=vec2(i,j,kz2)
        vec3(i,j,kz2)=vec33(i,j,kz2)
       enddo
      enddo
c
      ENDIF
c
      RETURN
      END
C
C
C---Subroutine Bound260------------------------E. Balaras-19/1/99------
C
      SUBROUTINE BOUND260(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,urun,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
C
C       Freestream boundary condition for turbulent bl
C       with dp/dx=0
C
C        - dudz=0
C        - dvdz=0
C        - w computed from running averages 
C          forcing mass conservation
C
C         Running averages routine should be active
C         Calls routine 'fit'
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax),urun(imax,kmax)
      real    dsdx(mmx),deltastar(mmx),uav2d(mmx,mmz)
      real    fjp
C
C----------------------------------------------------------------------
C
c...predictor step
      IF(.not.Ltstep) THEN

c...calculate delta*(i) from power laws (with n=7)
c      rexinl = ((1.0/ru1)/0.02)**(7.0/6.0)
c      xxo    = - rexinl*ru1/uinfty
c      do i = 1, ix2+1
c        rex = (xf(i)-xxo)/ru1
c        delta = 0.37*(xf(i)-xxo)*rex**(-0.2)
c        deltastar(i) = delta/8.0
c      write(6,*) i,deltastar(i)
c      end do

c...calculate delta*(i)from the runing averages
        do i = ix1, ix2
c
          deltastar(i) = 0.0
          do k = 2, kz2 
            und = urun(i,k)/urun(i,kz2) 
            if (und.le.0.99) then
              deltastar(i) = deltastar(i) + 
     *        (1.0 - und)*(zf(k) - zf(k-1))
            end if
          end do
c
        end do
c
c...calculate d*/dx by linear regression
        call fit(xf,deltastar,imax,mmx,ix2,sig,
     *           0,a,b,siga,sigb,chi2,q     )
c        write(*,*) 'linear regr. coef. (exp)',a,b

c...calculate d*/dx numericaly
c      do i = ix1+1, ix2
c        dsdx(i)= (deltastar(i)-deltastar(i-1))/(xf(i)-xf(i-1))
c      end do
      
 
c...impose normal velocity
        do j = jy1+1,jy2
          do i = ix1+1,ix2
            dsdx(i) = b
            vec3(i,j,kz2  )= uinfty*dsdx(i)      
          end do
        end do
c
      ELSE

c...corrector step
      dzodx= DZ/DX
      dzody= DZ/DY
      do i=ix1+1,ix2
       do j=jy1+1,jy2

         vec1(i,j,kz2+1) = vec1(i,j,kz2) +
     &                     dzodx*(au(i)/cw(kz2))*
     &                     (vec3(i+1,j,kz2)-vec3(i,j,kz2))
c
         vec2(i,j,kz2+1) = vec2(i,j,kz2) +
     &                     dzody*(bv(j)/cw(kz2))*
     &                     (vec3(i,j+1,kz2)-vec3(i,j,kz2))
c
         vec3(i,j,kz2  )=vec33(i,j,kz2)
c         vec1(i,j,kz2+1)=vec1(i,j,kz2)
c         vec2(i,j,kz2+1)=vec2(i,j,kz2)
       end do
      end do
c
      ENDIF
c
      RETURN
      END
C
C
C---Subroutine Bound265-------------------------------------------------
C
      SUBROUTINE BOUND265(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
C
C       Freestream boundary condition for turbulent bl
C       with zero pressure gradient
C
C        - u from omega2=0
C        - v from omega1=0
C        - w computed from power laws (Dirichlet) 
C
C
C------------------------------------------------------------------------
C
      include 'headers/common.h'

      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)

c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------
      IF(.not.Ltstep) THEN
c   Power law assumption (with n=7)
       rexinl = ((1.0/ru1)/0.02)**(7.0/6.0)
       xori   = - rexinl*ru1/uinfty
       do i=ix1+1,ix2
         xw = 0.5*(xf(i)+xf(i-1))
         rex=(xw-xori)*uinfty/ru1
         wtop = 0.017*uinfty*rex**(-1./7.)
        do j=jy1,jy2+1
         vec3(i,j,kz2)= wtop
        enddo
c       write(6,*) i,wtop
       enddo

c-------------------------------------------------------------------
c                                                 corrector step      
c-------------------------------------------------------------------
      ELSE
c..u^(n+1)=u* for normal component
      dzodx= DZ/DX
      dzody= DZ/DY
      do i = ix1+1, ix2
        do j = jy1+1, jy2
         vec1(i,j,kz2+1) = vec1(i,j,kz2) +
     &                     dzodx*(au(i)/cw(kz2))*
     &                     (vec3(i+1,j,kz2)-vec3(i,j,kz2))
c
         vec2(i,j,kz2+1) = vec2(i,j,kz2) +
     &                     dzody*(bv(j)/cw(kz2))*
     &                     (vec3(i,j+1,kz2)-vec3(i,j,kz2))
c
         vec3(i,j,kz2)=vec33(i,j,kz2)
        enddo
      enddo 
c
      ENDIF
c
      RETURN
      END
C
C
C---Subroutine Bound270------------------------------------------------
C
      SUBROUTINE BOUND270(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,urun,
     &                    imax,jmax,kmax,ibnd,Ltstep,uinf,duinf)
C
C       Freestream boundary condition for turbulent bl
C       with pressure gradient
C
C        - dudz=0
C        - dvdz=0
C        - w computed from running averages 
C          forcing mass conservation
C
C         Running averages routine should be active
C         requires Uinf(i)=uinf(i) and the dUinf/dx(i)=duinf(i)
C
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax),urun(imax,kmax)
      real    uinf(imax), duinf(imax)
      real    fjp
      real    dsdx(mmx), deltastar(mmx),uav2d(mmx,mmz)
      real    uinfc(mmx), dstarc(mmx)
C
C----------------------------------------------------------------------
C

c...predictor step
      IF(.not.Ltstep) THEN
c
        do i = ix1, ix2
          deltastar(i) = 0.0
          do k = 2, kz2 
            und = urun(i,k)/urun(i,kz2)
           if (und.le.0.99) then
              deltastar(i) = deltastar(i) + 
     *        (1.0 - und)*(zf(k) - zf(k-1))
           end if
          end do
        end do
c
c...calculate uinfty & d* at center cell
c
      do i = ix1+1, ix2
        dstarc(i) = 0.5*(deltastar(i-1) + deltastar(i))
        uinfc(i)  = 0.5*(uinf(i-1)      + uinf(i))
      end do
c
c...calculate d*/dx numericaly
c
      do i = ix1+1, ix2
        dsdx(i)= (deltastar(i)-deltastar(i-1))/(xf(i)-xf(i-1))
      end do
c
c...impose normal velocity
        do j = jy1+1,jy2
          do i = ix1+1,ix2
            vec3(i,j,kz2  ) =  uinfc(i)*dsdx(i)
     *                      + (dstarc(i)-zf(kz2))*duinf(i)       
          end do
        end do
c
      ELSE

c...corrector step
      dzodx= DZ/DX
      dzody= DZ/DY
      do i=ix1+1,ix2
       do j=jy1+1,jy2
         vec1(i,j,kz2+1) = vec1(i,j,kz2) +
     &                     dzodx*(au(i)/cw(kz2))*
     &                     (vec3(i+1,j,kz2)-vec3(i,j,kz2))
c
         vec2(i,j,kz2+1) = vec2(i,j,kz2) +
     &                     dzody*(bv(j)/cw(kz2))*
     &                     (vec3(i,j+1,kz2)-vec3(i,j,kz2))
c
         vec3(i,j,kz2)=vec33(i,j,kz2)
c
c         vec1(i,j,kz2+1)=vec1(i,j,kz2)
c         vec2(i,j,kz2+1)=vec2(i,j,kz2)
       end do
      end do
c
      ENDIF
      RETURN
c
      END

C
C
C---Subroutine Bound270------------------------------------------------
C
      SUBROUTINE BOUND275(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,urun,
     &                    imax,jmax,kmax,ibnd,Ltstep,uinf,duinf,t)
C
C       Freestream boundary condition for turbulent bl
C       with pressure gradient
C
C        - dudz=0
C        - dvdz=0
C        - w computed from running averages 
C          forcing mass conservation
C
C         Running averages routine should be active
C         requires Uinf(i)=uinf(i) and the dUinf/dx(i)=duinf(i)
C
C         In this routine Uinfty is varing gradually from a uniform
C         distribution until a steady state has been reached
C
C         Time over which weight is applied is specified in the 
C         routine
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C----------------------------------------------------------------------
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax),urun(imax,kmax)
      real    uinf(imax), duinf(imax)
      real    fjp
      real    dsdx(mmx), deltastar(mmx),uav2d(mmx,mmz)
      real    uinfc(mmx), dstarc(mmx)
      real    timeint
C
C----------------------------------------------------------------------
C

c      timeint = 400.
      timeint = -1.0

c...specify wight
      if (t.lt.timeint) then
        weight = - t/timeint + 1.0
      else 
        weight = 0.0
      end if

c...predictor step
      IF(.not.Ltstep) THEN
c
        do i = ix1, ix2
          deltastar(i) = 0.0
          do k = 2, kz2 
            und = urun(i,k)/urun(i,kz2)
           if (und.le.0.99) then
              deltastar(i) = deltastar(i) + 
     *        (1.0 - und)*(zf(k) - zf(k-1))
           end if
          end do
        end do
c
c...calculate uinfty & d* at center cell
c
      do i = ix1+1, ix2
        dstarc(i) = 0.5*(deltastar(i-1) + deltastar(i))
        uinfc(i)  = 0.5*(uinf(i-1)      + uinf(i))
      end do
c
c...calculate d*/dx numericaly
c
      do i = ix1+1, ix2
        dsdx(i)= (deltastar(i)-deltastar(i-1))/(xf(i)-xf(i-1))
      end do

c
c...impose normal velocity
c
        do j = jy1+1,jy2
          do i = ix1+1,ix2
            vec3(i,j,kz2  ) =     weight * ( uinfty*dsdx(i) ) + 
     *                       (1.0-weight)*
     *                       ( uinfc(i)*dsdx(i)
     *                     + (dstarc(i)-zf(kz2))*duinf(i) )
cc           vec3(i,j,kz2  ) =      weight*( uinfty*dsdx(i) ) +
cc    *                       (1.0-weight)* winf(i)
cc           vec3(i,j,kz2  ) =  winf(i)
          end do
        end do
c
      ELSE

c...corrector step
      dzodx= DZ/DX
      dzody= DZ/DY
      do i=ix1+1,ix2
       do j=jy1+1,jy2
         vec1(i,j,kz2+1) = vec1(i,j,kz2) +
     &                     dzodx*(au(i)/cw(kz2))*
     &                     (vec3(i+1,j,kz2)-vec3(i,j,kz2))
c
         vec2(i,j,kz2+1) = vec2(i,j,kz2) +
     &                     dzody*(bv(j)/cw(kz2))*
     &                     (vec3(i,j+1,kz2)-vec3(i,j,kz2))
c
         vec3(i,j,kz2)=vec33(i,j,kz2)
c
c         vec1(i,j,kz2+1)=vec1(i,j,kz2)
c         vec2(i,j,kz2+1)=vec2(i,j,kz2)
       end do
      end do
c
      ENDIF
      RETURN
c
      END

C
C
C---Subroutine Bound250------------------------------------------------
C
      SUBROUTINE BOUND280(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
C
C                         
C       Freestream boundary condition for
C       Wake-Boundary-Layer computations
C 
C       GRESHO and SANI approach
C       Imposes du*/dn=0  + correction for mass balance.
C       This correction is applied by MasWake routine
C       Other componets are obtained by applying omega=0 
C       ----- the routine works only for k=kz2!! ----    
C
C----------------------------------------------------------------------
C
      include 'headers/common.h'
C
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)

c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------

      IF(.not.Ltstep) THEN
c..impose homogeneous Newmann to u*_n  
      do i=ix1,ix2+1
       do j=jy1,jy2+1
        vec3(i,j,kz2)=vec3(i,j,kz2-1)
       enddo
      enddo
c..corect with uniform velocity to conserve mass
c      q1=0.0
c      k=kz2
c      do i=ix1+1,ix2
c       do j=jy1+1,jy2
c        q1=q1+vec3(i,j,k)*(xf(i)-xf(i-1))*(yf(j)-yf(j-1))
c       enddo
c      enddo
c
c      write(6,*) 'qtop=',q1
c
c      surf=(xf(ix2)-xf(ix1))*(yf(jy2)-yf(jy1))
c      wuni=-1.*(q1/surf)
c
c      write(6,*) 'uni w=',wuni
c
c      do i=ix1+1,ix2
c       do j=jy1+1,jy2
c        vec3(i,j,k)=vec3(i,j,k)+wuni
cc        vec3(i,j,k)=0.0
c       enddo
c      enddo
c-------------------------------------------------------------------
c                                                 corrector step      
c-------------------------------------------------------------------
      ELSE
c..u^(n+1)=u* for normal component
      do i=ix1,ix2+1
       do j=jy1,jy2+1
        vec3(i,j,kz2)=vec33(i,j,kz2)
       enddo
      enddo
c..derive bc for other two components from vorticity
      dzodx= DZ/DX
      dzody= DZ/DY
      do i = ix1+1, ix2
       do j = jy1+1, jy2
         vec1(i,j,kz2+1) = vec1(i,j,kz2) +
     &                     dzodx*(au(i)/cw(kz2))*
     &                     (vec3(i+1,j,kz2)-vec3(i,j,kz2))

c         vec2(i,j,kz2+1) = vec2(i,j,kz2) +
c     &                     dzody*(bv(j)/cw(kz2))*
c     &                     (vec3(i,j+1,kz2)-vec3(i,j,kz2))
c          vec1(i,j,kz2+1) = vec1(i,j,kz2)  
          vec2(i,j,kz2+1) = vec2(i,j,kz2) 
       enddo
      enddo 
c
      ENDIF
c
      RETURN
      END

*//
*//
*////////// BOUNDARY CONDITIONS FOR ANALYTICAL SOLUTION  ////////////
*//
*//
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BOUND901(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,tlevel,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
*                                                             *
*         left boundary for analytical solution
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'
c     ...External Functions...
      real     uex, wex
      external uex, wex
c
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)
      real   a, b, lamda, pi, u0, w0
      real   xlb, zu
c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------
      pi = 4.0*atan(1.0)
      u0 = 1.0
      w0 = 1.0
c      a  = 0.5*pi
c      b  = 0.5*pi
      a  = pi
      b  = pi
      lamda = ru1*(a**2+b**2)
      xlb = 0.0
c
      do k=2,kz2
       zu = 0.5*(zf(k)+zf(k-1))
       do j=2,jy2
        vec1(1,j,k) = uex(xlb,zu,tlevel,u0,lamda,a,b) ! it should be 0
        vec2(1,j,k) = - vec2(2,j,k)
        vec3(1,j,k) = 2.0*wex(xlb,zf(k),tlevel,w0,lamda,a,b)
     *               -vec3(2,j,k)
       enddo
      enddo
c
      RETURN
      END
c
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BOUND902(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,tlevel,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
*                                                             *
*          right boundary for analytical solution
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'
c     ...External Functions...
      real     uex, wex
      external uex, wex
c
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)
      real   a, b, lamda, pi, u0, w0
      real   xrb, zu
c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------
      pi = 4.0*atan(1.0)
      u0 = 1.0
      w0 = 1.0
c      a  = 0.5*pi
c      b  = 0.5*pi
      a  = pi
      b  = pi
      lamda = ru1*(a**2+b**2)
      xrb = xf(ix2)
c
      do k=2,kz2
       zu = 0.5*(zf(k)+zf(k-1))
       do j=2,jy2
        vec1(ix2,j,k)   = uex(xrb,zu,tlevel,u0,lamda,a,b) 
        vec2(ix2+1,j,k) = -vec2(ix2,j,k)
        vec3(ix2+1,j,k) = 2.0*wex(xrb,zf(k),tlevel,w0,lamda,a,b)
     *                    -vec3(ix2,j,k)
       enddo
      enddo
c
      RETURN
      END

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BOUND905(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,tlevel,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
*                                                             *
*         bottom boundary for analytical solution
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'
c     ...External Functions...
      real     uex, wex
      external uex, wex
c
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)
      real   a, b, lamda, pi, u0, w0
      real   zbb, xw
c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------
      pi = 4.0*atan(1.0)
      u0 = 1.0
      w0 = 1.0
c      a  = 0.5*pi
c      b  = 0.5*pi
      a  = pi
      b  = pi
      lamda = ru1*(a**2+b**2)
      zbb = 0.0
c
      do j=2,jy2
       do i=2,ix2
        xw = 0.5*(xf(i)+xf(i-1))
c        vec1(i,j,1) = 2.0*uex(xf(i),zbb,tlevel,u0,lamda,a,b)
c     *                - vec1(i,j,2)
        vec2(i,j,1) = - vec2(i,j,2)
        vec3(i,j,1) =   wex(xw,zbb,tlevel,w0,lamda,a,b)
       enddo
      enddo
c
      do j=2,jy2
       do i=1,ix2
c        xw = 0.5*(xf(i)+xf(i-1))
        vec1(i,j,1) = 2.0*uex(xf(i),zbb,tlevel,u0,lamda,a,b)
     *                - vec1(i,j,2)
c        vec2(i,j,1) = - vec2(i,j,2)
c        vec3(i,j,1) =   wex(xw,zbb,tlevel,w0,lamda,a,b)
       enddo
      enddo
c
      RETURN
      END

* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      SUBROUTINE BOUND906(vec1,vec2,vec3,vec11,vec22,vec33,
     &                    xf,yf,zf,tlevel,
     &                    imax,jmax,kmax,ibnd,Ltstep        )
*                                                             *
*         top boundary for analytical solution
*                                                             *
* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

      include 'headers/common.h'
c     ...External Functions...
      real     uex, wex
      external uex, wex
c
      logical Ltstep
      real    vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &        vec3(imax,jmax,kmax),
     &        vec11(imax,jmax,kmax),vec22(imax,jmax,kmax),
     &        vec33(imax,jmax,kmax)
      real    xf(imax),yf(jmax),zf(kmax)
      real   a, b, lamda, pi, u0, w0
      real   ztb, xw
c-------------------------------------------------------------------
c                                                 predictor step  
c-------------------------------------------------------------------
      pi = 4.0*atan(1.0)
      u0 = 1.0
      w0 = 1.0
c      a  = 0.5*pi
c      b  = 0.5*pi
      a  = pi
      b  = pi
      lamda = ru1*(a**2+b**2)
      ztb = zf(kz2)
c
      do j=2,jy2
       do i=2,ix2
        xw = 0.5*(xf(i)+xf(i-1))
c        vec1(i,j,kz2+1) = 2.0*uex(xf(i),ztb,tlevel,u0,lamda,a,b)
c     *                    - vec1(i,j,kz2)
        vec2(i,j,kz2+1) = -vec2(i,j,kz2)
        vec3(i,j,kz2  ) = wex(xw,ztb,tlevel,w0,lamda,a,b)
       enddo
      enddo
c
      do j=2,jy2
       do i=1,ix2
c        xw = 0.5*(xf(i)+xf(i-1))
        vec1(i,j,kz2+1) = 2.0*uex(xf(i),ztb,tlevel,u0,lamda,a,b)
     *                    - vec1(i,j,kz2)
c        vec2(i,j,kz2+1) = -vec2(i,j,kz2)
c        vec3(i,j,kz2  ) = wex(xw,ztb,tlevel,w0,lamda,a,b)
       enddo
      enddo
c
      RETURN
      END

*//
*//
*/////////// PERIODIC BOUNDARY CONDITIONS  //////////////////////////
*//
*//
C
C
C---Subroutine Bound500----------------------I. Balaras -19/2/99-----
C
      SUBROUTINE BOUND500(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
C
C     Imposes periodic boundary conditions to any set of     
C     vec1,vec2,vec3 vectors on any boundary        
C     In this routine the indexes on the surface             
C     the conditions are imposed are shifted      
C     i.e for isx=1 : jby=1 and jey=jy2+1                
C                   : kbz=1 and kez=kz2+1
C
C-----------------------------------------------------------------------
C
      include 'headers/common.h'
C
C-----------------------------------------------------------------------
C
      REAL   vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),
     &       vec3(imax,jmax,kmax)
C
C-----------------------------------------------------------------------
C
*
****  Starting and ending indexes for ibnd 
*
      if (mysize.gt.1.and.(ibnd.eq.5.or.ibnd.eq.6)) then 
       call bound5000(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
       return
      endif  
      ibx=ibegin(ibnd)
      iex=iend(ibnd)
      jby=jbegin(ibnd)
      jey=jend(ibnd)
      kbz=kbegin(ibnd)
      kez=kend(ibnd)
c      write(6,*) ibx, iex, jby, jey, kbz, kez
*
****  Find orientation of boundary : 
*      i=const (isx=1,jsy=0,ksz=0)
*      j=const (isx=0,jsy=1,ksz=0)
*      k=const (isx=0,jsy=0,ksz=1)
*     
      isx=ibx/iex
      jsy=jby/jey
      ksz=kbz/kez     
ccccccccccccccccccccccccccccccccccccccccccccccccc
cccccc  HARDWIRING NECESSARY IN CASE OF 2D FLOW:    VO 1/2006
      if (ibnd.eq.1.or.ibnd.eq.2.or.ibnd.eq.5.or.ibnd.eq.6.and.jey.eq.jby) jsy=0 !i.e. east/west/top/bottom bdrs need fixing
ccccccccccccccccccccccccccccccccccccccccccccccccc      
c      write(6,*) isx,jsy,ksz
*
****  Indices pointing the coresponding inner values
*
      icx=isx*( (ix1+1-ix2)*(ibx-ix1)/ix2 + ix2 - ibx )
      jcy=jsy*( (jy1+1-jy2)*(jby-jy1)/jy2 + jy2 - jby )
      kcz=ksz*( (kz1+1-kz2)*(kbz-kz1)/kz2 + kz2 - kbz )
c      write(6,*) icx, jcy, kcz
*
****  Shift the indices
*
       ibx=ibx-1*(jsy+ksz)
       iex=iex+1*(jsy+ksz)
       jby=jby-1*(isx+ksz)
       jey=jey+1*(isx+ksz)
       kbz=kbz-1*(isx+jsy)
       kez=kez+1*(isx+jsy)
c       write(6,*) ibx,iex,jby,jey,kbz,kez
*
****  Main loop
*
      DO 10 i=ibx,iex
         ii=i+icx
       DO 20 j=jby,jey
          jj=j+jcy
        DO 30 k=kbz,kez
           kk=k+kcz
         vec1(i,j,k) = vec1(ii,jj,kk)
         vec2(i,j,k) = vec2(ii,jj,kk)
         vec3(i,j,k) = vec3(ii,jj,kk)
30      CONTINUE
20     CONTINUE
10    CONTINUE

      RETURN
      END
c  
      subroutine bound5000(vec1,vec2,vec3,imax,jmax,kmax,ibnd)
c     this is the parallel version
c      V. Ovchinnikov 8/10/03. 
      implicit none  
      include 'headers/common.h'
      include 'mpif.h'
      integer imax,jmax,kmax,ibnd
      integer type,comm,ierr,status(MPI_STATUS_SIZE)
      real vec1(imax,jmax,kmax),vec2(imax,jmax,kmax),vec3(imax,jmax,kmax)

      comm=MPI_COMM_HYDRO
c      
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
#ifndef DOUBLE 
      type=MPI_REAL
#endif
      if (ibnd.eq.5) then
        if (myrank.eq.0) then            
          call MPI_RECV(vec1(1,1,1),imax*jmax,type,mysize-1,1,comm,status,ierr)
          call MPI_RECV(vec2(1,1,1),imax*jmax,type,mysize-1,2,comm,status,ierr)
          call MPI_RECV(vec3(1,1,1),imax*jmax,type,mysize-1,3,comm,status,ierr)
        else if (myrank.eq.mysize-1) then
          call MPI_SEND(vec1(1,1,kz2),imax*jmax,type,0,1,comm,ierr)
          call MPI_SEND(vec2(1,1,kz2),imax*jmax,type,0,2,comm,ierr)
          call MPI_SEND(vec3(1,1,kz2),imax*jmax,type,0,3,comm,ierr)
	endif
      elseif (ibnd.eq.6) then 
        if (myrank.eq.mysize-1) then            
          call MPI_RECV(vec1(1,1,kz2+1),imax*jmax,type,0,4,comm,status,ierr)
          call MPI_RECV(vec2(1,1,kz2+1),imax*jmax,type,0,5,comm,status,ierr)
          call MPI_RECV(vec3(1,1,kz2+1),imax*jmax,type,0,6,comm,status,ierr)
	else if (myrank.eq.0) then 
          call MPI_SEND(vec1(1,1,2),imax*jmax,type,mysize-1,4,comm,ierr)
          call MPI_SEND(vec2(1,1,2),imax*jmax,type,mysize-1,5,comm,ierr)
          call MPI_SEND(vec3(1,1,2),imax*jmax,type,mysize-1,6,comm,ierr)
	endif
      endif
      return
      end
	  
C
C///////////////////////////////////////////////////////////////////////
C
      SUBROUTINE LININT(X1, X2, F1, F2, NX1, NX1P, NX2, NX2P, IDUMP)
C
C///////////////////////////////////////////////////////////////////////
C
C    Author        : Andrea Pascarelli
C    Date          : February 23, 1999
C    Last update   : March 04, 1999 
C    Info          :       ***Purpose***
C                   
C                     Linear interpolation
C                         (x1,f1) ----> (x2,f2)
C
C    Input Variables:
C
C    Output variables:
C
C    Local variables:
C
C    External Subroutines:
C  
C///////////////////////////////////////////////////////////////////////
C
C.......................................................................
C     Declaration of formal arguments and local variables
C.......................................................................
C     ...Scalar Arguments...
      INTEGER IDUMP, NX1, NX2, NX1P, NX2P
C     ...Array Arguments...
      REAL X1(1:NX1P), X2(1:NX2P), F1(1:NX1P), F2(1:NX2P)
C     ...Local Scalars...
      INTEGER I, I1, ICOUNT, ISTAT
      REAL RATIO, FI, FIP, XI, XIP, XSTAT
C     ...Intrinsic Functions...
      INTRINSIC REAL
C.......................................................................
C     Executable Statements
C.......................................................................
C     Interpolate profile to mesh x2(j)
C.......................................................................
C.......................................................................
C     Break iteration to locate the right wall-normal indeces 
C.......................................................................
      ICOUNT = 1
      I1     = 1
      DO 10 ISTAT = 1, NX2
        XSTAT = X2(ISTAT)
        DO 20 I = I1, NX1-1
          XI  = X1(I)
          XIP = X1(I+1)
          IF ( (XSTAT.GE.XI) .AND. (XSTAT.LT.XIP) )  THEN
            FI  = F1(I)
            FIP = F1(I+1)
            RATIO = (XSTAT - XI)/(XIP - XI)
C.......................................................................
C           Interpolate f1
C.......................................................................
            F2(ISTAT) = FI + RATIO*(FIP - FI)
c           print*, ICOUNT, I, XI, XSTAT, XIP
            ICOUNT = ICOUNT + 1
            I1     = I
            GOTO 10
          END IF
   20   CONTINUE
        I1 = NX1
        F2(ISTAT) = (REAL(1-IDUMP))*F2(ISTAT-1)
   10 CONTINUE
C
      RETURN
      END
C
C///////////////////////////////////////////////////////////////////////
C///////////////////////////////////////////////////////////////////////
C
      SUBROUTINE BLASIUS(Y, U, V, NY, NPP, NP, NU, UINF, X, XHAT, DSHAT)
C
C///////////////////////////////////////////////////////////////////////
C
C    Author        : Andrea Pascarelli
C    Date          : January 12, 1999
C    Last update   : 
C    Info          :       ***Purpose***
C                   
C                  Boundary layer velocity profiles: Blasius solution.
C    Input Variables:
C
C    Output variables:
C
C    Local variables:
C
C    External Subroutines:
C  
C///////////////////////////////////////////////////////////////////////
C
C.......................................................................
C     Declaration of formal arguments and local variables
C.......................................................................
C     ...Parameters...
      INTEGER ITEMAX, NETA
      REAL EPS, ETAMAX
      PARAMETER (ITEMAX=50, EPS=1.0E-6, ETAMAX = 10.0, NETA = 2000)
C     ...Scalar Arguments...
      INTEGER    NY
      REAL DSHAT, NU, UINF, X, XHAT
C     ...Array Arguments...
      REAL   Y(1:NPP), U(1:NP),  V(1:NP)
C     ...Intrinsic Functions...
      INTRINSIC REAL, SQRT
C     ...Local Scalars...
      INTEGER COUNTER, I, J, J1, JCOUNT, JSTAT
      REAL ALPHA, DETA, DF1N, DSTAR, ETAMIN, HSHAP
      REAL FGUESS, F0, G0, H0, I0, L, DL
      REAL RATIO, THETA, UJ, UJP, X0, YJ, YJP, YSTAT, VJ, VJP
C     ...Local Arrays...
      REAL ETA(1:NETA), F(1:NETA), F1(1:NETA), F2(1:NETA)
C.......................................................................
C     Executable Statements
C.......................................................................
      ETAMIN = 0.0
      DETA   = (ETAMAX - ETAMIN)/ REAL(NETA-1)
C.......................................................................
C     Boundary conditions @ eta = 0
C.......................................................................
      F (1)   = 0.0
      F1(1)   = 0.0
      ETA(1)  = ETAMIN
C.......................................................................
C     Guess the value of f''(0) (greather than 0)
C.......................................................................
      FGUESS = 0.46960
      COUNTER = 0
C.......................................................................
C     Repeat until construct: 
C                           integrate Blasius equation from 0 to etamax 
C.......................................................................
  100 CONTINUE
c        WRITE(UNIT=*,FMT=1) ' Integration #', COUNTER
c    1   FORMAT(/,1X, A, I2)
C
        F2(1) = FGUESS
        DO 10 I = 2, NETA
          ETA(I) = ETA(I-1)+ DETA
          F(I)   = F(I-1)  + F1(I-1)*DETA
          F1(I)  = F1(I-1) + F2(I-1)*DETA
          F2(I)  = F2(I-1) - F(I-1)*F2(I-1)*DETA 
   10   CONTINUE
C.......................................................................
C       Check to see whether f'(infty) = 1 is fulfilled, and if not,
C       improve the guess for f''(0).
C.......................................................................
        ALPHA  = 1.0/SQRT(F1(NETA))
        DF1N   = ABS(1.0-F1(NETA))
        FGUESS = FGUESS*ALPHA
        COUNTER = COUNTER + 1
      IF (.NOT.((COUNTER.GT.ITEMAX).OR.(DF1N.LT.EPS))) GOTO 100
C     Writing   eta, f , f' & f''
c     WRITE(UNIT=*,FMT=*) (I, ETA(I), F(I), F1(I), F2(I), I = 1, NETA) 
C.......................................................................
C     Compute delta*, theta and shape factor
C.......................................................................
      DSTAR = 0.0
      THETA = 0.0
      DO 20 I = 2, NETA
        F0 = 1.0 - F1(I-1)
        G0 = 1.0 - F1(I)
        H0 = (1.0 - F1(I-1))*F1(I-1)
        I0 = (1.0 - F1(I-1))*F1(I-1) 
        DSTAR = DSTAR  + 0.5*(ETA(I)-ETA(I-1))*(F0 + G0)
        THETA = THETA  + 0.5*(ETA(I)-ETA(I-1))*(H0 + I0)
   20 CONTINUE
c      print*,dstar,theta
      X0   = XHAT - 0.5*UINF*((DSHAT/DSTAR)**2)/NU
c      write(6,*) 'x0=',x0
      L    = SQRT(2.0*NU*(X-X0)/UINF)
      DL   = 0.5*SQRT(2.0*NU/(UINF*(X-X0)))
      DSTAR = L * DSTAR 
      THETA = L * THETA 
      HSHAP = DSTAR/THETA 
C
c      print*, ' '
c      print*, ' dstar    ', DSTAR 
c      print*, ' '
c      print*, ' theta    ', THETA 
c      print*, ' '
c      print*, ' Hshape   ', HSHAP
c      print*, ' '
c      print*, ' Re_dstar ', DSTAR*UINF/NU
c      print*, ' '
c      print*, ' Re_theta ', THETA*UINF/NU
c      print*, ' '
C.......................................................................
C     Interpolate selfsimilar profiles to mesh y(j)
C.......................................................................
C.......................................................................
C     Break iteration to locate the right wall-normal indeces 
C.......................................................................
      JCOUNT = 1
      J1     = 1
      DO 30 JSTAT = 2, NY
        DO 40 J = J1, NETA-1
C.......................................................................
C         Convert eta ---> y
C.......................................................................
          YJ  = L*ETA(J)
          YJP = L*ETA(J+1)
          YSTAT = 0.5*(Y(JSTAT)+Y(JSTAT-1))
          IF ( (YSTAT.GE.YJ) .AND. (YSTAT.LT.YJP) )  THEN
C.......................................................................
C           Convert normalized quantities  ---> dimensional quantities
C.......................................................................
            UJ  = UINF*F1(J)
            UJP = UINF*F1(J+1)
            RATIO = (YSTAT - YJ)/(YJP - YJ)
C.......................................................................
C           Interpolate u velocities
C.......................................................................
            U(JSTAT) = min((UJ + RATIO*(UJP - UJ)),uinf)
c           print*, JCOUNT, J, YJ, YSTAT, YJP
            JCOUNT = JCOUNT + 1
            J1     = J
            GOTO 30
          END IF
   40   CONTINUE
        J1 = NETA
         U(jstat)=uinf
   30 CONTINUE
C
      JCOUNT = 1
      J1     = 1
      DO 31 JSTAT = 1, NY
        DO 41 J = J1, NETA-1
C.......................................................................
C         Convert eta ---> y
C.......................................................................
          YJ  = L*ETA(J)
          YJP = L*ETA(J+1)
          YSTAT = Y(JSTAT)
          IF ( (YSTAT.GE.YJ) .AND. (YSTAT.LT.YJP) )  THEN
C.......................................................................
C           Convert normalized quantities  ---> dimensional quantities
C.......................................................................
            VJ  = UINF*DL* (ETA(J)*F1(J)-F(J))
            VJP = UINF*DL* (ETA(J+1)*F1(J+1)-F(J+1))
            RATIO = (YSTAT - YJ)/(YJP - YJ)
C.......................................................................
C           Interpolate v velocities
C.......................................................................
            V(JSTAT) = VJ + RATIO*(VJP - VJ)
c           print*, JCOUNT, J, YJ, YSTAT, YJP
            JCOUNT = JCOUNT + 1
            J1     = J
            GOTO 31
          END IF
 41    CONTINUE
        J1 = NETA
         if (jstat.eq.1) then 
	   V(jstat)=uinf*dl*(eta(neta)*f1(neta)-f(neta))
         else 
	   V(JSTAT) = V(JSTAT-1)
	 endif
 31   CONTINUE
c 
c      write(6,*) u(ny),v(ny)
      RETURN
      END
C
C///////////////////////////////////////////////////////////////////////
C///////////////////////////////////////////////////////////////////////
C
      REAL FUNCTION G(Z,ZMIN,ZMAX,ALPHA)
C
C///////////////////////////////////////////////////////////////////////
C
C    Author        : Andrea Pascarelli
C    Date          : April 5, 1999
C    Last update   : 
C    Info          :       ***Purpose***
C                   
C    Input Variables:
C    Output variables:
C    Local variables:
C    External Subroutines:
C  
C///////////////////////////////////////////////////////////////////////
C     ...Scalar Arguments...    
      REAL ALPHA, Z, ZMIN, ZMAX
C     ...Intrinsic Functions...
      INTRINSIC TANH
C.......................................................................
C     Executable Statements
C.......................................................................
      IF (Z .LT. ZMIN) THEN
        G = 1.0
      ELSE
        G = 1.0 - TANH(ALPHA*(Z-ZMIN)/(ZMAX-ZMIN))/TANH(ALPHA)           
      END IF
C
      RETURN
      END
C
C///////////////////////////////////////////////////////////////////////
