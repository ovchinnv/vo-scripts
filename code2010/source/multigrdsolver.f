c
c
c----Subroutine multig3d--------------------------- 14/12/98 -----
c
      subroutine multig3d(jxc,jyc,jzc,
     &                    in_dx1,in_dx2,in_dx3,in_dx4,
     &                    in_sn1,in_sn2,in_sn3,in_sn4,
     &                    in_sp1,in_sp2,in_sp3,in_sp4,
     &                    in_st1,in_st2,in_st3,in_st4,
     &                    in_av1,in_av2,in_av3,in_av4,
     &                    in_in1,in_in2,in_in3,in_in4,
     &                    f,fi,rhsn,rhs2v,rhs2n,rhs3v,rhs3n,rhs4v,
     &                    pr2,pr3,pr4,
     &                    g11,g22,g33,g11_2,g22_2,g33_2,
     &                    g11_3,g22_3,g33_3,g11_4,g22_4,g33_4,
     &                    giac,
     &                    DIV,DPP,PPO,XF,YF,ZF,
     &                    IM,JM,KM,dtm1,nlevel,                      
     &                    n1,n2,n3,n12,n22,n32,n13,n23,n33,
     &                    n14,n24,n34,icycle                          )
c
c       Driver routine of the 3D  Multigrid Solver
c         --> Arranges the divergence 
c         --> Calls the multigrid solver  (multi)
c         --> Arranges the solution
c         --> Imposes bc on the solution
c
c-----------------------------------------------------------------
c
      include 'headers/common.h'
c
c-----------------------------------------------------------------
c
      INTEGER  jxc(0:4),jyc(0:4),jzc(0:4)
c
      INTEGER  in_dx1(n1,n2,n3),
     &         in_dx2(n12,n22,n32),
     &         in_dx3(n13,n23,n33),
     &         in_dx4(n14,n24,n34)
c
      INTEGER  in_sn1(n1,n2,n3),
     &         in_sn2(n12,n22,n32),
     &         in_sn3(n13,n23,n33),
     &         in_sn4(n14,n24,n34)
c
      INTEGER  in_sp1(n1,n2,n3),
     &         in_sp2(n12,n22,n32),
     &         in_sp3(n13,n23,n33),
     &         in_sp4(n14,n24,n34)
c
      INTEGER  in_st1(n1,n2,n3),
     &         in_st2(n12,n22,n32),
     &         in_st3(n13,n23,n33),
     &         in_st4(n14,n24,n34)
c
      INTEGER  in_av1(n1,n2,n3),
     &         in_av2(n12,n22,n32),
     &         in_av3(n13,n23,n33),
     &         in_av4(n14,n24,n34)
c
      INTEGER  in_in1(n1,n2,n3),
     &         in_in2(n12,n22,n32),
     &         in_in3(n13,n23,n33),
     &         in_in4(n14,n24,n34)
c
      REAL     rhsn(n1,n2,n3),
     &         rhs2v(n12,n22,n32),
     &         rhs2n(n12,n22,n32),
     &         rhs3v(n13,n23,n33),
     &         rhs3n(n13,n23,n33),
     &         rhs4v(n14,n24,n34)
c
      REAL     pr2(0:n12+1,0:n22+1,0:n32+1),
     &         pr3(0:n13+1,0:n23+1,0:n33+1),
     &         pr4(0:n14+1,0:n24+1,0:n34+1) 
c
      REAL     g11(0:n1,n2,n3),g22(n1,0:n2,n3),
     &         g33(n1,n2,0:n3),giac(n1,n2,n3)
c
      REAL     g11_2(0:n12,n22,n32),
     &         g22_2(n12,0:n22,n32),
     &         g33_2(n12,n22,0:n32),
     &         g11_3(0:n13,n23,n33),
     &         g22_3(n13,0:n23,n33),
     &         g33_3(n13,n23,0:n33),
     &         g11_4(0:n14,n24,n34),
     &         g22_4(n14,0:n24,n34),
     &         g33_4(n14,n24,0:n34)  
c
      REAL     fi(0:n1+1,0:n2+1,0:n3+1),f(n1,n2,n3)

      REAL     div(im,jm,km),DPP(im,jm,km),PPO(im,jm,km)
      REAL     zf(km),xf(im),yf(jm)
c
      INTEGER  nxu,nyu,nzu,i,j,k,ibnd
cc      REAL     omega,eps
      REAL     omega(4),eps
c
c-----------------------------------------------------------------
c                                     Set parameters for multigrid
c-----------------------------------------------------------------
c
cc      omega=wmega

cc      pi=4.*atan(1.)
cc      dx_1=xf(2)-xf(1)
cc      dy_1=yf(2)-yf(1)
cc      dz_1=zf(2)-zf(1)

cc      dx_2=2.*dx_1
cc      dy_2=2.*dy_1
cc      dz_2=2.*dz_1

cc      dx_3=2.*dx_2
cc      dy_3=2.*dy_2
cc      dz_3=2.*dz_2

cc      dx_4=2.*dx_3
cc      dy_4=2.*dy_3
cc      dz_4=2.*dz_3

cc      r_1 = ( cos(pi/real(n1)) / dx_1**2 +
cc     &        cos(2.*pi/real(n2)) / dy_1**2 +        
cc     &        cos(pi/real(n3)) / dz_1**2  ) /
cc     &      ( 1./ dx_1**2 +1./ dy_1**2 +1./ dz_1**2 )
c
cc      omega(1)=2./(1.+sqrt(1.-r_1**2))

cc      r_2 = ( cos(pi/real(n12)) / dx_2**2 +
cc     &        cos(2.*pi/real(n22)) / dy_2**2 +        
cc     &        cos(pi/real(n32)) / dz_2**2  ) /
cc     &      ( 1./ dx_2**2 +1./ dy_2**2 +1./ dz_2**2 )
c
cc      omega(2)=2./(1.+sqrt(1.-r_2**2))      

cc      r_3 = ( cos(pi/real(n13)) / dx_3**2 +
cc     &        cos(2.*pi/real(n23)) / dy_3**2 +        
cc     &        cos(pi/real(n33)) / dz_3**2  ) /
cc     &      ( 1./ dx_3**2 +1./ dy_3**2 +1./ dz_3**2 )
c
cc      omega(3)=2./(1.+sqrt(1.-r_3**2)) 

cc      r_4 = ( cos(pi/real(n14)) / dx_4**2 +
cc     &        cos(2.*pi/real(n24)) / dy_4**2 +        
cc     &        cos(pi/real(n34)) / dz_4**2  ) /
cc     &      ( 1./ dx_4**2 +1./ dy_4**2 +1./ dz_4**2 )
c
cc      omega(4)=2./(1.+sqrt(1.-r_4**2)) 

cc      do k=1,4
cc       write(6,*) k,omega(k)
cc      enddo

c..bl (128x64x64)
      omega(1)=1.65
      omega(2)=1.65
      omega(3)=1.65    
      omega(4)=wmega

c..bl (64x48x48)
c      omega(1)=1.8
c      omega(2)=1.75
c      omega(3)=1.70
c      omega(4)=wmega


      eps  =exitmulti
      nxu  =ix2-ix1
      nyu  =jy2-jy1
      nzu  =kz2-kz1
c
c-----------------------------------------------------------------
c                                                Initialise arrays
c-----------------------------------------------------------------
c
      do k=1,nzu
       do j=1,nyu
        do i=1,nxu
         f(i,j,k)=0.0
        enddo
       enddo
      enddo
c
      do k=0,nzu+1
       do j=0,nyu+1
        do i=0,nxu+1
         fi(i,j,k)=0.0
        end do
       end do
      end do
c
c-----------------------------------------------------------------
c                           rearrange divergence past to multigrid
c-----------------------------------------------------------------
c
      do k=kz1+1,kz2
       do j=jy1+1,jy2
        do i=ix1+1,ix2
         f(i-1,j-1,k-1)=div(i,j,k)*giac(i-1,j-1,k-1)
        enddo
       enddo
      enddo
c
c-----------------------------------------------------------------
c                                                 CALL Main Solver
c-----------------------------------------------------------------
c
      call multi(nxu,nyu,nzu,eps,omega,
     &           jxc,jyc,jzc,nlevel,
     &           in_dx1,in_dx2,in_dx3,in_dx4,
     &           in_sn1,in_sn2,in_sn3,in_sn4,
     &           in_sp1,in_sp2,in_sp3,in_sp4,
     &           in_st1,in_st2,in_st3,in_st4,
     &           in_av1,in_av2,in_av3,in_av4,
     &           in_in1,in_in2,in_in3,in_in4,
     &           f,fi,rhsn,rhs2v,rhs2n,rhs3v,rhs3n,rhs4v,
     &           pr2,pr3,pr4,
     &           g11,g22,g33,g11_2,g22_2,g33_2,
     &           g11_3,g22_3,g33_3,g11_4,g22_4,g33_4,
     &           giac,                   
     &           n1,n2,n3,n12,n22,n32,n13,n23,n33,
     &           n14,n24,n34,icycle,infnum                )
c
c-----------------------------------------------------------------
c                          rearrange solution coming from multigrd
c-----------------------------------------------------------------
c
      do k=1,nzu
       do j=1,nyu
        do i=1,nxu
         dpp(ix1+i,jy1+j,kz1+k)=fi(i,j,k)
        enddo
       enddo
      enddo
c
      do k=kz1+1,kz2
       do j=jy1+1,jy2
        do i=ix1+1,ix2
         ppo(i,j,k)=ppo(i,j,k)+dpp(i,j,k)
        enddo
       enddo
      enddo
c
c-----------------------------------------------------------------
c                          set boundary conditions for p and dp
c-----------------------------------------------------------------
c
      do ibnd=1,6
       if( itype(ibnd).eq.500 ) then
         CALL prbcpr(ppo,dpp,im,jm,km,ibnd)
       else
         CALL prbcnw(ppo,dpp,im,jm,km,ibnd)
       endif
      enddo

c
      return
      end
c
c
c---Subroutine PrBcNw-------------------------------------------
c
      SUBROUTINE prbcnw(vec1,vec2,imax,jmax,kmax,ibnd)
c
c       Hmogeneous Nwemann boundary conditions for
c       pressure on any boundary
c
c---------------------------------------------------------------
c
      include 'headers/common.h'
c
c---------------------------------------------------------------
c
      DIMENSION vec1(imax,jmax,kmax),vec2(imax,jmax,kmax)
c
c---------------------------------------------------------------
c
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
         vec1(i,j,k) = vec1(ii,jj,kk)
         vec2(i,j,k) = vec2(ii,jj,kk)
30      CONTINUE
20     CONTINUE
10    CONTINUE
c
      RETURN
      END
c
c
c---Subroutine PrBcPr-------------------------------------------
c
      SUBROUTINE prbcpr(vec1,vec2,imax,jmax,kmax,ibnd)
c
c       Periodic boundary conditions for
c       pressure on any boundary
c
c---------------------------------------------------------------
c
      include 'headers/common.h'
c
c---------------------------------------------------------------
c
      DIMENSION vec1(imax,jmax,kmax),vec2(imax,jmax,kmax)
c
c---------------------------------------------------------------
c
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
*
****  Indexes pointing the coresponding inner values
*
      icx=isx*( (ix1+1-ix2)*(ibx-ix1)/ix2 + ix2 - ibx )
      jcy=jsy*( (jy1+1-jy2)*(jby-jy1)/jy2 + jy2 - jby )
      kcz=ksz*( (kz1+1-kz2)*(kbz-kz1)/kz2 + kz2 - kbz )
*
****  Shift the indexes
*
       ibx=ibx-1*(jsy+ksz)
       iex=iex+1*(jsy+ksz)
       jby=jby-1*(isx+ksz)
       jey=jey+1*(isx+ksz)
       kbz=kbz-1*(isx+jsy)
       kez=kez+1*(isx+jsy)
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
30      CONTINUE
20     CONTINUE
10    CONTINUE
c
      RETURN
      END

c
c
c--------------------------------------------------- 14/12/98 -----
c
      subroutine multi(jx,jy,jz,eps,omega,
     &                 jxc,jyc,jzc,nlevel,
     &                 in_dx1,in_dx2,in_dx3,in_dx4,
     &                 in_sn1,in_sn2,in_sn3,in_sn4,
     &                 in_sp1,in_sp2,in_sp3,in_sp4,
     &                 in_st1,in_st2,in_st3,in_st4,
     &                 in_av1,in_av2,in_av3,in_av4,
     &                 in_in1,in_in2,in_in3,in_in4,
     &                 rhs,fi,rhsn,rhs2v,rhs2n,rhs3v,rhs3n,rhs4v,
     &                 pr2,pr3,pr4,
     &                 g11,g22,g33,g11_2,g22_2,g33_2,
     &                 g11_3,g22_3,g33_3,g11_4,g22_4,g33_4,
     &                 giac,            
     &                 n1,n2,n3,n12,n22,n32,n13,n23,n33,
     &                 n14,n24,n34,icycle,infnum                 )
c
c       3D MultiGrid Solver
c
c------------------------------------------------------------------
c
c
      INTEGER   kss(4) 
c
      INTEGER  jxc(0:4),jyc(0:4),jzc(0:4)
c
      INTEGER  in_dx1(n1,n2,n3),
     &         in_dx2(n12,n22,n32),
     &         in_dx3(n13,n23,n33),
     &         in_dx4(n14,n24,n34)
c
      INTEGER  in_sn1(n1,n2,n3),
     &         in_sn2(n12,n22,n32),
     &         in_sn3(n13,n23,n33),
     &         in_sn4(n14,n24,n34)
c
      INTEGER  in_sp1(n1,n2,n3),
     &         in_sp2(n12,n22,n32),
     &         in_sp3(n13,n23,n33),
     &         in_sp4(n14,n24,n34)
c
      INTEGER  in_st1(n1,n2,n3),
     &         in_st2(n12,n22,n32),
     &         in_st3(n13,n23,n33),
     &         in_st4(n14,n24,n34)
c
      INTEGER  in_av1(n1,n2,n3),
     &         in_av2(n12,n22,n32),
     &         in_av3(n13,n23,n33),
     &         in_av4(n14,n24,n34)
c
      INTEGER  in_in1(n1,n2,n3),
     &         in_in2(n12,n22,n32),
     &         in_in3(n13,n23,n33),
     &         in_in4(n14,n24,n34)
c
      REAL     rhs(n1,n2,n3),
     &         rhsn(n1,n2,n3),
     &         rhs2v(n12,n22,n32),
     &         rhs2n(n12,n22,n32),
     &         rhs3v(n13,n23,n33),
     &         rhs3n(n13,n23,n33),
     &         rhs4v(n14,n24,n34)
c
      REAL     pr2(0:n12+1,0:n22+1,0:n32+1),
     &         pr3(0:n13+1,0:n23+1,0:n33+1),
     &         pr4(0:n14+1,0:n24+1,0:n34+1) 
c
      REAL     g11(0:n1,n2,n3),g22(n1,0:n2,n3),
     &         g33(n1,n2,0:n3),giac(n1,n2,n3)
c
      REAL     g11_2(0:n12,n22,n32),
     &         g22_2(n12,0:n22,n32),
     &         g33_2(n12,n22,0:n32),
     &         g11_3(0:n13,n23,n33),
     &         g22_3(n13,0:n23,n33),
     &         g33_3(n13,n23,0:n33),
     &         g11_4(0:n14,n24,n34),
     &         g22_4(n14,0:n24,n34),
     &         g33_4(n14,n24,0:n34)  
c
      REAL     fi(0:n1+1,0:n2+1,0:n3+1)

      REAL     omega(4)

c
c initialize matrixes
c
      call mul_ini(n1,n2,n3,n12,n22,n32,n13,n23,n33,n14,n24,n34,
     >             rhsn,rhs2v,rhs3v,rhs4v,pr2,pr3,pr4)
c
c definizione dei colpi di SOR per ogni livello
c
      do n=1,nlevel
      kss(n)=10
      if (n.gt.1) kss(n)=20
      end do

      kss(1)=4
      kss(2)=4
      kss(3)=4
      kss(4)=20

c
c---------------------------ciclo iterativo-temporale-------------------------
c
      ktime=0
      resmax=1.
      do while (ktime.lt.50.and.resmax.ge.eps)
      ktime=ktime+1
      do k=0,jzc(2)+1
      do j=0,jyc(2)+1
      do i=0,jxc(2)+1
      pr2(i,j,k)=0.
      end do
      end do
      end do
      do k=0,jzc(3)+1
      do j=0,jyc(3)+1
      do i=0,jxc(3)+1
      pr3(i,j,k)=0.
      end do
      end do
      end do
      do k=0,jzc(4)+1
      do j=0,jyc(4)+1
      do i=0,jxc(4)+1
      pr4(i,j,k)=0.
      end do
      end do
      end do
c
c comincia a calcolare sui livelli a scendere (da fine a coarse)
c
      do 1 n=1,nlevel
c
c
      if (n.eq.1) then        !prima griglia
c
      call solut(n,n1,n2,n3,omega(1),kss,jxc,jyc,jzc,fi,rhs,
     >           g11,g22,g33,
     >           in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)
c
      call residmg(n,n1,n2,n3,jxc,jyc,jzc,fi,rhsn,rhs,
     >             g11,g22,g33,
     >             in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)
c
      call restrict(n,n1,n2,n3,n12,n22,n32,jxc,jyc,jzc,rhsn,rhs2v)
      sm1=0.
      do id1=1,n1
         do jd1=1,n2
            do  kd1=1,n3
               sm1=sm1+rhsn(id1,jd1,kd1)
            end do
         end do
      end do
    
      sm2=0.
      do id1=1,n12
         do jd1=1,n22
            do  kd1=1,n32
               sm2=sm2+rhs2v(id1,jd1,kd1)
            end do
         end do
      end do
c     
c
      else if (n.eq.2) then           !seconda griglia
c
      call solut(n,n12,n22,n32,omega(2),kss,jxc,jyc,jzc,pr2,rhs2v,
     >           g11_2,g22_2,g33_2,
     >           in_dx2,in_sn2,in_sp2,in_st2,in_av2,in_in2)
c
c
               if (n.lt.nlevel) then
c
      call residmg(n,n12,n22,n32,jxc,jyc,jzc,pr2,rhs2n,rhs2v,
     >             g11_2,g22_2,g33_2,
     >             in_dx2,in_sn2,in_sp2,in_st2,in_av2,in_in2)
c
      call restrict(n,n12,n22,n32,n13,n23,n33,jxc,jyc,jzc,rhs2n,rhs3v)
c
               end if
c
c
      else if (n.eq.3) then           !terza griglia
c
c
      call solut(n,n13,n23,n33,omega(3),kss,jxc,jyc,jzc,pr3,rhs3v,
     >           g11_3,g22_3,g33_3,
     >           in_dx3,in_sn3,in_sp3,in_st3,in_av3,in_in3)
c
c
                 if (n.lt.nlevel) then
c
      call residmg(n,n13,n23,n33,jxc,jyc,jzc,pr3,rhs3n,rhs3v,
     >             g11_3,g22_3,g33_3,
     >             in_dx3,in_sn3,in_sp3,in_st3,in_av3,in_in3)
c
      call restrict(n,n13,n23,n33,n14,n24,n34,jxc,jyc,jzc,rhs3n,rhs4v)
c

      sm1=0.
      do id1=1,n13
         do jd1=1,n23
            do  kd1=1,n33
               sm1=sm1+rhs3n(id1,jd1,kd1)
            end do
         end do
      end do
      sm2=0.
      do id1=1,n14
         do jd1=1,n24
            do  kd1=1,n34
               sm2=sm2+rhs4v(id1,jd1,kd1)
            end do
         end do
      end do

                 end if
c
c
      else if (n.eq.4) then           !quarta griglia
c
      call solut(n,n14,n24,n34,omega(4),kss,jxc,jyc,jzc,pr4,rhs4v,
     >           g11_4,g22_4,g33_4,
     >           in_dx4,in_sn4,in_sp4,in_st4,in_av4,in_in4)
c
      end if
c
  1   continue
c-------------------------------------------------------------------
c------------ inizio dei cicli a salire ----------------------------
c-------------------------------------------------------------------
c
      do n=nlevel-1,1,-1
c
c
      if (n.eq.3)      then           !terza griglia
c
c
      call prolong(n,n13,n23,n33,n14,n24,n34,jxc,jyc,jzc,pr3,pr4)
c
      call solut(n,n13,n23,n33,omega(3),kss,jxc,jyc,jzc,pr3,rhs3v,
     >           g11_3,g22_3,g33_3,
     >           in_dx3,in_sn3,in_sp3,in_st3,in_av3,in_in3)
c
c
      else if (n.eq.2) then           !seconda griglia
c
c
      call prolong(n,n12,n22,n32,n13,n23,n33,jxc,jyc,jzc,pr2,pr3)
c
      call solut(n,n12,n22,n32,omega(2),kss,jxc,jyc,jzc,pr2,rhs2v,
     >           g11_2,g22_2,g33_2,
     >           in_dx2,in_sn2,in_sp2,in_st2,in_av2,in_in2)
c
c
      else if (n.eq.1) then           !prima griglia
c
c
      call prolong(n,n1,n2,n3,n12,n22,n32,jxc,jyc,jzc,fi,pr2)
c
      call solut(n,n1,n2,n3,omega(1),kss,jxc,jyc,jzc,fi,rhs,
     >           g11,g22,g33,
     >           in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)
c
c
      end if
c
      end do
c
c calcolo del residuo sulla griglia base (prima griglia)
c
      n=1
      call residmg(n,n1,n2,n3,jxc,jyc,jzc,fi,rhsn,rhs,
     >             g11,g22,g33,
     >             in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)

      call mul_boun(n,n1,n2,n3,jxc,jyc,jzc,g11,g22,g33,fi)

c
c calcolo di resmax
c
      resmax=0.
      irsmax=1
      jrsmax=1
      krsmax=1
      do 6 k=1,jzc(n)
      do 6 j=1,jyc(n)
      do 6 i=1,jxc(n)
c
      aresi=abs(rhsn(i,j,k))
c      resmax=max(resmax,aresi)
       if(aresi .gt. resmax) then
         resmax=aresi
          irsmax =i
          jrsmax =j
          krsmax =k
         endif

c
 6    continue
c
      if(mod(icycle,infnum).eq.0) then
       write(*,*)'ktime, resmax',ktime,resmax,irsmax,jrsmax,krsmax
      endif
      end do
c            
 1234 format(6f12.5)
c
      return
      end
c
c
c----------------------------------------------------------------------
c
      subroutine mul_ini(n1,n2,n3,n12,n22,n32,n13,n23,n33,n14,n24,n34,
     >                   rhs1,rhs2,rhs3,rhs4,pr2,pr3,pr4              )
c
c     Initialize arrays for Multigrid
c
c----------------------------------------------------------------------
c     
      REAL    rhs1(n1,n2,n3),
     &        rhs2(n12,n22,n32),
     &        rhs3(n13,n23,n33),
     &        rhs4(n14,n24,n34)
c
      REAL    pr2(0:n12+1,0:n22+1,0:n32+1),
     &        pr3(0:n13+1,0:n23+1,0:n33+1),
     &        pr4(0:n14+1,0:n24+1,0:n34+1)
c
c inizializzazione delle matrici
c
      do 1 k=1,n3
      do 1 j=1,n2
      do 1 i=1,n1
c
      rhs1(i,j,k)=0.
c
  1   continue
c
      do 2 k=1,n32
      do 2 j=1,n22
      do 2 i=1,n12
c
      rhs2(i,j,k)=0.
c
  2   continue
c
      do 20 k=0,n32+1
      do 20 j=0,n22+1
      do 20 i=0,n12+1
c
      pr2(i,j,k)=0.
c
  20  continue
c
      do 3 k=1,n33
      do 3 j=1,n23
      do 3 i=1,n13
c
      rhs3(i,j,k)=0.
c
  3   continue
c
      do 30 k=0,n33+1
      do 30 j=0,n23+1
      do 30 i=0,n13+1
c
      pr3(i,j,k)=0.
c
  30  continue
c
      do 4 k=1,n34
      do 4 j=1,n24
      do 4 i=1,n14
c
      rhs4(i,j,k)=0.
c
  4   continue
c
      do 40 k=0,n34+1
      do 40 j=0,n24+1
      do 40 i=0,n14+1
c
      pr4(i,j,k)=0.
c
  40  continue
c
      return
      end
c
c
c----------------------------------------------------------------------
c
      subroutine solut(n,i1,j1,k1,omega,kss,jxc,jyc,jzc,pr,rh,
     >                 r11,r22,r33,
     >                 i_dx,i_sn,i_sp,i_st,i_av,i_in)
c
c         smoothing su tutti i livelli con SOR
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      integer  kss(4)
      integer  jxc(0:4),jyc(0:4),jzc(0:4)
      integer  i_dx(i1,j1,k1)
      integer  i_sn(i1,j1,k1)
      integer  i_sp(i1,j1,k1)
      integer  i_st(i1,j1,k1)
      integer  i_av(i1,j1,k1)
      integer  i_in(i1,j1,k1)
c
      REAL     r11(0:i1,j1,k1)
      REAL     r22(i1,0:j1,k1)
      REAL     r33(i1,j1,0:k1)
      REAL     pr(0:i1+1,0:j1+1,0:k1+1)
      REAL     rh(i1,j1,k1)
c
c----------------------------------------------------------------------
c
      ipot=2**(n-1)
      pot=float(ipot)
      ppot1=1./pot
      ppot2=1./pot/pot

c       write(6,*) n,omega


c    
      call mul_boun(n,i1,j1,k1,jxc,jyc,jzc,r11,r22,r33,pr)
c
      do 2 kk=1,kss(n)
c
c smoothing con SOR ( a otto colori)
c
c
      do 1 k=1,jzc(n)
      do 1 j=1,jyc(n)
      do 1 i=1,jxc(n)
c
c.....set periodic bc for pressure if needed
c
c..x-direction
       pr(0,j,k)=real(1-iper)*pr(jxc(n),j,k)+
     &           real(iper)*pr(0,j,k)
c
       pr(jxc(n)+1,j,k)=real(1-iper)*pr(1,j,k)+
     &                  real(iper)*pr(jxc(n)+1,j,k)
c..y-direction
       pr(i,0,k)=real(1-jper)*pr(i,jyc(n),k)+
     &           real(jper)*pr(i,0,k)
c
       pr(i,jyc(n)+1,k)=real(1-jper)*pr(i,1,k)+
     &                  real(jper)*pr(i,jyc(n)+1,k)
c..z-direction
       pr(i,j,0)=real(1-kper)*pr(i,j,jzc(n))+
     &           real(kper)*pr(i,j,0)
c
       pr(i,j,jzc(n)+1)=real(1-kper)*pr(i,j,1)+
     &                  real(kper)*pr(i,j,jzc(n)+1)
c
c-------------------- parte destra------------------------
c
       res_dx=(r11(i,j,k)*(pr(i+1,j,k)-pr(i,j,k))
     >  )*i_dx(i,j,k)*ppot1
c
       den_dx=i_dx(i,j,k)*r11(i,j,k)*ppot2
c
c-------------------- parte sinistra-------------------------
c
       res_sn=(r11(i-1,j,k)*(pr(i,j,k)-pr(i-1,j,k))
     > )*i_sn(i,j,k)*ppot1
c
      den_sn=i_sn(i,j,k)*r11(i-1,j,k)*ppot2
c
c-------------------- parte sopra-----------------------------
c
       res_sop=(r22(i,j,k)*(pr(i,j+1,k)-pr(i,j,k))
     > )*i_sp(i,j,k)*ppot1
c
       den_sop=i_sp(i,j,k)*r22(i,j,k)*ppot2
c
c-------------------- parte sotto-----------------------------
c
      res_sot=(r22(i,j-1,k)*(pr(i,j,k)-pr(i,j-1,k))
     >  )*i_st(i,j,k)*ppot1
c
      den_sot=i_st(i,j,k)*r22(i,j-1,k)*ppot2
c
c-------------------- parte avanti  -------------------------
c
      res_av=(r33(i,j,k)*(pr(i,j,k+1)-pr(i,j,k))
     > )*i_av(i,j,k)*ppot1
c
      den_av=i_av(i,j,k)*r33(i,j,k)*ppot2
c
c-------------------- parte indietro -------------------------
c
      res_ind=(r33(i,j,k-1)*(pr(i,j,k)-pr(i,j,k-1))
     > )*i_in(i,j,k)*ppot1
c
      den_ind=i_in(i,j,k)*r33(i,j,k-1)*ppot2
c
c------------------- calcolo di pr nuovo con SOR----------------
c
      resi=
     > (res_dx  - res_sn  + 
     >  res_sop - res_sot +
     >  res_av  - res_ind )*ppot1  -rh(i,j,k)
c
c      write(78,*)i,j,k,resi
      den= den_dx + den_sn + den_sop + den_sot + den_av + den_ind
c
      pr(i,j,k)=pr(i,j,k)+omega*resi/den
c     pr(i,j,k)=
c
 1    continue   ! done one complete sweep 

      call mul_boun(n,i1,j1,k1,jxc,jyc,jzc,r11,r22,r33,pr)

 2    continue   ! end of a pseudo time iteration 
c
      return
      end      
c
c
c----------------------------------------------------------------------
c
      subroutine mul_boun(n,i1,j1,k1,jxc,jyc,jzc,r11,r22,r33,pr)
c
c       aggiornamento condizioni al contorno 
c       per multigrid
c
c
c----------------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------------
c
      integer    jxc(0:4),jyc(0:4),jzc(0:4)
      REAL       r11(0:i1,j1,k1)
      REAL       r22(i1,0:j1,k1)
      REAL       r33(i1,j1,0:k1)
      REAL       pr(0:i1+1,0:j1+1,0:k1+1)
c
c----------------------------------------------------------------------
c
      ip=iper
      jp=jper
      kp=kper
c
c calcola le condizioni al contorno per la pressione
c
      if (n.eq.1) then    ! griglia iniziale 
       an=1.
      else                ! altre griglie
       an=0.
      end if

c
c calcolo della pressione nei punti fantasma sinistra e destra
c
      if(ip.eq.1) then
       do k=1,jzc(n)
        do j=1,jyc(n)
        pr(0,j,k)=pr(1,j,k) 
        pr(jxc(n)+1,j,k)=pr(jxc(n),j,k)
        end do
      end do
      else
       do k=1,jzc(n)
        do j=1,jyc(n)
         pr(0,j,k)=pr(jxc(n),j,k)
         pr(jxc(n)+1,j,k)=pr(1,j,k)
        end do
       end do
      end if
c     
c calcolo della pressione nei punti fantasma sotto e sopra
c
      if (jp.eq.1) then
       do i=1,jxc(n)
        do k=1,jzc(n)
         pr(i,0,k)=pr(i,1,k) 
         pr(i,jyc(n)+1,k)=pr(i,jyc(n),k) 
        end do
       end do
      else
       do i=1,jxc(n)
        do k=1,jzc(n)
         pr(i,       0,k)=pr(i,jyc(n),k)
         pr(i,jyc(n)+1,k)=pr(i,1     ,k)
        end do
       end do
      end if
c     
c calcolo della pressione nei punti fantasma indietro e avanti
c
      if (kp.eq.1) then
       do i=1,jxc(n)
        do j=1,jyc(n)
         pr(i,j,0)=pr(i,j,1) 
         pr(i,j,jzc(n)+1)=pr(i,j,jzc(n))
        end do
       end do
      else   
       do i=1,jxc(n)
        do j=1,jyc(n)
         pr(i,j,       0)=pr(i,j,jzc(n))
         pr(i,j,jzc(n)+1)=pr(i,j,     1)
        end do
       end do
      end if
c
      return
      end
c
c
c----------------------------------------------------------------------
c
      subroutine residmg(n,i1,j1,k1,jxc,jyc,jzc,pr,rh1,rh,
     >                   r11,r22,r33,
     >                   i_dx,i_sn,i_sp,i_st,i_av,i_in     )
c
c calcolo dei residui su griglia attuale 
c
c----------------------------------------------------------------------
c
      integer  jxc(0:4),jyc(0:4),jzc(0:4)
      integer  i_dx(i1,j1,k1)
      integer  i_sn(i1,j1,k1)
      integer  i_sp(i1,j1,k1)
      integer  i_st(i1,j1,k1)
      integer  i_av(i1,j1,k1)
      integer  i_in(i1,j1,k1)
c
      REAL     r11(0:i1,j1,k1)
      REAL     r22(i1,0:j1,k1)
      REAL     r33(i1,j1,0:k1)
      REAL     pr(0:i1+1,0:j1+1,0:k1+1)
      REAL     rh(i1,j1,k1)
      REAL     rh1(i1,j1,k1)
c
      ipot=2**(n-1)
      pot=float(ipot)
c
      do 1 k=1,jzc(n)
      do 1 j=1,jyc(n)
      do 1 i=1,jxc(n)
c
c-------------------- parte destra------------------------
c
       res_dx=(r11(i,j,k)*(pr(i+1,j,k)-pr(i,j,k))
     & )*i_dx(i,j,k)/pot
c
c-------------------- parte sinistra-------------------------
c
       res_sn=(r11(i-1,j,k)*(pr(i,j,k)-pr(i-1,j,k))
     & )*i_sn(i,j,k)/pot
c
c-------------------- parte sopra-----------------------------
c
       res_sop=(r22(i,j,k)*(pr(i,j+1,k)-pr(i,j,k))
     & )*i_sp(i,j,k)/pot
c
c-------------------- parte sotto-----------------------------
c
      res_sot=(r22(i,j-1,k)*(pr(i,j,k)-pr(i,j-1,k))
     &  )*i_st(i,j,k)/pot
c
c-------------------- parte avanti  -------------------------
c
      res_av=(r33(i,j,k)*(pr(i,j,k+1)-pr(i,j,k))
     & )*i_av(i,j,k)/pot
c
c-------------------- parte indietro -------------------------
c
      res_ind=(r33(i,j,k-1)*(pr(i,j,k)-pr(i,j,k-1))
     & )*i_in(i,j,k)/pot
c
c------------------- residuo finale---------------------------
c
       rh1(i,j,k)=
     > (res_dx  - res_sn  + 
     >  res_sop - res_sot +
     >  res_av  - res_ind )/pot - rh(i,j,k)
c
 1    continue
c
      return
      end
c
c
c----------------------------------------------------------------------
c
      subroutine restrict(n,i1,j1,k1,i2,j2,k2,jxc,jyc,jzc,rh,rh1)
c
c       calcolo di residuo su griglia piu' bassa
c       operazione di restriction con valore medio
c
c----------------------------------------------------------------------
c
      integer  jxc(0:4),jyc(0:4),jzc(0:4)
      real     rh(i1,j1,k1) 
      real     rh1(i2,j2,k2)
c
      do 1 k=1,jzc(n+1)
      do 1 j=1,jyc(n+1)
      do 1 i=1,jxc(n+1)
c
      id=2*i
      jd=2*j
      kd=2*k
c
      rh1(i,j,k)=rh(id,jd,kd)+rh(id-1,jd,kd)+rh(id,jd-1,kd)+
     >           rh(id-1,jd-1,kd)+rh(id,jd,kd-1)+rh(id-1,jd,kd-1)
     >           +rh(id,jd-1,kd-1)+rh(id-1,jd-1,kd-1)
c
      rh1(i,j,k)=-rh1(i,j,k)/8.
c
  1   continue
c
      return
      end
c
c
c----------------------------------------------------------------------
c
      subroutine prolong(n,i1,j1,k1,i2,j2,k2,jxc,jyc,jzc,prf,pr)
c
c     operazione di prolongation da griglia rada a griglia fine
c     con interpolazione bilineare
c
c----------------------------------------------------------------------
c
      integer  jxc(0:4),jyc(0:4),jzc(0:4)
      real     prf(0:i1+1,0:j1+1,0:k1+1)
      real     pr(0:i2+1,0:j2+1,0:k2+1)
c
c----------------------------------------------------------------------
c
c     prf pressione su griglia fine
c     pr pressione su griglia rada
c
      do 1 k=0,jzc(n+1)
      do 1 j=0,jyc(n+1)
      do 1 i=0,jxc(n+1)
c
      igf=2*i
      jf=2*j
      kf=2*k
c
      a1=27.
      a2=9.
      a3=9.
      a4=9.
      a5=3.
      a6=3.
      a7=3.
      a8=1.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
c
      igf=2*i+1
      jf=2*j
      kf=2*k
c
      a1=9.
      a2=27.
      a3=3.
      a4=3.
      a5=9.
      a6=9.
      a7=1.
      a8=3.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
c
      igf=2*i
      jf=2*j+1
      kf=2*k
c
      a1=9.
      a2=3.
      a3=27.
      a4=3.
      a5=9.
      a6=1.
      a7=9.
      a8=3.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
c
      igf=2*i
      jf=2*j
      kf=2*k+1
c
      a1=9.
      a2=3.
      a3=3.
      a4=27.
      a5=1.
      a6=9.
      a7=9.
      a8=3.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
c
      igf=2*i+1
      jf=2*j+1
      kf=2*k
c
      a1=3.
      a2=9.
      a3=9.
      a4=1.
      a5=27.
      a6=3.
      a7=3.
      a8=9.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
c
      igf=2*i+1
      jf=2*j
      kf=2*k+1
c
      a1=3.
      a2=9.
      a3=1.
      a4=9.
      a5=3.
      a6=27.
      a7=3.
      a8=9.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
c
      igf=2*i+1
      jf=2*j+1
      kf=2*k+1
c
      a1=1.
      a2=3.
      a3=3.
      a4=3.
      a5=9.
      a6=9.
      a7=9.
      a8=27.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
c
      igf=2*i
      jf=2*j+1
      kf=2*k+1
c
      a1=3.
      a2=1.
      a3=9.
      a4=9.
      a5=3.
      a6=3.
      a7=27.
      a8=9.
c
      prf(igf,jf,kf)=prf(igf,jf,kf)+
     >  (  a1*pr(i  ,j  ,k  )+
     >     a2*pr(i+1,j  ,k  )+
     >     a3*pr(i  ,j+1,k  )+
     >     a4*pr(i  ,j  ,k+1)+
     >     a5*pr(i+1,j+1,k  )+
     >     a6*pr(i+1,j  ,k+1)+
     >     a7*pr(i  ,j+1,k+1)+
     >     a8*pr(i+1,j+1,k+1) )/64.
c
 1    continue
c
      return
      end
