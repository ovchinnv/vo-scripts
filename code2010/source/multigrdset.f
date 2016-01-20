c
c
c---Subroutine MultiSetup-------------------E. Balaras 16/12/98----
c
      subroutine MultiSetup(nlevel,jxc,jyc,jzc,
     &                      xu,yv,zw,xmg,ymg,zmg,
     &                      g11,g22,g33,giac,
     &                      g11_2,g22_2,g33_2,             
     &                      g11_3,g22_3,g33_3, 
     &                      g11_4,g22_4,g33_4,
     &                      in_dx1,in_dx2,in_dx3,in_dx4,
     &                      in_sn1,in_sn2,in_sn3,in_sn4,
     &                      in_sp1,in_sp2,in_sp3,in_sp4,
     &                      in_st1,in_st2,in_st3,in_st4,
     &                      in_av1,in_av2,in_av3,in_av4,
     &                      in_in1,in_in2,in_in3,in_in4,   
     &                      n1,n2,n3,
     &                      n12,n22,n32,
     &                      n13,n23,n33,
     &                      n14,n24,n34,   
     &                      im,jm,km    )
c
c     Sets up the multigrid solver
c
c------------------------------------------------------------------
c
      include 'headers/common.h'
c
c------------------------------------------------------------------
c
      INTEGER jxc(0:4),jyc(0:4),jzc(0:4)
c
      REAL    xu(im),yv(jm),zw(km)
c
      REAL    xmg(-8:n1+8),ymg(-8:n2+8),zmg(-8:n3+8),
     &        g11(0:n1,n2,n3), g22(n1,0:n2,n3), g33(n1,n2,0:n3),
     &        giac(n1,n2,n3)
c
      REAL    g11_2(0:n12,n22,n32),
     &        g22_2(n12,0:n22,n32),
     &        g33_2(n12,n22,0:n32),
     &        g11_3(0:n13,n23,n33),
     &        g22_3(n13,0:n23,n33),
     &        g33_3(n13,n23,0:n33),
     &        g11_4(0:n14,n24,n34),
     &        g22_4(n14,0:n24,n34),
     &        g33_4(n14,n24,0:n34) 
c
      INTEGER in_dx1(n1,n2,n3),in_dx2(n12,n22,n32),
     &        in_dx3(n13,n23,n33),in_dx4(n14,n24,n34)
c
      INTEGER in_sn1(n1,n2,n3),in_sn2(n12,n22,n32),
     &        in_sn3(n13,n23,n33),in_sn4(n14,n24,n34)
c
      INTEGER in_sp1(n1,n2,n3),in_sp2(n12,n22,n32),
     &        in_sp3(n13,n23,n33),in_sp4(n14,n24,n34)
c
      INTEGER in_st1(n1,n2,n3),in_st2(n12,n22,n32),
     &        in_st3(n13,n23,n33),in_st4(n14,n24,n34)
c
      INTEGER in_av1(n1,n2,n3),in_av2(n12,n22,n32),
     &        in_av3(n13,n23,n33),in_av4(n14,n24,n34)
c
      INTEGER in_in1(n1,n2,n3),in_in2(n12,n22,n32),
     &        in_in3(n13,n23,n33),in_in4(n14,n24,n34)
c
      INTEGER nlevel
      INTEGER nxu,nyu,nzu
      REAL    lenx,leny,lenz
c
c------------------------------------------------------------------
c
      WRITE(6,*) '*..Multigrid Initialization...'

c..number of unknowns
      nxu=ix2-ix1  !n1
      nyu=jy2-jy1  !n2
      nzu=kz2-kz1  !n3
c      write(6,*) 'n1=',n1,'nxu=',nxu
c...coordinates of cell faces for multigrid
      do i=0,nxu
       xmg(i)=xu(i+1)
c       write(6,*) i, xmg(i)
      end do
      do j=0,nyu
       ymg(j)=yv(j+1)
      end do
      do k=0,nzu
       zmg(k)=zw(k+1)
      end do
c...add extra cells in i-dir
      lenx=xu(ix2)-xu(ix1)
c      write(6,*) ix2,ix1,xu(ix2),xu(ix1),lenx
      do i=1,8
       xmg(-i   )=xu(ix2-i)-lenx
       xmg(nxu+i)=xu(i+1  )+lenx
c       write(6,*) -i,xmg(-i),nxu+i,xmg(nxu+i)
      end do
c...add extra cells in j-dir
      leny=yv(jy2)-yv(jy1)
      do j=1,8
       ymg(   -j)=yv(jy2-j)-leny
       ymg(nyu+j)=yv(j+1  )+leny
      enddo
c...add extra cells in k-dir
      lenz=zw(kz2)-zw(kz1)
      do k=1,8
       zmg(   -k)=zw(kz2-k)-lenz 
       zmg(nzu+k)=zw(k+1  )+lenz 
      enddo
c      write(6,*) xmg
c
      call metrica(xmg,ymg,zmg,g11,g22,g33,giac,n1,n2,n3)
c
      call indy(nlevel,jxc,jyc,jzc) 
c
      call mul_met(xmg,ymg,zmg,nlevel,jxc,jyc,jzc,
     &             g11_2,g22_2,g33_2,             
     &             g11_3,g22_3,g33_3, 
     &             g11_4,g22_4,g33_4, 
     &             n1,n2,n3,
     &             n12,n22,n32,
     &             n13,n23,n33,
     &             n14,n24,n34                     )
c
      call wall(nlevel,jxc,jyc,jzc,
     &          in_dx1,in_dx2,in_dx3,in_dx4,
     &          in_sn1,in_sn2,in_sn3,in_sn4,
     &          in_sp1,in_sp2,in_sp3,in_sp4,
     &          in_st1,in_st2,in_st3,in_st4,
     &          in_av1,in_av2,in_av3,in_av4,
     &          in_in1,in_in2,in_in3,in_in4,   
     &          n1,n2,n3,                   
     &          n12,n22,n32,n13,n23,n33,
     &          n14,n24,n34                  )
c
      write(6,*) 'number of levels=',nlevel
      do n=1,nlevel
      write(*,*)jxc(n),jyc(n),jzc(n)
      end do
c
      return
      end
c
c
c----------------------------------------------------------------- 
c   
      SUBROUTINE metrica(X, Y, Z, G11, G22, G33, GIAC, N1, N2, N3)
c
c        calculates metric coefficients    
c
c        J(-1)ztz -----> ztz(k)   
c
c        tensore metrico controvariante * J(-1) 
c
c        J(-1)*g11  : g11(i,j,k)
c        J(-1)*g22  : g22(i,j,k)
c        J(-1)*g33  : g33(i,j,k) 
c
c----------------------------------------------------------------
c
      include 'headers/common.h'
c
c----------------------------------------------------------------
c     ...Scalar Arguments...
      INTEGER N1, N2, N3
C     ...Array Arguments...
      REAL X(-8:N1+8), Y(-8:N2+8), Z(-8:N3+8)
      REAL G11(0:N1,N2,N3), G22(N1,0:N2,N3), G33(N1,N2,0:N3)
      REAL GIAC(N1, N2, N3)
C     ...Local Scalars...
      INTEGER IIP, JP, KP
      INTEGER JX, JY, JZ

c...pressure bc flags
      iip = iper
      jp = jper
      kp = kper
c...beging and end indexes
      jx = ix2-ix1
      jy = jy2-jy1
      jz = kz2-kz1
c
c piani a csi costante
c
      do 1 i=0,jx
      do 1 j=1,jy
      do 1 k=1,jz
c
      if (i.eq.0.and.iip.eq.1) then
c
      x2 = x(i+2)
      x1 = x(i+1)
      x0 = x(i)
c
      y2 = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
      y1 = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
      y0 = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
c
      z2 = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
      z1 = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
      z0 = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
c
      xcsi = 0.5 * ( -3.*x0 + 4.*x1 - x2 )
      ycsi = 0.5 * ( -3.*y0 + 4.*y1 - y2 )
      zcsi = 0.5 * ( -3.*z0 + 4.*z1 - z2 )
c
      else if (i.eq.jx.and.iip.eq.1) then
c
      x2 = x(i-2)
      x1 = x(i-1)
      x0 = x(i)
c
      y2 = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
      y1 = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
      y0 = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
c
      z2 = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
      z1 = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
      z0 = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
c
      xcsi = 0.5 * ( 3.*x0 - 4.*x1 + x2 )
      ycsi = 0.5 * ( 3.*y0 - 4.*y1 + y2 )
      zcsi = 0.5 * ( 3.*z0 - 4.*z1 + z2 )
c
      else
c
      xdx = x(i+1)
      xsn = x(i-1)
c
      ydx = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
      ysn = 0.25*(y(j)+y(j)+y(j-1)+y(j-1))
c
      zdx = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
      zsn = 0.25*(z(k)+z(k-1)+z(k)+z(k-1))
c
      xcsi = 0.5*(xdx-xsn)
      ycsi = 0.5*(ydx-ysn)
      zcsi = 0.5*(zdx-zsn)
c
      end if
c
      xeta = 0.0
c
      yeta = y(j) - y(j-1) 
c
      zeta = 0.0
c
      xzet = 0.0
c
      yzet = 0.0
c
      zzet = z(k) - z(k-1)
c
      a_giac=     xcsi*(yeta*zzet-yzet*zeta)-
     >            xeta*(ycsi*zzet-yzet*zcsi)+
     >            xzet*(ycsi*zeta-yeta*zcsi)
c
      a_etx = yzet*zcsi - ycsi*zzet 
      a_ety = xcsi*zzet - xzet*zcsi
      a_etz = xzet*ycsi - xcsi*yzet
c
      a_ztx=ycsi*zeta-yeta*zcsi
      a_zty=xeta*zcsi-xcsi*zeta
      a_ztz=xcsi*yeta-xeta*ycsi
c
      g11(i,j,k)=((yeta*zzet - yzet*zeta)**2)/a_giac

 1    continue
c
c---------------------------------------------------------
c
c piani a eta costante
c
      do 2 i=1,jx
      do 2 j=0,jy
      do 2 k=1,jz
c
      if      (j.eq.0.and.jp.eq.1)  then
c
      x2 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x1 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x0 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
c
      y2 = y(j+2)
      y1 = y(j+1)
      y0 = y(j)
c
      z2 = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
      z1 = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
      z0 = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
c
      xeta = 0.5 * ( -3.*x0 + 4.*x1 - x2 )
      yeta = 0.5 * ( -3.*y0 + 4.*y1 - y2 )
      zeta = 0.5 * ( -3.*z0 + 4.*z1 - z2 )
c
      else if (j.eq.jy.and.jp.eq.1) then
c
      x2 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x1 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x0 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
c
      y2 = y(j-2)
      y1 = y(j-1)
      y0 = y(j)
c
      z2 = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
      z1 = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
      z0 = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
c
      xeta = 0.5 * ( 3.*x0 - 4.*x1 + x2 )
      yeta = 0.5 * ( 3.*y0 - 4.*y1 + y2 )
      zeta = 0.5 * ( 3.*z0 - 4.*z1 + z2 )
c
      else
c
      xsop = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      xsot = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
c
      ysop = y(j+1)
      ysot = y(j-1)
c
      zsop = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
      zsot = 0.25*(z(k)+z(k-1)+z(k-1)+z(k))
c
      xeta = 0.5*(xsop-xsot)
      yeta = 0.5*(ysop-ysot)
      zeta = 0.5*(zsop-zsot)
c     
      end if
c
      yzet = 0.0
c
      zcsi = 0.0
c
      ycsi = 0.0
c
      zzet = z(k) - z(k-1) 
c
      xcsi = x(i) - x(i-1) 
c
      xzet = 0.0
c
      a_csx = yeta*zzet - yzet*zeta
      a_csy = xzet*zeta - xeta*zzet
      a_csz = xeta*yzet - xzet*yeta
c
      a_ztx=ycsi*zeta-yeta*zcsi
      a_zty=xeta*zcsi-xcsi*zeta
      a_ztz=xcsi*yeta-xeta*ycsi
c
      a_giac=xcsi*(yeta*zzet-yzet*zeta)-
     >       xeta*(ycsi*zzet-yzet*zcsi)+
     >       xzet*(ycsi*zeta-yeta*zcsi)

      g22(i,j,k)=((xcsi*zzet - xzet*zcsi)**2)/a_giac

 2    continue
c
c-----------------------------------------------------------
c
c piani a zita costante
c
      do 3 i=1,jx
      do 3 j=1,jy
      do 3 k=0,jz
c
      if (k.eq.0.and.kp.eq.1) then
c
      x2 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x1 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x0 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
c
      y2 = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
      y1 = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
      y0 = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
c
      z2 = z(k+2)
      z1 = z(k+1)
      z0 = z(k)
c
      xzet = 0.5 * ( -3.*x0 + 4.*x1 - x2 )
      yzet = 0.5 * ( -3.*y0 + 4.*y1 - y2 )
      zzet = 0.5 * ( -3.*z0 + 4.*z1 - z2 )
c
      else if (k.eq.jz.and.kp.eq.1) then
c
      x2 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x1 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      x0 = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
c
      y2 = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
      y1 = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
      y0 = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
c
      z2 = z(k-2)
      z1 = z(k-1)
      z0 = z(k)
c
      xzet = 0.5 * ( 3.*x0 - 4.*x1 + x2 )
      yzet = 0.5 * ( 3.*y0 - 4.*y1 + y2 )
      zzet = 0.5 * ( 3.*z0 - 4.*z1 + z2 )
c
      else
c
      xav = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
      xdt = 0.25*(x(i)+x(i)+x(i-1)+x(i-1))
c
      yav = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
      ydt = 0.25*(y(j)+y(j-1)+y(j-1)+y(j))
c
      zav = z(k+1)
      zdt = z(k-1)
c
      xzet = 0.5*(xav-xdt)
      yzet = 0.5*(yav-ydt)
      zzet = 0.5*(zav-zdt)
c
      end if
c
      ycsi = 0.0
c
      zeta = 0.0
c
      yeta = y(j) - y(j-1) 
c
      zcsi = 0.0
c     
      xcsi = x(i) - x(i-1) 
c
      xeta = 0.0
c
      a_csx = yeta*zzet - yzet*zeta
      a_csy = xzet*zeta - xeta*zzet
      a_csz = xeta*yzet - xzet*yeta
c
      a_etx = yzet*zcsi - ycsi*zzet 
      a_ety = xcsi*zzet - xzet*zcsi
      a_etz = xzet*ycsi - xcsi*yzet
c
      a_giac=     xcsi*(yeta*zzet-yzet*zeta)-
     >            xeta*(ycsi*zzet-yzet*zcsi)+
     >            xzet*(ycsi*zeta-yeta*zcsi)

      g33(i,j,k)=((xcsi*yeta-xeta*ycsi)**2)/a_giac
c
  3   continue
c
c-------------------------------------------------------------
c
c giacobiano J(-1) : giac(i,j,k) definito a centro cella
c
C.......................................................................
C calcolo del jacobiano a centro cella
C
      do 4 i=1,jx
      do 4 j=1,jy
      do 4 k=1,jz
c
      xcsi = x(i) - x(i-1)
c
      ycsi = 0.0
c
      zcsi = 0.0
c
      xeta = 0.0
c
      yeta = y(j) - y(j-1)
c
      zeta = 0.0
c
      xzet = 0.0
c
      yzet = 0.0
c
      zzet = z(k) - z(k-1)
c
      giac(i,j,k)=xcsi*(yeta*zzet-yzet*zeta)-
     >            xeta*(ycsi*zzet-yzet*zcsi)+
     >            xzet*(ycsi*zeta-yeta*zcsi)
c
 4    continue
c
      return
      end
c
c
c-------------------------------------------------------------
c
      subroutine indy(nlevel,jxc,jyc,jzc)
c
c      finds the nember of grid levels and the number of
c      grid points for each one
c
c-------------------------------------------------------------
c
       include 'headers/common.h'
c
c-------------------------------------------------------------
c
       INTEGER  nlevel
       INTEGER  jxc(0:4),jyc(0:4),jzc(0:4)
c
c-------------------------------------------------------------
c...beging and end indexes
      jx=ix2-ix1
      jy=jy2-jy1
      jz=kz2-kz1

c...number of levels in x
      mez1=jx
      mez3=mez1
      mez2=mez1
c
      k=0
      do while (mez3.eq.mez1.and.k.lt.20)
c
      mez1=mez2
      k=k+1
      mez2=mez1/2
      mez3=2*mez2
c      write(6,*) mez1, mez3, k
c
      end do
      mezx=mez1
      kx=k
      if (mezx.eq.1) kx=kx-1
c      write(6,*) k

c...number of levels in y
      mez1=jy
      mez3=mez1
      mez2=mez1
c
      k=0
      do while (mez3.eq.mez1.and.k.lt.4)
c
      mez1=mez2
      k=k+1
      mez2=mez1/2
      mez3=2*mez2
c      write(6,*) mez1, mez2, mez3, k
c
      end do
      mezy=mez1
      ky=k
      if (mezy.eq.1) ky=ky-1
      

c...number of levels in z
      mez1=jz
      mez3=mez1
      mez2=mez1
c
      k=0
      do while (mez3.eq.mez1.and.k.lt.4)
c
      mez1=mez2
      k=k+1
      mez2=mez1/2
      mez3=2*mez2
c
      end do
      mezy=mez1
      kz=k
      if (mezy.eq.1) kz=kz-1
c
      nlevel=min(kx,ky,kz)

c...number of cells for each level
      jxc(0)=0
      jyc(0)=0
      jzc(0)=0
c
      jzc(1)=jz
      do j=2,nlevel
      jzc(j)=jzc(j-1)/2
      end do
c
      jyc(1)=jy
      do j=2,nlevel
      jyc(j)=jyc(j-1)/2
      end do
c
      jxc(1)=jx
      do i=2,nlevel
      jxc(i)=jxc(i-1)/2
      end do
c
      return
      end      
c
c
c-------------------------------------------------------------
c
      subroutine mul_met(x,y,z,nlevel,jxc,jyc,jzc,
     &                   g11_2,g22_2,g33_2,             
     &                   g11_3,g22_3,g33_3, 
     &                   g11_4,g22_4,g33_4, 
     &                   n1,n2,n3,           
     &                   n12,n22,n32,
     &                   n13,n23,n33,
     &                   n14,n24,n34                 )
c
c calcola i coefficienti per il multigrid e li mette 
c in matrice coe (dimensionata per 4 livelli) 
c
c-------------------------------------------------------------
c
      include 'headers/common.h'
c
c-------------------------------------------------------------
c     ...Scalar Arguments...
      INTEGER N1, N2, N3
C     ...Array Arguments...
      INTEGER JXC(0:4),JYC(0:4),JZC(0:4)
      REAL X(-8:N1+8), Y(-8:N2+8), Z(-8:N3+8)
c
      REAL    g11_2(0:n12,n22,n32),
     >        g22_2(n12,0:n22,n32),
     >        g33_2(n12,n22,0:n32),
     >        g11_3(0:n13,n23,n33),
     >        g22_3(n13,0:n23,n33),
     >        g33_3(n13,n23,0:n33),
     >        g11_4(0:n14,n24,n34),
     >        g22_4(n14,0:n24,n34),
     >        g33_4(n14,n24,0:n34)
c
c...pressure bc flags
      iip=iper
      jp=jper
      kp=kper
c...beging and end indexes
      jx=ix2-ix1
      jy=jy2-jy1
      jz=kz2-kz1
c
c opera su altri livelli per g11  g12 g13
c
      do 1 n=2,nlevel
c
      do 2 i=0,jxc(n)
      do 2 j=1,jyc(n)
      do 2 k=1,jzc(n)
c 
c piani a csi costante
c
      j2=2* j   *2**(n-2)
      j1=2*(j-1)*2**(n-2)
c
      k2=2* k   *2**(n-2)
      k1=2*(k-1)*2**(n-2)
c
      ic=2* i   *2**(n-2)
c
      if (i.eq.0.and.iip.eq.1) then
c
         if(n.eq.4) then
        icx=222

      end if
      i1=2*(i+1)*2**(n-2)
      i2=2*(i+2)*2**(n-2)
c
      x2 = x(i2)
      x1 = x(i1)
      x0 = x(0)
c
      y2 = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))
      y1 = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))
      y0 = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))

      z2 = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
      z1 = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
      z0 = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
c
      xcsi = 0.5 * ( -3.*x0 + 4.*x1 - x2 )
      ycsi = 0.5 * ( -3.*y0 + 4.*y1 - y2 )
      zcsi = 0.5 * ( -3.*z0 + 4.*z1 - z2 )
c
      else if (i.eq.jxc(n).and.iip.eq.1) then
c
         if(n.eq.4) then
        icx=222

      end if
      i1=2*(i-1)*2**(n-2)
      i2=2*(i-2)*2**(n-2)
c
      x2 = x(i2)
      x1 = x(i1)
      x0 = x(jx)
c
      y2 = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))
      y1 = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))
      y0 = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))
c
      z2 = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
      z1 = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
      z0 = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
c
      xcsi = 0.5 * ( 3.*x0 - 4.*x1 + x2 )
      ycsi = 0.5 * ( 3.*y0 - 4.*y1 + y2 )
      zcsi = 0.5 * ( 3.*z0 - 4.*z1 + z2 )
c
      else
c
      i2=2*(i+1)*2**(n-2)
      i1=2*(i-1)*2**(n-2)
c
      xdx = x(i2)
      xsn = x(i1)
c
      ydx = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))
      ysn = 0.25*(y(j2)+y(j2)+y(j1)+y(j1))
c
      zdx = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
      zsn = 0.25*(z(k2)+z(k1)+z(k2)+z(k1))
c
      xcsi = 0.5*(xdx-xsn)
      ycsi = 0.5*(ydx-ysn)
      zcsi = 0.5*(zdx-zsn)
c
      end if
c
      xeta = 0.0
c
      yeta = y(j2) -  y(j1) 
c
      zeta = 0.0
c
      xzet = 0.0
c
      yzet = 0.0
c
      zzet = z(k2) - z(k1)
c
c calcolo dei termini del tensore metrico controvariante
c su altre griglie (g11 g12 g13)
c
      if(i.eq.4.and.n.eq.4) then
      call gprima(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     >     t11,t22,t33)
      elseif (n.eq.4) then
      call gprima(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     >     t11,t22,t33)
      else
      call gprima(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     >     t11,t22,t33)
      end if
c
      if (n.eq.2) then
c
      pot=2.
      g11_2(i,j,k)=t11/pot
c
      else if (n.eq.3) then
c
      pot=4. 
      g11_3(i,j,k)=t11/pot
c
      else if (n.eq.4) then
c
      pot=8. 
      if(i.eq.4) then
       g11_4(i,j,k)=t11/pot
      else
       g11_4(i,j,k)=t11/pot
      end if
c
      end if
c
 2    continue
c
 1    continue

c
c ----------opera su altri livelli per g21  g22 g23-------------------
c

      do 3 n=2,nlevel
c
c piani a eta costante
c
      do 4 i=1,jxc(n)
      do 4 j=0,jyc(n)
      do 4 k=1,jzc(n)
c
      i2=2* i   *2**(n-2)
      i1=2*(i-1)*2**(n-2)
c
      k2=2* k   *2**(n-2)
      k1=2*(k-1)*2**(n-2)
c
      jc=2* j   *2**(n-2)
c
      if      (j.eq.0.and.jp.eq.1)  then
c
      j1=2*(j+1)*2**(n-2)
      j2=2*(j+2)*2**(n-2)
c
      x2 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x1 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x0 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
c
      y2 = y(j2)
      y1 = y(j1)
      y0 = y(j)
c
      z2 = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
      z1 = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
      z0 = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
c
      xeta = 0.5 * ( -3.*x0 + 4.*x1 - x2 )
      yeta = 0.5 * ( -3.*y0 + 4.*y1 - y2 )
      zeta = 0.5 * ( -3.*z0 + 4.*z1 - z2 )
c
      else if (j.eq.jyc(n).and.jp.eq.1) then
c
      j1=2*(j-1)*2**(n-2)
      j2=2*(j-2)*2**(n-2)
c
      x2 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x1 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x0 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
c
      y2 = y(j2)
      y1 = y(j1)
      y0 = y(jy)
c
      z2 = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
      z1 = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
      z0 = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
c
      xeta = 0.5 * ( 3.*x0 - 4.*x1 + x2 )
      yeta = 0.5 * ( 3.*y0 - 4.*y1 + y2 )
      zeta = 0.5 * ( 3.*z0 - 4.*z1 + z2 )
c
      else
c
      j2=2*(j+1)*2**(n-2)
      j1=2*(j-1)*2**(n-2)
c
      xsop = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      xsot = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
c
      ysop = y(j2)
      ysot = y(j1)
c
      zsop = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
      zsot = 0.25*(z(k2)+z(k1)+z(k1)+z(k2))
c
      xeta = 0.5*(xsop-xsot)
      yeta = 0.5*(ysop-ysot)
      zeta = 0.5*(zsop-zsot)
c     
      end if
c
      yzet = 0.0
c
      zcsi = 0.0
c
      ycsi = 0.0
c
      zzet = z(k2) - z(k1) 
c
      xcsi = x(i2) - x(i1)
c
      xzet = 0.0
c

      gia =xcsi*(yeta*zzet-yzet*zeta)-
     >     xeta*(ycsi*zzet-yzet*zcsi)+
     >     xzet*(ycsi*zeta-yeta*zcsi)

       if (gia.eq.0.) then
         eee=2323.
      end if
      
      call gseconda(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     >                 t11,t22,t33)
      if (n.eq.2) then
c
      pot=2.

      g22_2(i,j,k)=t22/pot
c
      else if (n.eq.3) then
c
      pot=4.
      g22_3(i,j,k)=t22/pot
c
      else if (n.eq.4) then
c
      pot=8.
      g22_4(i,j,k)=t22/pot
c
      end if
c
 4    continue
 3    continue
c
c
c ----------opera su altri livelli per g31  g32 g33-------------------
c
      do 6 n=2,nlevel
c
c piani a zita costante
c
      do 5 i=1,jxc(n)
      do 5 j=1,jyc(n)
      do 5 k=0,jzc(n)
c
      i2=2* i   *2**(n-2)
      i1=2*(i-1)*2**(n-2)
c
      j2=2* j   *2**(n-2)
      j1=2*(j-1)*2**(n-2)
c
      kc=2* k   *2**(n-2)
c
      if (k.eq.0.and.kp.eq.1) then
c
      k1=2*(k+1)*2**(n-2)
      k2=2*(k+2)*2**(n-2)
c
      x2 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x1 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x0 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
c
      y2 = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
      y1 = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
      y0 = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
c
      z2 = z(k2)
      z1 = z(k1)
      z0 = z(k)
c
      xzet = 0.5 * ( -3.*x0 + 4.*x1 - x2 )
      yzet = 0.5 * ( -3.*y0 + 4.*y1 - y2 )
      zzet = 0.5 * ( -3.*z0 + 4.*z1 - z2 )
c
      else if (k.eq.jzc(n).and.kp.eq.1) then
c
      k1=2*(k-1)*2**(n-2)
      k2=2*(k-2)*2**(n-2)
c
      x2 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x1 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      x0 = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
c
      y2 = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
      y1 = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
      y0 = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
c
      z2 = z(k2)
      z1 = z(k1)
      z0 = z(jz)
c
      xzet = 0.5 * ( 3.*x0 - 4.*x1 + x2 )
      yzet = 0.5 * ( 3.*y0 - 4.*y1 + y2 )
      zzet = 0.5 * ( 3.*z0 - 4.*z1 + z2 )
c
      else
c
      k2=2*(k+1)*2**(n-2)
      k1=2*(k-1)*2**(n-2)
c
      xav = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
      xdt = 0.25*(x(i2)+x(i2)+x(i1)+x(i1))
c
      yav = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
      ydt = 0.25*(y(j2)+y(j1)+y(j1)+y(j2))
c
c
      zav = z(k2)
      zdt = z(k1)
c
      xzet = 0.5*(xav-xdt)
      yzet = 0.5*(yav-ydt)
      zzet = 0.5*(zav-zdt)
c
      end if
c
      ycsi = 0.0
c
      zeta = 0.0
c
      yeta = y(j2) - y(j1) 
c
      zcsi = 0.0
c     
      xcsi = x(i2) - x(i1)
c
      xeta = 0.0
c
      call gterza(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     >                 t11,t22,t33)
c      
      if (n.eq.2) then
c
      pot=2.
      g33_2(i,j,k)=t33/pot
c
      else if (n.eq.3) then
c
      pot=4.
      g33_3(i,j,k)=t33/pot
c
      else if (n.eq.4) then
c
      pot=8.
      g33_4(i,j,k)=t33/pot
c
      end if
c
  5   continue
  6   continue    
c     
      return
      end
c
c
c--------------------------------------------------------------------
c
      subroutine gprima(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     &                  t11,t22,t33                                  )
c
c
c calcola gli elementi della prima linea del tensore metrico 
c controvariante su griglie alte per ciclo multigrid
c
c--------------------------------------------------------------------
c
       IMPLICIT REAL (a-h,o-z)
c
c--------------------------------------------------------------------
c
      giac=     xcsi*(yeta*zzet-yzet*zeta)-
     >          xeta*(ycsi*zzet-yzet*zcsi)+
     >          xzet*(ycsi*zeta-yeta*zcsi)
c
      csx = yeta*zzet - yzet*zeta
      csy = xzet*zeta - xeta*zzet
      csz = xeta*yzet - xzet*yeta
c
      etx = yzet*zcsi - ycsi*zzet 
      ety = xcsi*zzet - xzet*zcsi
      etz = xzet*ycsi - xcsi*yzet
c
      ztx=ycsi*zeta-yeta*zcsi
      zty=xeta*zcsi-xcsi*zeta
      ztz=xcsi*yeta-xeta*ycsi
c
      t11=(csx**2+csy**2+csz**2)/giac
c
      t22=(csx*etx+csy*ety+csz*etz)/giac
c
      t33=(csx*ztx+csy*zty+csz*ztz)/giac
c
      return
      end
c
c
c--------------------------------------------------------------------
c
      subroutine gseconda(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     >                 t11,t22,t33)
c
c calcola gli elementi della seconda linea del tensore metrico 
c controvariante su griglie alte per ciclo multigrid
c
c--------------------------------------------------------------------
c
       IMPLICIT REAL (a-h,o-z)
c
c--------------------------------------------------------------------
c
c
      etx = yzet*zcsi - ycsi*zzet 
      ety = xcsi*zzet - xzet*zcsi
      etz = xzet*ycsi - xcsi*yzet
c
      csx = yeta*zzet - yzet*zeta
      csy = xzet*zeta - xeta*zzet
      csz = xeta*yzet - xzet*yeta
c
      ztx=ycsi*zeta-yeta*zcsi
      zty=xeta*zcsi-xcsi*zeta
      ztz=xcsi*yeta-xeta*ycsi
c
      giac=xcsi*(yeta*zzet-yzet*zeta)-
     >     xeta*(ycsi*zzet-yzet*zcsi)+
     >     xzet*(ycsi*zeta-yeta*zcsi)
c
      t11=(etx*csx+ety*csy+etz*csz)/giac
c
      t22=(etx**2+ety**2+etz**2)/giac
c
      t33=(etx*ztx+ety*zty+etz*ztz)/giac
c
      return
      end
c
c
c--------------------------------------------------------------------
c
      subroutine gterza(xcsi,ycsi,zcsi,xeta,yeta,zeta,xzet,yzet,zzet,
     &                  t11,t22,t33                                  )
c
c calcola gli elementi della terza linea del tensore metrico 
c controvariant esu griglie alte per ciclo multigrid
c
c--------------------------------------------------------------------
c
       IMPLICIT REAL (a-h,o-z)
c
c--------------------------------------------------------------------
c
      csx = yeta*zzet - yzet*zeta
      csy = xzet*zeta - xeta*zzet
      csz = xeta*yzet - xzet*yeta
c
      etx = yzet*zcsi - ycsi*zzet 
      ety = xcsi*zzet - xzet*zcsi
      etz = xzet*ycsi - xcsi*yzet
c
      ztx=  ycsi*zeta - yeta*zcsi
      zty=  xeta*zcsi - xcsi*zeta
      ztz=  xcsi*yeta - xeta*ycsi
c
      giac=     xcsi*(yeta*zzet-yzet*zeta)-
     >          xeta*(ycsi*zzet-yzet*zcsi)+
     >          xzet*(ycsi*zeta-yeta*zcsi)
c
      t11=(ztx*csx+zty*csy+ztz*csz)/giac
c
      t22=(ztx*etx+zty*ety+ztz*etz)/giac
c
      t33=(ztx**2+zty**2+ztz**2)/giac
c
      return
      end
c
c
c-------------------------------------------------------------------
c
      subroutine wall(nlevel,jxc,jyc,jzc,
     &                in_dx1,in_dx2,in_dx3,in_dx4,
     &                in_sn1,in_sn2,in_sn3,in_sn4,
     &                in_sp1,in_sp2,in_sp3,in_sp4,
     &                in_st1,in_st2,in_st3,in_st4,
     &                in_av1,in_av2,in_av3,in_av4,
     &                in_in1,in_in2,in_in3,in_in4,   
     &                n1,n2,n3,                    
     &                n12,n22,n32,n13,n23,n33,
     &                n14,n24,n34                  )
c
c sono settati indici per il calcolo degli indici
c di pressione sulle pareti su tutte le griglie
c
c-------------------------------------------------------------------
c
      INTEGER    in_dx1(n1,n2,n3)
      INTEGER    in_dx2(n12,n22,n32)
      INTEGER    in_dx3(n13,n23,n33)
      INTEGER    in_dx4(n14,n24,n34)
c
      INTEGER    in_sn1(n1,n2,n3)
      INTEGER    in_sn2(n12,n22,n32)
      INTEGER    in_sn3(n13,n23,n33)
      INTEGER    in_sn4(n14,n24,n34)
c
      INTEGER    in_sp1(n1,n2,n3)
      INTEGER    in_sp2(n12,n22,n32)
      INTEGER    in_sp3(n13,n23,n33)
      INTEGER    in_sp4(n14,n24,n34)
c
      INTEGER    in_st1(n1,n2,n3)
      INTEGER    in_st2(n12,n22,n32)
      INTEGER    in_st3(n13,n23,n33)
      INTEGER    in_st4(n14,n24,n34)
c
      INTEGER    in_av1(n1,n2,n3)
      INTEGER    in_av2(n12,n22,n32)
      INTEGER    in_av3(n13,n23,n33)
      INTEGER    in_av4(n14,n24,n34)
c
      INTEGER    in_in1(n1,n2,n3)
      INTEGER    in_in2(n12,n22,n32)
      INTEGER    in_in3(n13,n23,n33)
      INTEGER    in_in4(n14,n24,n34)
c
      INTEGER    jxc(0:4),jyc(0:4),jzc(0:4)
c
c-------------------------------------------------------------------
c
c su tutte le griglie
c
      do 1 n=1,nlevel
c
      ix=jxc(n)
      iy=jyc(n)
      iz=jzc(n)
c
      if (n.eq.1) then
c
      call sett(ix,iy,iz,n1,n2,n3,
     >          in_dx1,in_sn1,in_sp1,in_st1,in_av1,in_in1)
c
      else if (n.eq.2) then 
c
      call sett(ix,iy,iz,n12,n22,n32,
     >          in_dx2,in_sn2,in_sp2,in_st2,in_av2,in_in2)
c
      else if (n.eq.3) then
c
      call sett(ix,iy,iz,n13,n23,n33,
     >          in_dx3,in_sn3,in_sp3,in_st3,in_av3,in_in3)
c
      else if (n.eq.4) then
c
      call sett(ix,iy,iz,n14,n24,n34,
     >          in_dx4,in_sn4,in_sp4,in_st4,in_av4,in_in4)
c
      end if
c
 1    continue
c
      return
      end

c
c
c-------------------------------------------------------------------
c
      subroutine sett(ix,iy,iz,nx,ny,nz,
     >          i_dx,i_sn,i_sp,i_st,i_av,i_in)
c
c
c settaggio indici
c
c-------------------------------------------------------------------
c
      include 'headers/common.h'
c
c-------------------------------------------------------------------
c
      integer  i_dx(nx,ny,nz)
      integer  i_sn(nx,ny,nz)
      integer  i_sp(nx,ny,nz)
      integer  i_st(nx,ny,nz)
      integer  i_av(nx,ny,nz)
      integer  i_in(nx,ny,nz)
c
c-------------------------------------------------------------------
c
      iip=iper
      jp=jper
      kp=kper
c
c in tutto il campo
c
      do k=1,iz
       do i=1,ix
        do j=1,iy
         i_dx(i,j,k) =1
         i_sn(i,j,k) =1
         i_sp(i,j,k) =1
         i_st(i,j,k) =1
         i_av(i,j,k) =1
         i_in(i,j,k) =1
        end do
       end do
      end do
c
c indici pareti sn e dx: (periodicita generalizzata)
c
      do 1 j=1,iy
      do 1 k=1,iz
c
      i_sn(1 ,j,k) =1-iip
      i_dx(ix,j,k) =1-iip
c
 1    continue
c
c parete sotto e sopra. (periodicita generalizzata)
c
      do 2 i=1,ix
      do 2 k=1,iz
c
      i_st(i,1 ,k) =1-jp
      i_sp(i,iy,k) =1-jp
c
 2    continue
c
c indici parete  indietro e avanti: (periodicita generalizzata)
c
      do 3 j=1,iy
      do 3 i=1,ix
c
      i_in(i,j,1 ) =1-kp
      i_av(i,j,iz) =1-kp    
c
 3    continue
c
      return
      end

