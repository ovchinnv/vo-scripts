c     5/04 spanwise corelations added
c     NOTE: this code is the wbl-box2 version (slightly modified from the standard)
c     10/04 surface pressure & cf planes added.
c     THIS CODE MAY HAVE HARDWIRINGS!
C-----------------------------------------------------------------------
C     ***************************                         
C     *         E D D Y         *                        
C     ***************************                       
C-----------------------------------------------------------------------
C     
C     THREEDIMENSIONAL INCOMPRESSIBLE NAVIER-STOKES SOLVER
C     (LARGE EDDY SIMULATION)                                              
C     ----------------------------
C     
C     PURPOSE:            - TIME ACCURATE NAVIER-STOKES SOLVER FOR 
C     TURBULENT FLOW FIELDS
C     - FINITE DIFFERENCE SCHEMES
C     
C     SUBGRID MODELS:     - SMAGORINSKI
C     - DYNAMIC MODEL
C     
C     
C     GRID:               - CARTESIAN RECTANGULAR GRID (SINGLE BLOCK)         
C     
C     SCHEMES:            - ADAMS-BASHFORTH (FRACTIONAL STEP)
C     - DIRECT POISSON SOLVER
C     - MULTIGRID POISSON SOLVER
C     
C     NON-DIMENSIONAL:    - BOUNDARY LAYER: 
C     U_INFTY, DISPLACEMENT THICKNESS
C     DENSITY = 1.
C     - CHANNEL:
C     U_TAU, HALF DIAMETER
C     - DUCT:
C     U_TAU, HALF DIAMETER
C     
C     INCLUDES:           - common.h
C     
C     INPUT:              - <dir>:    ../calc/
C     - <dir>:    ../init/
C     - NAMELIST
C     
C     OUTPUT:             - <dir>:    ../calc/
C     - RESTART-FILE
C     - STATISTICS
C     - INST. FLOW FIELD
C     
C     
C     COMPUTING TIME:     - in a typical run the single steps take 
C     the following comp. time (Adams Bashforth):
C     
C     - RHS                  1.707024       22.2 %
C     - velocity update      0.7407840       9.8 %
C     - b.c.                 4.2943999E-02   0.7 %
C     - divergence           0.2957280       3.9 %
C     - Poisson              0.5641280       7.5 %
C     - velocity correct.    0.6744159       8.8 %
C     - space average        0.1864160       2.6 %
C     - twall                2.9279999E-03   0.1 %
C     - b.c.                 4.1967999E-02   0.6 %
C     - strain               0.7076000       3.7 %
C     - turvis (Smagor.)     1.554768       20.2 %
C     - space time average   0.4792160       6.4 %
C     - screen info          0.7446880       9.8 %
C     - total (sec)          7.7426078
C     
C-----------------------------------------------------------------------
C     
C     (C) I. BALARAS
C     VON KARMAN INSTITUTE FOR FLUID DYNAMICS
C     CHAUSSEE DE WATERLOO, 72
C     B-1640 RHODE-ST-GENESE
C     
C     DATE OF CREATION :  20/10/1993                                
C     LAST UPDATE      :  05/10/2000
C     
C-----------------------------------------------------------------------
C     
      SUBROUTINE HYDRO_IB(koutput,communicator)
C     
      implicit none
C-----------------------------------------------------------------------
C     
      include 'mpif.h'
      include 'headers/common.h'
      include 'headers/dimension.h'
      include 'headers/ose.h'

#ifdef immbnd ! conditional compilation     
      include 'headers/immbd.h'
#endif
C     
C-----------------------------------------------------------------------
C     declarations 
C-----------------------------------------------------------------------
C     
      INTEGER nt
      INTEGER n1,n2,n3
c...  the input parameter
      integer koutput,communicator

c.... NT number of timesteps read each time as slice data
      PARAMETER (NT = 2)

c.... NT number of timesteps for second set of slice data
      integer ntnd
      PARAMETER (NTND = 2)

C     ...External Functions...
      REAL     PEX, UEX, WEX
      EXTERNAL PEX, UEX, WEX
C     
C.....arrays for the solver
C     
C.....xu, yv, wz:		cell faces
C.....xc, yc, zc:               cell centers
C.....uo, vo, wo:		inst. velocity field
C.....ustar,vstar,wstar:	predicted velocities
C.....p, dp, divv:		pressure, correction, divergence
C.....uab1, vab1, wab1:		RHS predictor step
C.....uab2, vab2, wab2:		RHS corrector step
C.....tv:                       turbulent viscosity
C.....sxx,sxy,sxz,syy,syz,szz:	strain rates
C.....dummy:                    dummy array for Poisson solver
C.....uoc,voc,woc:		velocities at cell centers
C.....upr, vpr, wpr:		fluctuating components
C.....wfrg:                     fringe velocity
C.....StrkPos,StrfIndx:         streaks
C     
C     allocate arrays
      real xu(nx),xc(nx),xcf(nx)
      real yv(ny),yc(ny),ycf(ny)
      real zw(nz),zc(nz),zcf(nz)
      real zwg(nzg),zcg(nzg)
c     
      real uo(nx,ny,nz),vo(nx,ny,nz),wo(nx,ny,nz)
      real ustar(nx,ny,nz),vstar(nx,ny,nz),wstar(nx,ny,nz)
      real p(nx,ny,nz),dp(nx,ny,nz),divv(nx,ny,nz)
      real uab1(nx,ny,nz),vab1(nx,ny,nz),wab1(nx,ny,nz)
      real uab2(nx,ny,nz),vab2(nx,ny,nz),wab2(nx,ny,nz)
      real upr(nx,ny,nz),vpr(nx,ny,nz),wpr(nx,ny,nz)
      real ucn(nx,ny,nz),vcn(nx,ny,nz),wcn(nx,ny,nz)
c
c     
      real tv(nx,ny,nz)
      real dummy(nx,ny,nz) !artifice for the doubly periodic poisson solver

#ifdef LES ! conditional compilation
c     LES DATA 
      real dxdydz(nx,ny,nz),g(nx,ny,nz)
      real sxx(nx,ny,nz),syy(nx,ny,nz),szz(nx,ny,nz)
      real sxy(nx,ny,nz),sxz(nx,ny,nz),syz(nx,ny,nz)
      real uoc(nx,ny,nz),voc(nx,ny,nz),woc(nx,ny,nz)
      real gbar(nx,ny,nz),lm(nx,ny,nz),mm(nx,ny,nz)
      real Ilm(nx,ny,nz),Imm(nx,ny,nz)
      real wrk1(nx,ny,nz),wrk2(nx,ny,nz),wrk3(nx,ny,nz)
      real wrk4(nx,ny,nz)
c      stats
      real tvavg(nx,nz),tvavgg(nx,nzg),tvi
#endif

      real wfrg(nx,ny,nz) ! probably for the "fringe region" boundary layer sim 
c     
      real u2d(nt,ny,nz),v2d(nt,ny,nz),w2d(nt,ny,nz)
      real u2dnd(ntnd,ny,nz),v2dnd(ntnd,ny,nz),w2dnd(ntnd,ny,nz)
c     
      real urun(nx,nz),vrun(nx,nz),wrun(nx,nz)
      real uinf(nx),duinf(nx)
      integer iuinf
c     
      real cutcoeffs(1:nx), u1(1:nx),w1(1:nx)
      real u2(1:nx),w2(1:nx)
cccccccc variables for 2D planes output
      logical planes
      integer nyplanes,nzplanes,nxplanes
      parameter (nyplanes=1, nxplanes=3, nzplanes=3)
      real nyploc(nyplanes), nxploc(nxplanes), nzploc(nzplanes)
      integer cpuxy(nzplanes),kxy(nzplanes)
      real u2dxz(nx,nz,2),v2dxz(nx,nz,2),w2dxz(nx,nz,2),p2dxz(nx,nz,2)
      real u2dxy(nx,ny),v2dxy(nx,ny),w2dxy(nx,ny),p2dxy(nx,ny)
      real u2dyz(ny,nz,2),v2dyz(ny,nz,2),w2dyz(ny,nz,2),p2dyz(ny,nz,2)
cccccccc
      real uinft(nz),uinftg(nzg),ubulk(nz),cf(nz),dst(nz),theta(nz),winst(nz),tauwg(nzg),utaug(nzg),ubulkg(nzg),dstg(nzg)
      real thetag(nzg),cfg(nzg),kpg(nzg),winstg(nzg)
      real uavg(nx,nz),vavg(nx,nz),wavg(nx,nz),uuavg(nx,nz),vvavg(nx,nz),wwavg(nx,nz),uwavg(nx,nz)
c      real omxavg(nx,nz),omyavg(nx,nz),omzavg(nx,nz)
      real wtemp(ny,nz),pavg(nx,nz),ppavg(nx,nz)
      real uuav(nz),vvav(nz),wwav(nz) ! real stat arrays
      real uav(nz),vav(nz),wav(nz) ! real stat arrays
      real guav(nzg),gvav(nzg),gwav(nzg) ! real stat arrays
      real guuav(nzg),gvvav(nzg),gwwav(nzg) ! real stat arrays
      real usave(nx,ny),vsave(nx,ny),wsave(nx,ny)
c 
      real uosm,vosm,wosm ! for rescaling inflow noise on the fly
c 
      integer statnum,samples,psample
c     arrays for square spectrum inflow
      real k1(1500),k2(1500),k3(1500),phase1u(3000),phase1v(3000),phase1w(3000),
     &        phase2u(1500),phase2v(1500),phase2w(1500),phase3u(1500),
     &        phase3v(1500),phase3w(1500),ampl(1500)
c     equivalent arrays for isotropic inflow     
      real kx(1500),ky(1500),kz(1500)
      complex uh(1500),vh(1500),wh(1500)
      integer count(1500)
      real uose(ny,nx),vose(ny,nx),wose(ny,nx)
      integer jcounter,ipos
c     other arrays for isotropic inflow
      integer is(1500,kose),kkmax ! kose defined in ose.h
c     scalars for isotropic inflow
      real x_u,x_w,y_v,y_u
      complex zh1,zh2,xh_u,xh_w,yh_v,yh_u,v1,v2,v3,temp
c 
      real uavgg(nx,nzg),vavgg(nx,nzg),wavgg(nx,0:nzg),uuavgg(nx,nzg),vvavgg(nx,nzg),
     &           wwavgg(nx,nzg),uwavgg(nx,nzg),pavgg(nx,nzg),ppavgg(nx,nzg)
c      real omxavgg(nx,ngz),omyavgg(nx,nzg),omzavgg(nx,nzg)
c     start of inflow        
      integer infstart,nper
c
c
      REAL fjp
C     
      CHARACTER*(30) InputFile
C     
      REAL TPRESSURE(3)
      REAL    qo,qold,dpdx,tke,divmax,dpdz
      INTEGER i, j, k, ii, jj, kk, iii, jjj, kkk
C     
      LOGICAL Ltstep            ! .f.-predictor .t.-corrector
      REAL    a1,a2,tmeth
C     
      REAL    delt, yplusmin, erm
      INTEGER icycle,istart,icount1,icount2,icinf,Niter,ibnd
      REAL    dtm1,dtold,adiv,zz
C     
      INTEGER ipart
c...  slice data
      INTEGER ic
      REAL    tlevel,t_old
      INTEGER iSliceSize,ksub
c...  series of restart files
      character*3 index(0:999)
      character*5 index2(0:99999)
      logical  LcreateR

c...  second series of slice data
      logical  Lndslice
c...  used to determine when to output
      logical scrivi
c     
      real   a, b, lamda, pi, u0, w0
c     
      INTEGER IERR
      real*8 cpue,cpub
      real*8 cputurvis
      real*8 cpucorr
      real*8 cpupres,cpupredict,cpurhs
c     
      integer nimu,nimv,nimw,nimp
      integer nbodyu,nbodyv,nbodyw,nbodyp
      integer nb
c
      LOGICAL wrim         !ugo's artifice (write immersed boundary data?)
      character*80 filestr
c
      integer type !double or single (for MPI)
c     
c     added for boundary layer simulation VO
      integer iistart,imid,imidf,ifin,nload,nsave,n,nn,status(MPI_STATUS_SIZE),counter,cpusave
      integer np,ksave,savenum,loadnum
      logical inflow, save,  rstat, spec
      real zstart, xloc, wnoise,  uinfin, uscale ,pscale, plevel, ux,vx,saveloc,scx
      real dummy1,dummy2,dummy3
c-----for box2 -- wbl-specific
      integer nb2
      parameter (nb2=232)
      logical box2
      integer nx_b2, ny_b2, i_cen, i_cor
      real dt_b2, x_b2(nb2), xc_b2(nb2), y_b2(ny), t_b2, u_b2(ny,nb2),v_b2(ny,nb2),w_b2(ny,nb2)
c-----------------------------
c--------------statistics variables-------------------------------------
      real uact,uint,u,v,w,uu,vv,ww,pp,ppp,uw,uactg,fac
c-----------------------------------------------------------------------
c     OSE FUNCTIONALITY FLIP CASE
ccccccccccccccccccccccccccccccccccccccccccccccccccc
c     GRID DATA
      real eta(nose),detacen(nose),detacor(nose)
c     OSE INTERPOLATION DATA
      real a_w(nx),b_w(nx),a_u(nx),b_u(nx)
      integer ps_w(nx),ps_u(nx),nxloc
ccccc BLASIUS VARIABLES      
      complex U_(nose),d2U(nose)
      real    Ur(2*nose),d2Ur(2*nose)
ccccc OSE VARIABLES
      complex phi(nose),phip(nose)
      real phir(2*nose),phipr(2*nose)
ccccc 
ccccc VARIABLES FOR STREAMWISE SPECTRA
      integer nspec,specfreq,skip
      parameter (nspec=14,specfreq=3,skip=3) !spectra to be computed at nspec locations; data sampled once in every specfreq timesteps at every "skip"'th z-point
      real specpos(nspec) 
      integer spec_cpu(nspec), spec_k(nspec)
ccccccccccccccccccccccccccccccccccccccccccccccccccc      
c
c     misc variables
      real topbc(3,ny,nz), botbc(3,ny,nz) !these variables are needed when the code is run implicitly
      real toinit
      real ran1 !from numerical recipes
      integer idum
      integer real1d, real2d, stat2d, cross2d !new MPI types
      real dstar ! boundary layer variable
      integer idir ! IO.F directive 
      integer m, m2
      real oldtime, znewtime ! variables for keeping track of time per timestep
      real amp,env !for perturbing laminar wake-bl interaction

c        unknown variables
      integer itstamp
      real dpdy ! it seems that someone wanted to introduce a spanwise or tranverse pressure gradient
      integer icnd, iplanend, nsubnd ! seems to be needed by slice data routines
c
c      temporary debugging variables
      real divmax_,divmin_,divavg_
      integer idmax, jdmax, kdmax, idmin, jdmin, kdmin
      real rsu, rsv, rsw ! currently exist for backward compatibility with "correct"        
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef immbnd ! conditional compilation
c     IMMERSED BOUNDARY DATA
      real sb(nnmax),yb(nnmax),xb(nnmax)
      real coefx(nnmax,3),coefy(nnmax,3)
      real uim(nnmax),vim(nnmax),wim(nnmax),tvim(nnmax),Ilmim(nnmax)
      integer nb
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc VARIABLES FOR THE PARALLELIZATION OF IMMERSED BOUNDARY cccccccccccccccc
c U
      integer uniml(0:ncpu-1),uiiml(0:ncpu-1,nnmax),ujiml(0:ncpu-1,nnmax),untop(0:ncpu-1),utk(0:ncpu-1,nimb),
     &        ustk1(0:ncpu-1,nnmax),ustk2(0:ncpu-1,nnmax),unbot(0:ncpu-1),ubk(0:ncpu-1,nimb),usti(0:ncpu-1,nnmax)
      real    utib_ind(0:ncpu-1,nimb),utib_c1(0:ncpu-1,nimb),utib_c2(0:ncpu-1,nimb),ubib_ind(0:ncpu-1,nimb),
     &        ubib_c1(0:ncpu-1,nimb),ubib_c2(0:ncpu-1,nimb),uc4intl(0:ncpu-1,nnmax,5)
c V
      integer vniml(0:ncpu-1),viiml(0:ncpu-1,nnmax),vjiml(0:ncpu-1,nnmax),vntop(0:ncpu-1),vtk(0:ncpu-1,nimb),
     &        vstk1(0:ncpu-1,nnmax),vstk2(0:ncpu-1,nnmax),vnbot(0:ncpu-1),vbk(0:ncpu-1,nimb),vsti(0:ncpu-1,nnmax)
      real    vtib_ind(0:ncpu-1,nimb),vtib_c1(0:ncpu-1,nimb),vtib_c2(0:ncpu-1,nimb),vbib_ind(0:ncpu-1,nimb),
     &        vbib_c1(0:ncpu-1,nimb),vbib_c2(0:ncpu-1,nimb),vc4intl(0:ncpu-1,nnmax,5)
c W
      integer wniml(0:ncpu-1),wiiml(0:ncpu-1,nnmax),wjiml(0:ncpu-1,nnmax),wntop(0:ncpu-1),wtk(0:ncpu-1,nimb),
     &        wstk1(0:ncpu-1,nnmax),wstk2(0:ncpu-1,nnmax),wnbot(0:ncpu-1),wbk(0:ncpu-1,nimb),wsti(0:ncpu-1,nnmax)
      real    wtib_ind(0:ncpu-1,nimb),wtib_c1(0:ncpu-1,nimb),wtib_c2(0:ncpu-1,nimb),wbib_ind(0:ncpu-1,nimb),
     &        wbib_c1(0:ncpu-1,nimb),wbib_c2(0:ncpu-1,nimb),wc4intl(0:ncpu-1,nnmax,5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif

c
cccccccccVARIABLES FOR REYNOLDS STRESS BUDGETS -- SEE RSBUD1.F routine
c        NOTE:  TO SAVE MEMORY rsblk_1 overwrites uab2 (1st 4*7 nx*nz planes) This works if ny>=28!!!
c                if LES is on rsblk_3 overwrites vab2  (1st 4*3 nx*nz planes)
      real rsblk_1(nx,nz,ny)
      real rsblk_2(nx,0:nz+1,33)
      real rsblkg_1(nx,nzg,7) 
#ifdef LES
      real rsblk_3 (nx,nz,ny)     !temporary LES block
      real rsblk_4 (nx,0:nz+1,25) !permanent LES block
      real rsblkg_3(nx,nzg,4)     !temporary global LES block
#endif

cccccccccc variables for two-pt correlations ccccccccccc           
      integer nxl,nzl,nz_act,nx_act
      integer zfreq,xfreq
      parameter (xfreq=5,zfreq=20)
      parameter (nxl=(nx-2+xfreq-1)/xfreq, nzl=(nz-2+zfreq-1)/zfreq) !   
      real xn(nxl), zn((nzg-2+zfreq-1)/zfreq)        
      integer ks(0:ncpu-1),ke(0:ncpu-1),npt(0:ncpu-1) ! these arrays are defined on all processors
      real um2pt(nxl,nzl), vm2pt(nxl,nzl), wm2pt(nxl,nzl)
      real u2pt(nxl,ny,nzl),v2pt(nxl,ny,nzl),w2pt(nxl,ny,nzl)
      real u2tmp(nxl,ny,nzl),v2tmp(nxl,ny,nzl),w2tmp(nxl,ny,nzl) ! temporary arrays
      logical two_pt_corr

c
      character*4 rsnames(0:3)
c
      equivalence (uab2,rsblk_1)
#ifdef LES
      equivalence (vab2,rsblk_3)
#endif
      equivalence (rsblkg_1(1,1,1),uavgg(1,1)),(rsblkg_1(1,1,2),vavgg(1,1)),
     &            (rsblkg_1(1,1,3),wavgg(1,1)),(rsblkg_1(1,1,4),uuavgg(1,1)),
     &            (rsblkg_1(1,1,5),vvavgg(1,1)),(rsblkg_1(1,1,6),wwavgg(1,1)),
     &            (rsblkg_1(1,1,7),pavgg(1,1))
#ifdef LES
      equivalence (rsblkg_3(1,1,1),uwavgg(1,1)),(rsblkg_3(1,1,2),tvavgg(1,1))     
#endif     
c      
      equivalence (U_,Ur),(d2U,d2Ur),(phi,phir),(phip,phipr)
c     
      equivalence (kx,k1),(ky,k2),(kz,k3),(phase1u,uh),(phase1v,vh),(phase1w,wh)
      equivalence (uose,usave),(vose,vsave),(wose,wsave)
c
#ifdef LES
#ifdef immbnd
      equivalence (tvim,Ilmim)
#endif
#endif
      equivalence (dp,divv)
c      equivalence (uo,ustar),(vo,vstar),(wo,wstar)
      equivalence (upr,ucn),(vpr,vcn),(wpr,wcn)
c      equivalence (u2tmp,ustar),(v2tmp,vstar),(w2tmp,wstar)
c      equivalence (u2dxy,uab2),(v2dxy,vab2),(w2dxy,wab2)
c      equivalence (u2dxz,uab2),(v2dxz,vab2),(w2dxz,wab2)
c      equivalence (u2dyz,uab2),(v2dyz,vab2),(w2dyz,wab2)
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     input data and setup
c-----------------------------------------------------------------------
#ifdef SINGLE
      TYPE=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif
c...  start mpi, get rank and number of processors
      MPI_COMM_HYDRO=communicator
      CALL MPI_COMM_SIZE(MPI_COMM_HYDRO,MYSIZE,IERR)
      CALL MPI_COMM_RANK(MPI_COMM_HYDRO,MYRANK,IERR)
c
      if (mysize.ne.ncpu) then 
        write(0,'("incorrect number of cpus requested (",I3,"): need", I3)') mysize,ncpu
	return 
      endif
c
c     
      toinit = ran1(-1)
      idum=10
      
c...  fill index array for various output purposes
      m=0
      m2=0
      do i=0,9
         do j=0,9
            do k=0,9
               index(m)=CHAR(i+48)//CHAR(j+48)//CHAR(k+48)
               m=m+1
      	       do jj=0,9
      	         do kk=0,9
      		   index2(m2)=char(i+48)//char(j+48)//char(k+48)//char(jj+48)//char(kk+48)
      		   m2=m2+1
                 enddo
      	       enddo	 
            enddo
         enddo
      enddo
c     
c...  set flag for immersed boundary conditions
c     * true or false
c     
      ibswitch=.false. 
#ifdef immbnd
      ibswitch=.true.
#endif       
c     
c...  input parameters (namelist-statements)
c     
      InputFile = 'eddywbl.input'
      IF(MYRANK.EQ.0) WRITE(6,*) InputFile,myrank,mysize
c     
      CALL INPALL(InputFile)
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '*..data read'
         write(6,*) 'ru1=', ru1
      ENDIF
c     
c...  grid and matrix coefficients
c     
      CALL Grid(str2,xu,yv,zw,zwg,nx,ny,nz,nzg)
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '*..grid read'
         write(6,*) ' dx=',dx,'dy=',dy,'dz=',dz
         write(6,*) ' ix1 =',ix1,' ix2 =',ix2
         write(6,*) ' jy1 =',jy1,' jy2 =',jy2
         write(6,*) ' kz1 =',kz1,' kz2 =',kz2
      ENDIF
c.....solver setup
      CALL SETUP

      CALL centercord(xu,yv,zw,xc,yc,zc,nx,ny,nz)
c     get global cell center coordinates
      zcg(1)=1.5*zwg(1)-0.5*zwg(2)
      do k=2,nzg
        zcg(k)=0.5*(zwg(k-1)+zwg(k))
      enddo
c             	       
      CALL centerfit (xu,yv,zw,xcf,ycf,zcf,nx,ny,nz)
#ifdef LES
      CALL filtersize(xu,yv,zw,dxdydz,nx,ny,nz)
      CALL REFRESHBC(dxdydz,nx*ny,(nz-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      if (myrank.eq.0) then 
        do i=2,ix2
	 do j=2,jy2
	  dxdydz(i,j,1)=dxdydz(i,j,2)
	 enddo
	enddo
      endif	  
      CALL weights(xc,yc,zc,nx,ny,nz)
#endif
c
      if (myrank.eq.0) then
         WRITE(6,*) '*...solver setup'
      endif
c     inflow noise setup VO -- obsolescent
c     
      inflow = .false.
      box2 = .true.
      save = .true.
      spec =.true.
      rstat = .true.
      two_pt_corr = .true.
      planes=.true.
c      save=.true.
      nper=0
c
c     initiaize statistics variables
      samples=0
      statnum=0
      do k=kz1,kz2+1
        uinft(k)=0.
	uuav(k)=0.
	vvav(k)=0.
	wwav(k)=0.
        uav(k)=0.
	vav(k)=0.
	wav(k)=0.
        do i=ix1,ix2+1
          uavg(i,k)=0.
	  vavg(i,k)=0.
	  wavg(i,k)=0.
	  pavg(i,k)=0.
	  uuavg(i,k)=0.
	  vvavg(i,k)=0.
	  wwavg(i,k)=0.
	  uwavg(i,k)=0.
c	  omxavg(i,k)=0.
c	  omyavg(i,k)=0.
c	  omzavg(i,k)=0.
#ifdef LES
	  tvavg(i,k)=0.
#endif
	enddo
      enddo
c	initialize reynolds stress variables
      rsnames(0)='rsuu'
      rsnames(1)='rsvv'
      rsnames(2)='rsww'
      rsnames(3)='rsuw'
      
      do j=1,33
        do k=1,nz
	  do i=1,nx
	    rsblk_2(i,k,j)=0.
	  enddo    
        enddo
      enddo	
#ifdef LES
      do j=1,25
        do k=1,nz
	  do i=1,nx
	    rsblk_4(i,k,j)=0.
	  enddo    
        enddo
      enddo	
#endif      
c      		
      call MPI_TYPE_CONTIGUOUS(nz-2,type,real1d,ierr)
      call MPI_TYPE_CONTIGUOUS((nz-2)*nx,type,stat2d,ierr)
      call MPI_TYPE_COMMIT(real1d,ierr)
      call MPI_TYPE_COMMIT(stat2d,ierr)
      call MPI_TYPE_VECTOR(nz,nx,ny*nx,type,real2d,ierr)
      call MPI_TYPE_VECTOR(nz*ny,1,nx,type,cross2d,ierr)
      call MPI_TYPE_COMMIT(cross2d,ierr)
      
      call MPI_TYPE_COMMIT(real2d,ierr)
c     
      imid=nx/2  !default values
      imidf=nx/2 !default values
c

      if (two_pt_corr) then 
c      compute 2-pt correlation grids (center)      
       ii=0
       do i=2,nx-1,xfreq
        ii=ii+1
        xn(ii)=xcf(i)
       enddo
       nx_act=ii
c
       kk=0
       do k=2,nzg-1,zfreq
        kk=kk+1
        zn(kk)=zcg(k)
       enddo
       
c  compute local z-indices
       nz_act=0
       do j=0,ncpu-1
        ks(j)=(nz-2)*j+2
	ke(j)=ks(j)+nz-3
	ks(j)=ks(j)+mod(zfreq-mod(ks(j)-2,zfreq),zfreq)! still global at this point
	ke(j)=ke(j)-mod(ke(j)-2,zfreq)      ! still global
	npt(j)=(ke(j)-ks(j))/zfreq+1        ! number of entries
        nz_act=nz_act+npt(j)
cccccc debug
c      if (myrank.eq.0) then 
c       write(35,'(5I)') j,ks(j),ke(j),ks(j)-(nz-2)*j,ke(j)-(nz-2)*j
c      endif
c      
cccccc debug
c       convert to local form
        ks(j)=ks(j)-(nz-2)*j	
	ke(j)=ke(j)-(nz-2)*j
       enddo
c      close(35)
   
c      do j=1,ncpu-1
c       write(36,'(6I)') j,nxl,nx_act,nzl,nz_act,npt(j)
c      enddo
c      close(36) 

      endif ! two-pt-corr 	      
c
      if (jdose) then
c     
         pi=4.*atan(1.0)
	 scx=2.*pi/zlen
         zstart=0.0
         wnoise=30.
         imid=ix1
         i=ix1
         do while (xu(i).lt.(zstart+5.)) ! where the stats are taken
            i=i+1
         enddo
	 imid=i
c
         do while (xu(i).lt.(zstart+wnoise))
            i=i+1
         enddo
	 imidf=i
c      target turbulence level
         plevel=u_rms ! u_rms is in the input file
c
c     read wavenumbers
         open(101,file=NOISEFILE,form='formatted',status='old') ! NOISEFILE is in the input
         read(101,'(I5)') jcounter,kkmax
	 if (kkmax.ne.kose) stop 'number of modes allocated is different from that specified (jd-ose)'
         if (jcounter.gt.1500) stop 'increase space for inflow wavenumbers (jd-ose)'
         do k=1,jcounter 
          read(101,'(2 I5,9G25.15)') i,count(k),kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
          do ipos=1,count(k)
           read(101,'(I5)') is(k,ipos)
          enddo	  	       
         enddo
         close(101)
         open(101,file='pscale.dat',form='formatted',status='old')
	   read(101,'(G25.15)') pscale
	   read(101,'(I6)') psample
	 close(101)  
c	 pscale=1.
c
	 uscale=sqrt(scx)*1. ! disturbance scaling parameter, do _not_ remove sqrt(scx)
c
c     set up ose solver and interpolation parameters
         nxloc=(nx-2)/mysize
	 if (mod(nx-2,mysize).ne.0) nxloc=nxloc+1
c	 
         if (-x_.ne.x_ref) stop 'xref mismatch between input file and ose.h' ! x_ obsolescent
	 
         L_=(ru1*(x_-xref_)/uinfty)**0.5 ! in ose.h ! x_ obsolescent
         Re_=uinfty*L_/ru1
c     set up grid
         if (abs(zlen-yv(ny-1)).gt.1E7) stop 'OSE: spanwise domain size not correct (see ose.h)'
         call grid_ose(eta,detacen,detacor,nose,zlen,stretch,L_) ! stretch and nose in ose.h
c
         call setinterp_ose(eta,xu(2),nose,nx-1,a_w(2),b_w(2),ps_w(2),L_) 
         call setinterp_ose(eta,xcf(2),nose,nx-1,a_u(2),b_u(2),ps_u(2),L_) ! because the 1st point is ghost
         call setup_ose(eta,detacor,d2U,U_,Ur,nose,ru1)

      endif
c   *** wbl box2-specific
        if (box2) then
	 if (myrank.eq.0) then 
	  loadnum=48
	  open(63, file='./box2/box2.fstle.6-'//index(loadnum),form='binary',status='old')
	  read(63) t_b2,nx_b2,ny_b2,(x_b2(i),i=1,nx_b2),(y_b2(j),j=1,ny_b2)
          tlevel=t_b2
	  if (lcfl) write(0,*) 'WARNING: box2 with cfl ON'
	  if (ny_b2.ne.ny) write(0,*) 'WARNING: box2 and box1 spanwise domains do not match'
	  if (nx_b2.ne.nb2) write(0,*) 'WARNING: box2 and box1 transverse domains do not match'
	  do i=2,nx_b2
	   xc_b2(i)=0.5*(x_b2(i)+x_b2(i-1))
	  enddo
	  xc_b2(1)=-xc_b2(2)

c         offset grids -- throw away near-wall region in the disturbance
          do i=1,nx_b2
	   x_b2(i)=x_b2(i)-0.75
	   xc_b2(i)=xc_b2(i)-0.75
	  enddo  

	  call setinterp_ose( x_b2, xu(2),nx_b2,nx-1,a_w(2),b_w(2),ps_w(2),1.) !w-interpolation (streamwise) 
	  call setinterp_ose(xc_b2,xcf(2),nx_b2,nx-1,a_u(2),b_u(2),ps_u(2),1.) !u,v-interpolation (span & transverse)
	 endif !myrank
	 call MPI_BCAST(tlevel,1,type,0,MPI_COMM_HYDRO,ierr)
ccccccccccccc skip until start index
         istart=12500

         if (myrank.eq.0) then 
          do k=1,mod(istart,500)  
          read(63) dt_b2
	   do i=2,nx_b2-1
	    read(63) (u_b2(j,i), v_b2(j,i), w_b2(j,i),j=2,ny-1)
           enddo !i
          enddo !k
	 endif !myrank 	   
	 
cccccccccccc
	endif !box2         

c   *** wbl box2-specific
       
      
c
      saveloc=215. ! where slices are to be saved
      if (save) then
c     
         nsave=0 
         savenum=25
         if (myrank.eq.0) then 
      	   if (lcfl) write(0,*) 'WARNING: saving slices with cfl ON' 
         endif
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
         if (myrank.eq.cpusave) open(59,file='./box3/box3.fstle-'//index(savenum),form='unformatted',status='unknown')
      endif
c
      if (spec) then
c      
         if (myrank.eq.0) then 
      	   if (lcfl) write(0,*) 'WARNING: spectra sampled with cfl ON!' 
         endif
c	 
         specpos(1)=1.
         specpos(2)=10.
         specpos(3)=20.
	 specpos(4)=30.
	 specpos(5)=40.
	 specpos(6)=50.
	 specpos(7)=60.
         specpos(8)=75.
	 specpos(9)=100.
	 specpos(10)=125.
	 specpos(11)=150.
	 specpos(12)=175.
	 specpos(13)=200.	 
         specpos(14)=230.
c        	 
         do i=1,nspec !same code as for 'save'
	 
          spec_k(i)=1
	  spec_cpu(i)=0 
          k=1
          do while (zwg(k).lt.(specpos(i)))
            k=k+1
          enddo
	  if (k.gt.1) then 
	    spec_cpu(i)=(k-2)/(nz-2)
	    spec_k(i)=mod(k-2,nz-2)+1
          endif
	  if (spec_cpu(i).gt.mysize-1) then 
	    spec_cpu(i)=spec_cpu(i)-1
	    spec_k(i)=nz
	  endif  
          if (myrank.eq.spec_cpu(i)) then 
      	     open(600+i,file='spec.fstle2.'//index(int(specpos(i))),form='unformatted',status='unknown')
             write(600+i) specpos(i),nx,(xu(k),k=1,nx,skip)
       	  endif      
         enddo
      endif
c
      if (inflow) then          !if shirani inflow is used 
       if (myrank.eq.0) then
         open(19,file='inflow20.ini', status='old',form='formatted')
         read(19,*) counter
c     write(20,*) counter
         do k=1,counter
	    read(19,*) k1(k),k2(k),k3(k),phase1u(k),phase1v(k),phase1w(k),
     &           phase2u(k),phase2v(k),phase2w(k),phase3u(k),phase3v(k),
     &           phase3w(k),ampl(k)
         enddo	   
         close(19)
c     uscale=1.5
         uscale=2.5	
         zstart=1.
         wnoise=20.
         iistart=ix1
         imid=ix1
         ifin=ix1
c     write(6,*) zw(kz2),2.*pi
         if ((zstart+wnoise).le.xu(nx)) then
            i=ix1
            do while (xu(i).lt.zstart)
               i=i+1
            enddo
            iistart=i
c     
            do while (xu(i).lt.(zstart+wnoise/4))
               i=i+1
            enddo
            imid=i
	    imidf=imid
c     
            do while (xu(i).lt.(zstart+wnoise))
               i=i+1
            enddo
            ifin=i
c                
         else 
            if (myrank.eq.0) write(6,*) 'error: domain size too small for noise inflow'
            stop
         endif
         if (myrank.eq.0) write(6,*) 'wavenumbers read'
c     write(6,*) cpustart,cpufin,myrank,kstart,kfin
       endif
      endif
ccccccccccccc OUTPUT PLANES      
      if (planes) then 
       nxploc(1)=0.2
       nxploc(2)=0.5
       nxploc(3)=1.
       nzploc(1)=45.
       nzploc(2)=75.
       nzploc(3)=100.
       nyploc(1)=25.
c      find save locations for xy planes       
       do n=1,nzplanes
         kxy(n)=1
	 cpuxy(n)=0 
         k=1
         do while (zwg(k).lt.(nzploc(n)))
            k=k+1
         enddo
	 if (k.gt.1) then 
	   cpuxy(n)=(k-2)/(nz-2)
	   kxy(n)=mod(k-2,nz-2)+1
         endif
	 if (cpuxy(n).gt.mysize-1) then 
	   cpuxy(n)=cpuxy(n)-1
	   kxy(n)=nz
	 endif  
	enddo !n
      endif !planes	

cccccccccccccccccccc PLANES      
c       
      IF(ifield.EQ.0) THEN
c     
      pi=4.0*atan(1.)
         do i=1,ix2+1
	    dummy1=xcf(i)/xu(ix2)
	    dummy2=4.*1./ru1/2.*dummy1*(1-dummy1)
            do j=1,jy2+1
               do k=1,kz2+1
                  uo(i,j,k)=0.0!+0.1/ru1*sin(2.*pi/2*xu(i))
                  vo(i,j,k)=0.0+0.1/ru1*sin(2.*pi/2*ycf(j))
                  wo(i,j,k)=dummy2!+sin(2*2.*pi/2*ycf(j))
c                  wo(i,j,k)=1./ru1/8*.66!uinfty+0.2*uinfty*(ran1(idum)-0.5)
               enddo
            enddo
         enddo
         tlevel=0.
c
      ELSE IF (IFIELD.EQ.4) THEN !initialize to Blasius profile
c     
         dstar=0.0              !displacement thickness at inflow 
                                !x_ref--location of inflow (in the input file)
c--------begin bl generation
c     1st point  
            xloc=zw(1)          !u-location
            call blasius(xu, u1, w1, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar) 
            xloc=0.5*(3*zw(1)-zw(2)) !w-location
c     
            u1(1)=-u1(2)
c     
            call blasius(xu, u2, w2, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar)  
            do i=ix1, ix2+1
               do j=jy1+1, jy2
                  wo(i,j,1)=u1(i)
                  uo(i,j,1)=w2(i)
                  vo(i,j,1)=0.
               enddo
            enddo
c     the rest	    
         do k=kz1+1, kz2+1
            xloc=zw(k)          !u-location
            call blasius(xu, u1, w1, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar) 
            xloc=0.5*(zw(k)+zw(k-1)) !w-location
c     
            u1(1)=-u1(2)
c     
            call blasius(xu, u2, w2, nx, nx, nx, ru1, uinfty, xloc, x_ref, dstar)  
            do i=ix1, ix2+1
               do j=jy1+1, jy2
                  wo(i,j,k)=u1(i)
                  uo(i,j,k)=w2(i)
                  vo(i,j,k)=0.
               enddo
            enddo
         enddo
c     --------Boundary conditions
            do j=jy1+1,jy2
               do k=kz1+1,kz2
                  uo(ix1,j,k)=0. !wall
                  wo(ix1,j,k)=-wo(ix1+1,j,k) !wall
                  vo(ix1,j,k)=-vo(ix1+1,j,k) !wall
               enddo
            enddo
           do i=ix1+1, ix2
              do k=kz1+1, kz2
                 uo(i,jy1,k)=uo(i,jy2,k)
                 uo(i,jy2+1,k)=uo(i,jy1+1,k)
                 vo(i,jy1,k)=vo(i,jy2,k)
                 wo(i,jy1,k)=wo(i,jy2,k)
                 wo(i,jy2+1,k)=wo(i,jy1+1,k)
              enddo
           enddo
c         
      ELSE if (ifield.eq.2) then
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOSCALAR(STR4,P,NX,NY,NZ,NZG,idir,TLEVEL)
         call IOVECTOR(STR3,UAB1,VAB1,WAB1,NX,NY,NZ,NZG,idir,TLEVEL)
#ifdef LES
         if (read_ilm_imm.and.Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif

      else if (ifield.eq.22) then	 
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOSCALAR(STR4,P,NX,NY,NZ,NZG,idir,TLEVEL)
#ifdef LES
         if (read_ilm_imm) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
c
      else if (ifield.eq.33) then	 
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
c        set pressure to 0 since it's not available
         do k=1,nz
	   do j=1,ny
	     do i=1,nx
	       p(i,j,k)=0.
	     enddo 
	   enddo
	 enddo          
#ifdef LES
         if (read_ilm_imm) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
c
      else if (ifield.eq.44) then
         idir  = -1
         call IOVECTOR(STR1,uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOVECTOR(STR3,UAB1,VAB1,WAB1,NX,NY,NZ,NZG,idir,TLEVEL)
c        set pressure to 0 since it's not available
         do k=1,nz
	   do j=1,ny
	     do i=1,nx
	       p(i,j,k)=0.
             enddo  
	   enddo
	 enddo   
#ifdef LES
         if (read_ilm_imm.and.Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         elseif (Isgsmodel.eq.30) then
	  if (read_vt_tv) then  
           call ioscalar(vt_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(tv_file,tv,NX,NY,NZ,NZG,idir,TLEVEL)
	  else 
	   do k=1,nz
	    do j=1,ny
	     do i=1,nx
	      Ilm(i,j,k)=ru1
	      tv(i,j,k)=ru1
	     enddo
	    enddo
	   enddo
	  endif  
         endif
#endif
      endif
c
      if (save) then 
        if (myrank.eq.cpusave) write(59) tlevel,tstep,nx,ny,(xu(i),i=1,nx),(yv(j),j=1,ny)
      endif
c     
      IF(MYRANK.EQ.0) THEN
         write(*,*) 'starting time=',tlevel
         write(*,*) 'Initial field'
      ENDIF
      CALL SCRUVWP(uo,vo,wo,p,nx,ny,nz)
      
c     CALL CALCTKE(tke,uo,vo,wo,urun,vrun,wrun,zw,xu,nx,ny,nz)
      
c...  set source terms if needed
      if (myrank.eq.0) then 
        dpdz=0.
        dpdy=0.
        IF (ftype .EQ. 1) then 
         if (itype(5).eq.500.and.itype(6).eq.500) then
          dpdx=0.
	  dpdz=-1. 
         else
          dpdx = -1.
	  dpdz=0.
         endif
        endif
      endif	
      call MPI_BCAST(dpdx,1,type,0,MPI_COMM_HYDRO,ierr)
      call MPI_BCAST(dpdy,1,type,0,MPI_COMM_HYDRO,ierr)
      call MPI_BCAST(dpdz,1,type,0,MPI_COMM_HYDRO,ierr)
c       	  
c
      IF (ftype .EQ. 2) dpdx =  0.0
      IF (ftype .EQ. 3) dpdx =  0.0
      IF (ftype .EQ. 4) dpdx =  0.0
c     i.e. assume running vertically if top and bottom are periodic
       
c     set nwrite
c     
      if(MYRANK .EQ. 0) then

         open(file='nwrite.dat',unit=67,form='formatted',
     $        status='old',err=1999)
         read(67,*) nwrite0
         rewind(67)
         close(67)
         go to 2000
 1999    continue
         open(file='nwrite.dat',unit=67,form='formatted',
     $        status='new',err=1999)
         nwrite0=0
         write(67,*) nwrite0
         close(67)
 2000    continue

      endif
      CALL MPI_BCAST(nwrite0,1,MPI_INTEGER,0,MPI_COMM_HYDRO,IERR)
c     
      Lndslice=.false.
#ifndef OUT
c...  fringe method
      IF (Lfrg) CALL FRGINIT(WFRG,xu,zw,nx,ny,nz)
c...  open file for writing space averages
      OPEN (166,status='unknown',form='formatted',
     &     file=str6                          )
c...  open file for sgs model stats
      OPEN (167,status='unknown',form='formatted',
     &     file=str8                          )
      
c.... open slice-data file...when slice data BC is selected  
      IF(Lslice) then
         ic=0
         t_old=-5.
         ksub_old=0
         CALL SetUpSlice
         OPEN(111,file=str10,status='unknown',form='unformatted')
         IF(iread.eq.1) then
            CALL AbsTime(tlevel,t_old,ksub)
            do i=1,ksub
               READ(111) u2d,v2d,w2d
            enddo
            write(6,*) 'read...ksub=',ksub
         ENDIF  
      ENDIF
      
c.... open second slice-data file...
      Lndslice = .false.
      nsubnd   =  200
      iplanend =  90
      IF(Lndslice) then
         OPEN(112,file=str9,status='unknown',form='unformatted')       
      ENDIF
#endif
c------------------------------------------------------------------------
c     Check Divergence of Starting Field
c-----------------------------------------------------------------------
      CALL CALCDIV(UO,VO,WO,NX,NY,NZ)
c-----------------------------------------------------------------------
c     setup external forcing
c-----------------------------------------------------------------------
#ifdef immbnd ! conditional compilation
      if (ibswitch) then
c.... obstacle location on the grid
               CALL ObjectGrdCyl(sb,xb,yb,coefx,coefy,nb,nnmax)
c.... locate boundary points and interpolation coeff. for u
               CALL bpoints(xu,zc,zcg,sb,xb,yb,nb,coefx,coefy,
     &                   uniml,uiiml,ujiml,untop,utk,ustk1,ustk2,usti,unbot,ubk,utib_ind,utib_c1,utib_c2,ubib_ind,
     &                   ubib_c1,ubib_c2,
     &                   uc4intl)
c.... locate boundary points and interpolation coeff. for v
               CALL bpoints(xc,zc,zcg,sb,xb,yb,nb,coefx,coefy,
     &                   vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vsti,vnbot,vbk,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                   vbib_c1,vbib_c2,
     &                   vc4intl)
c.... locate boundary points and interpolation coeff. for w
               CALL bpoints(xc,zw,zwg,sb,xb,yb,nb,coefx,coefy,
     &                   wniml,wiiml,wjiml,wntop,wtk,wstk1,wstk2,wsti,wnbot,wbk,wtib_ind,wtib_c1,wtib_c2,wbib_ind,
     &                   wbib_c1,wbib_c2,
     &                   wc4intl)
     
c
c      debug
c      write(myrank+200,'(I,6G20.10)') (i,(uc4intl(myrank,i,k),k=1,5),uc4intl(myrank,i,1)+uc4intl(myrank,i,2)+uc4intl(myrank
c     &     ,i,3)+uc4intl(myrank,i,4)+uc4intl(myrank,i,5)
c     &     ,i=1,uniml(myrank))
c      close(myrank+200)
c
c
c      write(myrank+100,'(2I3,7G20.10)') (utib_ind(myrank,i),utk(myrank,i),(uc4intl(myrank,utib_ind(myrank,i),k),k=1,5),
c     & utib_c1(myrank,i),utib_c2(myrank,i),i=1,untop(myrank))
c      write(myrank+100,*) '********************************************************'
c      write(myrank+100,'(2I3,7G20.10)') (ubib_ind(myrank,i),ubk(myrank,i),(uc4intl(myrank,ubib_ind(myrank,i),k),k=1,5),
c     &  ubib_c1(myrank,i),ubib_c2(myrank,i),i=1,unbot(myrank))
c      close(myrank+100)
c      return
c      debug
      
c...  number of points to be forced
               write(0,'("ib points: (me)",I2,": u-grid=",2 I5)') myrank,uniml(myrank)
               write(0,'("ib points: (me)",I2,": v-grid=",2 I5)') myrank,vniml(myrank)
               write(0,'("ib points: (me)",I2,": w-grid=",2 I5)') myrank,wniml(myrank)
c...  set up immersed boundary velocities
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
#endif      
c-----------------------------------------------------------------------
c................initialize counters...................
C.....istart:	contains number of previous iterations
C.....icycle:	iteration counter
C.....itmax:	number of iterations
C.....itini:	start space-time sampling
C.....infnum:	information counter
C.....ipart:	particle tracking counter

      ICOUNT1=0
      ICOUNT2=0
      ICINF  =0
      if (.not.box2) istart=0
      ICYCLE = istart
      ipart  =0
      icnd   =0
c-----------------------------------------------------------------------
c     compute coefficients
c-----------------------------------------------------------------------
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '*..compute coefficients'
      ENDIF
c
#ifdef LES
c
      dtm1=1. ! should be ineffectual in this first call , but included to avoid an undefined variable error
      CALL turvis(uo,vo,wo,dxdydz,zw,urun,
     &     sxx,syy,szz,sxy,sxz,syz,g,
     &     uoc,voc,woc,upr,vpr,wpr,gbar,
     &     lm,mm,Ilm,Imm,xcf,ycf,zcf,
     &     wrk1,wrk2,wrk3,wrk4,
     &     tv,dtm1,icycle,nx,ny,nz       )
#ifdef immbnd
      if(ibswitch.and.Isgsmodel.ne.0) then
c
             CALL Forcing(tv,tvim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
c...  numerator for largangian dynamic model
           if(Isgsmodel.eq.20) then
c
             CALL Forcing(Ilm,Ilmim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
             CALL REFRESHBC(Ilm,nx*ny,(nz-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
           endif
         CALL REFRESHBC(tv, nx*ny,(nz-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      endif
#endif      
#endif
#ifndef LES
c      DNS ONLY!!!!!!!
      do i=1,nx
        do j=1,ny
	  do k=1,nz
	    tv(i,j,k)=0.
	  enddo
	enddo
      enddo
#endif      	            

C-----------------------------------------------------------------------
      IF(MYRANK.EQ.0) THEN
         WRITE(6,*) '*..start flow'
      ENDIF
C-----------------------------------------------------------------------
C     begin main loop
C-----------------------------------------------------------------------
C     
      oldtime=0.
 1000 icycle = icycle + 1
c
      
      znewtime=mpi_wtime()
      if (myrank.eq.0) then
       open(666,file='fort.666',form='formatted',status='unknown')
       write(666,'(I7,G25.15)') icycle, znewtime-oldtime
       close(666)
      endif
      oldtime=znewtime
c
      cpurhs=0
      cpupredict=0
      cpupres=0
      cpucorr=0
      cputurvis=0
c.....time step
      if(Lcfl) then
         CALL TSTEPCEN(XU,YV,ZW,UO,VO,WO,TV,nx,ny,nz,DTM1)
      else
         dtm1=tstep
      endif
c     for phase reshuffling
      if (jdose) then 
       if (myrank.eq.0) then 
        dummy1=mod(tlevel+zw(1)/uinfty,zlen)
        dummy2=mod(tlevel+dtm1+zw(1)/uinfty,zlen)
c	write(0,*) dummy1,dummy2
	if (dummy2.lt.dummy1) then 
c        reshuffle phases
         write(0,*) 'reshuffling phases'
         do kk=1,jcounter
	  if (count(kk).ne.0) then 
c u
	    i=int(ran1(idum)*count(kk))+1
	    j=is(kk,i)
	    temp=uh(kk)
	    uh(kk)=uh(j)
	    uh(j)=temp	
c v
	    i=int(ran1(idum)*count(kk))+1
	    j=is(kk,i)
	    temp=vh(kk)
	    vh(kk)=vh(j)
	    vh(j)=temp	
c w
	    i=int(ran1(idum)*count(kk))+1
	    j=is(kk,i)
	    temp=wh(kk)
	    wh(kk)=wh(j)
	    wh(j)=temp
	  endif
         enddo
        endif
       endif
       call MPI_BCAST(uh(1),3000,type,0,MPI_COMM_HYDRO,ierr)		     	
       call MPI_BCAST(vh(1),3000,type,0,MPI_COMM_HYDRO,ierr)		     	
       call MPI_BCAST(wh(1),3000,type,0,MPI_COMM_HYDRO,ierr)		     	
      endif
c
      if (box2) then 
        if (myrank.eq.0) then 
          if ((icycle-istart).gt.1.and.mod(icycle,500).eq.1) then 
 	   loadnum=loadnum+1
	   open(63, file='./box2/box2.fstle.6-'//index(loadnum),form='binary',status='old')
	   read(63) t_b2,nx_b2,ny_b2,(x_b2(i),i=1,nx_b2),(y_b2(j),j=1,ny_b2)
	   if (lcfl) write(0,*) 'WARNING: box2 with cfl ON'
	   if (ny_b2.ne.ny) write(0,*) 'WARNING: box2 and box1 spanwise domains do not match'
	   if (nx_b2.ne.nb2) write(0,*) 'WARNING: box2 and box1 transverse domains do not match'
          endif !icycle
	  read(63) dt_b2
	  dtm1=dt_b2
        endif !myrank
	call MPI_BCAST(dtm1,1,type,0,MPI_COMM_HYDRO,ierr)
      endif !box2 

      TLEVEL=TLEVEL+DTM1
C-----------------------------------------------------------------------
C     predictor-corrector step
C-----------------------------------------------------------------------
      ITSTAMP = ITMAX
      ICOUNT1 = MOD(ICYCLE,100000)
      
      IF ((ICOUNT1.EQ.0).OR.((IFIELD.NE.2).AND.(ICYCLE.EQ.1))) THEN
c...  Euler step
         a1   =0.
         a2   =1.
         tmeth=1.
         CALL Rhs(uo,vo,wo,ucn,vcn,wcn,tv,wfrg,uab2,vab2,wab2,nx,ny,nz )        
      ELSE
c...  Adams-Bashforth step
         a1   =0.5
         a2   =1.5
         tmeth=1.
         cpub=mpi_wtime()
         CALL Rhs(uo,vo,wo,ucn,vcn,wcn,tv,wfrg,uab2,vab2,wab2,nx,ny,nz )   
         cpurhs=mpi_wtime()-cpub
      ENDIF
      
C.....Predictor step + boundary conditions for Ustar
C     
      if(impfla.eq.1.and.icycle.eq.1) then
c...  boundary condition setup for implicit part

         CALL BOUN_IMP(botbc,topbc,nx,ny,nz)
      endif
c
      cpub=mpi_wtime()
      CALL predictor(UO,VO,WO,P,UAB1,VAB1,WAB1,
     &     UAB2,VAB2,WAB2,ucn,vcn,wcn,tv,
     &     USTAR,VSTAR,WSTAR,xu,
     &     dpdx,dpdy,dpdz,dtm1,a1,a2,botbc,topbc,
     &     nx,ny,nz                          )
      Ltstep = .FALSE.
c
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbc'
      CALL REFRESHBC(USTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(VSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(WSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'    
ccccccIMMERSED BOUNDARY
#ifdef immbnd
c....modify viscosity at the immersed boundary points
c    DEAL WITH THIS LATER
c      if(ibswitch.and.ftype.eq.1) then !remove pressure from boundary points
c         do i=1,nbodyu
c            ustar(ibodyu(i),j,kbodyu(i)) = ustar(ibodyu(i),j,kbodyu(i))
c     &           +dpdx*dtm1
c         enddo 
c         do i=1,nbodyv
c            vstar(ibodyv(i),j,kbodyv(i)) = vstar(ibodyv(i),j,kbodyv(i))
c     &           +dpdy*dtm1
c         enddo 
c      endif
c
c....compute forcing in x momentum
      if(ibswitch) then
              CALL Forcing(ustar,uim,
     &                    uniml,uiiml,ujiml,untop,utk,ustk1,ustk2,unbot,ubk,usti,utib_ind,utib_c1,utib_c2,ubib_ind,
     &                    ubib_c1,ubib_c2,
     &                    uc4intl) 
c.... compute forcing in y momentum
              CALL Forcing(vstar,vim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
c.... compute forcing in z momentum
              CALL Forcing(wstar,wim,
     &                    wniml,wiiml,wjiml,wntop,wtk,wstk1,wstk2,wnbot,wbk,wsti,wtib_ind,wtib_c1,wtib_c2,wbib_ind,
     &                    wbib_c1,wbib_c2,
     &                    wc4intl) 
      endif
c
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbc after immersed boundary'
      CALL REFRESHBC(USTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(VSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(WSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'    
#endif  
ccccccIMMERSED BOUNDARY
c     need to call boundary because refreshbc changes wstar at 1st point!!!
      CALL Boundary(USTAR,VSTAR,WSTAR,UO,VO,WO,
     &     u2d,v2d,w2d,tlevel,t_old,ksub,nt,
     &     urun,wrun,
     &     XU,YV,ZW,nx,ny,nz,Ltstep,uinf,duinf,dtm1,x_ref)
c
      goto 12345
c     CODE BELOW IS A HARDWIRING!!!!
c     random forcing if necessary ....! ugo's artifice
c
      if (icycle .le. 2) then
         do i=2,ix2
            amp = 0.1
            env = 0.25*amp*(1+TANH(4.*(xu(i)-1.)))
     $           *(1.-tanh(xu(i)-10.))
c            write(6,*) k,zw(k),env
            do j=2,jy2
               do k=2,kz2
	          if (zw(k).gt.12) then ! put noise behind the cylinder only 
                   ustar(i,j,k) = ustar(i,j,k) + env*(ran1(idum)-0.5)
                   vstar(i,j,k) = vstar(i,j,k) + env*(ran1(idum)-0.5)
                   wstar(i,j,k) = wstar(i,j,k) + env*(ran1(idum)-0.5)
		  endif 
               enddo
            enddo      
         enddo
      endif
c  CODE ABOVE IS A HARDWIRING!!!!!!
12345  continue
c     
      if (inflow) then
       if (myrank.eq.0) then  
c     generate shirani inflow (call 8020 for uniform inflow)
         ux=zw(1)+tlevel*uinfty
         vx=zcf(1)+tlevel*uinfty
c
         do j=2,jy2
            do i=iistart,ifin
               ustar(i,j,1)=ustar(i,j,1)/uscale    
               wstar(i,j,1)=wstar(i,j,1)/uscale
               vstar(i,j,1)=vstar(i,j,1)/uscale
               do kk=1,counter
                  wstar(i,j,1)=wstar(i,j,1)+ampl(kk)*sin(ux*k1(kk)+phase1u(kk))*sin(ycf(j)*k2(kk)+phase2u(kk))*
     &                 sin((xcf(k)-zstart)*k3(kk)+phase3u(kk))
                  vstar(i,j,1)=vstar(i,j,1)+ampl(kk)*sin(vx*k1(kk)+phase1v(kk))*sin(yv(j)*k2(kk)+phase2v(kk))*
     &                 sin((xcf(k)-zstart)*k3(kk)+phase3v(kk))
                  ustar(i,j,1)=ustar(i,j,1)+ampl(kk)*sin(vx*k1(kk)+phase1w(kk))*sin(ycf(j)*k2(kk)+phase2w(kk))*
     &                 sin((xu(k)-zstart)*k3(kk)+phase3w(kk))
               enddo
               ustar(i,j,1)=ustar(i,j,1)*uscale
               vstar(i,j,1)=vstar(i,j,1)*uscale
               wstar(i,j,1)=wstar(i,j,1)*uscale 
            enddo
         enddo	   
        endif
      endif
c     
      if (jdose) then 
c     generate inflow with os wall-normal basis
      	 do i=2,nx-1
          do j=2,ny-1
           uose(j,i)=0.    
           vose(j,i)=0.    
           wose(j,i)=0.
          enddo
         enddo
c       	  
       x_u=tlevel+zwg(1)/uinfty
       x_w=tlevel+0.5*(3.*zwg(1)-zwg(2))/uinfty

       do kk=1,jcounter

         xh_u=exp((0.,1.)*kx(kk)*x_u)
         xh_w=exp((0.,1.)*kx(kk)*x_w)
	        
         call ose_(Re_,kx(kk)*L_/uinfty,kz(kk)*L_,U_,d2U,phi,phip,eta,detacen,detacor,nose)
c
c     NOTE:FLIPPED CASE -- distribute load between cpus
         do iii=nxloc*myrank+2,min(nxloc*(myrank+1)+1,nx-1)
	  i=ps_u(iii) ! u and v coords
	  zh1=(a_u(iii)*phip(i)+b_u(iii)*phip(i+1))/L_
	  i=ps_w(iii) ! w coords
	  zh2=(a_w(iii)* phi(i)+b_w(iii)* phi(i+1))
             do jjj=2,ny-1
                y_v=yv(jjj)
		y_u=ycf(jjj)
                yh_v=exp((0.,1.)*ky(kk)*y_v)
                yh_u=exp((0.,1.)*ky(kk)*y_u)
	        v1=xh_u*yh_u*zh1
		v2=xh_w*yh_v*zh1
		v3=xh_w*yh_u*zh2
                wose(jjj,iii)=wose(jjj,iii)+2.*real(uh(kk)*v1)
                vose(jjj,iii)=vose(jjj,iii)+2.*real(vh(kk)*v2)
                uose(jjj,iii)=uose(jjj,iii)+2.*real(wh(kk)*v3)
             enddo
         enddo
       enddo  
c
         if (myrank.eq.0) then 
           do k=1,mysize-1
	     call MPI_RECV(uose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,0,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(vose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,1,MPI_COMM_HYDRO,status,ierr)
	     call MPI_RECV(wose(1,nxloc*k+2),ny*min(nxloc,nx-1-(nxloc*k+1)),type,k,2,MPI_COMM_HYDRO,status,ierr)
	   enddo
	 else
	     call MPI_SEND(uose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,0,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(vose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,1,MPI_COMM_HYDRO,ierr)
	     call MPI_SEND(wose(1,nxloc*myrank+2),ny*min(nxloc,nx-1-(nxloc*myrank+1)),type,0,2,MPI_COMM_HYDRO,ierr)
         endif         
c
         if (myrank.eq.0) then 
c     add the perturbation
	   do i=2,nx-1 
            do j=2,ny-1
             ustar(i,j,1)=ustar(i,j,1)+uose(j,i)*uscale*pscale    
             vstar(i,j,1)=vstar(i,j,1)+vose(j,i)*uscale*pscale    
             wstar(i,j,1)=wstar(i,j,1)+wose(j,i)*uscale*pscale
	    enddo
	   enddo    
        endif ! myrank
      endif !jdose

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          

        if (box2) then 
	 if (myrank.eq.0) then 
	  do i=2,nx_b2-1
	   read(63) (u_b2(j,i), v_b2(j,i), w_b2(j,i),j=2,ny-1)
          enddo !i
c   bcs   
          do j=2,ny-1 
           u_b2(j,1)=-u_b2(j,2)
	   v_b2(j,1)=-v_b2(j,2)	       
	   w_b2(j,1)=0.	       
           u_b2(j,nb2)=u_b2(j,nb2-1)
	   v_b2(j,nb2)=v_b2(j,nb2-1)	       
	   w_b2(j,nb2)=w_b2(j,nb2-1)
c   interpolate
	  enddo !j
	  do i=2,nx-1
	   do j=2,ny-1
	    i_cen=ps_u(i)
	    i_cor=ps_w(i)
	    ustar(i,j,1)=a_w(i)*u_b2(j,i_cor)+b_w(i)*u_b2(j,i_cor+1)
	    vstar(i,j,1)=a_u(i)*v_b2(j,i_cen)+b_u(i)*v_b2(j,i_cen+1)
	    wstar(i,j,1)=a_u(i)*w_b2(j,i_cen)+b_u(i)*w_b2(j,i_cen+1)
	   enddo !j
	  enddo !i 	       
	 endif !myrank
	endif !box2         

c   *** wbl box2-specific

c
      cpupredict=mpi_wtime()-cpub
c      write(0,*) myrank,'called predictor'
c     
      call MasBlas(USTAR,VSTAR,WSTAR,xu,yv,zw,icycle,nx,ny,nz)
c     
      cpub=mpi_wtime()
      CALL Divergence(USTAR,VSTAR,WSTAR,nx,ny,nz,DTM1,dp)
c
      if(IPresMethod .eq. 1) then
         CALL PRESDP2P(dp,dp,P,ZWG,YV,XU,DUMMY,
     &        MYRANK,MPI_COMM_HYDRO,TPRESSURE)
      elseif(IPresMethod.eq.11) then 
         call pfft_yz(p,dp,dp,dummy,dummy,au,ap,bp(1),cwg(1),myrank,MPI_COMM_HYDRO,type)
      elseif(IPresMethod .eq. 22) then
c         write(0,*) myrank,'calling pressure'

         CALL presdirect(dp,dp,p,xu,yv,zwg,au,cwg,icycle-istart,
     &        myrank,mpi_comm_hydro,tpressure)
      else
         write(0,*) ' Undefined pressure solver '
         call mpi_finalize
         return
      endif
c
      cpupres=mpi_wtime()-cpub
C     
C.....Update velocity values , compute residuals and store RHS
C     
      cpub=mpi_wtime()
      CALL Correct(USTAR,VSTAR,WSTAR,DP,DTM1,nx,ny,nz,UO,VO,WO,
     &     ICYCLE,RSU,RSV,RSW)
C     
c      write(0,*) myrank,'called corrector'
C-----------------------------------------------------------------------
C     boundary conditions to the updated values
C-----------------------------------------------------------------------
C     
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbs'
      CALL REFRESHBC(UO,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(VO,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(WO,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
c     .....divergence
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
      CALL Storen(UAB1,VAB1,WAB1,UAB2,VAB2,WAB2,tmeth,nx,ny,nz)
c     
      Ltstep = .TRUE. 
      CALL Boundary(UO,VO,WO,USTAR,VSTAR,WSTAR,
     &     u2d,v2d,w2d,tlevel,t_old,ksub,nt,
     &     urun,wrun,
     &     XU,YV,ZW,nx,ny,nz,Ltstep,uinf,duinf,dtm1,x_ref)
c     
c     wbl - box2 specific -- boundary condition
      if (myrank.eq.0) then 
       do i=1,nx
        do j=1,ny
	 uo(i,j,1)=ustar(i,j,1)
	 vo(i,j,1)=vstar(i,j,1)
	 wo(i,j,1)=wstar(i,j,1)
        enddo
       enddo
      endif 	
c     wbl - box2 specific   
c
      cpucorr=mpi_wtime()-cpub
c     
      if (save) then 
         nsave=nsave+1
c	 call MPI_BARRIER(MPI_COMM_HYDRO,ierr)
         if (myrank.eq.cpusave) then 
      	    if ((icycle-istart).gt.1.and.mod(icycle-istart,500).eq.1) then
	      savenum=savenum+1 
      	      close(59)
      	      open(59,file='./box3/box3.fstle-'//index(savenum),form='unformatted',status='unknown')
              write(59) tlevel,tstep,nx,ny,(xu(i),i=1,nx),(yv(j),j=1,ny)
	    endif
c	       
c            open(59, file='inflow/144x180inf6%.'//index(nsave), form='unformatted', status='unknown')
            do i=2,nx-1
               write(59) (uo(i,j,ksave), vo(i,j,ksave), wo(i,j,ksave),j=2,ny-1)
            enddo
         endif
      endif   
      

C-----------------------------------------------------------------------
C     eddy viscosity  and forcing
C-----------------------------------------------------------------------
c.... compute runing average of velocities for inflow bc
c     CALL RunAve(uo,vo,wo,urun,vrun,wrun,xu,yv,zw,
c     &            tlevel,dtm1,icycle,nx,ny,nz       )

C.....turbulent viscosity
      cpub=mpi_wtime()
#ifdef LES
c      write(0,*) myrank,' calling turvis'
      CALL turvis(uo,vo,wo,dxdydz,zw,urun,
     &     sxx,syy,szz,sxy,sxz,syz,g,
     &     uoc,voc,woc,upr,vpr,wpr,gbar,
     &     lm,mm,Ilm,Imm,xcf,ycf,zcf,
     &     wrk1,wrk2,wrk3,wrk4,
     &     tv,dtm1,icycle,nx,ny,nz       )
c      write(0,*) myrank,' back from turvis'
#endif
      cputurvis=mpi_wtime()-cpub
c      
#ifdef immbnd
#ifdef LES
c
c....modify viscosity at the immersed boundary points
      if(ibswitch.and.Isgsmodel.ne.0) then

             CALL Forcing(tv,tvim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
c
c...  numerator and denominator for largangian dynamic model
           if(Isgsmodel.eq.20) then
c
             CALL Forcing(Ilm,Ilmim,
     &                    vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vnbot,vbk,vsti,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                    vbib_c1,vbib_c2,
     &                    vc4intl) 
             CALL REFRESHBC(Ilm,nx*ny,(nz-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
           endif
         CALL REFRESHBC(tv, nx*ny,(nz-2),MYSIZE,MYRANK,MPI_COMM_HYDRO)
      endif
c      write(0,*) myrank,' back from les-immbnd loop'
#endif 
#endif     
c     output for streamwise spectra
      if (spec) then
       if (mod(icycle,specfreq).eq.0) then 
        j=ny/2  
        do n=1,nspec
	 if (myrank.eq.spec_cpu(n)) then
	   write(600+n) tlevel,(uo(i,j,spec_k(n)),i=2,nx-1,skip),
     &	                       (vo(i,j,spec_k(n)),i=2,nx-1,skip), 
     &                         (wo(i,j,spec_k(n)),i=2,nx-1,skip),
     &	                       (p (i,j,spec_k(n)),i=2,nx-1,skip)
	 endif
	enddo
       endif	
      endif	     
c     
C     ... Particle trajectory
c     if(koutput.gt.0) then     !send the field to the receiving process (rank=0 in comm_world)
c     if(myrank.eq.0) then   !send time and time step to receiving proc.
c     call MPI_SEND(dtm1,1,TYPE,0,0,MPI_COMM_WORLD,STATUS,IERR)
c     call MPI_SEND(tlevel,1,TYPE,0,0,MPI_COMM_WORLD,STATUS,IERR)
c     endif
c     CALL MPI_SEND(UO(1,1,1),NX*NY*NZ,TYPE,0,0,
c     &        MPI_COMM_WORLD,STATUS,IERR)
c     CALL MPI_SEND(VO(1,1,1),NX*NY*NZ,TYPE,0,0,
c     &        MPI_COMM_WORLD,STATUS,IERR)
c     CALL MPI_SEND(WO(1,1,1),NX*NY*NZ,TYPE,0,0,
c     &        MPI_COMM_WORLD,STATUS,IERR)
c     endif
C.....forcing terms
C     CALL Forterm(UAV,YV,ZW,QO,QOLD,DPDX,DTM1,NY,NZ)
C     CALL FORTERMNEW(UAV,YV,ZW,QO,QOLD,DPDX,DTM1,nx,NY,NZ)
c     endif
C     
C---------------------------------------------------------------------
C.... write slice-data file..
C---------------------------------------------------------------------
C     
      if(Lslice .and. iread .eq. 0) THEN 
         ic=ic+1
         if(ic.eq.nsub) then
            CALL StoreSlice(uo,vo,wo,u2d,v2d,w2d,ic,nt,nx,ny,nz)
            ic=0
            ksub=ksub+1
            write(111) u2d,v2d,w2d
         else
            CALL StoreSlice(uo,vo,wo,u2d,v2d,w2d,ic,nt,nx,ny,nz)
         endif
      endif
      
c.....write a 2nd slice data set while reading the first one
      IF(Lndslice) THEN
         icnd=icnd+1
         if(icnd.eq.nsubnd) then
            do j=jy1,jy2+1
               do k=kz1,kz2+1
                  u2dnd(icnd,j,k)=uo(iplanend,j,k)
                  v2dnd(icnd,j,k)=vo(iplanend,j,k)
                  w2dnd(icnd,j,k)=wo(iplanend,j,k)
               enddo
            enddo       
            icnd=0
            write(112) u2dnd,v2dnd,w2dnd
         else
            do j=jy1,jy2+1
               do k=kz1,kz2+1
                  u2dnd(icnd,j,k)=uo(iplanend,j,k)
                  v2dnd(icnd,j,k)=vo(iplanend,j,k)
                  w2dnd(icnd,j,k)=wo(iplanend,j,k)
               enddo
            enddo     
         endif       
      ENDIF
c     
c     
 555  continue
      if(.not.CalcStat) then
cccccccccccccccccccccccccccc OUTPUT one xy PLANE cccccccccccccccccccccccccccccccccccccccccccc
         if(mod(ICYCLE-istart-1,infnum).eq.0) then
 1234            j=ny/2
	    if (myrank.eq.0) then 
	    write(66+mod(icycle,5),*) myrank,': matlab output'
               open(62,file=STR5,status='unknown',form='unformatted')
c     write(60,*) 'variables="x","z","u","v","w"'
c     write(60,*) 'zone i=',nx-2,',j=',nz-2,'f=point'
               write(62) nx-2,nzg-2 
               write(62) (xcf(i),i=2,ix2)		 
               write(62) (zcg(k),k=2,nzg-1)
               do k=2,kz2
                  write(62) (   !xcf(i),zcf(k),
     &                 0.5*(uo(i,j,k)+uo(i-1,j,k)),
     &                 vo(i,j,k),
     &                 0.5*(wo(i,j,k)+wo(i,j,k-1)),
     &                 p(i,j,k),
     &                          !char(10),
     &                 i=2,ix2)
               enddo
               do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(u2dxz(1,1,1),nx*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dxz(1,1,1),nx*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dxz(1,1,1),nx*nz,type,n,nn+2,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(p2dxz(1,1,1),nx*nz,type,n,nn+3,mpi_comm_hydro,status,ierr)
c     
                  do k=2,nz-1
                     write(62) (0.5*(u2dxz(i,k,1)+u2dxz(i-1,k,1)),
     &                          v2dxz(i,k,1),
     &                          0.5*(w2dxz(i,k,1)+w2dxz(i,k-1,1)),
     &                          p2dxz(i,k,1),
     &                          i=2,nx-1)
                  enddo    
               enddo		           	       
               write(62) ycf(j)
               close(62)
	    write(66+mod(icycle,5),*) myrank,': done'
            else
	       nn=myrank*10
               call MPI_Send(uo(1,j,1),1,real2d,0,nn,  MPI_COMM_HYDRO,ierr)
               call MPI_Send(vo(1,j,1),1,real2d,0,nn+1,MPI_COMM_HYDRO,ierr)
               call MPI_Send(wo(1,j,1),1,real2d,0,nn+2,MPI_COMM_HYDRO,ierr)
               call MPI_Send(p (1,j,1),1,real2d,0,nn+3,MPI_COMM_HYDRO,ierr)	      
            endif !myrank
         endif !mod
      endif ! .not.CalcStat
C     
C-----------------------------------------------------------------------
C     create a series of restart files
C-----------------------------------------------------------------------
C     
      
      LcreateR=.true.
      if(MOD(ICYCLE,ievery).EQ.0) then
c.... write restart file
         idir  = 1
	 
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': restart file'
         call IOVECTOR('res.int.vel_f',uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
         call IOSCALAR('res.int.p_f',p,nx,ny,nz,nzg,idir,tlevel)
         call IOVECTOR('res.int.ab_f',uab1,vab1,wab1,nx,ny,nz,nzg,idir,tlevel)
#ifdef LES
         if (Isgsmodel.eq.20) then 
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
	 
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
c     
c     write latest inflow phase information(note: this file overwrites the original)
         if (jdose) then 
          if (myrank.eq.0) then        
           open(101,file=NOISEFILE,form='formatted',status='unknown')
           write(101,'(I5)') jcounter,kkmax
             do k=1,jcounter 
               write(101,'(2 I5,9G25.15)') k,count(k),kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
               do ipos=1,count(k)
                 write(101,'(I5)') is(k,ipos)
               enddo	  	       
             enddo
            close(101)
	  endif
            open(101,file='pscale.dat',form='formatted',status='old')
              write(101,'(G25.15)') pscale
	      write(101,'(I6)') psample
	    close(101)  
	 endif     
c      write 3D velocity .........................
c     
         goto 999 !!OBSOLESCENT -- RECOMMENDED FOR DELETION
         if(LcreateR) then
            nwrite=nwrite0+(icycle/ievery)
c$$$  
c$$$  write(6,1995)  icycle,icycle/ievery,myrank,nwrite0,nwrite
c$$$  1995       format(' into LcreateR if; icycle, myrank, nwrite =',
c$$$  $           10i5)

            if(myrank .eq. 0) then 
               open(file='nwrite.dat',unit=67,form='formatted',
     $              status='unknown',err=1999)
               rewind(67)
               write(67,*) nwrite
               close(67)
            endif
            call MPI_BARRIER(MPI_COMM_HYDRO,IERR)
c$$$  write(6,*) myrank,nwrite
            call IOVSS('eddywbl3D.'//index(nwrite),
     $           uo,vo,wo,p,tv,nx,ny,nz,nzg,idir,tlevel)
         endif
         
c     
c     write 2D file
c     
c$$$  call statistiche(uo,vo,wo,p,'eddywbl2D.'//index(nwrite),
c$$$  &           nx,ny,nz,tlevel)
c     
 999  endif
      continue 
C     
C-----------------------------------------------------------------------
C     write intermediate restart files
C-----------------------------------------------------------------------
c     
      ICOUNT2 = MOD(ICYCLE,INFNUM)
      IF((ICOUNT2.EQ.0)) THEN
         IF(MYRANK.EQ.0) THEN
C-----------------------------------------------------------------------
C     screen and file information
C-----------------------------------------------------------------------
            write(6,*) '*********************************************'
            write(6,*) '*...Iteration number= ',icycle
            write(6,*) '*...Time Step       = ',dtm1
            write(6,*) '*...Total Time      = ',tlevel
            write(6,*) '*...CPU RHS         = ',cpurhs
            write(6,*) '*...CPU PRED        = ',cpupredict
            write(6,*) '*...CPU PRESS       = ',cpupres
            write(6,*) '*...CPU CORR        = ',cpucorr
            write(6,*) '*...CPU TURVIS      = ',cputurvis
         ENDIF
         CALL CALCDIV(uo,vo,wo,nx,ny,nz)
         CALL CALCME(uo,vo,wo,xu,yv,zw,dtm1,nx,ny,nz)
         CALL SCRUVWP(uo,vo,wo,p,nx,ny,nz)
      ENDIF
C     
C-----------------------------------------------------------------------
C     end main loop
C-----------------------------------------------------------------------
C     
C.....check counters
c
      if (CalcStat.eq..true..and.mod(icycle,ISTATS).eq.1) then 
         samples =samples+1      
         if (samples.gt.statsample) then ! write new stats files
           statnum=statnum+1
           samples=1
	 endif  	  

c     goto 333
c     compute uinfinity -- all processors, then gather on all CPUs
            uact=0.
            iuinf=ix2-3    
            do k=kz1,kz2
               uint=0. 
               do i=iuinf,ix2
                  do j=jy1+1,jy2
                     uint=uint+wo(i,j,k)/ap(i)
                  enddo
               enddo
c	       
               uint=uint/(xu(ix2)-xu(iuinf-1))/(jy2-jy1)
	       uinft(k)=uinft(k)+(uint-uinft(k))/samples
               uact=uact+uinft(k)/cw(k) ! for integration over x
            enddo
c     gather
            if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': compute uinft,uact'
            call MPI_GATHER(uinft(2),1,real1d,uinftg(2),1,real1d,0,MPI_COMM_HYDRO,ierr)
c     include first point (optional)
            if (myrank.eq.0) uinftg(1)=uinft(1)
	    call MPI_BCAST(uinftg(1),nzg,type,0,MPI_COMM_HYDRO,ierr)	     

            call MPI_ALLREDUCE(uact,uactg,1,TYPE,MPI_SUM,MPI_COMM_HYDRO,ierr)
c     get actual uinfinity
            fac=1./(1.5*zwg(nzg-1)-0.5*zwg(nzg-2)-(1.5*zwg(1)-0.5*zwg(2)))
            uactg=uactg*fac
      
            if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'	
c	    
c     compute statistics
c     
            do k=kz1,kz2
               u=0.
               v=0.
               w=0.
               uu=0.
               vv=0.
               ww=0.
c               wwd=0.
c               vvd=0.
c     
               do i=imid,imidf
                  do j=jy1+1,jy2
c     spatial averaging	  
                     uu=uu+(uo(i,j,k))**2 
                     u=u+uo(i,j,k) 
                     vv=vv+(vo(i,j,k))**2 
                     v=v+vo(i,j,k) 
                     ww=ww+(wo(i,j,k))**2 
                     w=w+wo(i,j,k) 
                  enddo
               enddo
               np=(jy2-jy1)*(imidf-imid+1)
               uu=uu/np
               u=u/np
               v=v/np
               vv=vv/np
               w=w/np
               ww=ww/np
c     
               uav(k)=uav(k)+(u-uav(k))/samples
               vav(k)=vav(k)+(v-vav(k))/samples
               wav(k)=wav(k)+(w-wav(k))/samples
               uuav(k)=uuav(k)+(uu-uuav(k))/samples
               vvav(k)=vvav(k)+(vv-vvav(k))/samples
               wwav(k)=wwav(k)+(ww-wwav(k))/samples
c     
            enddo 
	    call MPI_GATHER(uav(1),1,real1d,guav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
	    call MPI_GATHER(vav(1),1,real1d,gvav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
	    call MPI_GATHER(wav(1),1,real1d,gwav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
	    call MPI_GATHER(uuav(1),1,real1d,guuav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
	    call MPI_GATHER(vvav(1),1,real1d,gvvav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   
	    call MPI_GATHER(wwav(1),1,real1d,gwwav(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)   

            if (myrank.eq.0) then 
c     recomupute scaling parameter
               uosm=sqrt(uuav(1)-uav(1)**2)
               vosm=sqrt(vvav(1)-vav(1)**2)
               wosm=sqrt(wwav(1)-wav(1)**2)
	       	    
             if (jdose) then 
	       pscale=(pscale*psample+pscale*plevel*3./(uosm+vosm+wosm))
	       psample=psample+1
	       pscale=pscale/psample
               write(0,'("inlet noise rms: u,v,w, scale ",4 G20.10,I)') uosm,vosm,wosm,pscale,psample
             endif
 
	    endif   
	       call MPI_BCAST(pscale,1,type,0,MPI_COMM_HYDRO,ierr)
	       
c	       
            if (myrank.eq.mysize-1) then ! first processor writes this stat file only
c     include last point
	       guav(nzg-1)=uav(kz2)
	       gvav(nzg-1)=vav(kz2)
	       gwav(nzg-1)=wav(kz2)
	       guuav(nzg-1)=uuav(kz2)
	       gvvav(nzg-1)=vvav(kz2)
	       gwwav(nzg-1)=wwav(kz2)
c	       
               open(21,file='h'//index(statnum)//'.gnu',status='unknown',form='formatted')
               write(21,'(4 F20.15)') (zwg(k)-x_ref, gwwav(k)-gwav(k)**2, gvvav(k)-gvav(k)**2, 
     &                         	       guuav(k)-guav(k)**2,k=1,nzg-1)
               close(21)
c	       	       
            endif
c  
c      goto 333
 222     continue
c     compute turbulent boundary layer statistics
c     
         if (myrank.eq.mysize-1) then 
            open(22,file='shape'//index(statnum)//'.gnu',status='unknown',form='formatted')
            open(23,file='stress'//index(statnum)//'.gnu',status='unknown',form='formatted')
            open(27,file='cf'//index2(icycle)//'.gnu',status='unknown',form='formatted')

         endif
         np=jy2-jy1
         do k=kz1,kz2
            do i=ix1,ix2
               u=0.
               v=0.
               w=0.
	       pp=0.
               uu=0.
               ww=0.
               vv=0.
               ppp=0.
               uw=0.
#ifdef LES	       
	       tvi=0.
#endif
               do j=jy1+1,jy2
                  u=u+uo(i,j,k)
                  v=v+vo(i,j,k)
                  w=w+wo(i,j,k)
		  pp=pp+p(i,j,k)
                  uu=uu+uo(i,j,k)**2
                  vv=vv+vo(i,j,k)**2
                  ww=ww+wo(i,j,k)**2
		  ppp=ppp+p(i,j,k)**2
                  uw=uw+0.25*(uo(i+1,j,k+1)+uo(i,j,k+1))*(wo(i+1,j,k+1)+wo(i+1,j,k))! at k+1 point
#ifdef LES
                  tvi=tvi+tv(i,j,k+1)
#endif
               enddo  
               u=u/np
               v=v/np
               w=w/np
	       pp=pp/np
               uu=uu/np
               vv=vv/np
               ww=ww/np
      	       ppp=ppp/np
               uw=uw/np
               uavg(i,k)=uavg(i,k)+(u-uavg(i,k))/samples
               vavg(i,k)=vavg(i,k)+(v-vavg(i,k))/samples
               wavg(i,k)=wavg(i,k)+(w-wavg(i,k))/samples
	       pavg(i,k)=pavg(i,k)+(pp-pavg(i,k))/samples
               uuavg(i,k)=uuavg(i,k)+(uu-uuavg(i,k))/samples
               vvavg(i,k)=vvavg(i,k)+(vv-vvavg(i,k))/samples
               wwavg(i,k)=wwavg(i,k)+(ww-wwavg(i,k))/samples
               ppavg(i,k)=ppavg(i,k)+(ppp-ppavg(i,k))/samples       
               uwavg(i+1,k+1)=uwavg(i+1,k+1)+(uw-uwavg(i+1,k+1))/samples
#ifdef LES
	       tvi=tvi/np
	       tvavg(i,k+1)=tvavg(i,k+1)+(tvi-tvavg(i,k+1))/samples
#endif
            enddo
            
c     compute displacement thickness, momentum thickness, viscous units, skin friction coeff.
c     
            ubulk(k)=0.
            theta(k)=0.
            dst(k)=0.
	    winst(k)=0.
            do i=ix1+1,ix2
               ubulk(k)=ubulk(k)+wavg(i,k)/ap(i)
               dst(k)=dst(k)+(uinft(k)-wavg(i,k))/ap(i)
               theta(k)=theta(k)+wavg(i,k)*(uinft(k)-wavg(i,k))/ap(i)
            enddo
            ubulk(k)=ubulk(k)/(xu(ix2)-xu(ix1)) 
            dst(k)=dst(k)/uinft(k)
            theta(k)=theta(k)/(uinft(k)**2)
            do j=jy1+1,jy2
		winst(k)=winst(k)+wo(2,j,k)
	    enddo
	    winst(k)=winst(k)/np	
         enddo
c     
c     SUM CONTRIBUTIONS FROM DIFFERENT PROCESSORS
          if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': reducing stats'
         call MPI_GATHER(ubulk(1),1,real1d,ubulkg(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(dst(1),1,real1d,dstg(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(theta(1),1,real1d,thetag(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
c        for instantaneous cf curve
         call MPI_GATHER(winst(1),1,real1d,winstg(1),1,real1d,mysize-1,MPI_COMM_HYDRO,ierr)
c        
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': gathering stats'
c     gather data on the last processor for load balancing
         call MPI_GATHER(uavg(1,1),1,stat2d,uavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(vavg(1,1),1,stat2d,vavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(wavg(1,1),1,stat2d,wavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
	 call MPI_GATHER(pavg(1,1),1,stat2d,pavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(uuavg(1,1),1,stat2d,uuavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(vvavg(1,1),1,stat2d,vvavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(wwavg(1,1),1,stat2d,wwavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
	 call MPI_GATHER(ppavg(1,1),1,stat2d,ppavgg(1,1),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
         call MPI_GATHER(uwavg(1,2),1,stat2d,uwavgg(1,2),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr) ! from 2nd point
#ifdef LES
         call MPI_GATHER(tvavg(1,2),1,stat2d,tvavgg(1,2),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
#endif
c       
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
c     
c *************     -- only if one needs to know te stats at k=kz2, otherwise comment
c     include the last point (k=kz2)
       if (myrank.eq.mysize-1) then
          ubulkg(nzg-1)=ubulk(kz2)
          dstg(nzg-1)  =dst(kz2)
          thetag(nzg-1)=theta(kz2)
	  winstg(nzg-1)=winst(kz2)
          do i=ix1+1,ix2 
            uavgg(i,nzg-1) =uavg(i,kz2)
            vavgg(i,nzg-1) =vavg(i,kz2)
            wavgg(i,nzg-1) =wavg(i,kz2)
            uuavgg(i,nzg-1)=uuavg(i,kz2)
            vvavgg(i,nzg-1)=vvavg(i,kz2)
            wwavgg(i,nzg-1)=wwavg(i,kz2)
            ppavgg(i,nzg-1)=ppavg(i,kz2)
            uwavgg(i,1)=0. ! because stats here cannot be computed (we start at k=2)
            wavgg(i,0)=0. ! for consistency with output wavgg(i,0) has to be defined
#ifdef LES
            tvavgg(i,1)=0. ! tv is unknown at first point
#endif
            pavgg (i,nzg-1)=pavg (i,kz2)
       	 enddo
       endif	 
c ************
c	 
         if (myrank .eq.mysize-1) then
c           write number of streamwise stations
	    if (mod(nzg-1,statxfreq).gt.0) then 
              write(23,'(2I)') (nzg-1)/statxfreq+1,nx-2
	    else
	      write(23,'(2I)') (nzg-1)/statxfreq,nx-2
	    endif    
        
            do k=kz1,nzg-1 
               tauwg(k)=ru1*wavgg(2,k)*au(1)*2. !second order
               utaug(k)=sqrt(tauwg(k))
               cfg(k)=tauwg(k)*2./uactg**2
ccc 	for instantaneous cf curve
	       winstg(k)=ru1*winstg(k)*au(1)*2. !second order tau_wall
	       winstg(k)=winstg(k)*2./uactg**2  !uactg is NOT instant, but should be ok
ccc	        
	       kpg(k)=-ru1/(0.5*(uinftg(k)+uinftg(k+1)))**2*(uinftg(k+1)-uinftg(k))/(zwg(k+1)-zwg(k))
c     
               write(22,'(6 G20.12)') (zwg(k)-x_ref)*uactg/ru1,cfg(k),dstg(k),thetag(k),kpg(k),zwg(k)-x_ref
               write(27,'(2 G20.12)') (zwg(k)-x_ref)*uactg/ru1,winstg(k)
ccccccccccccccc the data below _can_ be used but is obsolete
               if (mod(k,statxfreq).eq.1) then ! NOTE: make sure (above) stats at position k=1 are _DEFINED_
                  write(23,*) char(10),char(10)
                  write(23,*) k, zwg(k)
                  do i=ix1+1,ix2
                     write(23,'(19 G20.12)') xu(i)*utaug(k)/ru1,xu(i)/dstg(k),0.5*(xu(i)+xu(i-1))*utaug(k)/ru1,
     &                    uavgg(i,k)/utaug(k),vavgg(i,k)/utaug(k),wavgg(i,k)/utaug(k),
     &                    (uuavgg(i,k)-uavgg(i,k)**2)/tauwg(k),(vvavgg(i,k)-vavgg(i,k)**2)/tauwg(k),
     &                    (wwavgg(i,k)-wavgg(i,k)**2)/tauwg(k),
     &                    (uwavgg(i,k)-uavgg(i,k)*wavgg(i,k))/tauwg(k),  
c                          normalize by macro scales
     &                    (uavgg(i,k)/uactg),(vavgg(i,k)/uactg),(wavgg(i,k)/uactg),
     &                    (uuavgg(i,k)-uavgg(i,k)**2)/uactg**2,(vvavgg(i,k)-vavgg(i,k)**2)/uactg**2,
     &                    (wwavgg(i,k)-wavgg(i,k)**2)/uactg**2,
     &                    (uwavgg(i,k)-0.25*(uavgg(i,k)+uavgg(i-1,k))*(wavgg(i,k)+wavgg(i,k-1)))/uactg**2 !data invalid at 1st point
#ifdef LES 
     &                    ,tvavgg(i,k)/ru1
#endif     
                  enddo
               endif ! mod		
            enddo !k
c	        	      
            close(22)
            close(23)
	    close(27)
c     write averaged field
               open(62,file='avg'//index(statnum)//'.mat',status='unknown',form='unformatted')
c              SAME OUTPUT -- ROTATE GRAPH	       
               write(62) nx-2,nzg-2 
               write(62) (xcf(i),i=2,ix2)		 
               write(62) (zcg(k),k=2,nzg-1)
               do k=2,nzg-1
                  write(62) (
     &                 0.5*(uavgg(i,k)+uavgg(i-1,k)),
     &                 vavgg(i,k),
     &                 0.5*(wavgg(i,k)+wavgg(i,k-1)),
     &                 pavgg(i,k),
     &                          !char(10),
     &                 i=2,ix2)
               enddo
               close(62)
c     write averaged _stress_ field
               open(62,file='avg_str'//index(statnum)//'.mat',status='unknown',form='unformatted')
c              SAME OUTPUT -- ROTATE GRAPH	       
               write(62) nx-2,nzg-2 
               write(62) (xcf(i),i=2,ix2)		 
               write(62) (zcg(k),k=2,nzg-1)
               do k=2,nzg-1
                  write(62) (
     &                 0.5*(uuavgg(i,k)+uuavgg(i-1,k)-uavgg(i,k)**2-uavgg(i-1,k)**2),
     &                 vvavgg(i,k)-vavgg(i,k)**2,
     &                 0.5*(wwavgg(i,k)+wwavgg(i,k-1)-wavgg(i,k)**2-wavgg(i,k-1)**2),
     &                 uwavgg(i,k)-0.25*(uavgg(i,k)+uavgg(i-1,k))*(wavgg(i,k)+wavgg(i,k-1)),
     &                 ppavgg(i,k)-pavgg(i,k)**2,
     &                          !char(10),
     &                 i=2,ix2)
               enddo
               close(62)
c     write _averaged_ vorticity field
               open(62,file='vortavg'//index(statnum)//'.mat',status='unknown',form='unformatted')
c              SAME OUTPUT -- ROTATE GRAPH	       
               write(62) nx-4,nzg-3 !x-dir lacks 1st and last pt; z-dir lacks last pt 
               write(62) (xcf(i),i=3,ix2-1)		 
               write(62) (zcg(k),k=2,nzg-2)
               do k=2,nzg-2
                  write(62) (
     &                  0.5*((vavgg(i,k)-vavgg(i,k-1))*cwg(k-1)
     &                      +(vavgg(i,k+1)-vavgg(i,k))*cwg(k)),  !omx
     &                -0.25*((uavgg(i,k)-uavgg(i,k-1)+uavgg(i-1,k)-uavgg(i-1,k-1))*cwg(k-1)+
     &                       (uavgg(i,k+1)-uavgg(i,k)+uavgg(i-1,k+1)-uavgg(i-1,k))*cwg(k)-
     &			     (wavgg(i,k)-wavgg(i-1,k)+wavgg(i,k-1)-wavgg(i-1,k-1))*au(i-1)-
     &			     (wavgg(i+1,k)-wavgg(i,k)+wavgg(i+1,k-1)-wavgg(i,k-1))*au(i)),  !omy
     &                 -0.5*((vavgg(i,k)-vavgg(i-1,k))*au(i-1)+
     &			     (vavgg(i+1,k)-vavgg(i,k))*au(i)),
#ifdef LES
     &                        tvavgg(i,k),
#endif
     &                 i=3,ix2-1)
               enddo
               close(62)
         endif	!myrank
         	       
cccccccccccccccccccccccccccc PERIODIC OUTPUT OF PLANES cccccccccccccccccccccccc
c     includes U V W P OMX OMY OMZ NUT
      if (planes) then
cccccccc xz
      do m=1,nyplanes
        j=1
	do while (yv(j).le.nyploc(m).and.j.le.jy2)
         j=j+1
	enddo
      	 
         do k=2,kz2
	  do i=2,ix2
           u2dxz(i,k,2)=-0.25*((wo(i,j,k)-wo(i,j-1,k)+wo(i,j,k-1)-wo(i,j-1,k-1))*bv(j-1)+
     &                    (wo(i,j+1,k)-wo(i,j,k)+wo(i,j+1,k-1)-wo(i,j,k-1))*bv(j)-
     &                    (vo(i,j,k)-vo(i,j,k-1)+vo(i,j-1,k)-vo(i,j-1,k-1))*cw(k-1)-
     &                    (vo(i,j,k+1)-vo(i,j,k)+vo(i,j-1,k+1)-vo(i,j-1,k))*cw(k))  !omx
           v2dxz(i,k,2)=-0.25*((uo(i,j,k)-uo(i,j,k-1)+uo(i-1,j,k)-uo(i-1,j,k-1))*cw(k-1)+
     &                    (uo(i,j,k+1)-uo(i,j,k)+uo(i-1,j,k+1)-uo(i-1,j,k))*cw(k)-
     &                    (wo(i,j,k)-wo(i-1,j,k)+wo(i,j,k-1)-wo(i-1,j,k-1))*au(i-1)-
     &		          (wo(i+1,j,k)-wo(i,j,k)+wo(i+1,j,k-1)-wo(i,j,k-1))*au(i))  !omy
           w2dxz(i,k,2)=-0.25*((vo(i,j,k)-vo(i-1,j,k)+vo(i,j-1,k)-vo(i-1,j-1,k))*au(i-1)+
     &		          (vo(i+1,j,k)-vo(i,j,k)+vo(i+1,j-1,k)-vo(i,j-1,k))*au(i)-
     &                    (uo(i,j,k)-uo(i,j-1,k)+uo(i-1,j,k)-uo(i-1,j-1,k))*bv(j-1)-
     &                    (uo(i,j+1,k)-uo(i,j,k)+uo(i-1,j+1,k)-uo(i-1,j,k))*bv(j))  ! omz
          enddo
	 enddo 
ccccccccccccccccccccccc
ccccccccccccccccccccccc
         if (myrank.eq.0) then !1st processor writes datafile
           open(62,file='planes/fstle.b2.'//index(m)//'.xz-'//index2(icycle),status='new',form='unformatted')
           write(62) nx-2,nzg-2  
           write(62) (xcf(i),i=2,ix2)		 
           write(62) (zcg(k),k=2,nzg-1)
               do k=2,kz2
                  write(62) (0.5*(uo(i,j,k)+uo(i-1,j,k)),vo(i,j,k),0.5*(wo(i,j,k)+wo(i,j,k-1)),p(i,j,k),                   
     &		             u2dxz(i,k,2),v2dxz(i,k,2),w2dxz(i,k,2),                
#ifdef LES
     &                       tv(i,j,k),                 
#endif
     &                                i=2,ix2)
              enddo
c              receive velocity from other processors
           do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(u2dxz(1,1,1),nx*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dxz(1,1,1),nx*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dxz(1,1,1),nx*nz,type,n,nn+2,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(p2dxz(1,1,1),nx*nz,type,n,nn+3,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(u2dxz(1,1,2),nx*nz,type,n,nn+4,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dxz(1,1,2),nx*nz,type,n,nn+5,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dxz(1,1,2),nx*nz,type,n,nn+6,mpi_comm_hydro,status,ierr)
#ifdef LES
                  call MPI_Recv(p2dxz(1,1,2),nx*nz,type,n,nn+7,mpi_comm_hydro,status,ierr)
#endif
c     
                  do k=2,nz-1
                    write(62) (0.5*(u2dxz(i,k,1)+u2dxz(i-1,k,1)),v2dxz(i,k,1),0.5*(w2dxz(i,k,1)+w2dxz(i,k-1,1)),p2dxz(i,k,1),
     &                         u2dxz(i,k,2),v2dxz(i,k,2),w2dxz(i,k,2),
#ifdef LES
     &		               p2dxz(i,k,2),
#endif     
     &		    i=2,ix2)
                  enddo
            enddo		           	       
	    write(62) ycf(j) ! write position last		           	       
            close(62)
         else
           nn=myrank*10
           call MPI_Send(uo(1,j,1),1,real2d,0,nn,  MPI_COMM_HYDRO,ierr)
           call MPI_Send(vo(1,j,1),1,real2d,0,nn+1,MPI_COMM_HYDRO,ierr)
           call MPI_Send(wo(1,j,1),1,real2d,0,nn+2,MPI_COMM_HYDRO,ierr)
           call MPI_Send(p (1,j,1),1,real2d,0,nn+3,MPI_COMM_HYDRO,ierr)	      
           call MPI_Send(u2dxz(1,1,2),nx*nz,type,0,nn+4,MPI_COMM_HYDRO,ierr)
           call MPI_Send(v2dxz(1,1,2),nx*nz,type,0,nn+5,MPI_COMM_HYDRO,ierr)
           call MPI_Send(w2dxz(1,1,2),nx*nz,type,0,nn+6,MPI_COMM_HYDRO,ierr)
#ifdef LES
           call MPI_Send(tv (1,j,1),1,real2d,0,nn+7,MPI_COMM_HYDRO,ierr)	      
#endif
         endif !(myrank.eq.0) 
      enddo !m
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
cccccccc yz
      do m=1,nxplanes
        i=1
	do while (xu(i).le.nxploc(m).and.i.le.ix2)
         i=i+1
	enddo
      	 
         do k=2,kz2
	  do j=2,jy2
           u2dyz(j,k,2)=-0.25*((wo(i,j,k)-wo(i,j-1,k)+wo(i,j,k-1)-wo(i,j-1,k-1))*bv(j-1)+
     &                    (wo(i,j+1,k)-wo(i,j,k)+wo(i,j+1,k-1)-wo(i,j,k-1))*bv(j)-
     &                    (vo(i,j,k)-vo(i,j,k-1)+vo(i,j-1,k)-vo(i,j-1,k-1))*cw(k-1)-
     &                    (vo(i,j,k+1)-vo(i,j,k)+vo(i,j-1,k+1)-vo(i,j-1,k))*cw(k))  !omx
           v2dyz(j,k,2)=-0.25*((uo(i,j,k)-uo(i,j,k-1)+uo(i-1,j,k)-uo(i-1,j,k-1))*cw(k-1)+
     &                    (uo(i,j,k+1)-uo(i,j,k)+uo(i-1,j,k+1)-uo(i-1,j,k))*cw(k)-
     &                    (wo(i,j,k)-wo(i-1,j,k)+wo(i,j,k-1)-wo(i-1,j,k-1))*au(i-1)-
     &		          (wo(i+1,j,k)-wo(i,j,k)+wo(i+1,j,k-1)-wo(i,j,k-1))*au(i))  !omy
           w2dyz(j,k,2)=-0.25*((vo(i,j,k)-vo(i-1,j,k)+vo(i,j-1,k)-vo(i-1,j-1,k))*au(i-1)+
     &		          (vo(i+1,j,k)-vo(i,j,k)+vo(i+1,j-1,k)-vo(i,j-1,k))*au(i)-
     &                    (uo(i,j,k)-uo(i,j-1,k)+uo(i-1,j,k)-uo(i-1,j-1,k))*bv(j-1)-
     &                    (uo(i,j+1,k)-uo(i,j,k)+uo(i-1,j+1,k)-uo(i-1,j,k))*bv(j))  ! omz
          enddo
	 enddo 
ccccccccccccccccccccccc
ccccccccccccccccccccccc
         if (myrank.eq.0) then !1st processor writes datafile
           open(62,file='planes/fstle.b2.'//index(m)//'.yz-'//index2(icycle),status='new',form='unformatted')
           write(62) ny-2,nzg-2  
           write(62) (ycf(j),j=2,jy2)		 
           write(62) (zcg(k),k=2,nzg-1)
               do k=2,kz2
                  write(62) (uo(i,j,k),0.5*(vo(i,j,k)+vo(i,j-1,k)),0.5*(wo(i,j,k)+wo(i,j,k-1)),p(i,j,k),                   
     &		             u2dyz(j,k,2),v2dyz(j,k,2),w2dyz(j,k,2),                
#ifdef LES
     &                       tv(i,j,k),                 
#endif
     &                                j=2,jy2)
              enddo
c              receive velocity from other processors
           do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(u2dyz(1,1,1),ny*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dyz(1,1,1),ny*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dyz(1,1,1),ny*nz,type,n,nn+2,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(p2dyz(1,1,1),ny*nz,type,n,nn+3,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(u2dyz(1,1,2),ny*nz,type,n,nn+4,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2dyz(1,1,2),ny*nz,type,n,nn+5,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2dyz(1,1,2),ny*nz,type,n,nn+6,mpi_comm_hydro,status,ierr)
#ifdef LES
                  call MPI_Recv(p2dyz(1,1,2),ny*nz,type,n,nn+7,mpi_comm_hydro,status,ierr)
#endif
c     
                  do k=2,nz-1
                    write(62) (u2dyz(j,k,1),0.5*(v2dyz(j,k,1)+v2dyz(j-1,k,1)),0.5*(w2dyz(j,k,1)+w2dyz(j,k-1,1)),p2dyz(j,k,1),
     &                         u2dyz(j,k,2),v2dyz(j,k,2),w2dyz(j,k,2),
#ifdef LES
     &		    p2dxz(j,k,2),
#endif     
     &		    j=2,jy2)
                  enddo
            enddo
	    write(62) xcf(i) ! write position last		           	       
            close(62)
         else
           nn=myrank*10
           call MPI_Send(uo(i,1,1),1,cross2d,0,nn,  MPI_COMM_HYDRO,ierr)
           call MPI_Send(vo(i,1,1),1,cross2d,0,nn+1,MPI_COMM_HYDRO,ierr)
           call MPI_Send(wo(i,1,1),1,cross2d,0,nn+2,MPI_COMM_HYDRO,ierr)
           call MPI_Send(p (i,1,1),1,cross2d,0,nn+3,MPI_COMM_HYDRO,ierr)	      
           call MPI_Send(u2dyz(1,1,2),ny*nz,type,0,nn+4,MPI_COMM_HYDRO,ierr)
           call MPI_Send(v2dyz(1,1,2),ny*nz,type,0,nn+5,MPI_COMM_HYDRO,ierr)
           call MPI_Send(w2dyz(1,1,2),ny*nz,type,0,nn+6,MPI_COMM_HYDRO,ierr)
#ifdef LES
           call MPI_Send(tv (i,1,1),1,cross2d,0,nn+7,MPI_COMM_HYDRO,ierr)	      
#endif
         endif !(myrank.eq.0) 
      enddo !m
      
cccccccc surface pressure/skin friction
c
         i=2
ccccccccccccccccccccccc
         if (myrank.eq.0) then !1st processor writes datafile
           open(62,file='planes/fstle.b2.pcf0-'//index2(icycle),status='new',form='unformatted')
           write(62) ny-2,nzg-2  
           write(62) (ycf(j),j=2,jy2)		 
           write(62) (zcg(k),k=2,nzg-1)
               do k=2,kz2
                  write(62) (p(i,j,k),ru1*0.5*(wo(i,j,k)+wo(i,j,k-1))*au(1)*2.*2./uinfty**2,j=2,jy2) !p & cf
               enddo
c              receive velocity from other processors
               do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(w2dyz(1,1,1),ny*nz,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(p2dyz(1,1,1),ny*nz,type,n,nn+1,mpi_comm_hydro,status,ierr)
c     
                  do k=2,nz-1
                    write(62) (p2dyz(j,k,1),ru1*0.5*(w2dyz(j,k,1)+w2dyz(j,k-1,1))*au(1)*2.*2./uinfty**2,j=2,jy2)
                  enddo
            enddo		           	       
            write(62) xcf(i)
            close(62)
         else
           nn=myrank*10
           call MPI_Send(wo(i,1,1),1,cross2d,0,nn,  MPI_COMM_HYDRO,ierr)
           call MPI_Send(p (i,1,1),1,cross2d,0,nn+1,MPI_COMM_HYDRO,ierr)	      
         endif !(myrank.eq.0) 
ccccccccccccccccccccccccccccc
cccccccccccccc xy ccccccccccc
      do m=1,nzplanes
       if (myrank.eq.cpuxy(m)) then 
        k=kxy(m)
c
         do j=2,jy2
	  do i=2,ix2
           u2dxy(i,j)=-0.25*((wo(i,j,k)-wo(i,j-1,k)+wo(i,j,k-1)-wo(i,j-1,k-1))*bv(j-1)+
     &                    (wo(i,j+1,k)-wo(i,j,k)+wo(i,j+1,k-1)-wo(i,j,k-1))*bv(j)-
     &                    (vo(i,j,k)-vo(i,j,k-1)+vo(i,j-1,k)-vo(i,j-1,k-1))*cw(k-1)-
     &                    (vo(i,j,k+1)-vo(i,j,k)+vo(i,j-1,k+1)-vo(i,j-1,k))*cw(k))  !omx
           v2dxy(i,j)=-0.25*((uo(i,j,k)-uo(i,j,k-1)+uo(i-1,j,k)-uo(i-1,j,k-1))*cw(k-1)+
     &                    (uo(i,j,k+1)-uo(i,j,k)+uo(i-1,j,k+1)-uo(i-1,j,k))*cw(k)-
     &                    (wo(i,j,k)-wo(i-1,j,k)+wo(i,j,k-1)-wo(i-1,j,k-1))*au(i-1)-
     &		          (wo(i+1,j,k)-wo(i,j,k)+wo(i+1,j,k-1)-wo(i,j,k-1))*au(i))  !omy
           w2dxy(i,j)=-0.25*((vo(i,j,k)-vo(i-1,j,k)+vo(i,j-1,k)-vo(i-1,j-1,k))*au(i-1)+
     &		          (vo(i+1,j,k)-vo(i,j,k)+vo(i+1,j-1,k)-vo(i,j-1,k))*au(i)-
     &                    (uo(i,j,k)-uo(i,j-1,k)+uo(i-1,j,k)-uo(i-1,j-1,k))*bv(j-1)-
     &                    (uo(i,j+1,k)-uo(i,j,k)+uo(i-1,j+1,k)-uo(i-1,j,k))*bv(j))  ! omz
          enddo
	 enddo 
ccccccccccccccccccccccc
ccccccccccccccccccccccc
           open(62,file='planes/fstle.b2.'//index(m)//'.xy-'//index2(icycle),status='new',form='unformatted')
           write(62) nx-2,ny-2  
           write(62) (xcf(i),i=2,ix2)		 
           write(62) (ycf(j),j=2,ny-1)
               do j=2,jy2
                  write(62) (0.5*(uo(i,j,k)+uo(i-1,j,k)),0.5*(vo(i,j,k)+vo(i,j-1,k)),wo(i,j,k),p(i,j,k),                   
     &		             u2dxy(i,j),v2dxy(i,j),w2dxy(i,j),                
#ifdef LES
     &                       tv(i,j,k),                 
#endif
     &                                i=2,ix2)
               enddo
               write(62) zcf(kxy(m))
               close(62)
       endif !myrank
      enddo !m
      
      endif !planes  
ccccccccccccccccccccccccccccccccREYNOLDS STRESS BUDGETS -- SEE 'RSBUD1.F cccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (rstat) then
c
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': computing reynolds stress budgets'
c        compute stress budget data
         call rsbud(rsblk_1(1,1,1),rsblk_1(1,1,2),rsblk_1(1,1,3),rsblk_1(1,1,4),rsblk_1(1,1,5),rsblk_1(1,1,6),rsblk_1(1,1,7),
     &              rsblk_1(1,1,8),rsblk_1(1,1,9),rsblk_1(1,1,10),rsblk_1(1,1,11),rsblk_1(1,1,12),rsblk_1(1,1,13),rsblk_1(1,1,14),
     &              rsblk_1(1,1,15),rsblk_1(1,1,16),rsblk_1(1,1,17),rsblk_1(1,1,18),rsblk_1(1,1,19),rsblk_1(1,1,20),rsblk_1(1,1,21),
     &              rsblk_1(1,1,22),rsblk_1(1,1,23),rsblk_1(1,1,24),rsblk_1(1,1,25),rsblk_1(1,1,26),rsblk_1(1,1,27),rsblk_1(1,1,28),
#ifdef LES
     &              rsblk_3(1,1,1),rsblk_3(1,1,2),rsblk_3(1,1,3),rsblk_3(1,1,4),rsblk_3(1,1,5),rsblk_3(1,1,6),rsblk_3(1,1,7),
     &              rsblk_3(1,1,8),rsblk_3(1,1,9),rsblk_3(1,1,10),rsblk_3(1,1,11),rsblk_3(1,1,12),!12
#endif     
     &              rsblk_2(1,0,1),rsblk_2(1,0,2),rsblk_2(1,0,3),rsblk_2(1,0,4),rsblk_2(1,0,5),rsblk_2(1,0,6),rsblk_2(1,0,7),
     &              rsblk_2(1,0,8),rsblk_2(1,0,9),rsblk_2(1,0,10),rsblk_2(1,0,11),rsblk_2(1,0,12),rsblk_2(1,0,13),rsblk_2(1,0,14),    	           
     &              rsblk_2(1,0,15),rsblk_2(1,0,16),rsblk_2(1,0,17),rsblk_2(1,0,18),rsblk_2(1,0,19),rsblk_2(1,0,20),rsblk_2(1,0,21),
     &              rsblk_2(1,0,22),rsblk_2(1,0,23),rsblk_2(1,0,24),rsblk_2(1,0,25),rsblk_2(1,0,26),rsblk_2(1,0,27),rsblk_2(1,0,28),
     &              rsblk_2(1,0,29),rsblk_2(1,0,30),rsblk_2(1,0,31),rsblk_2(1,0,32),rsblk_2(1,0,33),!33
#ifdef LES
     &              rsblk_4(1,0,1),rsblk_4(1,0,2),rsblk_4(1,0,3),rsblk_4(1,0,4),rsblk_4(1,0,5),rsblk_4(1,0,6),rsblk_4(1,0,7),
     &              rsblk_4(1,0,8),rsblk_4(1,0,9),rsblk_4(1,0,10),rsblk_4(1,0,11),rsblk_4(1,0,12),rsblk_4(1,0,13),rsblk_4(1,0,14),    	           
     &              rsblk_4(1,0,15),rsblk_4(1,0,16),rsblk_4(1,0,17),rsblk_4(1,0,18),rsblk_4(1,0,19),rsblk_4(1,0,20),rsblk_4(1,0,21),
     &              rsblk_4(1,0,22),rsblk_4(1,0,23),rsblk_4(1,0,24),rsblk_4(1,0,25),
#endif
     &              uo,vo,wo,p,
#ifdef LES
     &              tv,
#endif          
     &              samples)
c
         do n=0,3     
	  do j=1,7
	    call MPI_GATHER(rsblk_1(1,2,n*7+j),1,stat2d,rsblkg_1(1,2,j),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
	  enddo
#ifdef LES
	  do j=1,3
	    call MPI_GATHER(rsblk_3(1,2,n*3+j),1,stat2d,rsblkg_3(1,2,j),1,stat2d,mysize-1,MPI_COMM_HYDRO,ierr)
	  enddo
#endif
	  
          goto 444 ! matlab output instead

	  if (myrank.eq.mysize-1) then 
c         write rs budget file
            open(26,file=rsnames(n),form='formatted',status='unknown')
	    do k=3,nzg-2,statxfreq
	        write(26,*) char(10),char(10)    
		write(26,*) k
	        do i=3,ix2-1
#ifndef LES
      	         write(26,'(8 G20.12)')  
#endif
#ifdef LES
      	         write(26,'(11 G20.12)')  
#endif
     &		                     xcf(i)*utaug(k)/ru1,-rsblkg_1(i,k,1)/tauwg(k)**2*ru1,-rsblkg_1(i,k,2)/tauwg(k)**2*ru1,
     &                                                   -rsblkg_1(i,k,3)/tauwg(k)**2*ru1,-rsblkg_1(i,k,4)/tauwg(k)**2*ru1,
     &                                                   -rsblkg_1(i,k,5)/tauwg(k)**2*ru1,-rsblkg_1(i,k,6)/tauwg(k)**2*ru1,
     &                                                   -rsblkg_1(i,k,7)/tauwg(k)**2*ru1
c     &          (-rsblkg_1(i,k,1)-rsblkg_1(i,k,2)-rsblkg_1(i,k,3)-rsblkg_1(i,k,4)-rsblkg_1(i,k,5)-rsblkg_1(i,k,6)-rsblkg_1(i,k,7))
c     &           /tauwg(k)**2*ru1,
#ifdef LES
     &                                                   ,-rsblkg_3(i,k,1)/tauwg(k)**2*ru1,rsblkg_3(i,k,2)/tauwg(k)**2*ru1
     &                                                   ,rsblkg_3(i,k,3)/tauwg(k)**2*ru1
#endif
                enddo !i
            enddo ! k
	    close(26)
c
          endif !myrank
 444      continue ! matlab output instead!!	  

            if (myrank.eq.mysize-1) then 
c	    
               open(62,file=rsnames(n)//index(statnum)//'.dat',status='unknown',form='unformatted')
               write(62) ru1
               write(62) nx-4,nzg-4 
               write(62) (xcf(i),i=3,ix2-1)		 
               write(62) (zcg(k),k=3,nzg-2)
               do k=3,nzg-2
                write(62) utaug(k),(
     &              -rsblkg_1(i,k,1)/tauwg(k)**2*ru1,-rsblkg_1(i,k,2)/tauwg(k)**2*ru1,
     &              -rsblkg_1(i,k,3)/tauwg(k)**2*ru1,-rsblkg_1(i,k,4)/tauwg(k)**2*ru1,
     &              -rsblkg_1(i,k,5)/tauwg(k)**2*ru1,-rsblkg_1(i,k,6)/tauwg(k)**2*ru1,
     &              -rsblkg_1(i,k,7)/tauwg(k)**2*ru1,
#ifdef LES
     &              -rsblkg_3(i,k,1)/tauwg(k)**2*ru1,rsblkg_3(i,k,2)/tauwg(k)**2*ru1,
     &               rsblkg_3(i,k,3)/tauwg(k)**2*ru1,
#endif
     &                 i=3,ix2-1)
               enddo
               close(62)
            endif !myrank
          enddo !n
	 endif !rstat
ccccccccccccccccccccccccccccccccREYNOLDS STRESS BUDGETS -- SEE 'RSBUD1.F cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccc two-point spanwise correlations cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (two_pt_corr) then 
	  call span_corr(uo,vo,wo,u2pt,v2pt,w2pt,um2pt,vm2pt,wm2pt,nxl,nzl,ks(myrank),ke(myrank),zfreq,xfreq,samples)
c      write output file       

c          subtract mean from correlations	    
	      do k=1,nzl
	       do i=1,nxl
	        do j=1,ny
		 u2tmp(i,j,k)=u2pt(i,j,k)-um2pt(i,k)**2
		 v2tmp(i,j,k)=v2pt(i,j,k)-vm2pt(i,k)**2
		 w2tmp(i,j,k)=w2pt(i,j,k)-wm2pt(i,k)**2
		enddo
		 u2tmp(i,1,k)=u2tmp(i,1,k)/u2tmp(i,2,k)
		 v2tmp(i,1,k)=v2tmp(i,1,k)/v2tmp(i,2,k)
		 w2tmp(i,1,k)=w2tmp(i,1,k)/w2tmp(i,2,k)
	        do j=ny,2,-1
		 u2tmp(i,j,k)=u2tmp(i,j,k)/u2tmp(i,2,k)
		 v2tmp(i,j,k)=v2tmp(i,j,k)/v2tmp(i,2,k)
		 w2tmp(i,j,k)=w2tmp(i,j,k)/w2tmp(i,2,k)
		enddo
	       enddo
	      enddo 	 

	    if (myrank.eq.0) then 
               open(62,file='two_pt_'//STR5,status='unknown',form='unformatted')
               write(62) nxl,ny-2,(nzg-3)/zfreq+1 
               write(62) (xn(i),i=1,nxl)
	       write(62) (ycf(j),j=2,ny-1)		 
               write(62) (zn(k),k=1,(nzg-3)/zfreq+1)
               do k=1,npt(myrank)
                  write(62) ((
     &                 u2tmp(i,j,k),
     &                 v2tmp(i,j,k),
     &                 w2tmp(i,j,k),
     &                 j=2,ny-1),i=1,nxl)
               enddo
               do n=1,mysize-1
	          nn=n*10
                  call MPI_Recv(u2tmp,nxl*nzl*ny,type,n,nn+0,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(v2tmp,nxl*nzl*ny,type,n,nn+1,mpi_comm_hydro,status,ierr)
                  call MPI_Recv(w2tmp,nxl*nzl*ny,type,n,nn+2,mpi_comm_hydro,status,ierr)
c     
                  do k=1,npt(n)
                    write(62) ((
     &                 u2tmp(i,j,k),
     &                 v2tmp(i,j,k),
     &                 w2tmp(i,j,k),
     &                 j=2,ny-1),i=1,nxl)
                  enddo
               enddo		           	       
               close(62)
	    write(66+mod(icycle,5),*) myrank,': done'
            else
c         now send	    
	       nn=myrank*10
               call MPI_Send(u2tmp,nxl*nzl*ny,type,0,nn,  MPI_COMM_HYDRO,ierr)
               call MPI_Send(v2tmp,nxl*nzl*ny,type,0,nn+1,MPI_COMM_HYDRO,ierr)
               call MPI_Send(w2tmp,nxl*nzl*ny,type,0,nn+2,MPI_COMM_HYDRO,ierr)
            endif 
         endif ! two-pt-corr	   
cccccccccccccccccccccccccccccccc two-point spanwise correlations cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	   
 333     continue 
      endif ! stat loop

      close(66+mod(icycle,5))      
      IF (icycle .LT. itmax) GO TO 1000
C     
      Niter = Itmax
      IF(MYRANK.eq.0) THEN
         WRITE(6,*) '*.........................'
         WRITE(6,*) '*..EDDY finished.'
      ENDIF
c     
c     write final restart file
c     
      ifield = 2
      idir   = 1
       if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': final restart'
      call IOVECTOR('res.fin.vel_f',uo,vo,wo,nx,ny,nz,nzg,idir,tlevel)
      call IOSCALAR('res.fin.p_f',p,nx,ny,nz,nzg,idir,tlevel)
      call IOVECTOR('res.fin.ab_f',uab1,vab1,wab1,nx,ny,nz,nzg,idir,
     &     tlevel)
#ifdef LES
         if (Isgsmodel.eq.20) then  
           call ioscalar(ilm_file,Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar(imm_file,Imm,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
c     write latest inflow phase information(note: this file overwrites the original)
         if (jdose) then 
          if (myrank.eq.0) then        
           open(101,file=NOISEFILE,form='formatted',status='unknown')
           write(101,'(I5)') jcounter,kkmax
             do k=1,jcounter 
               write(101,'(2 I5,9G25.15)') k,count(k),kx(k),kz(k),ky(k),uh(k),wh(k),vh(k)
               do ipos=1,count(k)
                 write(101,'(I5)') is(k,ipos)
               enddo	  	       
             enddo
            close(101)
c
            open(101,file='pscale.dat',form='formatted',status='old')
              write(101,'(G25.15)') pscale
	      write(101,'(I6)') psample
	    close(101)  
	  endif
	 endif     
c
       if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,' done'
      call MPI_TYPE_FREE(real1d,ierr)
      call MPI_TYPE_FREE(real2d,ierr)
      call MPI_TYPE_FREE(stat2d,ierr)
      call MPI_TYPE_FREE(cross2d,ierr)
c     
      if(myrank .eq. 0) then 
         close(67)
	 close(66+mod(icycle,5))
      endif
      if (save) then 
         if (myrank.eq.cpusave) close(59)
      endif
      if (spec) then 
	 do i=1,nspec
	   if (myrank.eq.spec_cpu(i)) close(600+i)
	 enddo  
      endif
c
      if (box2) then 
        if (myrank.eq.0) close(63)
      endif
c      	
      return
      end
C
C-----------------------------------------------------------------------
C
