c     6/05 
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
C-----------------------------------------------------------------------
C     
C     (C) I. BALARAS
C     VON KARMAN INSTITUTE FOR FLUID DYNAMICS
C     CHAUSSEE DE WATERLOO, 72
C     B-1640 RHODE-ST-GENESE
C     
C     DATE OF CREATION :  20/10/1993                                
C     LAST UPDATE      :  05/10/2000

C     PARALLELIZATION

C     (C) A. SCOTTI 2002

C     (C) V. OVCHINNIKOV 2005


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
      real tvavg(nx,nz)!,tvavgg(nx,nzg),tvi
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
cccccccc variables for 2D planes output
      include 'headers/planes.h'

      integer cpuxy(nzplanes),kxy(nzplanes)
c      real u2dxz(nx,nz,2),v2dxz(nx,nz,2),w2dxz(nx,nz,2),p2dxz(nx,nz,2)
c      real u2dxy(nx,ny),v2dxy(nx,ny),w2dxy(nx,ny),p2dxy(nx,ny)
c      real u2dyz(ny,nz,2),v2dyz(ny,nz,2),w2dyz(ny,nz,2),p2dyz(ny,nz,2)
cccccccc
      real uinft(nz),uinftg(nzg),ubulk(nz),cf(nz),dst(nz),theta(nz),winst(nz),tauwg(nzg),utaug(nzg),ubulkg(nzg),dstg(nzg)
      real thetag(nzg),cfg(nzg),kpg(nzg),winstg(nzg)
      real uavg(nx,nz),vavg(nx,nz),wavg(nx,nz),uuavg(nx,nz),vvavg(nx,nz),wwavg(nx,nz),uwavg(nx,nz),uuuavg(nx,nz),
     &     wwuavg(nx,nz) 
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
c      real uavgg(nx,nzg),vavgg(nx,nzg),wavgg(nx,0:nzg),uuavgg(nx,nzg),vvavgg(nx,nzg),
c     &           wwavgg(nx,nzg),uwavgg(nx,nzg),pavgg(nx,nzg),ppavgg(nx,nzg),uuuavgg(nx,nzg),wwuavg(nx,nzg)
c      real omxavgg(nx,ngz),omyavgg(nx,nzg),omzavgg(nx,nzg)
c     start of inflow        
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
      REAL    tlevel,t_prev,t_old
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
c
      integer type !double or single (for MPI)
c     
c     added for boundary layer simulation VO
      integer iistart,imid,imidf,ifin,nload,n,nn,status(MPI_STATUS_SIZE),counter,cpusave
      integer np,ksave
      real xloc, uscale, puscale, pvscale, pwscale, plevel, ux,vx,scx
      integer pind
      real dummy1,dummy2,dummy3
c--------------statistics variables-------------------------------------
      real uactg
c-----------------------------------------------------------------------
c     OSE/ROGALLO FUNCTIONALITY FLIP CASE
      include 'headers/ose.h'
c      arrays for isotropic inflow     
      real kx(nmodes),ky(nmodes),kz(nmodes)
      complex uh(nmodes),vh(nmodes),wh(nmodes)
      real uose(ny,nx),vose(ny,nx),wose(ny,nx)
      real cutoffu(nx),cutoffw(nx),shape
      integer jcounter, count(nmodes), is(nmodes,kose) 
c     other vars for isotropic inflow
      integer seed,kkmax ! kose defined in ose.h
c      
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
      include 'headers/spec.h'
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
      integer m, m2, mm
      real oldtime, znewtime ! variables for keeping track of time per timestep
      real amp,env !for perturbing laminar wake-bl interaction

c        unknown variables
      integer itstamp
      real dpdy ! it seems that someone wanted to introduce a spanwise or tranverse pressure gradient
c
c      temporary debugging variables
      real divmax_,divmin_,divavg_
      integer idmax, jdmax, kdmax, idmin, jdmin, kdmin
      real rsu, rsv, rsw ! currently exist for backward compatibility with "correct"        
      real umx,vmx,wmx
c
cccccccccVARIABLES FOR REYNOLDS STRESS BUDGETS -- SEE RSBUD1.F routine
c        NOTE:  TO SAVE MEMORY rsblk_1 overwrites uab2 (1st 4*7 nx*nz planes) This works if ny>=28!!!
c                if LES is on rsblk_3 overwrites vab2  (1st 4*4 nx*nz planes)
      real rsblk_1(nx,nz,ny)
      real rsblk_2(nx,0:nz+1,33)
      real rsblkg_1(nx,nzg,7) 
#ifdef LES
      real rsblk_3 (nx,nz,ny)     !temporary LES block
      real rsblk_4 (nx,0:nz+1,25) !permanent LES block
      real rsblkg_3(nx,nzg,4)     !temporary global LES block
#endif
c
#ifdef immbnd ! conditional compilation
c     IMMERSED BOUNDARY DATA
      real sb(0:nnmax),yb(0:nnmax),xb(0:nnmax)
      real coefx(nnmax,3),coefy(nnmax,3)
      real uim(nnmax),vim(nnmax),wim(nnmax),tvim(nnmax),Ilmim(nnmax)
      integer nb
      integer npres!, ipres(nimb), kpres(nimb) ! for surface pressure computation
      integer npresl, ipresl(nimb), kpresl(nimb), ! for surface pressure computation (local)
     &        npres_array(ncpu), npres_disp(ncpu)
      real     zpres(nimb), xpres(nimb), pres(nimb)              
      real     zpresl(nimb), xpresl(nimb), presl(nimb)              
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc VARIABLES FOR THE PARALLELIZATION OF IMMERSED BOUNDARY cccccccccccccccc
c U
      integer uniml(0:ncpu-1),uiiml(0:ncpu-1,nnmax),ujiml(0:ncpu-1,nnmax),untop(0:ncpu-1),utk(0:ncpu-1,nimb),
     &        ustk1(0:ncpu-1,nnmax),ustk2(0:ncpu-1,nnmax),unbot(0:ncpu-1),ubk(0:ncpu-1,nimb),usti(0:ncpu-1,nnmax)
      real    utib_ind(0:ncpu-1,nimb),utib_c1(0:ncpu-1,nimb),utib_c2(0:ncpu-1,nimb),ubib_ind(0:ncpu-1,nimb),
     &        ubib_c1(0:ncpu-1,nimb),ubib_c2(0:ncpu-1,nimb),uc4intl(0:ncpu-1,nnmax,5)
      integer uibod((nz-2)*(nx-2)),ujbod((nz-2)*(nx-2))
      integer unbod!(0:ncpu-1)
c V
      integer vniml(0:ncpu-1),viiml(0:ncpu-1,nnmax),vjiml(0:ncpu-1,nnmax),vntop(0:ncpu-1),vtk(0:ncpu-1,nimb),
     &        vstk1(0:ncpu-1,nnmax),vstk2(0:ncpu-1,nnmax),vnbot(0:ncpu-1),vbk(0:ncpu-1,nimb),vsti(0:ncpu-1,nnmax)
      real    vtib_ind(0:ncpu-1,nimb),vtib_c1(0:ncpu-1,nimb),vtib_c2(0:ncpu-1,nimb),vbib_ind(0:ncpu-1,nimb),
     &        vbib_c1(0:ncpu-1,nimb),vbib_c2(0:ncpu-1,nimb),vc4intl(0:ncpu-1,nnmax,5)
      integer vibod((nz-2)*(nx-2)),vjbod((nz-2)*(nx-2))
      integer vnbod!(0:ncpu-1)
c W
      integer wniml(0:ncpu-1),wiiml(0:ncpu-1,nnmax),wjiml(0:ncpu-1,nnmax),wntop(0:ncpu-1),wtk(0:ncpu-1,nimb),
     &        wstk1(0:ncpu-1,nnmax),wstk2(0:ncpu-1,nnmax),wnbot(0:ncpu-1),wbk(0:ncpu-1,nimb),wsti(0:ncpu-1,nnmax)
      real    wtib_ind(0:ncpu-1,nimb),wtib_c1(0:ncpu-1,nimb),wtib_c2(0:ncpu-1,nimb),wbib_ind(0:ncpu-1,nimb),
     &        wbib_c1(0:ncpu-1,nimb),wbib_c2(0:ncpu-1,nimb),wc4intl(0:ncpu-1,nnmax,5)
      integer wibod((nz-2)*(nx-2)),wjbod((nz-2)*(nx-2))
      integer wnbod!(0:ncpu-1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccc  VARIABLES FOR PRESSURE RECOVERY
      integer ptrashn,ptrashi(0:ncpu-1,nnmax), ptrashj(0:ncpu-1,nnmax)
      integer me
      logical found
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif
cccccccccc variables for two-pt correlations ccccccccccc           
      include 'headers/two_pt.h'

      integer nz_act,nx_act
      real xn(nxl), zn((nzg-2+zfreq-1)/zfreq)        
      integer ks(0:ncpu-1),ke(0:ncpu-1),npt(0:ncpu-1) ! these arrays are defined on all processors
      real um2pt(nxl,nzl), vm2pt(nxl,nzl), wm2pt(nxl,nzl)
      real u2pt(nxl,ny,nzl),v2pt(nxl,ny,nzl),w2pt(nxl,ny,nzl)
      real u2tmp(nxl,ny,nzl),v2tmp(nxl,ny,nzl),w2tmp(nxl,ny,nzl) ! temporary arrays
c
      character*4 rsnames(0:3)
c
      equivalence (uab2,rsblk_1)
#ifdef LES
      equivalence (vab2,rsblk_3)
#endif
c      equivalence (rsblkg_1(1,1,1),uavgg(1,1)),(rsblkg_1(1,1,2),vavgg(1,1)),
c     &            (rsblkg_1(1,1,3),wavgg(1,1)),(rsblkg_1(1,1,4),uuavgg(1,1)),
c     &            (rsblkg_1(1,1,5),vvavgg(1,1)),(rsblkg_1(1,1,6),wwavgg(1,1)),
c     &            (rsblkg_1(1,1,7),pavgg(1,1))
c#ifdef LES
c      equivalence (rsblkg_3(1,1,1),uwavgg(1,1)),(rsblkg_3(1,1,2),tvavgg(1,1))     
c#endif     
c      
      equivalence (U_,Ur),(d2U,d2Ur),(phi,phir),(phip,phipr)
c     
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
      equivalence (u2tmp,uab2),(v2tmp,vab2),(w2tmp,wab2)

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
      InputFile = 'eddy.input'
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
c    random number seed
      seed=2
c
c     initiaize statistics variables
      call init_stat(samples,statnum,uinft,uav,vav,wav,uuav,vvav,wwav,
     &  uavg,vavg,wavg,pavg,uuavg,vvavg,wwavg,ppavg,uwavg,uuuavg,wwuavg,
#ifdef LES
     & tvavg,
#endif
     & rsblk_2, 
#ifdef LES     
     & rsblk_4,
#endif
     & rsnames)
c
c     set up types      		
      call MPI_TYPE_CONTIGUOUS(nz-2,type,real1d,ierr)
      call MPI_TYPE_CONTIGUOUS((nz-2)*nx,type,stat2d,ierr)
      call MPI_TYPE_COMMIT(real1d,ierr)
      call MPI_TYPE_COMMIT(stat2d,ierr)
      call MPI_TYPE_VECTOR(nz,nx,ny*nx,type,real2d,ierr)
      call MPI_TYPE_VECTOR(nz*ny,1,nx,type,cross2d,ierr)
      call MPI_TYPE_COMMIT(cross2d,ierr)
      
      call MPI_TYPE_COMMIT(real2d,ierr)
c     
cccccccccccccc READ INITIAL DATA ccccccccccccccccccc      
      call readini(uo,vo,wo,p,uab1,vab1,wab1,tv,
#ifdef LES
     &             ilm,imm,
#endif LES
     &             xu,zw,xcf,ycf,tlevel)
ccccccccccccccccccccccccccccccccccccccccccccccccccccc       
c
      imid=nx/2-1  !default values
      imidf=nx/2+1 !default values
c      compute 2-pt correlation grids (center)
      if (two_pt_corr) call two_pt_corr_setup(nx_act,nz_act,xn,zn,ks,ke,npt,xcf,zcg)
c
      if (jdose.and.rogal) then 
       write(0,*) 'error: jd-ose and rogallo noise are both on!'
       return
      endif
      if (box2.and.(jdose.or.rogal)) then 
       write(0,*) 'error: box2 and jd-ose (or rogallo) noise are both on!'
       return
      endif
c       
      if (rogal) call rogal_new_setup(jcounter,kx,ky,kz,uh,vh,wh,cutoffu,cutoffw,uscale,xu,xcf,yv)
c      if (rogal) call rogal_old_setup(jcounter,kx,ky,kz,count,is,uh,vh,wh,cutoffu,cutoffw,uscale,xu,xcf,yv)
      if (jdose) call jdose_new_setup(jcounter,kx,ky,kz,uh,vh,wh,cutoffu,cutoffw,uscale,eta,detacor,d2U,U_,Ur,xu,xcf,yv,
     &                                a_w,b_w,a_u,b_u,ps_w,ps_u)
      if (jdose.or.rogal) then 
         open(101,file='pscale.dat',form='formatted',status='old')
	   read(101,*) puscale, pvscale, pwscale
	   read(101,*) psample
	 close(101)  
c      target turbulence level
         plevel=u_rms ! u_rms is in the input file
	 pind=1
	 do while (zwg(pind).le.u_rms_loc)
	   pind=pind+1
	 enddo
      endif	
c
      if (box2) call box2_setup(a_u,a_w,b_u,b_w,ps_u,ps_w,tlevel,xu,xcf,index,type)
c
      if (save) call save_setup(cpusave,ksave,xu,yv,zwg,tlevel,t_prev,index)
      if (spec) call spec_setup(spec_cpu,spec_k,xu,zwg,index)
c
ccccccccccccc setup for output planes cccccccccccccc      
      if (planes) call planes_setup(cpuxy,kxy,zwg)
cccccccccccccccccccccccccccccccccccccccccccccccccccc

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
      IF (ftype .EQ. 2) dpdx =  0.0
      IF (ftype .EQ. 3) dpdx =  0.0
      IF (ftype .EQ. 4) dpdx =  0.0
c     i.e. assume running vertically if top and bottom are periodic
       
c------------------------------------------------------------------------
c     Check Divergence of Starting Field
c-----------------------------------------------------------------------
      CALL CALCDIV(UO,VO,WO,NX,NY,NZ)
c-----------------------------------------------------------------------
c     setup immersed boundary
c-----------------------------------------------------------------------
#ifdef immbnd ! conditional compilation
      if (ibswitch) then
c.... obstacle location on the grid
               CALL ObjectGrdCyl(sb,xb,yb,coefx,coefy,nb,nnmax)
c.... locate boundary points and interpolation coeff. for u
               CALL bpoints(xu,zc,zcg,sb,xb,yb,nb,coefx,coefy,
     &                   uniml,uiiml,ujiml,untop,utk,ustk1,ustk2,usti,unbot,ubk,utib_ind,utib_c1,utib_c2,ubib_ind,
     &                   ubib_c1,ubib_c2,
     &                   uc4intl,unbod,uibod,ujbod,'u.dat')
c.... locate boundary points and interpolation coeff. for v
               CALL bpoints(xc,zc,zcg,sb,xb,yb,nb,coefx,coefy,
     &                   vniml,viiml,vjiml,vntop,vtk,vstk1,vstk2,vsti,vnbot,vbk,vtib_ind,vtib_c1,vtib_c2,vbib_ind,
     &                   vbib_c1,vbib_c2,
     &                   vc4intl,vnbod,vibod,vjbod,'p.dat')
c.... locate boundary points and interpolation coeff. for w
               CALL bpoints(xc,zw,zwg,sb,xb,yb,nb,coefx,coefy,
     &                   wniml,wiiml,wjiml,wntop,wtk,wstk1,wstk2,wsti,wnbot,wbk,wtib_ind,wtib_c1,wtib_c2,wbib_ind,
     &                   wbib_c1,wbib_c2,
     &                   wc4intl,wnbod,wibod,wjbod,'w.dat')
cccccccccccc experimental cccccccccccc surface pressure computation cccccccccccccccccc

c      call pressure_setup(nb,sb,xb,yb,coefx,coefy,xc,zc,zcg,npres,ipres,kpres,npresl,ipresl,kpresl,myrank)
      call pressure_setup(nb,sb,xb,yb,coefx,coefy,xc,zc,zcg,npresl,ipresl,kpresl,myrank)
       
ccccc debug; write out points
c     form local arrays
      do i=1,npresl
       xpresl(i)=xc(ipresl(i))
       zpresl(i)=zc(kpresl(i))
       presl(i)=p(ipresl(i), ny/2+1, kpresl(i)) ! pressure
      enddo
c     collect data              
      call mpi_reduce(npresl, npres, 1, MPI_INTEGER, MPI_SUM, ncpu-1, MPI_COMM_HYDRO, ierr)
      call mpi_gather(npresl, 1, MPI_INTEGER, npres_array, 1, MPI_INTEGER, ncpu-1, MPI_COMM_HYDRO, ierr)
c
      npres_disp(1)=0
      do i=2,ncpu
       npres_disp(i)=npres_array(i-1)
      enddo 
c      write(0,*) npresl, npres, npres_array, npres_disp
      call mpi_gatherv(zpresl, npresl, type, zpres, npres_array, npres_disp, type, ncpu-1, MPI_COMM_HYDRO, ierr)
      call mpi_gatherv(xpresl, npresl, type, xpres, npres_array, npres_disp, type, ncpu-1, MPI_COMM_HYDRO, ierr)
      call mpi_gatherv(presl, npresl, type, pres, npres_array, npres_disp, type, ncpu-1, MPI_COMM_HYDRO, ierr)
c
      if (myrank.eq.ncpu-1) then 
        open(99,file='pres.dat',form='unformatted', status='unknown')
	write(99) npres, (xpres(i), zpres(i), pres(i), i=1,npres)
        close(99)
      endif		

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       for pressure recovery -- experimental

c      me=myrank
c      do i=1,vnbod
c        write(300+myrank,'(3 I)') i,vibod(i),vjbod(i)
c      enddo
c      close (300+myrank)
ccccccccccccccccccccccccccccc EXPERIMENTAL AND NON_ESSENTIAL ccccccccccccccccccccccccccc
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
c     
c      do i=1, ptrashn
c        write(me+250,'(3 I)') i,ptrashi(me,i),ptrashj(me,i)  
c      enddo
c        close (me+250)
c
c      do i=1, uniml(me)
c        write(me+300,'(7 I)') i,uiiml(me,i),ujiml(me,i),viiml(me,i),vjiml(me,i),wiiml(me,i),wjiml(me,i)  
c      enddo
c        close (me+300)
cccccccccccccccccccccc      for pressure recovery experimental	
     
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
      
 9875  continue                  
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
c      if (ifield.ne.2) then 
c       me=myrank
c       do j=jy1,jy2+1
c        do i=1,vnbod
c         p (vibod(i),j,vjbod(i))=1.
c         vo(vibod(i),j,vjbod(i))=0.
c        enddo
c        do i=1,unbod
c         uo(uibod(i),j,ujbod(i))=0.
c        enddo
c        do i=1,wnbod
c         wo(wibod(i),j,wjbod(i))=0.
c        enddo
c       enddo
c      endif 
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
      ICYCLE =ITINI
      ipart  =0
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
      CALL turvis(uo,vo,wo,dxdydz,xu,yv,zw,urun,
     &     sxx,syy,szz,sxy,sxz,syz,g,
     &     uoc,voc,woc,upr,vpr,wpr,gbar,
     &     lm,mm,Ilm,Imm,xcf,ycf,zcf,
     &     wrk1,wrk2,wrk3,wrk4,
     &     tv,dtm1,icycle-itini,nx,ny,nz       )
c
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
      write(666,*) icycle, znewtime-oldtime
      oldtime=znewtime
      close(666)
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
c     for phase reshuffling -- new routine
      if (jdose.or.rogal) then 
       if (myrank.eq.0) then 
        dummy1=mod(uinfty*tlevel+zw(1),zlen)
        dummy2=mod(uinfty*tlevel+dtm1+zw(1),zlen)
c	write(0,'(5 G20.10)') dummy1,dummy2,tlevel,dtm1,zlen
	if (dummy2.lt.dummy1) call reshuffle(jcounter,uh,vh,wh,ky,kz,seed)
c	if (dummy2.lt.dummy1) call reshuffle_old(jcounter,uh,vh,wh,count,is,seed)
       endif !myrank

       call MPI_BCAST(uh(1),2*nmodes,type,0,MPI_COMM_HYDRO,ierr)		     	
       call MPI_BCAST(vh(1),2*nmodes,type,0,MPI_COMM_HYDRO,ierr)		     	
       call MPI_BCAST(wh(1),2*nmodes,type,0,MPI_COMM_HYDRO,ierr)		     	
       call MPI_BCAST(ky(1),nmodes,type,0,MPI_COMM_HYDRO,ierr)		     	
       call MPI_BCAST(kz(1),nmodes,type,0,MPI_COMM_HYDRO,ierr)		     	
      endif !jdose/rogal
c
      if (box2) call box2dt(dtm1,icycle,index,type)
c      
      TLEVEL=TLEVEL+DTM1
C-----------------------------------------------------------------------
C     predictor-corrector step
C-----------------------------------------------------------------------
      ITSTAMP = ITMAX
      ICOUNT1 = MOD(ICYCLE,100000)
      
      IF ((ICOUNT1.EQ.0).OR.((IFIELD.NE.2).AND.(ICYCLE-ITINI.EQ.1))) THEN
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
      if(impfla.eq.1.and.icycle-itini.eq.1) then
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
c     don't think this is needed anymore since new IB does _NOT_ use ghost points
c      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbc'
c      CALL REFRESHBC(USTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
c      CALL REFRESHBC(VSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
c      CALL REFRESHBC(WSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
c      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'    
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
c      me=myrank
c      do j=jy1,jy2+1
c       do i=1,vnbod
c        p (vibod(i),j,vjbod(i))=0.
c        vstar(vibod(i),j,vjbod(i))=0.
c       enddo
c       do i=1,unbod
c        ustar(uibod(i),j,ujbod(i))=0.
c       enddo
c       do i=1,wnbod
c        wstar(wibod(i),j,wjbod(i))=0.
c       enddo
c      enddo
c
#endif  
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbc after predictor (or IB)'
      CALL REFRESHBC(USTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(VSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(WSTAR,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'    
ccccccIMMERSED BOUNDARY
c     need to call boundary because refreshbc changes wstar at 1st point!!!
      CALL Boundary(USTAR,VSTAR,WSTAR,UO,VO,WO,
     &     u2d,v2d,w2d,tlevel,t_old,ksub,nt,
     &     urun,wrun,
     &     XU,YV,ZW,nx,ny,nz,Ltstep,uinf,duinf,dtm1,x_ref)
c
      if (jdose) then
       call jdorrsom(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,eta,detacen,detacor,U_,d2U,
     &                 a_w,b_w,a_u,b_u,ps_w,ps_u)
       if (myrank.eq.0) then 
c     add the perturbation
        do i=2,nx-1 
         do j=2,ny-1
          ustar(i,j,1)=ustar(i,j,1)+wose(j,i)*uscale*puscale    
          vstar(i,j,1)=vstar(i,j,1)+vose(j,i)*uscale*pvscale    
          wstar(i,j,1)=wstar(i,j,1)+uose(j,i)*uscale*pwscale
         enddo
        enddo    
       endif ! myrank
      endif !jdose
c          
      if (rogal) then 
c     generate inflow with fourier wall-normal basis cccccc
       call rogallo(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,type)
c       call rogallo_old(uose,vose,wose,jcounter,uh,vh,wh,kx,ky,kz,xu,xcf,yv,ycf,zwg,tlevel,type)
       if (myrank.eq.0) then
c     add the perturbation
        do i=2,nx-1 
         do j=2,ny-1
          ustar(i,j,1)=ustar(i,j,1)+wose(j,i)*uinfty*uscale*cutoffu(i)*puscale     
          vstar(i,j,1)=vstar(i,j,1)+vose(j,i)*uinfty*uscale*cutoffw(i)*pvscale   
          wstar(i,j,1)=wstar(i,j,1)+uose(j,i)*uinfty*uscale*cutoffw(i)*pwscale
	 enddo
        enddo    
c	   write(0,*) 'max(u,v,w):',umx,vmx,wmx
       endif ! myrank
      endif !rogal
c    
      if (box2) call box2read(ustar,vstar,wstar,a_u,a_w,b_u,b_w,ps_u,ps_w)
c      
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc          
      cpupredict=mpi_wtime()-cpub
c     
      call MasBlas(USTAR,VSTAR,WSTAR,xu,yv,zw,icycle-itini,nx,ny,nz)
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

         CALL presdirect(dp,dp,p,xu,yv,zwg,au,cwg,icycle-itini,
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
     &     ICYCLE-itini,RSU,RSV,RSW)
ccccc pressure recovery
c      me=myrank
c      do i=1,ptrashn
c       do j=jy1,jy2+1
c        p(ptrashi(me,i),j,ptrashj(me,i))=0.
c       enddo
c      enddo

#ifdef immbnd
      me=myrank
      do j=jy1,jy2+1
       do i=1,vnbod
        p (vibod(i),j,vjbod(i))=4.
c        vo(vibod(i),j,vjbod(i))=0.
       enddo
c       do i=1,unbod
c        uo(uibod(i),j,ujbod(i))=0.
c       enddo
c       do i=1,wnbod
c        wo(wibod(i),j,wjbod(i))=0.
c       enddo
      enddo
#endif
C-----------------------------------------------------------------------
C     boundary conditions to the updated values
C-----------------------------------------------------------------------
C     
      Ltstep = .TRUE. 
      CALL Boundary(UO,VO,WO,USTAR,VSTAR,WSTAR,
     &     u2d,v2d,w2d,tlevel,t_old,ksub,nt,
     &     urun,wrun,
     &     XU,YV,ZW,nx,ny,nz,Ltstep,uinf,duinf,dtm1,x_ref)
c     
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': refreshbs'
      CALL REFRESHBC(UO,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(VO,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
      CALL REFRESHBC(WO,NX*NY,NZ-2,MYSIZE,MYRANK,MPI_COMM_HYDRO)
c
      if (box2) call box2bc(ustar,vstar,wstar,uo,vo,wo)
c     .....divergence
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
      CALL Storen(UAB1,VAB1,WAB1,UAB2,VAB2,WAB2,tmeth,nx,ny,nz)
c     
      cpucorr=mpi_wtime()-cpub
c     ccccccccccccccccccccccccccccccccccc
      if (save) call save_slice(uo,vo,wo,cpusave,ksave,xu,yv,tlevel,t_prev,icycle,index)
c
c-----------------------------------------------------------------------
C     eddy viscosity  and forcing
C-----------------------------------------------------------------------
c.... compute runing average of velocities for inflow bc
c     CALL RunAve(uo,vo,wo,urun,vrun,wrun,xu,yv,zw,
c     &            tlevel,dtm1,icycle,nx,ny,nz       )

C.....turbulent viscosity
      cpub=mpi_wtime()
#ifdef LES
c      write(0,*) myrank,' calling turvis'
      CALL turvis(uo,vo,wo,dxdydz,xu,yv,zw,urun,
     &     sxx,syy,szz,sxy,sxz,syz,g,
     &     uoc,voc,woc,upr,vpr,wpr,gbar,
     &     lm,mm,Ilm,Imm,xcf,ycf,zcf,
     &     wrk1,wrk2,wrk3,wrk4,
     &     tv,dtm1,icycle-itini,nx,ny,nz       )
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
          call spec_write(uo,vo,wo,p,spec_cpu,spec_k,tlevel,j)
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
 555  continue
      if (.not.(CalcStat.and.planes)) then
cccccccccccccccccccccccccccc OUTPUT one xz PLANE cccccccccccccccccccccccccccccccccccccccccccc
         if(mod(ICYCLE-itini-1,infnum).eq.0) then
 1234            j=ny/2
                 write(66+mod(icycle,5),*) myrank,': matlab output'
                 call plane_xz(uo,vo,wo,p,xcf,ycf,zcg,STR5,j,real2d)
                 write(66+mod(icycle,5),*) myrank,': done'
         endif !mod
      endif ! .not.CalcStat
C     
C-----------------------------------------------------------------------
C     create a series of restart files
C-----------------------------------------------------------------------
      
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
         elseif (Isgsmodel.eq.30) then  
           call ioscalar('savt.res',Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar('satv.res',tv,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
	 
         if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'
c     
c     write latest inflow phase information(note: this file overwrites the original)
         if (jdose.or.rogal) then 
          call waves_write(jcounter,kx,ky,kz,uh,vh,wh)
c          call waves_write_old(jcounter,kx,ky,kz,count,is,uh,vh,wh)
c	  
          if (myrank.eq.0) then        
            open(101,file='pscale.dat',form='formatted',status='old')
              write(101,'(3 G25.15)') puscale, pvscale, pwscale
	      write(101,'(I6)') psample
	    close(101)  
	  endif     
 999  continue
         endif
      endif
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
  456 continue
      if (CalcStat.eq..true..and.mod(icycle,ISTATS).eq.1) then 
         samples =samples+1      
         if (samples.gt.statsample) then ! write new stats files
           statnum=statnum+1
           samples=1
	 endif  	  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc	 
c     compute uinfinity 
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': compute uinft,uact'
      call uinfinity(wo,xu,zwg,uinft,uinftg,uactg,samples,real1d)      
      if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,': done'	
c     compute free-stream statistics
      call stats_old(uo,vo,wo,zwg,uav,vav,wav,uuav,vvav,wwav,guav,gvav,gwav,guuav,gvvav,gwwav,samples,statnum,imid,imidf,
     &               index,real1d)
cccccccccc DYNAMIC NOISE ADJUSTMENT ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc      
        if (0) then!(jdose.or.rogal) then 
          if (myrank.eq.mysize-1) then 
c     recomupute scaling parameter
            uosm=sqrt(guuav(pind)-guav(pind)**2)/uinfty
            vosm=sqrt(gvvav(pind)-gvav(pind)**2)/uinfty
            wosm=sqrt(gwwav(pind)-gwav(pind)**2)/uinfty
	       	    
      	    puscale=(puscale*psample+puscale*plevel/(uosm))
      	    pvscale=(pvscale*psample+pvscale*plevel/(vosm))
      	    pwscale=(pwscale*psample+pwscale*plevel/(wosm))
      
            psample=psample+1
            puscale=puscale/psample
            pvscale=pvscale/psample
            pwscale=pwscale/psample
            write(0,'("inlet noise rms: u,v,w & scales ",6 G13.5,I)') uosm,vosm,wosm,puscale,pvscale,pwscale,psample
          endif
c
          call MPI_BCAST(puscale,1,type,mysize-1,MPI_COMM_HYDRO,ierr)
          call MPI_BCAST(pvscale,1,type,mysize-1,MPI_COMM_HYDRO,ierr)
          call MPI_BCAST(pwscale,1,type,mysize-1,MPI_COMM_HYDRO,ierr)
          call MPI_BCAST(psample,1,MPI_INTEGER,mysize-1,MPI_COMM_HYDRO,ierr)
        endif   

 222  continue
c     compute turbulent boundary layer statistics
      call stats_new(uo,vo,wo,p,tv,xcf,xu,zcg,zwg,uavg,vavg,wavg,pavg,uuavg,vvavg,wwavg,uwavg,ppavg,
     &                     uuuavg,wwuavg,utaug,tauwg,
     &                     uinft,uinftg,dstg,thetag,uactg,winstg,cfg,
     &                     kpg,ubulkg,samples,statnum,real1d,stat2d,index,icycle)
c
        if (myrank.eq.mysize-1) then 
            open(22,file='shape'//index(statnum)//'.gnu',status='unknown',form='formatted')
            open(27,file='cf'//index2(icycle)//'.gnu',status='unknown',form='formatted')
            do k=kz1,nzg-1 
               write(22,'(6 G20.12)') (zwg(k)-x_ref)*uactg/ru1,cfg(k),dstg(k),thetag(k),kpg(k),zwg(k)-x_ref
               write(27,'(3 G20.12)') (zwg(k)-x_ref)*uactg/ru1,winstg(k),zwg(k)-x_ref
            enddo !k
            close(22)
      	    close(27)
        endif	!myrank
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccc INSTANTANEOUS PLANES cccccccccccccccccccccccccccc
         if (planes) call write_planes(uo,vo,wo,p,
#ifdef LES
     &                        tv,
#endif
     &  xu,xcf,yv,ycf,zcf,zcg,cpuxy,kxy,type,real2d,cross2d,index,index2,icycle)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccREYNOLDS STRESS BUDGETS -- SEE 'RSBUD1.F cccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (rstat) call rstress(rsblk_1,rsblkg_1,rsblk_2,
#ifdef LES      
     &                           rsblk_3,rsblk_4,rsblkg_3,
#endif     
     &                           uo,vo,wo,p,tv,utaug,tauwg,samples,statnum,rsnames,xcf,zcg,stat2d,index)
ccccccccccccccccccccccccccccccccREYNOLDS STRESS BUDGETS -- SEE 'RSBUD1.F cccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccc two-point spanwise correlations cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         if (two_pt_corr) then 
	  call span_corr(uo,vo,wo,u2pt,v2pt,w2pt,um2pt,vm2pt,wm2pt,ks(myrank),ke(myrank),samples)
c      normalize & write output file       
          call write_corr(u2pt,v2pt,w2pt,um2pt,vm2pt,wm2pt,nx_act,nz_act,xn,ycf,zn,npt,type)
          if (myrank.eq.0) write(66+mod(icycle,5),*) myrank,' two_pt_corr: done'
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
         elseif (Isgsmodel.eq.30) then  
           call ioscalar('savt.res',Ilm,NX,NY,NZ,NZG,idir,TLEVEL)
           call ioscalar('satv.res',tv,NX,NY,NZ,NZG,idir,TLEVEL)
         endif
#endif
c     write latest inflow phase information(note: this file overwrites the original)
         if (jdose.or.rogal) then 
          call waves_write(jcounter,kx,ky,kz,uh,vh,wh)
c          call waves_write_old(jcounter,kx,ky,kz,count,is,uh,vh,wh)
c
          if (myrank.eq.0) then        
            open(101,file='pscale.dat',form='formatted',status='old')
              write(101,'(3 G25.15)') puscale, pvscale, pwscale
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
      if (save) call save_close(cpusave)
      if (spec) call spec_close(spec_cpu)
      if (box2) call box2close()
c
      return
      end
C
C-----------------------------------------------------------------------
C