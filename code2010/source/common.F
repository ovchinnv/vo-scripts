! 
!-----INCLUDE-FILE-- common.h -----------------I. Balaras-20/05/1998----
!
!-----------------------------------------------------------------------
!                                                                      
!         I. BALARAS                          
!         LAST UPDATE      :  29/11/1998
!
!-----------------------------------------------------------------------
!                                                             parameters
!-----------------------------------------------------------------------
!
!     
      CHARACTER*80  str1,str2,str3,str4,str5,str6,str7,str8,str9,str10

      LOGICAL       CalcStat,save,spec,rstat,planes,two_pt_corr, pres_out
      LOGICAL       Lcfl
      INTEGER       ftype, Ipresmethod, Isgsmodel

      INTEGER       ix1,ix2,jy1,jy2,kz1,kz2
      INTEGER       ibu,ieu,jbv,jev,kbw,kew
      REAL          dx,dy,dz

      INTEGER       itype(6),
     &     ibegin(6),iend(6),jbegin(6),jend(6),
     &     kbegin(6),kend(6)

      REAL          app(mmx),ap(mmx),au(mmx),av(mmx),aw(mmx),
     &     bpp(mmy),bp(mmy),bu(mmy),bv(mmy),bw(mmy),
     &     cpp(mmz),cp(mmz),cu(mmz),cv(mmz),cw(mmz),cwg(mmz)

      REAL          Xpos(10)
      INTEGER       ifield,itmax,itini,infnum,ievery,istats
      INTEGER       ibsp,iesp
      REAL          ru1,uinfty,css,angle,vspan
      REAL          tstep,tsub
      INTEGER       kst,isub,ksub_old
      INTEGER       ibsub(1000),iesub(1000)
      real          qinout,uconv
      REAL          whx(1:mmx,-1:1,1:2),why(1:mmy,-1:1,1:2),
     &              whz(1:mmz,-1:1,1:2)
      INTEGER       iper,jper,kper
      INTEGER       impfla
      REAL          ibc(3,6)
!
      real          x_ref, u_rms, u_rms_loc
      character*80  NOISEFILE
      logical       jdose,rogal,box2
!
      logical       read_ilm_imm,read_vt_tv
      character*80  ilm_file, imm_file, vt_file, tv_file
!
      integer       statsample
!
      LOGICAL ibswitch
      REAL    ibheight
!                                                                       
!-----------------------------------------------------------------------
!                                                                commons
!-----------------------------------------------------------------------
!
      COMMON /IND/    ix1,ix2,jy1,jy2,kz1,kz2
      COMMON /CMP/    dx,dy,dz
      COMMON /MOM/    ibu,ieu,jbv,jev,kbw,kew
      COMMON /BND/    itype,ibegin,iend,jbegin,jend,kbegin,kend
      COMMON /MTCOF/  app,ap,au,av,aw,bpp,bp,bu,bv,bw,
     $     cpp,cp,cu,cv,cw,cwg
      COMMON /DAT1/   ifield,itmax,itini,infnum,ievery,istats,
     $     ftype,ibsp,iesp
      COMMON /DAT2/   ru1,uinfty,css,angle,vspan,xpos
      COMMON /SCHEME/ Isgsmodel,Ipresmethod,read_ilm_imm, ilm_file, imm_file,
     &      read_vt_tv, vt_file, tv_file
      COMMON /IOFILE1/ str1,str2,str3,str4,str5,
     &                 str6,str7,str8,str9,str10
      COMMON /timestep1/ Lcfl
      COMMON /timestep2/ tstep
      COMMON /boundary1/ qinout,uconv
      COMMON /boundary2/ ibc
      COMMON /weight/ whx,why,whz
      COMMON /multigrd1/ iper,jper,kper
      COMMON /implicit/ impfla
      COMMON /immerse1/ ibswitch
      COMMON /immerse2/ ibheight
      COMMON /parallel/ mysize,myrank,MPI_COMM_HYDRO
      COMMON /BNDLYR/ box2,jdose, rogal, x_ref, u_rms, u_rms_loc, NOISEFILE
      COMMON /RUNSTAT/ CalcStat, statsample, save,spec, rstat, planes, two_pt_corr, pres_out


      
