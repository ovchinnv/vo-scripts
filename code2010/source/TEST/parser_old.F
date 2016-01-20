! original 2006 code
!
!      
      
      

      CHARACTER*1   TEXT(72)
      int       i,type,ier,ierr
      
!
!-----------------------------------------------------------------------
!                                                    NAMELIST STATEMENTS
!-----------------------------------------------------------------------
!
      NAMELIST /STRPRM/ STR1,STR2,STR3,STR4,STR5,
     &                  STR6,STR7,STR8,STR9,STR10
      NAMELIST /BNDPRM/ ITYPE
      NAMELIST /SMEPRM/ Isgsmodel, Ipresmethod, read_ilm_imm, read_vt_tv,
     &                  ilm_file,imm_file, vt_file, tv_file
      NAMELIST /BGNPRM/ IFIELD
      NAMELIST /STATYESNO/ CalcStat, statsample, spec, save, rstat, planes, two_pt_corr, pres_out
      NAMELIST /IMP_EXP/ impfla
      NAMELIST /ITRPRM/ ITMAX,ITINI,INFNUM,IEVERY,ISTATS
      NAMELIST /FLOPRM/ Ftype,ru1,css,angle,uinfty,ibsp,iesp,vspan
      NAMELIST /BNDLYR/ box2, jdose, rogal, x_ref, u_rms, u_rms_loc, NOISEFILE
      NAMELIST /TIMEST/ Lcfl,tstep
      
#ifdef SINGLE
      type=MPI_REAL
#endif
#ifdef DOUBLE
      type=MPI_DOUBLE_PRECISION
#endif

!
!-----------------------------------------------------------------------
!                                                            READ UNIT 1
!-----------------------------------------------------------------------
!
!      FNAME = InputFile
!      OPEN (1,STATUS='OLD',FILE=FNAME)
      IF(MYRANK.EQ.0) THEN
         OPEN (1,STATUS='OLD',FILE=InputFile)
         OPEN (2,STATUS='SCRATCH') 
! 
!-----------------------------------------------------------------------
!                            COPY NAMELIST-DIRECTED INPUT DATA TO UNIT 2
!-----------------------------------------------------------------------
!
         REWIND 1 
         REWIND 2 
 10      READ (1,1000,END=20) (TEXT(I),I=1,72)
 1000    FORMAT (72(A1))
!     WRITE(6,2000) (TEXT(I),I=1,72)
!     6000 FORMAT (1X,72(A1))
         IF (TEXT(1) .EQ. '*') GOTO 10
         WRITE(2,2000) (TEXT(I),I=1,72)
!         WRITE(6,2000) (TEXT(I),I=1,72)
!     write(*.*) (TEXT(i),i=1,72)
 2000    FORMAT (1X,72(A1))
         GOTO 10
 20      REWIND 2
      ENDIF

! 
!-----------------------------------------------------------------------
!                                    SET DEFAULT-VALUES & READ NAMELISTS
!-----------------------------------------------------------------------
!
!....."STRPRM" 
! 
      STR1 = 'res.fin.vel'
      STR2 = 'cfd.grd'
      STR3 = 'res.fin.ab'
      STR4 = 'res.fin.p'
      STR5 = 'field.mat'
      STR6 = 'test.sta'
      STR7 = 'test.vis'
      STR8 = 'stilltobeused'
      STR9 = 'comp_hist.dat'
      STR10= 'slice.dat'
!      Xpos(1)  = 10.
!      Xpos(2)  = 20.
!      Xpos(3)  = 30.
!
      IF(MYRANK.EQ.0) THEN
         READ (2,STRPRM) 
      ENDIF

      CALL MPI_BCAST(STR1,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(STR2,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(STR3,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(STR4,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(STR5,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)

!
!....."BNDPRM" 
! 
      ITYPE(1)=500
      ITYPE(2)=500
      ITYPE(3)=500
      ITYPE(4)=500
      ITYPE(5)=110
      ITYPE(6)=110
!
      IF (MYRANK.EQ.0) THEN
         READ (2, BNDPRM)
      ENDIF

      CALL MPI_BCAST(ITYPE,6,MPI_int,0,MPI_COMM_HYDRO,IER)

       IF(ncpu.GT.1) THEN
         IF(me.EQ.0) THEN
            if (bound(top).ne.500) bound(top)=0 
         ELSEIF(me.EQ.ncpu-1) THEN
            if (bound(bottom).ne.500) bound(bottom)=0
         ELSE
            bottom(bottom)=0
            bottom(top)=0
         ENDIF
      ENDIF


!
!....."SMEPRM" 
!
      Isgsmodel = 1
      Ipresmethod = 1
      read_ilm_imm =.false.
      read_vt_tv = .false.
      ilm_file = 'ilm.res'
      imm_file = 'imm.res'
      vt_file = 'savt.res'
      tv_file = 'satv.res'
!

      IF(MYRANK.EQ.0) READ (2, SMEPRM)
      CALL MPI_BCAST(ISGSMODEL,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(IPRESMETHOD,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(read_ilm_imm,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(read_vt_tv,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(ilm_file,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(imm_file,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(vt_file,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(tv_file,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,IER)
! 
!....."BGNPRM" 
! 
      IFIELD=0
!
      IF(MYRANK.EQ.0) READ (2, BGNPRM) 
      CALL MPI_BCAST(IFIELD,1,MPI_int,0,MPI_COMM_HYDRO,IER) 

!.... "STATYESNO" VO 7.2.03

      CalcStat=.false.
      statsample=50
      save=.false.
      spec=.false.
      rstat=.false.
      planes=.false.
      two_pt_corr=.false.
      pres_out=.false.
 
      IF (MYRANK.EQ.0) read (2, STATYESNO)
      CALL MPI_BCAST(CalcStat,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(statsample,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(save,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(spec,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(rstat,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(planes,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(two_pt_corr,1,MPI_bool,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(pres_out,1,MPI_bool,0,MPI_COMM_HYDRO,IER)

!
!....."IMP_EXP"
      impfla=0.
      
      if (myrank.eq.0) read(2,IMP_EXP)
      call MPI_BCAST(impfla,1,MPI_int,0,MPI_COMM_HYDRO,ier)

!     
!....."ITRPRM" 
! 
      ITMAX=50
      ITINI=10
      INFNUM=10
      IEVERY=10
      ISTATS=20
!
      IF(MYRANK.EQ.0) READ (2,ITRPRM) 
      CALL MPI_BCAST(ITMAX,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(ITINI,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(INFNUM,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(IEVERY,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(ISTATS,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      
      
!
!....."FLOPRM"
!     
      Ftype=2
      ru1=1500
      css=0.11
      angle=0.0
      uinfty=1.
      ibsp=42
      iesp=94
      vspan=0.5
 
      IF(MYRANK.EQ.0) READ (2,FLOPRM) 
      ru1=1./ru1
      CALL MPI_BCAST(FTYPE,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(IBSP,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(IESP,1,MPI_int,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(RU1,1,TYPE,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(CSS,1,TYPE,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(ANGLE,1,TYPE,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(UINFTY,1,TYPE,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(VSPAN,1,TYPE,0,MPI_COMM_HYDRO,IER)

!....."BNDLYR"
      box2=.false.
      jdose=.false.
      rogal=.false.
      x_ref=-15.17
      u_rms=0.06
      u_rms_loc=0.0
      NOISEFILE='jdwaves6.dat'
      
      if (myrank.eq.0) read (2,BNDLYR)
      call MPI_BCAST(box2,1,MPI_bool,0,MPI_COMM_HYDRO,ier)
      call MPI_BCAST(jdose,1,MPI_bool,0,MPI_COMM_HYDRO,ier)
      call MPI_BCAST(rogal,1,MPI_bool,0,MPI_COMM_HYDRO,ier)
      call MPI_BCAST(x_ref,1,TYPE,0,MPI_COMM_HYDRO,ier)
      call MPI_BCAST(u_rms,1,TYPE,0,MPI_COMM_HYDRO,ier)
      call MPI_BCAST(u_rms_loc,1,TYPE,0,MPI_COMM_HYDRO,ier)
      call MPI_BCAST(x_ref,1,TYPE,0,MPI_COMM_HYDRO,ier)
      call MPI_BCAST(NOISEFILE,80,MPI_CHARACTER,0,MPI_COMM_HYDRO,ier)
!
!....."TIMEST"
!
      Lcfl=.true.
      tstep=0.01
      
      IF(MYRANK.EQ.0) READ (2,TIMEST)
      CALL MPI_BCAST(TSTEP,1,TYPE,0,MPI_COMM_HYDRO,IER)
      CALL MPI_BCAST(LCFL,1,MPI_bool,0,MPI_COMM_HYDRO,IER)

      IF(MYRANK.EQ.0) THEN
         IF((MOD(NY-2,MYSIZE).GT.0).OR.(MOD(NZG-2,MYSIZE).GT.0)) THEN
            PRINT*,'NY AND/OR NZ ARE NOT DIVISIBLE BY # PROC'
            PRINT*,'ABORTING NOW....'
            CALL MPI_FINALIZE(IERR)
            STOP
         ENDIF
      ENDIF

      IF(MYRANK.EQ.0) THEN
         CLOSE (1)
         CLOSE (2)
      ENDIF

!
      RETURN 
      END 
!
! 

