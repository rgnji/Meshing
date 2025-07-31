C
C======================================================================C
C                A TOOL FOR PARALLEL MULTI-BLOCK FDNSRFV               C
C                    BY HUAN-MIN SHANG, 12/02/1997                     C
C                      ENGINEERING SCIENCES, INC.                      C
C                  MODIFIED BY GARY CHENG, 11/14/2000                  C
C                              SECA, INC.                              C
C                      -------------------------                       C
C                           MASTER PROGRAM                             C
C======================================================================C
C
c-----------------------------------------------------------------------
c     lblk: maximum original blocks
c     mproc: maximum process; mblk: maximum fdns block per host
c-----------------------------------------------------------------------
      include 'mpif.h'
      PARAMETER (IIQMAX =  1000000, NSPM = 11, IVOFMX = 2)
      PARAMETER (MZON = 200)
c
      character nstcoef(50,2)*80
c
      DIMENSION X(IIQMAX),Y(IIQMAX),Z(IIQMAX),DEN(IIQMAX),U(IIQMAX),
     &  V(IIQMAX),W(IIQMAX),P(IIQMAX),TM(IIQMAX),DK(IIQMAX),
     &  DE(IIQMAX),FM(IIQMAX,NSPM),VOFS(IIQMAX,IVOFMX),WORK(IIQMAX),
     &  VOF(IIQMAX),AP(IIQMAX)  ! GC: array for quality
C-GC
      PARAMETER (MEXT=100) ! DO NOT CHANGE THIS DIMENSION
      DIMENSION AEXT(MEXT),IEXT(MEXT)
      parameter (mproc=1000, mblk=mzon, lblk=1000)
      dimension lproc(lblk),lizon(mproc),lzone(mblk,mproc)
      dimension izs(lblk),izt(lblk),jzt(lblk),kzt(lblk)
      character system_command*160,null*80,fname*80
      character myname*80,myarch*80,cnumber*8
C-GC1
      COMMON /COMMPI/NPROC,IPROC,MYID
      call MPI_INIT(IERR)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, isize, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, myid, ierr)
      nproc = isize
      iproc = myid+1
      IF(MYID .EQ. 0) THEN
        print*, 'My id : ', myid 
        print*, 'NO. of Processes Working : ' , isize 
        print *, '## hello from process ', myid
      ENDIF
      if (nproc .eq. 0) then
        print *, 'Are you out of your mind? Exit and run again'
c-gc    call exitnet
        CALL MPI_FINALIZE(IERR)
        STOP
      endif
      if (nproc .eq. 1) then
        print *, 'This is a serial computation'
        myid = 0
      endif
C-GC2
c
c-----------------------------------------------------------------------
      do i=1,80 ! assign an null array
        null(i:i)=' '
      enddo
c--------------------catch input file-----------------------------------
      call read_fort11(0,0,izon,izs,izt,jzt,kzt,lproc,nstcoef)
      IF(MYID .EQ. 0) print *, 'reading fort11 done'
c-----------------------------------------------------------------------
      do i=1,nproc
        lizon(i)=0
      enddo
      do ii=1,izon
        i=lproc(ii)
        lizon(i)=lizon(i)+1
        lzone(lizon(i),i)=ii
      enddo
c-----------------------------------------------------------------------
c
 500  continue
      IF(MYID .EQ. 0) THEN
        print*,' '
        print*,'ICONT = 0: eixt'              
        print*,'        1: distribute FDNS input file:         fort.11' 
        print*,'        2: distribute FDNS/PLOT3D grid file:   fort.12'
        print*,'        3: distribute FDNS restart file:       fort.13'
        print*,'        4: collect FDNS/PLOT3D grid file:      fort.22'
        print*,'        5: collect FDNS restart file:          fort.23'
        print*,'        6: collect PLOT3D solution files:      fort.9*'
ccc     print*,'        7: mv fort.22 fort.12 at current (master) dir'
ccc     print*,'        8: mv fort.23 fort.13 at current (master) dir'
ccc     print*,'        9: mv f(proc).22 f(proc).12 for restarting at ',
ccc  &                    'working dir'
ccc     print*,'       10: mv f(proc).23 f(proc).13 for restarting at ',
ccc  &                    'working dir'
ccc     print*,'       11: make   working dir'
ccc     print*,'       12: delete working dir'
ccc     print*,'Caution: Options 9,10,12 may not work properly'
        print*,'For WIN32, fort.11, fort.12, fort.13 must be put on ',
     &         'mpi_tmp dir'
        print*,' '
 510    continue
        print*,'Please input ICONT'
        read(*,*) ICONT
        if(.not.(icont.ge.0.and.icont.le.12)) then
          print*,'Wrong option. Try again.'
          goto 510
        endif
 520    continue
        if(icont.ge.2.and.icont.le.6) then
          print*,'Please input IFMT1,IFMT2 for input & output file',
     &           ' format:'
          print*,' = 1: unformatted; 2: formatted'              
          read(*,*) IFMT1,IFMT2
          if(ifmt1.lt.1.or.ifmt1.gt.2) then
            print*,'Wrong ifmt1 option. Try again.'
            goto 520
          endif
          if(ifmt2.lt.1.or.ifmt2.gt.2) then
            print*,'Wrong ifmt2 option. Try again.'
            goto 520
          endif
        endif
        if(icont.eq.6) then
          print*,'Please input NGAS (number of species) for fort.94-99'
          read(*,*) ngas ! output plot3d q file for mass fraction
        endif
      ENDIF
c-----------------------------------------------------------------------
      CALL MPI_BCAST(ICONT, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
      if(icont.eq.0) goto 900
c-gc  print*, 'My id, ICONT: ',myid ,ICONT
      if(icont.eq.0) goto 900
      if(icont.ge.2 .and. icont.le.6) then
        CALL MPI_BCAST(IFMT1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
        CALL MPI_BCAST(IFMT2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
        IF(ICONT .EQ. 6)
     &    CALL MPI_BCAST(NGAS,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      endif
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      if(icont.eq.1) call write_fort11_parallel(icont,
     &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt)
c
C-GC  if(icont.eq.2) call distribute_plot3dg(icont,ifmt1,ifmt2,
C-GC &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
C-GC &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
      if(icont.eq.2) call distribute_grid(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,work)
c
C-GC  if(icont.eq.3) call distribute_restart(icont,ifmt1,ifmt2,
C-GC &  lproc,lizon,lzone,mblk,mext,aext,iext,
C-GC &  izon,izs,izt,jzt,kzt,
C-GC &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
      if(icont.eq.3) call distribute_flow(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,mext,nspm,izon,
     &  izs,izt,jzt,kzt,iiqmax,den,u,v,w,p,tm,dk,de,fm,work,vof,ap)
c
C-GC  if(icont.eq.4) call collect_plot3dg(icont,ifmt1,ifmt2,
C-GC &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
C-GC &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
      if(icont.eq.4) call collect_grid(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,work)
c
C-GC  if(icont.eq.5) call collect_restart(icont,ifmt1,ifmt2,
C-GC &  lproc,lizon,lzone,mblk,mext,aext,iext,
C-GC &  izon,izs,izt,jzt,kzt,
C-GC &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
      if(icont.eq.5) call collect_flow(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,mext,nspm,izon,izs,izt,jzt,
     &  kzt,iiqmax,den,u,v,w,p,tm,dk,de,fm,work,vof,ap)
c
      if(icont.eq.6) then
        call     collect_plot3dq(icont,ifmt1,ifmt2,
     &    '92',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &    iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
        call     collect_plot3dq(icont,ifmt1,ifmt2,
     &    '93',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &    iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
c-gc    print*,'Please input NGAS (number of species) for fort.94-99'
c-gc    read(*,*) ngas ! output plot3d q file for mass fraction
        if(ngas.ge. 2)
     &    call   collect_plot3dq(icont,ifmt1,ifmt2,
     &      '94',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &      iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
        if(ngas.ge. 6)
     &    call   collect_plot3dq(icont,ifmt1,ifmt2,
     &      '95',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &      iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
        if(ngas.ge.11)
     &    call   collect_plot3dq(icont,ifmt1,ifmt2,
     &      '96',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &      iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
        if(ngas.ge.16)
     &    call   collect_plot3dq(icont,ifmt1,ifmt2,
     &      '97',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &      iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
        if(ngas.ge.21)
     &    call   collect_plot3dq(icont,ifmt1,ifmt2,
     &      '98',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &      iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
        if(ngas.ge.26.and.ngas.le.30)
     &    call   collect_plot3dq(icont,ifmt1,ifmt2,
     &      '99',lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &      iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
      endif
c
      if(icont.eq.7) call system('mv fort.22 fort.12')
c
      if(icont.eq.8) call system('mv fort.23 fort.13')
c
      if(icont.ge.9.and.icont.le.12) then
c-gc    call inumber_cnumber(my_rank+1,cnumber) ! process number
        call inumber_cnumber(iproc,cnumber) ! process number
        system_command=null//null
        fname(1:5)='f'//cnumber(6:8)//'.' 
        if(icont.eq. 9) system_command = 'mv ' 
     &    //fname(1:6)//'22' //fname(1:6)//'12'
        if(icont.eq.10) system_command = 'mv ' 
     &    //fname(1:6)//'23' //fname(1:6)//'13'
        print*,system_command
c-gc    call info_mpi(info,nbtid)
        mtag=10000
        if (myid .eq. 0) then 
           call system(system_command) 
c          call netcast_string(system_command,1)
           leng=len(system_command)
           CALL MPI_BCAST(LENG,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
           CALL MPI_BCAST(SYSTEM_COMMAND,LENG,MPI_CHARACTER,0,
     &                    MPI_COMM_WORLD,IERR)
        else
c          call netcast_string(system_command,1) 
           CALL MPI_BCAST(LENG,1,MPI_INTEGER,0,MPI_COMM_WORLD,IERR)
           CALL MPI_BCAST(SYSTEM_COMMAND,LENG,MPI_CHARACTER,0,
     &                    MPI_COMM_WORLD,IERR)
           write(*,'(a)') system_command
           call system(system_command)
           write(*,*) 'finished'
        endif        
      endif
c
      goto 500
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
 900  continue
c
      call MPI_FINALIZE(ierr)
c
      stop
      end
C
C-GC  subroutine write_fort11_parallel(icont,nproc,
      subroutine write_fort11_parallel(icont,
     &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt)
C
      PARAMETER (MZON = 200, MBIF = 200, MBIO = 100, MBWA = 200,
     &           MBSN =   6, MPNT =  50, MZPO =  20)
      PARAMETER (NSPM =  11, MSP  =NSPM, MST  = 100)
C
c MZON     : Maximum number of blocks (zones)    ( > IZON   )
c MBIF     : Maximum number of interfaces        ( > IZFACE )
c MBIO     : Maximum number of flow boundaries   ( > IBND   )
c MBWA     : Maximum number of wall segments     ( > ID     )
c MBSN     : Maximum number of singularity lines ( > ISNGL  )
c MPNT     : Maximum number of user print lines  ( > INPNT  )
c MZPO     : Maximum number of porosity inputs   ( > IDPORO )
c NSPM     : Maximum number of chemistry species ( > NGAS   )
C
c
      character nstcoef(50,2)*80,aa80(80)*1
c
      DIMENSION IZS(*),IZT(*),JZT(*),KZT(*),LPROC(*),LIZON(*),
     &  LZONE(MBLK,*)
c
      COMMON /FDNSCOM/ITITLE,IDIM,IZFACE,IBND,ID,ISNGL,INPNT,
     & CBG(mzon,3),CBV(mzon,3),THCYCX(mbif),THCYCY(mbif),THCYCZ(mbif),
     & IFCYC(MBIF),IZB1(mbif),IZB2(mbif),IZF1(mbif),IZF2(mbif),
     & IJZ11(mbif),IJZ12(mbif),JKZ11(mbif),JKZ12(mbif),INONUF(mbif),
     & IJZ21(mbif),IJZ22(mbif),JKZ21(mbif),JKZ22(mbif),
     & IBCZON(mbio),IDBC(mbio),ITYBC(mbio),IJBB(mbio),IJBS(mbio),
     & IJBT(mbio),IKBS(mbio),IKBT(mbio),IVFINT(mbio),
     & IWBZON(mbwa),L1(mbwa),L2(mbwa),M1(mbwa),M2(mbwa),
     & N1(mbwa),N2(mbwa),
     & IWTM(mbwa),HQDOX(mbwa),IWALL(mbwa),DENNX(mbwa),VISWX(mbwa),
     & ISNZON(mbsn),ISNBC(mbsn),ISNAX(mbsn),ISNBS(mbsn),ISNBT(mbsn),
     & IPRZON(50),LP1(50),LP2(50),MP1(50),MP2(50),NP1(50),NP2(50),
     & IDATA,IGEO,ITT,ITPNT,ICOUP,NLIMT,IAX,ICYC,
     & DTT,IREC,REC,THETA,BETAP,IEXX,PRAT,
     & IPC,JPC,IPEX,JPEX,IMN,JMN,GFORC1,GFORC2,GFORC3,
     & VISC,IG,ITURB,AMC,GAMA,CBE,CBH,EREXT,
     & ISWU,ISWP,ISWK,ISKEW,INSO(12),IUNIT,
     & DNREFX,UREFX,TREFX,XREFX,PREFX,IOFINN,IOFOUT,IGDINN,IOP3DOUT,
     & NGAS,GAS_NAME(nspm),CEC_DATA(nspm),WTMOLE(nspm),HF(7,2,nspm),
     & NREACT, ITHIRD(mst),IGLOB(mst),ARRHA(mst),ARRHN(mst),ARRHB(mst),
     & STCOEF(mst,msp),STCOEG(mst,msp),
     & PTFREE,TMFREE,XMFREE,ISPN2,ISPO2,
     & ISPVOF,DENVOF(mbio),IVFPTL(10),
     & IEMBED,IEMBLK(50),IEMIIS(50),IEMIIT(50),IEMINC(50),
     & JEMIIS(50),JEMIIT(50),JEMINC(50),KEMIIS(50),
     & KEMIIT(50),KEMINC(50)
C-GC & DTT,IREC,REC,THETA,BETAP,IEXX,PRAT(mbio),
C-GC & IPC,JPC,IPEX(mbio),JPEX(mbio),IMN,JMN,GFORC1,GFORC2,GFORC3,
C
      COMMON/WEDGE/IWEDGE,IWEDGE_IZ(10),IWEDGE_I1(10),IWEDGE_I2(10),
     &  IWEDGE_J1(10),IWEDGE_J2(10),IWEDGE_K1(10),IWEDGE_K2(10)
      COMMON/Q2NAME/P3DQ2_NAME
      CHARACTER P3DQ2_NAME(5)*12
      COMMON/LGROUP/LPTFREE,LISPVOF,LNGAS,LNREACT,LIEMBED,LPLOT3D,
     &              LFLUID  !GC: flag for real fluid model control cards
C-GC  Spark igniter group
      COMMON /SPARK/ LISPARK,ISPARK,ISPKMIN,ISPKMAX,ISPKON(10),
     &  ISPKZN(10),ISPKI1(10),ISPKIM(10),ISPKJ1(10),ISPKJM(10),
     &  ISPKK1(10),ISPKKM(10),TMSPK(10),ISPKDBG(10),ISPKCNT(10)
C
      CHARACTER GAS_NAME*24,CEC_DATA*80,ITITLE(72)*1
c
      common /zonmpi/iproc1(1000),iproc2(1000),idface(1000) ! mpi.04
c
C-GC1
      CHARACTER SPECIES*20
      COMMON /FLUIDS/IMIX,ISRK,ITERMAX,IPSPECS,IPGEN,IPCRIT,IPSAT,
     &               IPITER,TOLER0,SMALLEST,SPECIES(NSPM),ITHERMO(NSPM)
      COMMON /COMMPI/NPROC,IPROC,MYID
C-GC2
      parameter (mproc=1000,lblk=1000)
      dimension lzface(100,mproc),liproc1(100,mproc),liproc2(100,mproc),
     &  lzbnd(100,mproc),lzd(100,mproc),lzsngl(100,mproc),
     &  lznpnt(100,mproc),lprozon(lblk),lzonzon(lblk),
     &  lizface(mproc),libnd(mproc),lid(mproc),
     &  lisngl(mproc),linpnt(mproc)
      character*80 fname
      character null*80,cnumber*8,acha(1000)*160
      include 'mpif.h'
c
      do i=1,80
        null(i:i)=' '
      enddo
C-GC1
cc    IF(MYID .NE. 0) RETURN
      IWEDGE = 0
      LPTFREE = 0
      LISPVOF = 0
      LNGAS = 0
      LNREACT = 0
      LIEMBED = 0
      LPLOT3D = 0
      LFLUID = 0
C-GC2
c
      call read_fort11(icont,0,izon,izs,izt,jzt,kzt,lproc,nstcoef)
c
      izs(1)=0
      do iz=2,izon
        izs(iz)=izs(iz-1)+izt(iz-1)*jzt(iz-1)*kzt(iz-1)
      enddo
c
      if(icont.ne.1) return
c
C=======================================================================
C     in the current implement, the whole grid block must 
C     stay in one process. further implement will allow put
C     one grid block into different processes.
C=======================================================================
c
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C     determine lizon()
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-----lizon(i): number of zones in i-th process
c-----lprozon(iz): process number for iz-th zone 
c-----lzonzon(iz): local zone number for golbal iz-th zone 
c-----lzone(j,i): global zone number for zone j in i process
      do i=1,nproc
        lizon(i)=0
      enddo
      IF(MYID .EQ. 0) print*,'izon=',izon,' nproc=',nproc
      do ii=1,izon
        i=lproc(ii)
        lprozon(ii)=i 
        lizon(i)=lizon(i)+1
        lzonzon(ii)=lizon(i)
        lzone(lizon(i),i)=ii
      enddo
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-----determine lizface()
c     lizface(i): number of zonal faces in i-th process
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i=1,nproc
        lizface(i)=0
      enddo
      do ii=1,izface
        idface(ii)=ii
        ip1=lprozon(izb1(ii))
        ip2=lprozon(izb2(ii))
c
        jf1=lizface(ip1)+1
        lizface(ip1)=jf1
        liproc1(jf1,ip1)=ip1
        liproc2(jf1,ip1)=ip2
        lzface(jf1,ip1)=ii
c
c       if(izb1(ii).ne.izb2(ii)) then !such as cylclic bc in same zone
        if(ip1.ne.ip2) then 
          jf2=lizface(ip2)+1
          lizface(ip2)=jf2 
          liproc1(jf2,ip2)=ip1
          liproc2(jf2,ip2)=ip2
          lzface(jf2,ip2)=ii
        endif
      enddo
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-----determine libnd()
c     libnd(i): number of boundary conditions in i-th process
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i=1,nproc
        libnd(i)=0
      enddo
      do ii=1,ibnd
        i=lprozon(ibczon(ii))
        libnd(i)=libnd(i)+1
        lzbnd(libnd(i),i)=ii 
      enddo
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-----determine lid()
c     lid(i): number of wall bcs in i-th process
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i=1,nproc
        lid(i)=0
      enddo
      do ii=1,id
        i=lprozon(iwbzon(ii))
        lid(i)=lid(i)+1
        lzd(lid(i),i)=ii
      enddo
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-----determine lisngl()
c     lisngl(i): number of singular bcs in i-th process
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i=1,nproc
        lisngl(i)=0
      enddo
      do ii=1,isngl
        i=lprozon(isnzon(ii))
        lisngl(i)=lisngl(i)+1
        lzsngl(lisngl(i),i)=ii
      enddo
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c-----determine linpnt()
c     linpnt(i): number of printing blocks in i-th process
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      do i=1,nproc
        linpnt(i)=0
      enddo
      do ii=1,inpnt
        i=lprozon(iabs(iprzon(ii)))
        linpnt(i)=linpnt(i)+1
        lznpnt(linpnt(i),i)=ii
      enddo
c
C=======================================================================
C
      IW=40
c
c-----open output files fort.11
      i1=1 ! write 'f'//cnumber(6:8)//'.11' first in master dir
c-gc  DO I=1,NPROC
      I=MYID+1
        call inumber_cnumber(i,cnumber) ! process number
        fname(1:7)='f'//cnumber(6:8)//'.11' 
c-gc    print*,fname(1:7)
        open(IW,file=fname(1:8),form='formatted')
C
C=======================================================================
c
card0 title 
C-GC    WRITE(IW,'(80(1H#))')
c       WRITE(IW,'(7HTITLE: ,72A1)') (ITITLE(II),II=1,72)
        WRITE(IW,'(72A1)') (ITITLE(II),II=1,72)
C-GC    WRITE(IW,'(80(1H#))')
card1 dimension (2d or 3d)
        WRITE(IW,'("  IDIM")')
        WRITE(IW,'(I6)') IDIM
card2 zonal information and number bcs
C-GC    WRITE(IW,'("#IZON,IZFACE,  IBND,    ID, ISNGL, INPNT")')
        WRITE(IW,'("  IZON,IZFACE,  IBND,    ID, ISNGL")')
C-GC    WRITE(IW,'(I5,7I7)') LIZON(I),LIZFACE(I),LIBND(I),LID(I),
C-GC &    LISNGL(I),LINPNT(I)
        WRITE(IW,'(I6,7I7)') LIZON(I),LIZFACE(I),LIBND(I),LID(I),
     &    LISNGL(I)
card3 zonal grid size and rotational/translational speeds
C-GC    WRITE(IW,'("#  IZT,  JZT,  KZT,LPROC,      CBG1,      ",
C-GC &    "CBV2")')
        WRITE(IW,'("   IZT, JZT, KZT,LPROC,    CBGX,    ",
     &    "CBGY,    CBGZ,    CBVX,    CBVY,    CBVZ")')
        DO II=1,LIZON(I)
          IZ=LZONE(II,I) ! IZ/II is the global/local zone number
          WRITE(IW,'(I6,2I5,I6,1P6E9.1)') IZT(IZ),JZT(IZ),KZT(IZ),
     &      LPROC(IZ),CBG(IZ,1),CBG(IZ,2),CBG(IZ,3),CBV(IZ,1),
     &      CBV(IZ,2),CBV(IZ,3)
C-GC      WRITE(IW,'(4I6,1P6E11.3)') IZT(IZ),JZT(IZ),KZT(IZ),
C-GC &      LPROC(IZ),CBG(IZ,1),CBV(IZ,2)
        ENDDO
card4 zonal interface matching indices
C-GC    WRITE(IW,'("#THCYCX, IZB1, IZF1,IJZ11,IJZ12,JKZ11,JKZ12,",
        WRITE(IW,'(" IFCYC, IZB1, IZF1,IJZ11,IJZ12,JKZ11,JKZ12,",
     &    "INONUF,IPROC1")')
C-GC    WRITE(IW,'("#        IZB2, IZF2,IJZ21,IJZ22,JKZ21,JKZ22,",
        WRITE(IW,'("        IZB2, IZF2,IJZ21,IJZ22,JKZ21,JKZ22,",
     &    "IDFACE,IPROC2")')
        DO II=1,LIZFACE(I)
          IZ=LZFACE(II,I) ! IZ/II is the global/local face number
          IF(INONUF(IZ).EQ.0.AND.(LIPROC1(II,I).NE.LIPROC2(II,I))) 
     &      INONUF(IZ)=51
C-GC      WRITE(IW,'(F7.2,6I6,2I7)') THCYCX(IZ),LZONZON(IZB1(IZ)),
          WRITE(IW,'(7I6,2I7)') IFCYC(IZ),LZONZON(IZB1(IZ)),
     &      IZF1(IZ),IJZ11(IZ),IJZ12(IZ),JKZ11(IZ),
     &      JKZ12(IZ),INONUF(IZ),LIPROC1(II,I)
C-GC      WRITE(IW,'(7X,6I6,2I7)')              LZONZON(IZB2(IZ)),
          WRITE(IW,'(6X,6I6,2I7)')              LZONZON(IZB2(IZ)),
     &      IZF2(IZ),IJZ21(IZ),IJZ22(IZ),JKZ21(IZ),
     &      JKZ22(IZ),IDFACE(IZ),LIPROC2(II,I) 
        ENDDO
c
card5 flow boundary
C-GC    WRITE(IW,'("#IBCZON,IDBC,ITYBC,IJBB,IJBS,IJBT,JKBS,",
C-GC &    "JKBT,IVFINT,      PRAT,IPZ,IPI,IPJ,IPK")')
        WRITE(IW,'("IBCZON,IDBC,ITYBC,IJBB,IJBS,IJBT,JKBS,",
     &    "JKBT,IBG")')
        DO II=1,LIBND(I)
          IZ=LZBND(II,I) ! IZ/II is the global/local bc number
C-GC      CALL BINDEX(IPEX(IZ),IPZN,IPII,IPJJ,IPKK,IDIM,IZON,
C-GC &      IZS,IZT,JZT,KZT) ! ipzn=jpex(ii)
c---------determine jpex reference indices at each process  
c         negative number means the reference point is not in this process
C-GC      LIPZN=LPROZON(JPEX(IZ)) ! PROCESSOR NUMBER
C-GC      IF(LIPZN.EQ.I) THEN
C-GC        LIPZN=LZONZON(JPEX(IZ)) ! hms
C-GC        LIPII=IPII
C-GC        LIPJJ=IPJJ
C-GC        LIPKK=IPKK
C-GC      ELSE
C-GC        LIPZN=-LIPZN 
C-GC        LIPII=IZ
C-GC        LIPJJ=0
C-GC        LIPKK=0
C-GC      ENDIF
C---------IF LIPZN>0: IN THE CURRENT PROCESSOR,  LIPZN=THE # OF LOCAL ZONE
C---------IF LIPZN<0: IN THE OTHER   PROCESSOR, -LIPZN=THE # OF PROCESSOR,
C---------            LPII=GLOBAL IBCZON INDEX, LPJJ=LPKK=0   
C-GC      WRITE(IW,'(I7,I5,I6,5I5,I7,1PE11.3,5I4)') LZONZON(IBCZON(IZ)),
C-GC &      IDBC(IZ),ITYBC(IZ),IJBB(IZ),IJBS(IZ),IJBT(IZ),IKBS(IZ),
C-GC &      IKBT(IZ),IVFINT(IZ),PRAT(IZ),LIPZN,LIPII,LIPJJ,LIPKK,IZ
          WRITE(IW,'(I6,I5,I6,5I5,I4)') LZONZON(IBCZON(IZ)),
     &      IDBC(IZ),ITYBC(IZ),IJBB(IZ),IJBS(IZ),IJBT(IZ),IKBS(IZ),
     &      IKBT(IZ),IZ
        ENDDO
c
card6 wall block indices
C-GC    WRITE(IW,'("#IWBZON,  L1,  L2,  M1,  M2,  N1,  N2,IWTM,     ",
        WRITE(IW,'("IWBZON,  L1,  L2,  M1,  M2,  N1,  N2,IWTM,     ",
     &    "HQDOX,IWALL,    DENNX,    VISWX")')
        DO II=1,LID(I)
          IZ=LZD(II,I) ! IZ/II is the global/local wall number
C-GC      WRITE(IW,'(I7,6I5,I5,1PE11.3,I6,2E10.3)')  
          WRITE(IW,'(I6,6I5,I5,1PE11.3,I6,2E10.3)')  
     &      LZONZON(IWBZON(IZ)),L1(IZ),L2(IZ),M1(IZ),M2(IZ),N1(IZ),
     &      N2(IZ),IWTM(IZ),HQDOX(IZ),IWALL(IZ),DENNX(IZ),VISWX(IZ)
        ENDDO
c
card7 singularity lines
C-GC    WRITE(IW,'("#ISNZON, ISNBC, ISNAX, ISNBS, ISNBT")')
        WRITE(IW,'("ISNZON, ISNBC, ISNAX, ISNBS, ISNBT")')
        DO II=1,LISNGL(I)
          IZ=LZSNGL(II,I) ! IZ/II is the global/local singularity number
C-GC      WRITE(IW,'(9I7)') LZONZON(ISNZON(IZ)),ISNBC(IZ),ISNAX(IZ),
          WRITE(IW,'(I6,4I7)') LZONZON(ISNZON(IZ)),ISNBC(IZ),ISNAX(IZ),
     &      ISNBS(IZ),ISNBT(IZ)
        ENDDO
c
card8 zonal printing indices
C-GC    WRITE(IW,'("#IPRZON,  LP1,  LP2,  MP1,  MP2,  NP1,  NP2")')
C-GC    DO II=1,LINPNT(I)
C-GC      IZ=LZNPNT(II,I) ! IZ/II is the global/local zonal printing number
C-GC      JJ=IPRZON(IZ)
C-GC      WRITE(IW,'(I7,9I6)') ISIGN(LZONZON(IABS(JJ)),JJ),LP1(IZ),
C-GC &      LP2(IZ),MP1(IZ),MP2(IZ),NP1(IZ),NP2(IZ),IZ,IABS(JJ)
C-GC    ENDDO
c
card9 i/o and problem control parameters
C-GC    WRITE(IW,'("#IDATA, IGEO,  ITT,ITPNT,ICOUP,NLIMT,  IAX, ",
        WRITE(IW,'(" IDATA,  IGEO,   ITT, ITPNT, ICOUP, NLIMT,   IAX,",
     &    "  ICYC")')
        WRITE(IW,'(9(I6,1X))')IDATA,IGEO,ITT,ITPNT,ICOUP,NLIMT,IAX,ICYC
c
card10 time step soze and schemes
C-GC    WRITE(IW,'("#       DTT,IREC,   REC, THETA, BETAP,IEXX")')
        WRITE(IW,'("        DTT,IREC,   REC, THETA, BETAP,IEXX,"
     &             "  PRAT")')
C-GC    WRITE(IW,'(1PE11.3,I5,0P3F7.4,I5)') 
C-GC &    DTT,IREC,REC,THETA,BETAP,IEXX
        WRITE(IW,'(1PE11.3,I5,0P3F7.4,I5,F7.2)') 
     &    DTT,IREC,REC,THETA,BETAP,IEXX,PRAT
c
card11 reference and mointering points
c-------determine ljpc and ljmn reference indices at each process  
c       negative number means the reference point is not in this process
        ljpc=-lprozon(jpc)
        if(lprozon(jpc).eq.i) ljpc= lzonzon(jpc)
        ljmn=-lprozon(jmn)
        if(lprozon(jmn).eq.i) ljmn= lzonzon(jmn)
C-GC    WRITE(IW,'("#   IPC,   JPC,   IMN,   JMN,    ",
C-GC &    "GFORC1,    GFORC2,    GFORC3")')
C-GC    WRITE(IW,'(4I7,1P3E11.3)') IPC,LJPC,IMN,LJMN,GFORC1,GFORC2,
C-GC &    GFORC3
C-GC
        LJPEX=-LPROZON(JPEX)
        IF(LPROZON(JPEX).EQ.I) LJPEX= LZONZON(JPEX)
        WRITE(IW,'("    IPC,   JPC,  IPEX,  JPEX,   IMN,   JMN")')
        WRITE(IW,'(6I7)') IPC,LJPC,IPEX,LJPEX,IMN,LJMN
c
card12 reference viscosity, mach number and turbulence model
C-GC    WRITE(IW,'("#      VISC,IG,ITURB,   AMC,  GAMA,       ",
        WRITE(IW,'("       VISC,IG,ITURB,   AMC,  GAMA,       ",
     &    "CBE,       CBH,     EREXT")')
        WRITE(IW,'(1PE12.5,I3,I6,0PF7.4,F7.4,1P3E11.3)') 
     &    VISC,IG,ITURB,AMC,GAMA,CBE,CBH,EREXT
c
card13 martix solver and skewness term
C-GC    WRITE(IW,'("#ISWU,ISWP,ISWK,ISKEW")')
        WRITE(IW,'(" ISWU,ISWP,ISWK,ISKEW")')
        WRITE(IW,'(3I5,I6)') ISWU,ISWP,ISWK,ISKEW
c
card14 solution variables
C-GC    WRITE(IW,'("#INSO(IEQ):")')
C-GC    WRITE(IW,'("# U, V, W, T,DK,DE,07,08,09,VS,FM,SP")')
        WRITE(IW,'(" INSO(IEQ):")')
        WRITE(IW,'("  U, V, W, T,DK,DE,07,08,09,VS,FM,SP")')
        WRITE(IW,'(12I3)') (INSO(IEQ),IEQ=1,12)
c
card15 unit and rerefnces
C-GC    WRITE(IW,'("#IUNIT,    DENREF,      UREF,      TREF,      ",
C-GC &    "XREF")')
        WRITE(IW,'(" NGAS, NREACT, IUNIT,    DENREF,      UREF,",
     &             "      TREF,      XREF, PREF")')
        WRITE(IW,'(I5,I8,I7,1P5E11.3)')NGAS,NREACT,IUNIT,DNREFX,UREFX,
     &                                 TREFX,XREFX,PREFX
c
C-GC   Spark igniter group
        WRITE(IW,'("ISPARK, ISPKMIN, ISPKMAX")')
        IF(NREACT .GE. 1) THEN
          LSPARK = 0
          IF(ISPARK .GE. 1) THEN
            DO KK = 1,ISPARK
              IF(LPROZON(ISPKZN(KK)) .EQ. I) LSPARK = LSPARK+1
            ENDDO
          ENDIF
          WRITE(IW,'(I6,I9,I9)')LSPARK,ISPKMIN,ISPKMAX
        ENDIF
        WRITE(IW,'("ISPKON, ISPKZN, ISPKI1, ISPKIM, ISPKJ1,",
     &    " ISPKJM, ISPKK1, ISPKKM, TMSPK, ISPKDBG")')
        IF(LSPARK .GE. 1) THEN
          DO KK = 1,ISPARK
            IF(LPROZON(ISPKZN(KK)) .EQ. I) THEN
              LSPKZN = LZONZON(ISPKZN(KK))
              WRITE(IW,'(8(I6,2X),F6.1,2X,I7)')ISPKON(KK),LSPKZN,
     &          ISPKI1(KK),ISPKIM(KK),ISPKJ1(KK),ISPKJM(KK),
     &          ISPKK1(KK),ISPKKM(KK),TMSPK(KK),ISPKDBG(KK)
            ENDIF
          ENDDO
        ENDIF
c
card16 input/output restart file format and output plot3d format
C-GC    WRITE(IW,'("#IGDINN,IOFINN,IOFOUT,IOP3DOUT ",
        WRITE(IW,'(" IGDINN,IOFINN,IOFOUT,IOP3DOUT ",
     &    "(1:Unf, 2:Fmt) ! IO FORMAT CONTROL")')
        WRITE(IW,'(3I7,I9)') IGDINN,IOFINN,IOFOUT,IOP3DOUT
card17 
C-GC    WRITE(IW,'(80(1H*))')
C-GC1
c-------NGAS GROUP
        IF(NGAS.GE.1) THEN
          DO II=1,NGAS
C-GC        WRITE(IW,'(A80)'          ) CEC_DATA(II)
C-GC        WRITE(IW,'(A24,44X,F12.5)') GAS_NAME(II),WTMOLE(II)
            WRITE(IW,'(A24,42X,F12.5)') GAS_NAME(II),WTMOLE(II)
            WRITE(IW,'(1P5E15.8)') ((HF(LL,KK,II),LL=1,7),KK=1,2)
          ENDDO
C-GC      WRITE(IW,'(80(1H*))')
        ENDIF
C-------NREACT GROUP
        IF(NREACT.GE.1) THEN
          WRITE(IW,'("REACTION:",15A6)') (GAS_NAME(K)(1:6),K=1,NGAS)
          DO II=1,NREACT
            WRITE(IW,'(I5,1P3E12.4,2I5)') II,ARRHA(II),ARRHN(II),
     &                              ARRHB(II),ITHIRD(II),IGLOB(II)
C-GC        call cutblank(nstcoef(ii,1),aa80,ilast)
C-GC        write(iw,'(80a1)') (aa80(ik),ik=1,ilast)
C-GC        if(iglob(ii).ge.2) then
C-GC          call cutblank(nstcoef(ii,2),aa80,ilast)
C-GC          write(iw,'(80a1)') (aa80(ik),ik=1,ilast)
C-GC        endif
            WRITE(IW,'(5X,15F6.1)') (STCOEF(II,K),K=1,NGAS)
            IF(IGLOB(II).GE.2) WRITE(IW,'(5X,15F6.1)')
     &                (STCOEG(II,K),K=1,NGAS)
          ENDDO
C-GC      WRITE(IW,'(80(1H*))')
        ENDIF
C-GC2
c-------OUTLET_PREF POINT ! ONLY READ BY PARALLEL
C-GC    LOUTLET_PREF=1
C-GC    IF(LOUTLET_PREF.EQ.1) THEN
C-GC      WRITE(IW,'("OUTLET_PREF POINT")')
C-GC      WRITE(IW,*) IBND
C-GC      DO II=1,IBND
C-GC        LIPROC=LPROZON(JPEX(II)) ! PROCESSOR NUMBER
C-GC        LIPZN =LZONZON(JPEX(II)) ! LOCAL ZONE NUMBER
C-GC        LIPEX=0
C-GC        DO J=2,LIPZN
C-GC          IZ=LZONE(J-1,LIPROC)
C-GC          LIPEX=LIPEX+IZT(IZ)*JZT(IZ)*KZT(IZ)
C-GC        ENDDO
C-GC        LIPEX=LIPEX+IPEX(II)-IZS(JPEX(II)) ! LOCAL INDEX NUMBER
C-GC        WRITE(IW,*) IPEX(II),LIPEX,LIPROC,PRAT(II)
C-GC      ENDDO
C-GC      WRITE(IW,'(80(1H*))')
C-GC    ENDIF
c
c-------PLOT3D VARIABLES
        IF(LPLOT3D.EQ.1) THEN
          WRITE(IW,'("PLOT3D fort.93 VARIABLES:")')
          WRITE(IW,'(5(A12,1X))') (P3DQ2_NAME(II),II=1,5)
          WRITE(IW,'(80(1H*))')
        ENDIF
c-------PTFREE GROUP
        IF(LPTFREE.EQ.1) THEN
          WRITE(IW,'("PTFREE,         TMFREE,    XMFREE, ",
     &      "ISPN2, ISPO2")')
          WRITE(IW,'(1P3E11.3,2I7)')  PTFREE,TMFREE,XMFREE,ISPN2,ISPO2
          WRITE(IW,'(80(1H*))')
        ENDIF
c
c-------ISPVOF GROUP
        IF(LISPVOF.EQ.1) THEN
          WRITE(IW,'("ISPVOF")')
          WRITE(IW,*) ISPVOF
          WRITE(IW,'("#DENVOF(K),K=1,ISPVOF")')
          WRITE(IW,*) (DENVOF(K),K=1,ISPVOF)
          WRITE(IW,'("#IVFPTL(K),K=1,ISPVOF")')
          WRITE(IW,*) (IVFPTL(K),K=1,ISPVOF)
          WRITE(IW,'(80(1H*))')
        ENDIF
c
c-------NGAS GROUP
        IF(LNGAS.EQ.1) THEN
          WRITE(IW,'("NGAS")')
          WRITE(IW,*) NGAS
          DO II=1,NGAS
ccc         WRITE(IW,'(A24,44X,F12.5)') GAS_NAME(II),WTMOLE(II)
            WRITE(IW,'(A80)'          ) CEC_DATA(II)
            WRITE(IW,'(5E15.8)') ((HF(LL,KK,II),LL=1,7),KK=1,2)
          ENDDO
          WRITE(IW,'(80(1H*))')
        ENDIF
c
c-------NREACT GROUP
        IF(LNREACT.EQ.1) THEN
          WRITE(IW,'("NREACT")')
          WRITE(IW,*) NREACT
ccc       WRITE(IW,'("#CHEMICAL REACTIONS:")')
          DO II=1,NREACT
            WRITE(IW,'(I5,3E12.4,2I5)') II,ARRHA(II),ARRHN(II),ARRHB(II)
     &               ,ITHIRD(II),IGLOB(II)
c-ysc
            call cutblank(nstcoef(ii,1),aa80,ilast)
            write(iw,'(80a1)') (aa80(ik),ik=1,ilast)
            if(iglob(ii).ge.2) then
              call cutblank(nstcoef(ii,2),aa80,ilast)
              write(iw,'(80a1)') (aa80(ik),ik=1,ilast)
            endif
c-ysc
c-ysc       WRITE(IW,'(5X,50F8.3)') (STCOEF(II,K),K=1,NGAS)
c-ysc       IF(IGLOB(II).GE.2) WRITE(IW,'(5X,50F8.3)')
c-ysc&                (STCOEG(II,K),K=1,NGAS)
c-ysc
          ENDDO
          WRITE(IW,'(80(1H*))')
        ENDIF
c
c-------IEMBED GROUP
        IF(LIEMBED.EQ.1) THEN
c
          iembed_i=0 ! count number of embeded block in this process
          do kk=1,iembed
            do j=1,lizon(i)
              iz=lzone(j,i)
              if(iz.eq.iemblk(kk)) iembed_i=iembed_i+1
            enddo
          enddo
c
          WRITE(IW,'("IEMBED")')
ccc       WRITE(IW,*) IEMBED
          WRITE(IW,*) iembed_i
          WRITE(IW,'("#")')
          WRITE(IW,'("#(I = 1,IEMBED) for the following input sets ",
     &      "--- IEMBED sets")')
          DO KK=1,IEMBED
            WRITE(IW,'("#     IEMBLK")')
            do j=1,lizon(i)
              iz=lzone(j,i)
              if(iz.eq.iemblk(kk)) then
ccc             WRITE(IW,*) IEMBLK(KK)
                WRITE(IW,*) j
                WRITE(IW,'("#          IEMIIS  IEMIIT  IEMINC")')
                WRITE(IW,*) IEMIIS(KK),IEMIIT(KK),IEMINC(KK)
                WRITE(IW,'("#          JEMIIS  JEMIIT  JEMINC")')
                WRITE(IW,*) JEMIIS(KK),JEMIIT(KK),JEMINC(KK)
                WRITE(IW,'("#          KEMIIS  KEMIIT  KEMINC")')
                WRITE(IW,*) KEMIIS(KK),KEMIIT(KK),KEMINC(KK)
              endif
            enddo
          ENDDO
          WRITE(IW,'("#")')
          WRITE(IW,'(80(1H*))')
        ENDIF
C
c-------IWEDGE GROUP
        IF(IWEDGE.GE.1) THEN
c
          iwedge_i=0 ! count number of wedge block in this process
          do kk=1,iwedge
            do j=1,lizon(i)
              iz=lzone(j,i)
              if(iz.eq.iwedge_iz(kk)) iwedge_i=iwedge_i+1
            enddo
          enddo
c
          if(iwedge_i.ge.1) then
          WRITE(IW,'("IWEDGE")')
ccc       WRITE(IW,*) IWEDGE
          WRITE(IW,*) iwedge_i
          WRITE(IW,'("#")')
          WRITE(IW,'("#(I = 1,IWEDGE) for the following input sets ")')
          WRITE(IW,'("#   IZ    I1    I2    J1    J2    K1    K2")')
          DO KK=1,IWEDGE
            do j=1,lizon(i)
              iz=lzone(j,i)
              if(iz.eq.iwedge_iz(kk)) then
                WRITE(IW,'(7I6)')    LZONZON(IZ),IWEDGE_I1(KK),
     &               IWEDGE_I2(KK),IWEDGE_J1(KK),IWEDGE_J2(KK),
     &               IWEDGE_K1(KK),IWEDGE_K2(KK)
              endif
            enddo
          ENDDO
          WRITE(IW,'("#")')
          WRITE(IW,'(80(1H*))')
          endif
        ENDIF
C-GC1
C-------FLUID GROUP
        IF(LFLUID.GE.1) THEN
          WRITE(IW,'("FLUID")')
c         WRITE(IW,'("c")')
c         WRITE(IW,'("c     file fluid.inp - flags for real fluids ",
c    &               "model")')
c         WRITE(IW,'("c")')
c         WRITE(IW,'(" 4,      iunits (1= Eng,lbm/ft3; 2= Eng slugs/",
c    &               "ft3; 3= cgs; 4= mks)")')
c         WRITE(IW,'(I2,",      imix   (1= ideal gas mixture (not ",
c    &             "available yet); 2= ideal solution)")')IMIX
c         WRITE(IW,'(I2,",      isrk   (1= use SRK mixing rules for ",
c    &             "first guesses; 0= don,t)")') ISRK
c         WRITE(IW,'(I4,",    itermax (max iterations for fluids)")')
c    &             ITERMAX
c         WRITE(IW,'(I2,",      ipspecs(n= print intermediate ",
c    &             "results for species n, otherwise=0)")') IPSPECS
c         WRITE(IW,'("c     the following prints occur for species n",
c    &             " when ipspecs=n")')
c         WRITE(IW,'(I2,",      ipgen  (1= general fluids printout,",
c    &             " 2= to pause after printout)")') IPGEN
c         WRITE(IW,'(I2,",      ipcrit (1= print crit. props, 2= to",
c    &             " pause after printout)")') IPCRIT
c         WRITE(IW,'(I2,",      ipsat  (1= print saturation props, ",
c    &             "2= to pause after printout)")') IPSAT
c         WRITE(IW,'(I2,",      ipiter (1= print iteration results,",
c    &             " 2= to pause after printout)")') IPITER
c         WRITE(IW,'(1PE8.1,", toler (convergence tolerance, ",
c    &             "fractional difference allowed)")') TOLER0
c         WRITE(IW,'(1PE8.1,", smallest (smallest mole fraction ",
c    &             "considered)")') SMALLEST
          WRITE(IW,'("c     Species(a20)  , ideal gas(=0), real fluid",
     &             "(=1)")')
          DO KK=1,NGAS
            WRITE(IW,'(A20,I5)') SPECIES(KK),ITHERMO(KK)
          ENDDO
          WRITE(IW,'("DONE")')
c         WRITE(IW,'("c     interaction parameters for SRK mixing ",
c    &             "rules")')
c         WRITE(IW,'(" 0, 0, 0.0")')
        ENDIF
C-GC2
C
        CLOSE(IW)
c-gc  ENDDO
c-------------------------------------------------------------------
C
      RETURN
      END
c
      subroutine cutblank(nstcoef,aa80,ilast)
      character nstcoef*80,aa80(80)*1
c
        do i=1,80
          k=1+80-i
          if(nstcoef(k:k).ne.' ') then
            ilast = k
            go to 9901
          endif
        enddo
 9901   continue
        do i=1,ilast
          aa80(i) = nstcoef(i:i)
        enddo
c
      return
      end
c
      subroutine distribute_grid(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,work)
c
      DIMENSION X(*),Y(*),Z(*),WORK(*)
      dimension izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      COMMON /COMMPI/NPROC,IPROC,MYID
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----INPUT RESTART FILE
      iunit=30
      iuclose=1
      if(ifmt1.eq.1) open(iunit,file='fort.12',form='unformatted')
      if(ifmt1.eq.2) open(iunit,file='fort.12',form=  'formatted')
      CALL READ_GRID_RFV(IUNIT,IFMT1,IUCLOSE,IGDMAX,IIQMAX,IZON,
     &     IZT,JZT,KZT,X,Y,Z)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      igdmax = 0
      IU=40
      IUCLOSE=1
      I=IPROC
c
ccc   iproc = i
ccc   print*,'$$$ iproc=',i
C-GC  np = igdmax
ccc   print *,'np =',np
c
c-gc  izon = lizon(i)
      IZONL = lizon(i)
      do j=1,lizon(i)
        iz=lzone(j,i)
        igdmax=igdmax+izt(iz)*jzt(iz)*kzt(iz)
      enddo
c-gc  call check_izon_max(izon,mblk,info)
      call check_izon_max(IZONL,mblk,info)
      call check_igdmax_max(igdmax,iiqmax,info)
c
      call inumber_cnumber(iproc,cnumber) ! process number
c-------fort.12 formatted and unfrmatted
c
      fname(1:7)='f'//cnumber(6:8)//'.12'
ccc   print*,'fname=',fname(1:7)
      IF(IFMT2.EQ.1) THEN ! unformatted file
        open(IU,file=fname(1:7),form='unformatted')
        WRITE(IU) IZONL
        IF(MYID .EQ. 0)
     &    WRITE(6,'(" FDNSRFV GRID DATA FILE, IZON=",I5)') IZONL
        DO J=1,IZONL
          II=LZONE(J,I)
          IF(MYID .EQ. 0)
     &      WRITE(6,'("II=",4I8)') II,IZT(II),JZT(II),KZT(II)
          WRITE(IU) IZT(II),JZT(II),KZT(II)
        ENDDO
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          IF(II .GT. 1 .AND. L .LE. 1) 
     &      WRITE(6,*)'@@ Error in Distribute_Grid: IZ, IZS =',II,L
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU) (X(L+IJK), IJK = 1,LMN)
          WRITE(IU) (Y(L+IJK), IJK = 1,LMN)
          WRITE(IU) (Z(L+IJK), IJK = 1,LMN)
        ENDDO
      ENDIF
      IF(IFMT2.EQ.2) THEN ! formatted file
        open(IU,file=fname(1:7),form='formatted')
        WRITE(IU,*) IZONL
        IF(MYID .EQ. 0)
     &    WRITE(6,'(" FDNSRFV GRID DATA FILE, IZON=",I5)') IZONL
        DO J=1,IZONL
          II=LZONE(J,I)
          IF(MYID .EQ. 0)
     &      WRITE(6,'("II=",4I8)') II,IZT(II),JZT(II),KZT(II)
          WRITE(IU,*) IZT(II),JZT(II),KZT(II)
        ENDDO
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          IF(II .GT. 1 .AND. L .LE. 1) 
     &      WRITE(6,*)'@@ Error in Distribute_Grid: IZ, IZS =',II,L
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU,1001) (X(L+IJK), IJK = 1,LMN)
          WRITE(IU,1001) (Y(L+IJK), IJK = 1,LMN)
          WRITE(IU,1001) (Z(L+IJK), IJK = 1,LMN)
        ENDDO
      ENDIF
c
      REWIND(IU)
      IF(IUCLOSE.EQ.1) CLOSE(IU)
C-GC  CALL WRITE_GRID_RFV(IUNIT,IFMT2,IUCLOSE,IGDMAX,IIQMAX,IZON,
C-GC &  IZT(I),JZT(I),KZT(I),X,Y,Z,NP)
      IF(MYID .EQ. 0)
     &  WRITE(6,*) '*** END OF OUTPUT FDNSRFV GRID FILE, IPROC =',I
 1001 FORMAT(5(1P,E16.8))
c
      return
      end
c
      subroutine distribute_flow(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,mext,nspm,izon,izs,izt,jzt,kzt,
     &  iiqmax,den,u,v,w,p,tm,dk,de,fm,work,vof,ap)
c
      DIMENSION DEN(*),U(*),V(*),W(*),P(*),TM(*),DK(*),
     &  DE(*),FM(IIQMAX,*),WORK(*),VOF(*),AP(*)
      dimension izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      COMMON /COMMPI/NPROC,IPROC,MYID
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----INPUT RESTART FILE
      iunit=30
      iuclose=1
      if(ifmt1.eq.1) open(iunit  ,file='fort.13',form='unformatted')
      if(ifmt1.eq.2) open(iunit  ,file='fort.13',form=  'formatted')
cc    if(ifmt1.eq.1) open(iunit+1,file='fort.14',form='unformatted')
cc    if(ifmt1.eq.2) open(iunit+1,file='fort.14',form=  'formatted')
      CALL READ_FLOW_RFV(IVERSION,IUNIT,IFMT1,IUCLOSE,IGDMAX,
     &  NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,INSO1,INSO4,INSOKE,
     &  INSOFL,IZON,IZT,JZT,KZT,DEN,U,V,W,P,TM,DK,DE,FM,VOF,AP)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      IU=40
      IUCLOSE=1
      igdmax = 0
      I=IPROC
ccc   iproc = i
ccc   print*,'$$$ iproc=',i
c-gc  np = igdmax
ccc   print *, 'np =', np
c
c-gc  izon = lizon(i)
      IZONL = lizon(i)
      do j=1,lizon(i)
        iz=lzone(j,i)
        igdmax = igdmax+izt(iz)*jzt(iz)*kzt(iz)
      enddo
c
c-gc  call check_izon_max(izon,mblk,info)
      call check_izon_max(IZONL,mblk,info)
      call check_igdmax_max(igdmax,iiqmax,info)
c
      call inumber_cnumber(iproc,cnumber) ! process number
c------fort.13 formatted and unfrmatted
c
      fname(1:7)='f'//cnumber(6:8)//'.13'
ccc   print*,'fname=',fname(1:7)
      IF(IFMT2.EQ.1)THEN ! unformatted file
        open(IU,file=fname(1:7),form='unformatted')
        WRITE(IU) INSO1,INSO4,INSOKE,INSOFL,NGAS
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          IF(II .GT. 1 .AND. L .LE. 1) 
     &      WRITE(6,*)'@@ Error in Distribute_Flow: IZ, IZS =',II,L
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU) (DEN(L+IJK), IJK = 1,LMN)
          WRITE(IU) (  U(L+IJK), IJK = 1,LMN)
          WRITE(IU) (  V(L+IJK), IJK = 1,LMN)
          WRITE(IU) (  W(L+IJK), IJK = 1,LMN)
          WRITE(IU) (  P(L+IJK), IJK = 1,LMN)
          IF(INSO4.EQ.1) WRITE(IU) ( TM(L+IJK), IJK = 1,LMN)
          IF(INSOKE.EQ.1) THEN
            WRITE(IU) ( DK(L+IJK), IJK = 1,LMN)
            WRITE(IU) ( DE(L+IJK), IJK = 1,LMN)
          ENDIF
          WRITE(IU) ( AP(L+IJK), IJK = 1,LMN)
          IF(INSOFL.GE.1) WRITE(IU) ( VOF(L+IJK), IJK = 1,LMN)
          IF(NGAS.GT.0) THEN
            DO KK=1,NGAS
              DO IJK = L+1,L+LMN
                FM(IJK,KK) = AMAX1(1.E-20,AMIN1(1.0,FM(IJK,KK)))
              ENDDO
              WRITE(IU) ( FM(L+IJK,KK), IJK = 1,LMN)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF(IFMT2.EQ.2) THEN ! formatted file
        open(IU,file=fname(1:7),form='formatted')
        WRITE(IU,1000) INSO1,INSO4,INSOKE,INSOFL,NGAS
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          IF(II .GT. 1 .AND. L .LE. 1) 
     &      WRITE(6,*)'@@ Error in Distribute_Flow: IZ, IZS =',II,L
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU,1001) (DEN(L+IJK), IJK = 1,LMN)
          WRITE(IU,1001) (  U(L+IJK), IJK = 1,LMN)
          WRITE(IU,1001) (  V(L+IJK), IJK = 1,LMN)
          WRITE(IU,1001) (  W(L+IJK), IJK = 1,LMN)
          WRITE(IU,1001) (  P(L+IJK), IJK = 1,LMN)
          IF(INSO4.EQ.1) WRITE(IU,1001) ( TM(L+IJK), IJK = 1,LMN)
          IF(INSOKE.EQ.1) THEN
            WRITE(IU,1001) ( DK(L+IJK), IJK = 1,LMN)
            WRITE(IU,1001) ( DE(L+IJK), IJK = 1,LMN)
          ENDIF
          WRITE(IU,1001) ( AP(L+IJK), IJK = 1,LMN)
          IF(INSOFL.GE.1) WRITE(IU,1001) ( VOF(L+IJK), IJK = 1,LMN)
          IF(NGAS.GT.0) THEN
            DO KK=1,NGAS
              DO IJK = L+1,L+LMN
                FM(IJK,KK) = AMAX1(1.E-20,AMIN1(1.0,FM(IJK,KK)))
              ENDDO
              WRITE(IU,1001) ( FM(L+IJK,KK), IJK = 1,LMN)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
C
      REWIND(IU)
      IF(IUCLOSE.EQ.1) CLOSE(IU)
c
C-GC  CALL WRITE_FLOW_RFV(IVERSION,IUNIT,IFMT2,IUCLOSE,IGDMAX,
C-GC &  NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,INSO1,INSO4,
C-GC &  INSOKE,INSOFL,IZON,IZT(I),JZT(I),KZT(I),DEN,U,V,W,P,
C-GC &  TM,DK,DE,FM,VOF,AP,NP)
c
      IF(MYID .EQ. 0)
     &  WRITE(6,*) '*** END OF OUTPUT FDNSRFV FLOW FILE, IPROC =',I
 1000 FORMAT(8I5)
 1001 FORMAT(5(1P,E16.8))
c
      return
      end
c
      subroutine collect_grid(icont,ifmt1,ifmt2,lproc,
     &  lizon,lzone,mblk,izon,izs,izt,jzt,kzt,iiqmax,x,y,z,work)
c
      DIMENSION X(*),Y(*),Z(*),WORK(*)
      dimension izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMMON /COMMPI/NPROC,IPROC,MYID
      DIMENSION LIZT(1000),LJZT(1000),LKZT(1000)
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c
      IF(MYID .EQ. 0) print*,'iproc=',iproc,' ifmt =',ifmt1
c
      call inumber_cnumber(iproc,cnumber) ! process number
c-----------fort.12 formatted and unfrmatted
c-gc  ii=index(fname,' ')-1
      fname(1:7)='f'//cnumber(6:8)//'.22'
      print*,'fname=',fname(1:7)
      IUNIT=30
      IUCLOSE=1
      if(IFMT1.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
      if(IFMT1.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
C
      CALL READ_GRID_RFV(IUNIT,IFMT1,IUCLOSE,IGDMAX,IIQMAX,IZONL,
     &  LIZT,LJZT,LKZT,X,Y,Z)
C
      IF(MYID .EQ. 0) THEN
        DO IJK=1,IGDMAX
          WORK(IJK)=X(IJK)
        ENDDO
        call recv_and_swap(x,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=Y(IJK)
        ENDDO
        call recv_and_swap(y,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=Z(IJK)
        ENDDO
        call recv_and_swap(z,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
C
        DO I=2,NPROC
        DO J=1,3
          CALL MPI_RECV(work,IIQMAX,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,
     &              MPI_COMM_WORLD,status,ierr)
          IP=STATUS(MPI_SOURCE)+1
          IC=STATUS(MPI_TAG)
          print*,'Receiving package from IPROC =',ip,', IC =',ic
C
          IF(IC .EQ. 1) THEN
            call recv_and_swap(x,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 2) THEN
            call recv_and_swap(y,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 3) THEN
            call recv_and_swap(z,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE
            print *,'@@@ Error: Wrong tag no. =',IC
            print *,'           check SUBROUTINE COLLECT_GRID'
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF
C
        ENDDO
        ENDDO
C
        IUNIT=30
        IUCLOSE=1
        if(ifmt2.eq.1) open(iunit,file='fort.22',form='unformatted')
        if(ifmt2.eq.2) open(iunit,file='fort.22',form=  'formatted')
        CALL WRITE_GRID_RFV(IUNIT,IFMT2,IUCLOSE,IGDMAX,IIQMAX,IZON,
     &       IZT,JZT,KZT,X,Y,Z)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSE
        print*,'IPROC =',IPROC,', start send in FDNS-RFV GRID ...'
        CALL MPI_SEND(x,IIQMAX,MPI_REAL,0,1,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(y,IIQMAX,MPI_REAL,0,2,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(z,IIQMAX,MPI_REAL,0,3,MPI_COMM_WORLD,ierr)
      ENDIF
c
      return
      end
C
      subroutine collect_flow(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,mext,nspm,izon,izs,izt,jzt,kzt,
     &  iiqmax,den,u,v,w,p,tm,dk,de,fm,work,vof,ap)
c
      DIMENSION DEN(*),U(*),V(*),W(*),P(*),TM(*),DK(*),
     &  DE(*),FM(IIQMAX,*),VOF(*),AP(*),WORK(*)
      dimension izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMMON /COMMPI/NPROC,IPROC,MYID
      DIMENSION LIZT(1000),LJZT(1000),LKZT(1000)
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      IF(MYID .EQ. 0) print*,'iproc=',iproc,' ifmt =',ifmt1
c
      call inumber_cnumber(iproc,cnumber) ! process number
c-----------fort.13 formatted and unfrmatted
      IUNIT=40
      IUCLOSE=1
      fname(1:7)='f'//cnumber(6:8)//'.23'
      print*,'fname=',fname(1:7)
      if(IFMT1.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
      if(IFMT1.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
c
      IZONL=LIZON(IPROC)
      DO J=1,LIZON(IPROC)
        IZ = LZONE(J,IPROC)
        LIZT(J)=IZT(IZ)
        LJZT(J)=JZT(IZ)
        LKZT(J)=KZT(IZ)
      ENDDO
      CALL READ_FLOW_RFV(IVERSION,IUNIT,IFMT1,IUCLOSE,IGDMAX,
     &  NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,INSO1,INSO4,INSOKE,
     &  INSOFL,IZONL,LIZT,LJZT,LKZT,DEN,U,V,W,P,TM,DK,DE,FM,VOF,AP)
c
      IF(MYID .EQ. 0) THEN
        DO IJK=1,IGDMAX
          WORK(IJK)=DEN(IJK)
        ENDDO
        call recv_and_swap(den,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=U(IJK)
        ENDDO
        call recv_and_swap(u,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=V(IJK)
        ENDDO
        call recv_and_swap(v,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=W(IJK)
        ENDDO
        call recv_and_swap(w,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=P(IJK)
        ENDDO
        call recv_and_swap(p,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        IF(INSO4 .GE. 1) THEN
          DO IJK=1,IGDMAX
            WORK(IJK)=TM(IJK)
          ENDDO
          call recv_and_swap(tm,work,lizon(1),lzone(1,1),izs,izt,jzt,
     &                       kzt)
        ENDIF
        IF(INSOKE .GE. 1) THEN
          DO IJK=1,IGDMAX
            WORK(IJK)=DK(IJK)
          ENDDO
          call recv_and_swap(DK,work,lizon(1),lzone(1,1),izs,izt,jzt,
     &                       kzt)
          DO IJK=1,IGDMAX
            WORK(IJK)=DE(IJK)
          ENDDO
          call recv_and_swap(DE,work,lizon(1),lzone(1,1),izs,izt,jzt,
     &                       kzt)
        ENDIF
        DO IJK=1,IGDMAX
          WORK(IJK)=AP(IJK)
        ENDDO
        call recv_and_swap(ap,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        IF(INSOFL .GE. 1) THEN
          DO IJK=1,IGDMAX
            WORK(IJK)=VOF(IJK)
          ENDDO
          call recv_and_swap(vof,work,lizon(1),lzone(1,1),izs,izt,jzt,
     &                       kzt)
        ENDIF
        DO K=1,NGAS
          DO IJK=1,IGDMAX
            WORK(IJK)=FM(IJK,K)
          ENDDO
          call recv_and_swap(fm(1,k),work,lizon(1),lzone(1,1),
     &                       izs,izt,jzt,kzt)
        ENDDO
C
C       ICT=5+NGAS
        ICT=6+NGAS
        IF(INSO4 .GE. 1) ICT=ICT+1
        IF(INSOKE .GE. 1) ICT=ICT+2
        IF(INSOFL .GE. 1) ICT=ICT+1
        DO I=2,NPROC
        DO J=1,ICT
C
          CALL MPI_RECV(work,IIQMAX,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,
     &                  MPI_COMM_WORLD,status,ierr)
          IP=STATUS(MPI_SOURCE)+1
          IC=STATUS(MPI_TAG)
          print*,'Receiving package from IPROC =',ip,', IC =',ic
          IF(IC .EQ. 1) THEN
            call recv_and_swap(DEN,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 2) THEN
            call recv_and_swap(U,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 3) THEN
            call recv_and_swap(V,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 4) THEN
            call recv_and_swap(W,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 5) THEN
            call recv_and_swap(P,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 6) THEN
            call recv_and_swap(TM,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 7) THEN
            call recv_and_swap(DK,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 8) THEN
            call recv_and_swap(DE,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 9) THEN
            call recv_and_swap(AP,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 10) THEN
            call recv_and_swap(VOF,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
C-GC      ELSE IF(IC .GT. 10 .AND. IC .LE. ICT) THEN
          ELSE IF(IC .GT. 10) THEN
            K = IC-10
            call recv_and_swap(FM(1,K),work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE
            print *,'@@@ Error: Wrong tag no. =',IC
            print *,'           check SUBROUTINE COLLECT_FLOW'
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF
        ENDDO
        ENDDO
C
C-----  OUTPUT Collected flowfield data
        iunit=30
        iuclose=1
        if(ifmt2.eq.1) open(iunit  ,file='fort.23',form='unformatted')
        if(ifmt2.eq.2) open(iunit  ,file='fort.23',form=  'formatted')
ccc     if(ifmt2.eq.1) open(iunit+1,file='fort.24',form='unformatted')
ccc     if(ifmt2.eq.2) open(iunit+1,file='fort.24',form=  'formatted')
        CALL WRITE_FLOW_RFV(IVERSION,IUNIT,IFMT2,IUCLOSE,IGDMAX,
     &    NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,INSO1,INSO4,INSOKE,
     &    INSOFL,IZON,IZT,JZT,KZT,DEN,U,V,W,P,TM,DK,DE,FM,VOF,AP)
C
      ELSE
        print*,'IPROC =',IPROC,', start send in FDNS-RFV FLOW DATA ...'
        CALL MPI_SEND(den,IIQMAX,MPI_REAL,0,1,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(u,IIQMAX,MPI_REAL,0,2,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(v,IIQMAX,MPI_REAL,0,3,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(w,IIQMAX,MPI_REAL,0,4,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(p,IIQMAX,MPI_REAL,0,5,MPI_COMM_WORLD,ierr)
        IF(INSO4 .GE. 1)
     &    CALL MPI_SEND(tm,IIQMAX,MPI_REAL,0,6,MPI_COMM_WORLD,ierr)
        IF(INSOKE .GE. 1) THEN
          CALL MPI_SEND(dk,IIQMAX,MPI_REAL,0,7,MPI_COMM_WORLD,ierr)
          CALL MPI_SEND(de,IIQMAX,MPI_REAL,0,8,MPI_COMM_WORLD,ierr)
        ENDIF
        CALL MPI_SEND(ap,IIQMAX,MPI_REAL,0,9,MPI_COMM_WORLD,ierr)
        IF(INSOFL .GE. 1)
     &    CALL MPI_SEND(vof,IIQMAX,MPI_REAL,0,10,MPI_COMM_WORLD,ierr)
        do k=1,ngas
          IC=10+K
          do ijk=1,igdmax
            work(ijk)=fm(ijk,k)
          enddo
          CALL MPI_SEND(work,IIQMAX,MPI_REAL,0,IC,MPI_COMM_WORLD,ierr)
        enddo
      ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      return
      end
c
c**********************************************************************
      subroutine distribute_plot3dg(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
c**********************************************************************
      DIMENSION X(*),Y(*),Z(*),DEN(*),U(*),V(*),W(*),P(*),TM(*),DK(*),
     &  DE(*),FM(IIQMAX,*),VOFS(IIQMAX,*),WORK(*)
      dimension izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      COMMON /COMMPI/NPROC,IPROC,MYID
      DIMENSION LIZT(1000),LJZT(1000),LKZT(1000)
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----INPUT RESTART FILE
      iunit=30
      iuclose=1
      if(ifmt1.eq.1) open(iunit,file='fort.12',form='unformatted')
      if(ifmt1.eq.2) open(iunit,file='fort.12',form=  'formatted')
      CALL READ_PLOT3DG(IUNIT,IFMT1,IUCLOSE,IGDMAX,IIQMAX,IZON,
     &     IZT,JZT,KZT,X,Y,Z)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      igdmax = 0
c-gc  iproc = i
c-gc  np = igdmax
c-gc  print *,'np =', np
c-gc  print*,'$$$ iproc=',i
c-gc  izon = lizon(i)
C-GC1
      I=IPROC
      IZONL = lizon(i)
C-GC2
      do j=1,lizon(i)
        iz=lzone(j,i)
        igdmax=igdmax+izt(iz)*jzt(iz)*kzt(iz)
C-GC1
        LIZT(J)=IZT(IZ)
        LJZT(J)=JZT(IZ)
        LKZT(J)=KZT(IZ)
C-GC2
      enddo
c     print*,'iproc=',iproc,' ifmt  =',ifmt  
c-gc  call check_izon_max(izon,mblk,info)
      call check_izon_max(IZONL,mblk,info)
      call check_igdmax_max(igdmax,iiqmax,info)
c
      call inumber_cnumber(iproc,cnumber) ! process number
c-----------fort.12 formatted and unfrmatted
      IU=40
      IUCLOSE=1
c
      fname(1:7)='f'//cnumber(6:8)//'.12'
      print*,'fname=',fname(1:7)
C-GC1
      IF(IFMT2.eq.1) THEN
        open(IUNIT,file=fname(1:7),form='unformatted')
        WRITE(IU) IZONL
        IF(MYID .EQ. 0)
     &    WRITE(6,'(" PLOT3D GRID DATA FILE, IZON=",I5)') IZONL
        WRITE(IU) (LIZT(II),LJZT(II),LKZT(II),II=1,IZONL)
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU) (X(L+IJK), IJK = 1,LMN)
     &             ,(Y(L+IJK), IJK = 1,LMN)
     &             ,(Z(L+IJK), IJK = 1,LMN)
        ENDDO
      ENDIF
      IF(IFMT2.eq.2) THEN
        open(IUNIT,file=fname(1:7),form='formatted')
        WRITE(IU,*) IZONL
        IF(MYID .EQ. 0)
     &    WRITE(6,'(" PLOT3D GRID DATA FILE, IZON=",I5)') IZONL
        WRITE(IU,*) (LIZT(II),LJZT(II),LKZT(II),II=1,IZONL)
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU,'(1P5E16.8)') (X(L+IJK), IJK = 1,LMN)
     &                          ,(Y(L+IJK), IJK = 1,LMN)
     &                          ,(Z(L+IJK), IJK = 1,LMN)
        ENDDO
      ENDIF
      REWIND(IU)
      IF(IUCLOSE.EQ.1) CLOSE(IU)
C-GC2
c
c-gc  if(IFMT2.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
c-gc  if(IFMT2.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
c-gc  CALL WRITE_PLOT3DG(IUNIT,IFMT2,IUCLOSE,IGDMAX,IIQMAX,IZON,
c-gc &     IZT(i),JZT(i),KZT(i),X,Y,Z,np)
c
      return
      end
c
c***********************************************************************
      subroutine distribute_restart(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,mext,aext,iext,
     &  izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
c***********************************************************************
      DIMENSION X(*),Y(*),Z(*),DEN(*),U(*),V(*),W(*),P(*),TM(*),DK(*),
     &  DE(*),FM(IIQMAX,*),VOFS(IIQMAX,*),WORK(*)
      dimension aext(*),iext(*),izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      COMMON /COMMPI/NPROC,IPROC,MYID
      DIMENSION LIZT(1000),LJZT(1000),LKZT(1000)
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-----INPUT RESTART FILE
      iunit=30
      iuclose=1
      if(ifmt1.eq.1) open(iunit  ,file='fort.13',form='unformatted')
      if(ifmt1.eq.2) open(iunit  ,file='fort.13',form=  'formatted')
cc    if(ifmt1.eq.1) open(iunit+1,file='fort.14',form='unformatted')
cc    if(ifmt1.eq.2) open(iunit+1,file='fort.14',form=  'formatted')
      CALL READ_RESTART_0402(IVERSION,IUNIT,IFMT1,IUCLOSE,IGDMAX,
     &  NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,XMACH,REYNL,ALPHX,
     &  MEXT,AEXT,IEXT,IZON,IZT,JZT,KZT,DEN,U,V,W,P,TM,DK,DE,FM,VOFS)
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      igdmax = 0
c-gc  iproc = i
c-gc  print*,'$$$ iproc=',i
c-gc  np = igdmax
c-gc  print *, 'np =', np
c
c-gc  izon = lizon(i)
C-GC1
      I = IPROC
      IZONL = lizon(i)
C-GC2
      do j=1,lizon(i)
        iz=lzone(j,i)
        igdmax = igdmax+izt(iz)*jzt(iz)*kzt(iz)
C-GC1
        LIZT(J)=IZT(IZ)
        LJZT(J)=JZT(IZ)
        LKZT(J)=KZT(IZ)
C-GC2
      enddo
c
c-gc  print*,'finished sending process'
c-gc  call check_izon_max(izon,mblk,info)
      call check_izon_max(IZONL,mblk,info)
      call check_igdmax_max(igdmax,iiqmax,info)
c 
      call inumber_cnumber(iproc,cnumber) ! process number
c------------fort.13 formatted and unfrmatted
      IU=40
      IUCLOSE=1
c
      fname(1:7)='f'//cnumber(6:8)//'.13'
      print*,'fname=',fname(1:7)
C-GC1
      IF(IFMT2.eq.1) THEN
        open(IUNIT,file=fname(1:7),form='unformatted')
        WRITE(IU) IZONL,IVERSION,ITPRE,TIMET,NGAS,ISPVOF,
     &    (AEXT(I),IEXT(I),I=1,MEXT)
        WRITE(IU) (LIZT(II),LJZT(II),LKZT(II),II=1,IZONL)
        DO J=1,IZONL ! THIS PART CAN BE READED BY PLOT3D
          II=LZONE(J,I)
          L = IZS(II)
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU) XMACH,REYNL,ALPHX,TIMET
          WRITE(IU) (DEN(L+IJK), IJK = 1,LMN)
     &             ,(  U(L+IJK), IJK = 1,LMN)
     &             ,(  V(L+IJK), IJK = 1,LMN)
     &             ,(  W(L+IJK), IJK = 1,LMN)
     &             ,(  P(L+IJK), IJK = 1,LMN)
        ENDDO
C
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU) ( TM(L+IJK), IJK = 1,LMN)
     &             ,( DK(L+IJK), IJK = 1,LMN)
     &             ,( DE(L+IJK), IJK = 1,LMN)
          IF(NGAS.GT.0) THEN ! MASS FRACTION
            DO KK=1,NGAS
            DO IJK = L+1,L+LMN
              FM(IJK,KK) = AMAX1(1.0E-20,AMIN1(1.0,FM(IJK,KK)))
            ENDDO
            ENDDO
            WRITE(IU) ((FM(L+IJK,KK),IJK=1,LMN),KK=1,NGAS)
          ENDIF
        ENDDO
C
        IF(ISPVOF.GE.1) THEN ! VOF FUNCTION
          DO J=1,IZONL
            II=LZONE(J,I)
            L = IZS(II)
            LMN = IZT(II)*JZT(II)*KZT(II)
            WRITE(IU) ((VOFS(L+IJK,KK),IJK=1,LMN),KK=1,ISPVOF)
          ENDDO
        ENDIF
      ENDIF
      IF(IFMT2.EQ.2) THEN
        open(IUNIT,file=fname(1:7),form='formatted')
        WRITE(IU,*) IZONL,IVERSION,ITPRE,TIMET,NGAS,ISPVOF,
     &    (AEXT(I),IEXT(I),I=1,MEXT)
        WRITE(IU,*) (LIZT(II),LJZT(II),LKZT(II),II=1,IZONL)
        DO J=1,IZONL ! THIS PART CAN BE READED BY PLOT3D
          II=LZONE(J,I)
          L = IZS(II)
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU,*) XMACH,REYNL,ALPHX,TIMET
          WRITE(IU,*) (DEN(L+IJK), IJK = 1,LMN)
     &               ,(  U(L+IJK), IJK = 1,LMN)
     &               ,(  V(L+IJK), IJK = 1,LMN)
     &               ,(  W(L+IJK), IJK = 1,LMN)
     &               ,(  P(L+IJK), IJK = 1,LMN)
        ENDDO
C
        DO J=1,IZONL
          II=LZONE(J,I)
          L = IZS(II)
          LMN = IZT(II)*JZT(II)*KZT(II)
          WRITE(IU,*) ( TM(L+IJK), IJK = 1,LMN)
     &               ,( DK(L+IJK), IJK = 1,LMN)
     &               ,( DE(L+IJK), IJK = 1,LMN)
          IF(NGAS.GT.0) THEN ! MASS FRACTION
            DO KK=1,NGAS
            DO IJK = L+1,L+LMN
              FM(IJK,KK) = AMAX1(1.0E-20,AMIN1(1.0,FM(IJK,KK)))
            ENDDO
            ENDDO
            WRITE(IU,*) ((FM(L+IJK,KK),IJK=1,LMN),KK=1,NGAS)
          ENDIF
        ENDDO
C
        IF(ISPVOF.GE.1) THEN ! VOF FUNCTION
          DO J=1,IZONL
            II=LZONE(J,I)
            L = IZS(II)
            LMN = IZT(II)*JZT(II)*KZT(II)
            WRITE(IU,*) ((VOFS(L+IJK,KK),IJK=1,LMN),KK=1,ISPVOF)
          ENDDO
        ENDIF
      ENDIF
      REWIND(IU)
      IF(IUCLOSE.EQ.1) CLOSE(IU)
C-GC2
c 
c-gc  if(IFMT2.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
c-gc  if(IFMT2.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
c-gc  CALL WRITE_RESTART_0402(IVERSION,IUNIT,IFMT2,IUCLOSE,IGDMAX,
c-gc &  NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,XMACH,REYNL,
c-gc &  ALPHX,MEXT,AEXT,IEXT,IZON,IZT(i),JZT(i),KZT(i),DEN,U,V,W,P,
c-gc &  TM,DK,DE,FM,VOFS,np)
c
      return
      end
c
c***********************************************************************
      subroutine collect_plot3dg(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
c***********************************************************************
      DIMENSION X(*),Y(*),Z(*),DEN(*),U(*),V(*),W(*),P(*),TM(*),DK(*),
     &  DE(*),FM(IIQMAX,*),VOFS(IIQMAX,*),WORK(*)
      dimension izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMMON /COMMPI/NPROC,IPROC,MYID
      DIMENSION LIZT(1000),LJZT(1000),LKZT(1000)
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c-gc  mtag=10000
c-gc  do 500 i=1,nproc
c
c-gc    if(myid .eq. 0) then
c-gc      call sendit2(i,ifmt1,nbtid,mtag)
c-gc    else
c-gc      call recvit2(iproc,ifmt,nbtid,mtag)
c-gc      print*,'iproc=',iproc,' ifmt =',ifmt 
        IF(MYID .EQ. 0) print*,'iproc=',iproc,' ifmt =',ifmt1
c
          call inumber_cnumber(iproc,cnumber) ! process number
c-----------fort.12 formatted and unfrmatted
          ii=index(fname,' ')-1
          fname(1:7)='f'//cnumber(6:8)//'.22'
          print*,'fname=',fname(1:7)
          IUNIT=30
          IUCLOSE=1
C-GC1
          if(IFMT1.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
          if(IFMT1.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
          CALL READ_PLOT3DG(IUNIT,IFMT1,IUCLOSE,IGDMAX,IIQMAX,IZONL,
     &         LIZT,LJZT,LKZT,X,Y,Z)
C-GC2
c-gc      if(IFMT.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
c-gc      if(IFMT.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
c-gc      CALL READ_PLOT3DG(IUNIT,IFMT,IUCLOSE,IGDMAX,IIQMAX,IZON,
c-gc &         IZT,JZT,KZT,X,Y,Z)
c-gc      print*,'start send in plot3dg ....'
c-gc      call sendnb(x,1,igdmax,nbtid,mtag)
c-gc      call sendnb(y,1,igdmax,nbtid,mtag)
c-gc      call sendnb(z,1,igdmax,nbtid,mtag)
c-gc    endif
c-gc    if(myid .eq. 0) then
c-gc      print*,'start recv in colllection info=',info
c-gc      call recv_and_swap(x,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(y,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(z,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc    endif
c-gc 500  continue
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-GC1
      IF(MYID .EQ. 0) THEN
        DO IJK=1,IGDMAX
          WORK(IJK)=X(IJK)
        ENDDO
        call recv_and_swap(x,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=Y(IJK)
        ENDDO
        call recv_and_swap(y,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=Z(IJK)
        ENDDO
        call recv_and_swap(z,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
C
        DO I=2,NPROC
        DO J=1,3
          CALL MPI_RECV(work,IIQMAX,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,
     &              MPI_COMM_WORLD,status,ierr)
          IP=STATUS(MPI_SOURCE)+1
          IC=STATUS(MPI_TAG)
          print*,'Receiving package from IPROC =',ip,', IC =',ic
C
          IF(IC .EQ. 1) THEN
            call recv_and_swap(x,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 2) THEN
            call recv_and_swap(y,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 3) THEN
            call recv_and_swap(z,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE
            print *,'@@@ Error: Wrong tag no. =',IC
            print *,'           check SUBROUTINE COLLECT_PLOT3DG'
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF
        ENDDO
        ENDDO
C
        IUNIT=30
        IUCLOSE=1
        if(ifmt2.eq.1) open(iunit,file='fort.22',form='unformatted')
        if(ifmt2.eq.2) open(iunit,file='fort.22',form=  'formatted')
        CALL WRITE_PLOT3DG(IUNIT,IFMT2,IUCLOSE,IGDMAX,IIQMAX,IZON,
     &       IZT,JZT,KZT,X,Y,Z)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      ELSE
        print*,'IPROC =',IPROC,', start send in FDNS GRID ...'
        CALL MPI_SEND(x,IIQMAX,MPI_REAL,0,1,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(y,IIQMAX,MPI_REAL,0,2,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(z,IIQMAX,MPI_REAL,0,3,MPI_COMM_WORLD,ierr)
      ENDIF
c-gc  IUNIT=30
c-gc  IUCLOSE=1
c-gc  if(ifmt2.eq.1) open(iunit,file='fort.22',form='unformatted')
c-gc  if(ifmt2.eq.2) open(iunit,file='fort.22',form=  'formatted')
c-gc  CALL WRITE_PLOT3DG(IUNIT,IFMT2,IUCLOSE,IGDMAX,IIQMAX,IZON,
c-gc &     IZT,JZT,KZT,X,Y,Z)
c
      return
      end
C
c***********************************************************************
      subroutine collect_restart(icont,ifmt1,ifmt2,
     &  lproc,lizon,lzone,mblk,mext,aext,iext,
     &  izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
c***********************************************************************
      DIMENSION X(*),Y(*),Z(*),DEN(*),U(*),V(*),W(*),P(*),TM(*),DK(*),
     &  DE(*),FM(IIQMAX,*),VOFS(IIQMAX,*),WORK(*)
      dimension aext(*),iext(*),izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80,cnumber*8
      include 'mpif.h'
C-GC1
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMMON /COMMPI/NPROC,IPROC,MYID
      DIMENSION LIZT(1000),LJZT(1000),LKZT(1000)
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c-gc  mtag=10000
c-gc  do 500 i=1,nproc
c-gc    if(myid .eq. 0) then
c-gc      call sendit2(i,ifmt1,nbtid,mtag)
c-gc    else
c-gc      call recvit2(iproc,ifmt,nbtid,mtag)
c-gc      print*,'iproc=',iproc,' ifmt =',ifmt 
      IF(MYID .EQ. 0) print*,'iproc=',iproc,' ifmt =',ifmt1
c-----------fort.23 formatted and unfrmatted
      call inumber_cnumber(iproc,cnumber) ! process number
      IUNIT=40
      IUCLOSE=1
      fname(1:7)='f'//cnumber(6:8)//'.23'
      print*,'fname=',fname(1:7)
C-GC1
      if(IFMT1.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
      if(IFMT1.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
c
      CALL READ_RESTART_0402(IVERSION,IUNIT,IFMT,IUCLOSE,IGDMAX,
     &     NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,XMACH,REYNL,
     &     ALPHX,MEXT,AEXT,IEXT,IZONL,LIZT,LJZT,LKZT,DEN,U,V,W,P,TM,
     &     DK,DE,FM,VOFS)
      ICS=8
      ICG=ICS+NGAS
      ICT=ICG+ISPVOF
C-GC2
c-gc      if(IFMT.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
c-gc      if(IFMT.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
c
c-gc      CALL READ_RESTART_0402(IVERSION,IUNIT,IFMT,IUCLOSE,IGDMAX,
c-gc &         NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,XMACH,REYNL,
c-gc &         ALPHX,MEXT,AEXT,IEXT,IZON,IZT,JZT,KZT,DEN,U,V,W,P,TM,
c-gc &         DK,DE,FM,VOFS)
c
c-gc      print*,'start send in read_restart ....'
c-gc      call sendit4(ngas,ispvof,IVERSION,ITPRE,nbtid,mtag)
c-gc      call sendnb4(TIMET,XMACH,REYNL,ALPHX,nbtid,mtag)
c
c-gc      call sendnb(aext,1,mext,nbtid,mtag)
c-gc      call sendit(iext,1,mext,nbtid,mtag)
c
c-gc      call sendnb(den ,1,igdmax,nbtid,mtag)
c-gc      call sendnb(u   ,1,igdmax,nbtid,mtag)
c-gc      call sendnb(v   ,1,igdmax,nbtid,mtag)
c-gc      call sendnb(w   ,1,igdmax,nbtid,mtag)
c-gc      call sendnb(p   ,1,igdmax,nbtid,mtag)
c-gc      call sendnb(tm  ,1,igdmax,nbtid,mtag)
c-gc      call sendnb(dk  ,1,igdmax,nbtid,mtag)
c-gc      call sendnb(de  ,1,igdmax,nbtid,mtag)
c-gc      do k=1,ngas
c-gc        do ijk=1,igdmax
c-gc          work(ijk)=fm(ijk,k)
c-gc        enddo
c-gc        call sendnb(work,1,igdmax,nbtid,mtag)
c-gc      enddo
c-gc      do k=1,ispvof
c-gc        do ijk=1,igdmax
c-gc          work(ijk)=vofs(ijk,k)
c-gc        enddo
c-gc        call sendnb(work,1,igdmax,nbtid,mtag)
c-gc      enddo
c-gc    endif
c    
c-gc    if(myid .eq. 0) then
c-gc      print*,'start recv in colllection. info=',info
c-gc      call recvit4(ngas,ispvof,IVERSION,ITPRE,nbtid,mtag)
c-gc      call recvnb4(TIMET,XMACH,REYNL,ALPHX,nbtid,mtag)
c
c-gc      call recvnb(aext,1,mext1,nbtid,mtag)
c-gc      call recvit(iext,1,mext1,nbtid,mtag)
c
c-gc      call recv_and_swap(den ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(u   ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(v   ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(w   ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(p   ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(tm  ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(dk  ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      call recv_and_swap(de  ,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc      do k=1,ngas
c-gc        call recv_and_swap(fm(1,k),work,lizon(i),lzone(1,i),
c-gc &           izs,izt,jzt,kzt,nbtid,mtag)
c-gc      enddo
c-gc      do k=1,ispvof
c-gc        call recv_and_swap(vofs(1,k),work,lizon(i),lzone(1,i),
c-gc &           izs,izt,jzt,kzt,nbtid,mtag)
c-gc      enddo
c-gc    endif
c-gc 500  continue
C-GC1
      IF(MYID .EQ. 0) THEN
        DO IJK=1,IGDMAX
          WORK(IJK)=DEN(IJK)
        ENDDO
        call recv_and_swap(den,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=U(IJK)
        ENDDO
        call recv_and_swap(u,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=V(IJK)
        ENDDO
        call recv_and_swap(v,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=W(IJK)
        ENDDO
        call recv_and_swap(w,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=P(IJK)
        ENDDO
        call recv_and_swap(p,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=TM(IJK)
        ENDDO
        call recv_and_swap(tm,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=DK(IJK)
        ENDDO
        call recv_and_swap(DK,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=DE(IJK)
        ENDDO
        call recv_and_swap(DE,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO K=1,NGAS
          DO IJK=1,IGDMAX
            WORK(IJK)=FM(IJK,K)
          ENDDO
          call recv_and_swap(fm(1,k),work,lizon(1),lzone(1,1),
     &                       izs,izt,jzt,kzt)
        ENDDO
        DO K=1,ISPVOF
          DO IJK=1,IGDMAX
            WORK(IJK)=VOFS(IJK,KK)
          ENDDO
          call recv_and_swap(vofs(1,k),work,lizon(1),lzone(1,1),
     &                       izs,izt,jzt,kzt)
        ENDDO
C
        DO I=2,NPROC
        DO J=1,ICT
C
          CALL MPI_RECV(work,IIQMAX,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,
     &                  MPI_COMM_WORLD,status,ierr)
          IP=STATUS(MPI_SOURCE)+1
          IC=STATUS(MPI_TAG)
          print*,'Receiving package from IPROC =',ip,', IC =',ic
          IF(IC .EQ. 1) THEN
            call recv_and_swap(DEN,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 2) THEN
            call recv_and_swap(U,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 3) THEN
            call recv_and_swap(V,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 4) THEN
            call recv_and_swap(W,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 5) THEN
            call recv_and_swap(P,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 6) THEN
            call recv_and_swap(TM,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 7) THEN
            call recv_and_swap(DK,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 8) THEN
            call recv_and_swap(DE,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .GT. 8 .AND. IC .LE. ICG) THEN
            K = IC-9
            call recv_and_swap(FM(1,K),work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .GT. ICG .AND. IC .LE. ICT) THEN
            K = IC-ICG
            call recv_and_swap(VOFS(1,k),work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE
            print *,'@@@ Error: Wrong tag no. =',IC
            print *,'           check SUBROUTINE COLLECT_RESTART'
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF
        ENDDO
        ENDDO
C
C-----  OUTPUT Collected flowfield data
        iunit=30
        iuclose=1
        if(ifmt2.eq.1) open(iunit  ,file='fort.23',form='unformatted')
        if(ifmt2.eq.2) open(iunit  ,file='fort.23',form=  'formatted')
ccc     if(ifmt2.eq.1) open(iunit+1,file='fort.24',form='unformatted')
ccc     if(ifmt2.eq.2) open(iunit+1,file='fort.24',form=  'formatted')
        CALL WRITE_RESTART_0402(IVERSION,IUNIT,IFMT2,IUCLOSE,IGDMAX,
     &    NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,XMACH,REYNL,ALPHX,
     &    MEXT,AEXT,IEXT,IZON,IZT,JZT,KZT,DEN,U,V,W,P,TM,DK,DE,FM,VOFS)
C
      ELSE
        print*,'IPROC =',IPROC,', start sending in FDNS FLOW DATA ...'
        CALL MPI_SEND(den,IIQMAX,MPI_REAL,0,1,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(u,IIQMAX,MPI_REAL,0,2,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(v,IIQMAX,MPI_REAL,0,3,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(w,IIQMAX,MPI_REAL,0,4,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(p,IIQMAX,MPI_REAL,0,5,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(tm,IIQMAX,MPI_REAL,0,6,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(dk,IIQMAX,MPI_REAL,0,7,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(de,IIQMAX,MPI_REAL,0,8,MPI_COMM_WORLD,ierr)
        do k=1,ngas
          IC=ICS+K
          do ijk=1,igdmax
            work(ijk)=fm(ijk,k)
          enddo
          CALL MPI_SEND(work,IIQMAX,MPI_REAL,0,IC,MPI_COMM_WORLD,ierr)
        enddo
        do k=1,ispvof
          IC=ICG+K
          do ijk=1,igdmax
            work(ijk)=vofs(ijk,k)
          enddo
          CALL MPI_SEND(work,IIQMAX,MPI_REAL,0,IC,MPI_COMM_WORLD,ierr)
        enddo
      ENDIF
C-GC2
c
C-----OUTPUT RESTART FILE
c-gc  iunit=30
c-gc  iuclose=1
c-gc  if(ifmt2.eq.1) open(iunit  ,file='fort.23',form='unformatted')
c-gc  if(ifmt2.eq.2) open(iunit  ,file='fort.23',form=  'formatted')
ccc   if(ifmt2.eq.1) open(iunit+1,file='fort.24',form='unformatted')
ccc   if(ifmt2.eq.2) open(iunit+1,file='fort.24',form=  'formatted')
c-gc  CALL WRITE_RESTART_0402(IVERSION,IUNIT,IFMT2,IUCLOSE,IGDMAX,
c-gc &  NGAS,ISPVOF,IIQMAX,NSPM,IVOFMX,ITPRE,TIMET,XMACH,REYNL,ALPHX,
c-gc &  MEXT,AEXT,IEXT,IZON,IZT,JZT,KZT,DEN,U,V,W,P,TM,DK,DE,FM,VOFS)
c
      return
      end
c
c***********************************************************************
      subroutine collect_plot3dq(icont,ifmt1,ifmt2,
     &  qfile,lproc,lizon,lzone,mblk,izon,izs,izt,jzt,kzt,
     &  iiqmax,x,y,z,den,u,v,w,p,tm,dk,de,fm,vofs,work)
c***********************************************************************
      DIMENSION X(*),Y(*),Z(*),DEN(*),U(*),V(*),W(*),P(*),TM(*),DK(*),
     &  DE(*),FM(IIQMAX,*),VOFS(IIQMAX,*),WORK(*)
      dimension izs(*),izt(*),jzt(*),kzt(*)
      dimension lproc(*),lizon(*),lzone(mblk,*)
      character fname*80, cnumber*8
      character qfile*2
      include 'mpif.h'
C-GC1
      INTEGER STATUS(MPI_STATUS_SIZE)
      COMMON /COMMPI/NPROC,IPROC,MYID
      DIMENSION LIZT(1000),LJZT(1000),LKZT(1000)
C-GC2
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C
c-gc  mtag=10000
c-gc  do 500 i=1,nproc
c
c-gc    if(myid .eq. 0) then
c-gc       call sendit2(i,ifmt1,nbtid,mtag)
c-gc       call send_string(qfile,  1,nbtid,mtag)
c-gc    else
c-gc       call recvit2(iproc,ifmt,nbtid,mtag)
c-gc       print*,'iproc=',iproc,' ifmt =',ifmt 
c-gc       call recv_string(qfile,iaux,nbtid,mtag)
ccc        print*,'fname=',fname(1:ii)
c
c-----------fort.9* formatted and unfrmatted
c
       call inumber_cnumber(iproc,cnumber) ! process number
       fname(1:7)='f'//cnumber(6:8)//'.'//qfile
       print*,'fname=',fname(1:7)
       IUNIT=30
       IUCLOSE=1
C-GC1
       if(IFMT1.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
       if(IFMT1.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
       CALL READ_PLOT3DQ(IUNIT,IFMT,IUCLOSE,IGDMAX,IIQMAX,IZONL,
     &     LIZT,LJZT,LKZT,X,Y,Z,U,V,XMACH,ALPHA,RENU,TIME)
C-GC2
c-gc       if(IFMT.eq.1) open(IUNIT,file=fname(1:7),form='unformatted')
c-gc       if(IFMT.eq.2) open(IUNIT,file=fname(1:7),form='formatted')
c
c-gc       CALL READ_PLOT3DQ(IUNIT,IFMT,IUCLOSE,IGDMAX,IIQMAX,IZON,
c-gc &         IZT,JZT,KZT,X,Y,Z,U,V,XMACH,ALPHA,RENU,TIME)
c
c-gc       print*,'start send in plot3dq ....'
c-gc       call sendnb4(XMACH,ALPHA,RENU,TIME,nbtid,mtag)
c
c-gc       call sendnb(x,1,igdmax,nbtid,mtag)
c-gc       call sendnb(y,1,igdmax,nbtid,mtag)
c-gc       call sendnb(z,1,igdmax,nbtid,mtag)
c-gc       call sendnb(u,1,igdmax,nbtid,mtag)
c-gc       call sendnb(v,1,igdmax,nbtid,mtag)   
c-gc    endif    
c
c-gc    if(myid .eq. 0) then   
c-gc       print*,'start recv in colllection info=',info
c-gc       call recvnb4(XMACH,ALPHA,RENU,TIME,nbtid,mtag)
C
c-gc       call recv_and_swap(x,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc       call recv_and_swap(y,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc       call recv_and_swap(z,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc       call recv_and_swap(u,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc       call recv_and_swap(v,work,lizon(i),lzone(1,i),
c-gc &         izs,izt,jzt,kzt,nbtid,mtag)
c-gc    endif
c
c-gc 500  continue
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C-GC1
      IF(MYID .EQ. 0) THEN
        DO IJK=1,IGDMAX
          WORK(IJK)=X(IJK)
        ENDDO
        call recv_and_swap(x,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=Y(IJK)
        ENDDO
        call recv_and_swap(y,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=Z(IJK)
        ENDDO
        call recv_and_swap(z,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=U(IJK)
        ENDDO
        call recv_and_swap(u,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO IJK=1,IGDMAX
          WORK(IJK)=V(IJK)
        ENDDO
        call recv_and_swap(v,work,lizon(1),lzone(1,1),izs,izt,jzt,kzt)
        DO I=2,NPROC
        DO J=1,5
C
          CALL MPI_RECV(work,IIQMAX,MPI_REAL,MPI_ANY_SOURCE,MPI_ANY_TAG,
     &                  MPI_COMM_WORLD,status,ierr)
          IP=STATUS(MPI_SOURCE)+1
          IC=STATUS(MPI_TAG)
          print*,'Receiving package from IPROC =',ip,', IC =',ic
          IF(IC .EQ. 1) THEN
            call recv_and_swap(X,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 2) THEN
            call recv_and_swap(Y,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 3) THEN
            call recv_and_swap(Z,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 4) THEN
            call recv_and_swap(U,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE IF(IC .EQ. 5) THEN
            call recv_and_swap(V,work,lizon(ip),lzone(1,ip),
     &                         izs,izt,jzt,kzt)
          ELSE
            print *,'@@@ Error: Wrong tag no. =',IC
            print *,'           check SUBROUTINE COLLECT_PLOT3DQ'
            CALL MPI_FINALIZE(IERR)
            STOP
          ENDIF
        ENDDO
        ENDDO
C
C-----  OUTPUT Collected flowfield data
C
        IU=30
        ICLOSE=1
        if(ifmt2.eq.1) open(iu,file='fort.'//qfile,form='unformatted')
        if(ifmt2.eq.2) open(iu,file='fort.'//qfile,form=  'formatted')
c
        CALL WRITE_PLOT3DQ(iu,IFMT2,IUCLOSE,IGDMAX,IIQMAX,IZON,
     &       IZT,JZT,KZT,X,Y,Z,U,V,XMACH,ALPHA,RENU,TIME)
      ELSE
        print*,'IPROC =',IPROC,', start send in PLOT3D Q-DATA ...'
        CALL MPI_SEND(x,IIQMAX,MPI_REAL,0,1,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(y,IIQMAX,MPI_REAL,0,2,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(z,IIQMAX,MPI_REAL,0,3,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(u,IIQMAX,MPI_REAL,0,4,MPI_COMM_WORLD,ierr)
        CALL MPI_SEND(v,IIQMAX,MPI_REAL,0,5,MPI_COMM_WORLD,ierr)
      ENDIF
C-GC2
c-gc  IUNIT=30
c-gc  ICLOSE=1
c-gc  if(ifmt2.eq.1) open(iunit,file='fort.'//qfile,form='unformatted')
c-gc  if(ifmt2.eq.2) open(iunit,file='fort.'//qfile,form=  'formatted')
c
c-gc  CALL WRITE_PLOT3DQ(iunit,IFMT2,IUCLOSE,IGDMAX,IIQMAX,IZON,
c-gc &     IZT,JZT,KZT,X,Y,Z,U,V,XMACH,ALPHA,RENU,TIME)
c
      return
      end
c
c**************************************************************
      subroutine get_dirlen_from_exe(jj,exe)
c**************************************************************
      character exe*80,fslash*1,bslash*1
      fslash='h'
c-gc  fslash=1h/
ccc   bslash=1h\  ! complianed by power challange
      bslash=char(92) ! same as '\'
      ii=index(exe,' ')-1
      jj=0
      do k=1,ii
        if(exe(k:k).eq.fslash.or.exe(k:k).eq.bslash) jj=k
      enddo
      return
      end
c
c*********************************************************************
      subroutine swap_and_send(var,work,lizon,lzone,izs,izt,jzt,kzt,
     &  nbtid,mtag)
c**********************************************************************
      dimension var(*),work(*),lzone(*),izs(*),izt(*),jzt(*),kzt(*)
c
      ijk2=0
      do j=1,lizon
        iz=lzone(j)
        nxyz=izt(iz)*jzt(iz)*kzt(iz)
        do ii=1,nxyz
          ijk2=ijk2+1
          ijk1=izs(iz)+ii
          work(ijk2)=var(ijk1)
        enddo
      enddo
c-gc  call sendnb(work,1,ijk2,nbtid,mtag)
      return
      end
c
c***********************************************************************
c-gc  subroutine recv_and_swap(var,work,lizon,lzone,izs,izt,jzt,kzt,
c-gc &  nbtid,mtag)
      subroutine recv_and_swap(var,work,lizon,lzone,izs,izt,jzt,kzt)
c***********************************************************************
      dimension var(*),work(*),lzone(*),izs(*),izt(*),jzt(*),kzt(*)
c
c-gc  call recvnb(work,1,ijk2,nbtid,mtag)
      ijk2=0
      do j=1,lizon
        iz=lzone(j)
        nxyz=izt(iz)*jzt(iz)*kzt(iz)
        do ii=1,nxyz
          ijk2=ijk2+1
          ijk1=izs(iz)+ii
          var(ijk1)=work(ijk2)
        enddo
      enddo
      return
      end
c
