        program fmain
        real rinitmodel(120), residual, wp, wsh, wrayl, &
     &   ralpha, rcmt, rdc, rbdep, redep, rdepi, rbdist, redist, &
     &   rdisti, elo, ela, evp, evs
        integer nnstap, nnstash, insectp, insectsh, ndep, unknown, &
     &   nnstarayl, nndist, insectrayl
        parameter(nnsta=200, maxpts=500, maxnsub=9, maxdep=71, &
     &   maxndist=250, maxn3d=1, maxnloc=120)
        character stainfop_stname(64,nnsta),seisdatap_stname(64,nnsta),&
     &   stainfosh_stname(64,nnsta),seisdatash_stname(64,nnsta), &
     &   seisdatarayl_stname(64,maxnloc)
        real stainfop_stdist(nnsta), stainfop_staz(nnsta),&
     &   stainfop_stcp(nnsta),stainfop_stcs(nnsta),&
     &   seisdatap_btime(nnsta),seisdatap_t1(nnsta),&
     &   stainfosh_stdist(nnsta), stainfosh_staz(nnsta),&
     &   stainfosh_stcp(nnsta),stainfosh_stcs(nnsta),&
     &   seisdatash_btime(nnsta),seisdatash_t1(nnsta),&
     &   fltsnp(30),fltsdp(30),fltsnsh(30),fltsdsh(30),&
     &   stainforayl_stlo(maxnloc), stainforayl_stla(maxnloc),&
     &   seisdatarayl_btime(maxnloc),seisdatarayl_t1(maxnloc),&
     &   fltsnrayl(30),fltsdrayl(30),&
     &   CMTsolution(5),strayl_wt(3,maxnloc)
        integer nl,nlnz,ran,rand,j3
        real gridxbg,gridxnd,gridybg,gridynd,griddx,randx,randy
        real, dimension(:), allocatable :: gridx, gridy, gridz, &
     &   gridxnz, gridynz, gridznz
        real, allocatable :: seisdp_xdata(:,:), gfdp(:,:), &
     &   seisdsh_xdata(:,:), gfdsh(:,:), seisdrayl_xdata(:,:),&
     &   gfdrayl(:,:)

        integer :: burn_in, nsample, neq_min, neq_max
        character(len=20) :: param_name
        real :: param_value
        character(len=256) :: srvmode, cmdline
        integer :: ios2, isub, jj7

        allocate(seisdp_xdata(maxpts,nnsta),&
     &   gfdp(maxpts,6*maxdep*nnsta),seisdsh_xdata(maxpts,nnsta),&
     &   gfdsh(maxpts,2*maxdep*nnsta), &
     &   seisdrayl_xdata(maxpts,maxnloc*3),&
     &   gfdrayl(maxpts,14*maxdep*maxndist))

        open(unit=10, file='search_par.file', status='old')
          do
          read(10, *, iostat=ierror) param_name, param_value
          if (ierror /= 0) exit
          select case (trim(param_name))
          case ('burn_in')
              burn_in = int(param_value)
          case ('nsample')
              nsample = int(param_value)
          case ('neq_min')
              neq_min = int(param_value)
          case ('neq_max')
              neq_max = int(param_value)
          end select
        end do
        close(10)

        open (14, FILE='totalmt.dat')
        do ii=1,5
                read(14,*) CMTsolution(ii)
        end do
        close(14)

        call sub_init(rbtimep, retimep, rbtimesh, retimesh, &
     &   rbtimerayl, retimerayl, sbtime,&
     &   setime, nptsp, nptssh, nptsrayl, rdelta, instap,&
     &   instash, instarayl, insubev, &
     &   rinitmodel, stainfop_stname, stainfop_stdist, stainfop_staz,&
     &   stainfop_stcp, stainfop_stcs, stainfosh_stname, &
     &   stainfosh_stdist, stainfosh_staz, stainfosh_stcp, & 
     &   stainfosh_stcs, stainforayl_stlo, &
     &   stainforayl_stla, strayl_wt, seisdatap_stname, &
     &   seisdatap_btime, seisdatap_t1, seisdp_xdata,seisdatash_stname,&
     &   seisdatash_btime, seisdatash_t1, seisdsh_xdata,&
     &   seisdatarayl_stname, seisdatarayl_btime, seisdatarayl_t1, &
     &   seisdrayl_xdata, gfdp, gfdsh, gfdrayl, &
     &   fltsnp, fltsdp, fltsnsh, fltsdsh, fltsnrayl, &
     &   fltsdrayl, insectp, insectsh, &
     &   insectrayl, wp, wsh, wrayl, ralpha, rcmt, &
     &   rdc, elo, ela, evp, evs, rbdist, redist, rdisti, nndist, &
     &   rbdep, redep, rdepi, ndep, unknown)
         insubev = neq_min
! FFWD_SERVER=1: persistent forward server for ensemble drivers (EKI/EKS).
! Data and Green's functions are loaded ONCE by sub_init above; each request
! (any line on stdin; "QUIT" exits) re-reads Input.model and evaluates.
        call get_environment_variable("FFWD_SERVER", srvmode)
        if (trim(srvmode) == "1") then
          do
            read(*,'(A)',iostat=ios2) cmdline
            if (ios2 /= 0) exit
            if (trim(cmdline) == "QUIT") exit
! model rows arrive on stdin right after the GO line (no NFS round-trips)
            do isub = 1, insubev
              read(*,*) (rinitmodel((isub-1)*7+jj7), jj7=1,7)
            enddo
        call sub_forward(rbtimep, retimep, rbtimesh, retimesh, &
     &   rbtimerayl, retimerayl, sbtime,&
     &   setime, nptsp, nptssh, nptsrayl, rdelta, instap,&
     &   instash, instarayl, insubev, &
     &   rinitmodel, stainfop_stname, stainfop_stdist, stainfop_staz,&
     &   stainfop_stcp, stainfop_stcs, stainfosh_stname, &
     &   stainfosh_stdist, stainfosh_staz, stainfosh_stcp, &
     &   stainfosh_stcs, stainforayl_stlo, &
     &   stainforayl_stla, strayl_wt, seisdatap_stname, &
     &   seisdatap_btime, seisdatap_t1, seisdp_xdata,seisdatash_stname,&
     &   seisdatash_btime, seisdatash_t1, seisdsh_xdata,&
     &   seisdatarayl_stname, seisdatarayl_btime, seisdatarayl_t1, &
     &   seisdrayl_xdata, gfdp, gfdsh, gfdrayl, &
     &   residual, fltsnp, fltsdp, fltsnsh,fltsdsh,fltsnrayl,&
     &   fltsdrayl, insectp, insectsh, &
     &   insectrayl, wp, wsh, wrayl, CMTsolution, &
     &   ralpha, rcmt, rdc, elo, ela, evp, evs, rbdist, redist, &
     &   rdisti, nndist, rbdep, redep, rdepi, ndep)
            write(*,'(A,F14.6)') "DONE ", residual
            flush(6)
          enddo
        else
        call sub_forward(rbtimep, retimep, rbtimesh, retimesh, &
     &   rbtimerayl, retimerayl, sbtime,&
     &   setime, nptsp, nptssh, nptsrayl, rdelta, instap,&
     &   instash, instarayl, insubev, &
     &   rinitmodel, stainfop_stname, stainfop_stdist, stainfop_staz,&
     &   stainfop_stcp, stainfop_stcs, stainfosh_stname, &
     &   stainfosh_stdist, stainfosh_staz, stainfosh_stcp, &
     &   stainfosh_stcs, stainforayl_stlo, &
     &   stainforayl_stla, strayl_wt, seisdatap_stname, &
     &   seisdatap_btime, seisdatap_t1, seisdp_xdata,seisdatash_stname,&
     &   seisdatash_btime, seisdatash_t1, seisdsh_xdata,&
     &   seisdatarayl_stname, seisdatarayl_btime, seisdatarayl_t1, &
     &   seisdrayl_xdata, gfdp, gfdsh, gfdrayl, &
     &   residual, fltsnp, fltsdp, fltsnsh,fltsdsh,fltsnrayl,&
     &   fltsdrayl, insectp, insectsh, &
     &   insectrayl, wp, wsh, wrayl, CMTsolution, &
     &   ralpha, rcmt, rdc, elo, ela, evp, evs, rbdist, redist, &
     &   rdisti, nndist, rbdep, redep, rdepi, ndep)
        endif
        deallocate(seisdp_xdata, gfdp, seisdsh_xdata, gfdsh,&
     &   seisdrayl_xdata, gfdrayl)
        end program fmain
