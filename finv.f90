program fmain
	include 'mpif.h'
! Fix the number of subevents
! ######################################
! Each chain is run for 'burn_in + nsample' steps in total. The first burn-in samples are discarded as burn-in steps, only after which the sampling algorithm is assumed to have converged. To eliminate dependent samples in the ensemble solution, every thinn model visited in the second part is selected for the ensemble.
      ! Reading parameters from file search_par.file
      integer :: burn_in, nsample, neq_min, neq_max
      real :: x_min, x_max, y_min, y_max, z_min, z_max
      real :: dura_min, dura_max, vr_min, vr_max
      real :: theta_min, theta_max, cen_min, cen_max
      character(len=20) :: param_name
      real :: param_value
      real :: v1x, v1y, v2x, v2y, v3x, v3y
      real :: dot_product1, dot_product2, nm1, nm2, nm3
      real :: angle1, angle2
      real, parameter :: pi = 3.141592653589793
      real :: deg_to_rad = 0.0174533
      integer, parameter :: nsubmax = 7  ! maximum nsubevent for arrays
      integer, parameter :: thin = 1  ! Thinning of the chain 
      integer, parameter :: display = 250
      integer, parameter :: best = 50 ! Best solution "best <= nsample/thin"
      real, parameter  :: T_gap = 2.5 ! minimum gap of subevent times, using 1/2 period here
      real, parameter :: ppara1 = 2  ! proposal on change in x
      real, parameter :: ppara2 = 3  ! y
      real, parameter :: ppara3 = 1  ! z
      real, parameter :: ppara4 = 0.5  ! Dur
      real, parameter :: ppara5 = 0.1  ! Vr
      real, parameter :: ppara6 = 10.0  ! ThetaR
      real, parameter :: ppara7 = 0.5  ! cen
	integer, parameter :: disX = 200
	integer, parameter :: disY = 200
	integer, parameter :: disZ = 30
	integer, parameter :: disDura = 20
	integer, parameter :: disCen = 40


! DECLARATIONS OF VARIABLES
	real , EXTERNAL    ::    gasdev,ran3
	real log, sqrt
	integer ra,ra1,ran,rank, rank0, nbproc, ierror, tag, status(MPI_STATUS_SIZE)
	real para(nsubmax,7),para_prop(nsubmax,7),para_min(nsubmax,7),par1(nsubmax,7)
	integer :: neq_prop, i1, i2, i3, min_idx
        real :: dx, dy, dz, temp_distance, distances(nsubmax), temp_xyz(3)
	real :: fit,like,like_min,like_prop, apvr(6),apdt(6), tdiff(6)
	real :: u, post(disX,disY),posts(disX,disY)
	integer :: para1,para2,para3,para4,para5,para6,para7,birth,death,accept
	real  Pr1,Pr2,Pr3,Pr4,Pr5,PrB,PrD
	real  Ac1,Ac2,Ac3,Ac4,Ac5,AcB,AcE
	integer :: i,j,k,tes,nb,ount,sample,outsample,sampleO
	integer :: neq,ind
	real :: logprob,out
	real,allocatable :: fitall(:,:,:)
	real :: fittmp(nsubmax,8)
        integer :: rupTF
	integer :: fx,fy,fz
	integer :: histo(nsubmax),histos(nsubmax)
	character :: filename*13, filename2*12,number*4
        real tbtemp(nsubmax), tetemp(nsubmax), dct
        real rinitmodel(120), residual,res1
        integer nnstap, nnstash, insectp, insectsh, ndep, unknown
        parameter(nnsta=200, maxpts=500, maxdep=41, &
     &   maxndist=250, maxn3d=1, maxnloc=70)
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
     &   fltsnrayl(30),fltsdrayl(30),CMTsolution(5),strayl_wt(3,maxnloc)
        real, allocatable :: seisdp_xdata(:,:), gfdp(:,:), &
     &   seisdsh_xdata(:,:), gfdsh(:,:), seisdrayl_xdata(:,:),&
     &   gfdrayl(:,:)

        integer nl,nlnz,rand,j2,j3,j4
        real gridxbg,gridxnd,gridybg,gridynd,griddx,randx,randy,&
     &   pprop,porg
        real, dimension(:), allocatable :: gridx, gridy, gridz, &
     &   gridxnz, gridynz, gridznz

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
          case ('x_min')
              x_min = param_value
          case ('x_max')
              x_max = param_value
          case ('y_min')
              y_min = param_value
          case ('y_max')
              y_max = param_value
          case ('z_min')
              z_min = param_value
          case ('z_max')
              z_max = param_value
          case ('dura_min')
              dura_min = param_value
          case ('dura_max')
              dura_max = param_value
          case ('vr_min')
              vr_min = param_value
          case ('vr_max')
              vr_max = param_value
          case ('theta_min')
              theta_min = param_value
          case ('theta_max')
              theta_max = param_value
          case ('cen_min')
              cen_min = param_value
          case ('cen_max')
              cen_max = param_value
          case default
              print *, 'Warning: Unknown parameter ', trim(param_name)
          end select
        end do
        close(10)
        open (13, FILE='rank.dat')
                read(13,*) rank0
        close(13)
!       Read seismicity density file
        nl=0
        nlnz=0
        open (13, file='seisdens.dat')
        do
                read(13,*,iostat=io) tempx, tempy, tempz
                if (io/=0) exit
                nl=nl+1
                if (tempz>0) then
                        nlnz=nlnz+1
                endif
        end do
        close (13)
        allocate(gridx(nl),gridy(nl),gridz(nl))
        allocate(gridxnz(nlnz),gridynz(nlnz),gridznz(nlnz))
        jj=0
        open (12, file='seisdens.dat')
        do ii=1,nl
                read(12,*) gridx(ii),gridy(ii),gridz(ii)
                if (gridz(ii)>0) then
                        jj=jj+1
                        gridxnz(jj)=gridx(ii)
                        gridynz(jj)=gridy(ii)
                        gridznz(jj)=gridz(ii)
                endif
        end do
        close (12)

        open (11, file="seisgrids.dat")
                read(11,*) gridxbg,gridxnd,gridybg,gridynd,griddx
        close (11)
! Start Parralelization of the code. From now on, the code is run on each processor independently,
! with ra = the number of the proc.
	CALL cpu_time(t1) 
	call MPI_INIT(ierror)
 	call MPI_COMM_SIZE(MPI_COMM_WORLD, nbproc, ierror)
 	call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
 	nb=nbproc
        ra=rank+rank0
        ran=rank+rank0
! matrix used to store the soluitons.
        dct=((cen_max-cen_min)-(neq_max-1)*T_gap)/neq_max
        do i=1,neq_max
                tbtemp(i)=(i-1)*(dct+T_gap)+cen_min
                tetemp(i)=(i-1)*(dct+T_gap)+dct+cen_min
        enddo
	outsample = ANINT(DFLOAT(nsample)/DFLOAT(thin))
	allocate(fitall(outsample,neq_max,8))
        do i=1,neq_max
                do j=1,8
                        fittmp(i,j)=0
                        do k=1,outsample
                                fitall(k,i,j)=0
                        enddo
                enddo
        enddo
! Read the data and parameters
	insubev = neq_min ! just need a number here
	write(*,*) "Read all the parameters ... Start."

        open (14, FILE='totalmt.dat')
        do ii=1,5
                read(14,*) CMTsolution(ii)
        end do
        close(14)
	write(*,*) "Read all the parameters ...  End."
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
! Draw the first model randomly from the prior distribution
	j=0
	residual=999
        neq=neq_min
	do while(j<=100) !100
		j=j+1
                j2=-1
                do while(j2<0)
		do i = 1,neq
!                       Generate random initials
                        jj=0
                        j3=-1
                        do while (jj<100 .and. j3<0)
                        rand=int(ran3(ra)*nlnz)
                        randx=gridxnz(rand)+(ran3(ra)-0.5)*griddx
                        randy=gridynz(rand)+(ran3(ra)-0.5)*griddx
                        if (ran3(ra)<seisprior(randx,randy,griddx,nl,gridx,&
     &                  gridy,gridz,p1)) then
                        j3=0
                        endif
                        jj=jj+1
                        enddo
                        par1(i,1) = max(min(randx,x_max),x_min)
                        par1(i,2) = max(min(randy,y_max),y_min)
			par1(i,3) = z_min + ran3(ra)*(z_max-z_min) 
                        par1(i,4) = dura_min + ran3(ra)*(min(par1(i,7)*2,dura_max)-dura_min)
                        par1(i,5) = 0.01
                        par1(i,6) = 0
                        par1(i,7) = tbtemp(i) + ran3(ra)*(tetemp(i)-tbtemp(i))
                        if(i .eq. 1) then
                        par1(i,1) = 0.0+ran3(ra)*(0.0-(-0.0)) ! fix the location of first event
                        par1(i,2) = 0.0+ran3(ra)*(0.0-(-0.0))
                        par1(i,4) = dura_min + ran3(ra)*(min(par1(i,7)*2,dura_max)-dura_min)
                        endif
		enddo
                if (neq==1) then
                        j2=0
                else
                do i = 2, neq
                        dx = par1(i, 1) - par1(1, 1)
                        dy = par1(i, 2) - par1(1, 2)
                        dz = par1(i, 3) - par1(1, 3)
                        distances(i) = sqrt(dx**2 + dy**2 + dz**2)
                end do
                distances(1) = 0.0  ! Distance from the first event to itself
                do i = 1, neq-1
                        min_idx = i
                        do k = i+1, neq
                            if (distances(k) < distances(min_idx)) then
                                min_idx = k
                            end if
                        end do
                        if (min_idx /= i) then
                            temp_distance = distances(i)
                            distances(i) = distances(min_idx)
                            distances(min_idx) = temp_distance
                            temp_xyz = par1(i, 1:3)
                            par1(i, 1:3) = par1(min_idx, 1:3)
                            par1(min_idx, 1:3) = temp_xyz
                        end if
                end do
                if (neq==2) then
                        j2=0
                endif
                j4=0
                if (neq>=3) then
                        do i = 3, neq
                        ! Compute vectors
                        v1x = par1(i, 1) - par1(i-1, 1) ! Ei-1 to Ei
                        v1y = par1(i, 2) - par1(i-1, 2)
                        v2x = par1(i, 1) - par1(1, 1) ! E1 to Ei
                        v2y = par1(i, 2) - par1(1, 2) 
                        v3x = par1(i-1, 1) - par1(1, 1) ! E1 to Ei-1
                        v3y = par1(i-1, 2) - par1(1, 2)
                        dot_product1 = v1x * v3x + v1y * v3y
                        dot_product2 = -v3x * v2x - v3y * v2y
                        nm1 = sqrt(v1x**2 + v1y**2)
                        nm2 = sqrt(v2x**2 + v2y**2)
                        nm3 = sqrt(v3x**2 + v3y**2)

                        ! Compute angles in degrees
                        if (nm1 > 0.0 .and. nm3 > 0.0) then
                            angle1 = acos(dot_product1 / (nm1 * nm3)) / deg_to_rad
                        else
                            angle1 = 90.0  ! Treat as 90 degrees if vectors have zero length
                        endif
                        if (nm2 > 0.0 .and. nm3 > 0.0) then
                            angle2 = acos(dot_product2 / (nm2 * nm3)) / deg_to_rad
                        else
                            angle2 = 90.0  ! Treat as 90 degrees if vectors have zero length
                        endif
                        if (angle1 > 45.0 .and. angle2 > 45.0) then
                            j4 = -1  ! Indicate an invalid configuration
                        endif
                        end do
                        j2=j4
                endif
                endif
         enddo
	 call convertModel(neq,neq_max,nsubmax,par1,rinitmodel)
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
     &   res1, fltsnp, fltsdp, fltsnsh,fltsdsh,fltsnrayl,&
     &   fltsdrayl, insectp, insectsh, &
     &   insectrayl, wp, wsh, wrayl, CMTsolution, &
     &   ralpha, rcmt, rdc, elo, ela, evp, evs, rbdist, redist, &
     &   rdisti, nndist, rbdep, redep, rdepi, ndep)
		if (res1<999.0 ) then
                        if (res1<residual) then
                                residual=res1
                                para=par1
                        endif
		endif
                !print *, para(1,1:2), para(2,1:2), para(3,1:2), res1, residual, j2, j4, j, angle1, angle2
        enddo
        para_min=para
        write(*,*) "Done Initializing."
! Get initial likelihood or misfit
	like = residual
	write(*,*) "Rank ",ran,"Like_init",like


! Start the sampling of the posterior distribution
	like_min = like
	sample=0
	sampleO=0
	ount=0
	th=0
	Pr1=0.1
	Pr2=0.1
	Pr3=0.1
	Pr4=0.1
	Pr5=0.1
        Pr6=0.1
        Pr7=0.1
	PrB=0
	PrD=0
	Ac1=0
	Ac2=0
	Ac3=0
	Ac4=0
	Ac5=0
        Ac6=0
        Ac7=0
	AcB=0
	AcD=0

	do while (sample<nsample)
		para_prop=para
		neq_prop=neq
		like_prop=like
		u=ran3(ra)
		ind=1
		out=1
		para1=0
		para2=0
		para3=0
		para4=0
		para5=0
                para6=0
                para7=0
		birth=0
		death=0
		logprob=0
                ra1=ra
                like_prop=like
!               Propose a new model
		if(u<0.2) then ! change the location of the sub-event x
			para1=1
			! for the ind's subevent
			ind=2+ANINT(ran3(ra)*(neq-2)) 
			Pr1 = Pr1 + 1
                        para_prop(ind,1)=para(ind,1)+gasdev(ra)*ppara1
                        if(para_prop(ind,1)<x_min .or. para_prop(ind,1)>x_max) then
                                out = 0
                        endif

			! check if outside bounds of prior
		elseif(u<0.4) then !y
			para2=1
			ind=2+ANINT(ran3(ra)*(neq-2))
			Pr2 = Pr2 + 1
			para_prop(ind,2)=para(ind,2)+gasdev(ra)*ppara2
                        if(para_prop(ind,2)<y_min .or. para_prop(ind,2)>y_max) then
                                out = 0
                        endif
		elseif(u<0.6) then !z
	        	para3=1
			ind=1+ANINT(ran3(ra)*(neq-1))
			Pr3 = Pr3 + 1!			para_prop(ind,3)=para(ind,3)+gasdev(ra)*ppara3
			para_prop(ind,3)=para(ind,3)+gasdev(ra)*ppara3
			if(para_prop(ind,3)<z_min .or. para_prop(ind,3)>z_max) then
				out = 0
			endif
		elseif(u<0.8) then ! dura
			para4=1
			ind=1+ANINT(ran3(ra)*(neq-1))
			Pr4 = Pr4 + 1
			para_prop(ind,4)=para(ind,4)+gasdev(ra)*ppara4
                        if(para_prop(ind,4)<dura_min .or. para_prop(ind,4)>dura_max) then
                                out = 0
                        endif
		else ! cen
			para7=1
			ind=1+ANINT(ran3(ra)*(neq-1))
			Pr7 = Pr7 + 1
                        para_prop(ind,7)=para(ind,7)+gasdev(ra)*ppara7
                        if(para_prop(ind,7)<cen_min .or. para_prop(ind,7)>cen_max) then
                                out=0
                        endif
		endif
                u=ran3(ra)
                do i=1,neq-1
                        tdiff(i)=abs(para_prop(i+1,7)-para_prop(i,7))
                        if (tdiff(i)<T_gap-(1e-5)) then
                                out=0
                        endif
                        if (para_prop(i,7)<para_prop(i,4)/2) then
                                out=0
                        endif
                enddo

                if (neq>=3) then
                        do i = 3, neq
                        ! Compute vectors
                        v1x = para_prop(i, 1) - para_prop(i-1, 1) ! Ei-1 to Ei
                        v1y = para_prop(i, 2) - para_prop(i-1, 2)
                        v2x = para_prop(i, 1) - para_prop(1, 1) ! E1 to Ei
                        v2y = para_prop(i, 2) - para_prop(1, 2)
                        v3x = para_prop(i-1, 1) - para_prop(1, 1) ! E1 to Ei-1
                        v3y = para_prop(i-1, 2) - para_prop(1, 2)
                        dot_product1 = v1x * v3x + v1y * v3y
                        dot_product2 = -v3x * v2x - v3y * v2y
                        nm1 = sqrt(v1x**2 + v1y**2)
                        nm2 = sqrt(v2x**2 + v2y**2)
                        nm3 = sqrt(v3x**2 + v3y**2)

                        ! Compute angles in degrees
                        if (nm1 > 0.0 .and. nm3 > 0.0) then
                            angle1 = acos(dot_product1 / (nm1 * nm3)) / deg_to_rad
                        else
                            angle1 = 90.0  ! Treat as 90 degrees if vectors have zero length
                        endif
                        if (nm2 > 0.0 .and. nm3 > 0.0) then
                            angle2 = acos(dot_product2 / (nm2 * nm3)) / deg_to_rad
                        else
                            angle2 = 90.0  ! Treat as 90 degrees if vectors have zero length
                        endif
                        if (angle1 > 45.0 .and. angle2 > 45.0) then
                                out=0
                        endif
                        end do
                endif

                ra1=ra1+1


! After moving a cell, Get the misfit of the proposed model
		if(out==1) then
                        ount=ount+1
			call convertModel(neq_prop,neq_max,nsubmax,para_prop,rinitmodel)
			insubev = neq
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

			like_prop = residual
		

!  Now, depending on the type of move, compute acceptance term, and see if we accept or reject the model
		accept = 0
                pprop=1
                porg=1
                do i=1,neq
                pprop=pprop*seisprior(para_prop(i,1),para_prop(i,2),griddx,nl,gridx,gridy,gridz,p1)
                porg=porg*seisprior(para(i,1),para(i,2),griddx,nl,gridx,gridy,gridz,p1)
                enddo
			if (log(ran3(ra)*porg/pprop)<(-like_prop+like)/like_min/2/0.05**2)then
				accept=1
				if(para1==1) then
					Ac1=Ac1+1
				elseif(para2==1) then
					Ac2=Ac2+1
				elseif(para3==1) then
					Ac3=Ac3+1
				elseif(para4==1) then
					Ac4=Ac4+1
				elseif(para5==1) then
					Ac5=Ac5+1
                                elseif(para6==1) then
                                        Ac6=Ac6+1
                                elseif(para7==1) then
                                        Ac7=Ac7+1
				endif			
			endif
		
!   If we accept the proposed model, update the status of the Markov Chain
		if (accept == 1) then
			para = para_prop
			like = like_prop
			neq = neq_prop
!   Output models in each step.
!			do i=1,neq
!			write(*,*) i,para(i,5),para(i,1),para(i,2),para(i,4),para(i,3)
!			enddo
!			write(*,*) "*****************************",like
                        i1=0

			if(like<like_min .and. i1==0) then
				like_min = like
				para_min = para
			endif
		endif
! Sotre models for ensemble solution


		if (ount.GT.burn_in) then
			sample=sample+1

			if (mod(ount,thin)==0) then
! Sotre all the misfits	for evety thin's model
			sampleO = sampleO+1
                        if (sampleO .le. outsample) then
			do i=1,neq
                                fitall(sampleO,i,1) = like
                                fitall(sampleO,i,2) = para(i,7)
                                fitall(sampleO,i,3) = para(i,1)
                                fitall(sampleO,i,4) = para(i,2)
                                fitall(sampleO,i,5) = para(i,4)
                                fitall(sampleO,i,6) = para(i,5)
                                fitall(sampleO,i,7) = para(i,6)
                                fitall(sampleO,i,8) = para(i,3)
			enddo
                        endif

!			histo(neq)=histo(neq)+1
!			do i=1,neq
!			        fx=ceiling(para(i,1)-x_min)*disX/(x_max-x_min)
!			        fy=ceiling(para(i,2)-y_min)*disY/(y_max-y_min)
!			        fz=ceiling(para(i,3)-z_min)*disZ/(z_max-z_min)
!			        post(fx,fy)=post(fx,fy)+1
!			enddo
			end if		
		end if	
!		endif
	if(mod(ount,display) .EQ. 0) then
		write(*,*) 'sample',ount,'/',burn_in+nsample,rank
		write(*,*) 'number of sub-events:',neq
		write(*,*) 'like',like,'line_min',like_min
                do i=1,neq
                        write(*,'(7F9.3)') para_min(i,7),para_min(i,1),para_min(i,2),&
                        & para_min(i,4),para_min(i,5),para_min(i,6),para_min(i,3)

                enddo
!		if (ount.GT.burn_in) write(*,*)'Accept rate',100*AcB/PrB,100*AcD/PrD
                if (ount.GT.burn_in) write(*,*)'Accept rate',100*Ac1/Pr1,100*Ac2/Pr2,100*Ac3/Pr3,100*Ac4/Pr4,100*Ac5/Pr5,&
                        & 100*Ac6/Pr6,100*Ac7/Pr7

	endif
        endif

	enddo

!	call MPI_REDUCE(post,posts,disX*disY,MPI_Real,MPI_Sum,0,MPI_COMM_WORLD,ierror)
! 	call MPI_REDUCE(histo,histos,neq_max,MPI_Integer,MPI_Sum,0,MPI_COMM_WORLD,ierror)

        write(number,'(I04)') ran
        filename=trim(adjustl(number))//'chain.dat'
        open(ran,file=filename,status='replace')
        do i=1,outsample
                if (fitall(i,1,1) > 0.0) then
                write(ran,'(I2,$)') neq
                write(ran,'(A8,$)') " <neq "
                do j=1,neq
                write(ran,'(8F10.4,$)') fitall(i,j,1), &
                &fitall(i,j,2),fitall(i,j,3),&
                &fitall(i,j,4),fitall(i,j,5),fitall(i,j,6)&
                &,fitall(i,j,7),fitall(i,j,8)
                enddo
                write(ran,*)
                endif
        enddo
        close(ran)

! sort all the misfits
	do i=1,outsample-1
	do j=i+1,outsample
		if(fitall(i,1,1) .gt. fitall(j,1,1)) then
			do l=1,neq			
			do k=1,8
				fittmp(l,k) = fitall(i,l,k)
				fitall(i,l,k) = fitall(j,l,k)
				fitall(j,l,k) = fittmp(l,k)
			enddo
			enddo
		endif
	enddo
	enddo

! store the best solutions
	filename2=trim(adjustl(number))//'best.dat'
	open(ran+100,file=filename2,status='replace')
	do i=1,best
		if (fitall(i,1,1) > 0.0) then
		write(ran+100,'(I2,$)') neq
                write(ran+100,'(A8,$)') " <neq "
		do j=1,neq
		write(ran+100,'(8F10.4,$)') fitall(i,j,1), &
                &fitall(i,j,2),fitall(i,j,3),&
		&fitall(i,j,4),fitall(i,j,5),fitall(i,j,6),&
                &fitall(i,j,7),fitall(i,j,8)
		enddo
                write(ran+100,*)
		endif
	enddo
	close(ran+100)
	close(ran)
	CALL cpu_time(t2)
	if (ran==0) write(*,*)'time taken by the code was',t2-t1,'seconds'

        call MPI_FINALIZE(ierror)

        deallocate(seisdp_xdata, gfdp, seisdsh_xdata, gfdsh,&
     &   seisdrayl_xdata, gfdrayl)
	deallocate(fitall)

end program

subroutine convertModel(neq,neq_max,nsubmax,para,rinitmodel)
	implicit none
	integer :: neq,neq_max,nsubmax,i
	real :: para(nsubmax,7)
	real :: rinitmodel(120)
	do i=1,neq
		
                rinitmodel((i-1)*7+1) = para(i,7)
                rinitmodel((i-1)*7+2) = para(i,1)
                rinitmodel((i-1)*7+3) = para(i,2)
                rinitmodel((i-1)*7+4) = para(i,4)
                rinitmodel((i-1)*7+5) = para(i,5)
                rinitmodel((i-1)*7+6) = para(i,6)
                rinitmodel((i-1)*7+7) = para(i,3)
	enddo
        if (neq .lt. neq_max) then
        do i=neq+1,neq_max
                para(i,5)=0
                para(i,1)=0
                para(i,2)=0
                para(i,3)=0
                para(i,4)=0
                para(i,6)=0
                para(i,7)=0
        enddo
        endif
	return
end subroutine

!       Function for density calculation 
        function seisprior(xx,yy,dx,nl,gridx,gridy,gridz)
        real :: xx,yy,min_d,x1,y1,p1,x2,y2,p2,x3,y3,p3,x4,y4,p4,&
         dx,nx,ny,seisprior
        integer :: nl,ii,jj
        real :: gridx(nl),gridy(nl),gridz(nl)
        min_d=999
        do ii=1,nl
                if (min_d>sqrt((xx-gridx(ii))**2+(yy-gridy(ii))**2)) then
                        min_d=sqrt((xx-gridx(ii))**2+(yy-gridy(ii))**2)
                        x1=gridx(ii)
                        y1=gridy(ii)
                        p1=gridz(ii)
                endif
        enddo
        if (xx<x1) then
                if (yy<y1) then
                        x4=x1
                        y4=y1
                        x3=x1-dx
                        y3=y1
                        x2=x1
                        y2=y1-dx
                        x1=x1-dx
                        y1=y1-dx
                else
                        x4=x1
                        y4=y1+dx
                        x3=x1-dx
                        y3=y1+dx
                        x2=x1
                        y2=y1
                        x1=x1-dx
                endif
        elseif (yy<y1) then
                x4=x1+dx
                y4=y1
                x3=x1
                y3=y1
                x2=x1+dx
                y2=y1-dx
                y1=y1-dx
        else
                x4=x1+dx
                y4=y1+dx
                x3=x1
                y3=y1+dx
                x2=x1+dx
                y2=y1
        endif
        p1=0
        p2=0
        p3=0
        p4=0
        do ii=1,nl
                if ((x1-gridx(ii))**2+(y1-gridy(ii))**2<1e-4) then
                        p1=gridz(ii)
                endif
                if ((x2-gridx(ii))**2+(y2-gridy(ii))**2<1e-4) then
                        p2=gridz(ii)
                endif
                if ((x3-gridx(ii))**2+(y3-gridy(ii))**2<1e-4) then
                        p3=gridz(ii)
                endif
                if ((x4-gridx(ii))**2+(y4-gridy(ii))**2<1e-4) then
                        p4=gridz(ii)
                endif
        enddo
        nx=(xx-x1)/dx
        ny=(yy-y1)/dx
        seisprior=p1*(1-nx)*(1-ny)+p2*nx*(1-ny)+p3*(1-nx)*ny+p4*nx*ny
        return
        end

!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator for 
!       a Gaussian distribution
!
! ----------------------------------------------------------------------------


FUNCTION GASDEV(idum)

!     ..Arguments..
integer          idum
real GASDEV

!     ..Local..
real v1,v2,r,fac
real ran3

if (idum.lt.0) iset=0
10   v1=2*ran3(idum)-1
v2=2*ran3(idum)-1
r=v1**2+v2**2
if(r.ge.1.or.r.eq.0) GOTO 10
fac=sqrt(-2*log(r)/r)
GASDEV=v2*fac

RETURN
END


!-------------------------------------------------------------------
!						
!	Numerical Recipes random number generator
!
! ----------------------------------------------------------------------------

FUNCTION ran3(idum)
INTEGER idum
INTEGER MBIG,MSEED,MZ
!     REAL MBIG,MSEED,MZ
REAL ran3,FAC
PARAMETER (MBIG=1000000000,MSEED=161803478,MZ=0,FAC=1./MBIG)
!PARAMETER (MBIG=4000000.,MSEED=1618033.,MZ=0.,FAC=1./MBIG)
INTEGER i,iff,ii,inext,inextp,k
INTEGER mj,mk,ma(55)
!     REAL mj,mk,ma(55)
SAVE iff,inext,inextp,ma
DATA iff /0/
! write(*,*)' idum ',idum
if(idum.lt.0.or.iff.eq.0)then
	iff=1
	mj=MSEED-iabs(idum)
	mj=mod(mj,MBIG)
	ma(55)=mj
	mk=1
	do 11 i=1,54
	ii=mod(21*i,55)
	ma(ii)=mk
	mk=mj-mk
	if(mk.lt.MZ)mk=mk+MBIG
	mj=ma(ii)
!  write(*,*)' idum av',idum
11      continue
	do 13 k=1,4
	do 12 i=1,55
	ma(i)=ma(i)-ma(1+mod(i+30,55))
	if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12        continue
13      continue
! write(*,*)' idum ap',idum
	inext=0
	inextp=31
	idum=1
endif
! write(*,*)' idum app ',idum
inext=inext+1
if(inext.eq.56)inext=1
inextp=inextp+1
if(inextp.eq.56)inextp=1
mj=ma(inext)-ma(inextp)
if(mj.lt.MZ)mj=mj+MBIG
ma(inext)=mj
ran3=mj*FAC
!  write(*,*)' idum ',idum
	
return
END


