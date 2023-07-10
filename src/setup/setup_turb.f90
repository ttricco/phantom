!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2023 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module setup
!
! Sets up a calculation of supersonic turbulence in a periodic box.
!  Works for hydro, mhd, and dusty turbulence.
!
! :References:
!    Price & Federrath (2010), MNRAS
!    Tricco, Price & Federrath (2016), MNRAS
!    Tricco, Price & Laibe (2017), MNRAS Letters
!
! :Owner: Daniel Price
!
! :Runtime parameters: None
!
! :Dependencies: boundary, dim, dust, io, mpidomain, mpiutils, options,
!   part, physcon, prompting, set_dust, setup_params, table_utils,
!   timestep, unifdis, units
!
 implicit none
 public :: setpart

 integer, private :: npartx,ilattice
 real,    private :: polykset
 private

contains

!----------------------------------------------------------------
!+
!  setup for uniform particle distributions
!+
!----------------------------------------------------------------
subroutine setpart(id,npart,npartoftype,xyzh,massoftype,vxyzu,polyk,gamma,hfact,time,fileprefix)
 use dim,          only:use_dust,maxdustsmall,maxvxyzu,periodic
 use options,      only:nfulldump,beta,use_dustfrac
 use setup_params, only:rhozero,npart_total,ihavesetupB
 use io,           only:master
 use unifdis,      only:set_unifdis,latticetype
 use boundary,     only:set_boundary,xmin,ymin,zmin,xmax,ymax,zmax,dxbound,dybound,dzbound
 use mpiutils,     only:bcast_mpi
 use part,         only:Bxyz,mhd,dustfrac,grainsize,graindens,ndusttypes,igas,idust,set_particle_type,ndustsmall
 use physcon,      only:pi,solarm,pc,km
 use units,        only:set_units,udist,umass
 use prompting,    only:prompt
 use dust,         only:grainsizecgs,graindenscgs
 use set_dust_options, only:set_dust_default_options,set_dust_interactively,dust_method,dust_to_gas,ndustsmallinp
 use set_dust,     only:set_dustfrac,set_dustbinfrac
 use timestep,     only:dtmax,tmax
 use table_utils,  only:logspace
 use mpidomain,    only:i_belong
 integer,           intent(in)    :: id
 integer,           intent(inout) :: npart
 integer,           intent(out)   :: npartoftype(:)
 real,              intent(out)   :: xyzh(:,:)
 real,              intent(out)   :: massoftype(:)
 real,              intent(out)   :: vxyzu(:,:)
 real,              intent(out)   :: polyk,gamma,hfact
 real,              intent(inout) :: time
 character(len=20), intent(in)    :: fileprefix
 character(len=26)                :: filename
 integer :: i, j
 logical :: iexist
 real :: totmass,deltax
 real :: Bz_0
 real :: smincgs,smaxcgs,sindex,dustbinfrac(maxdustsmall)

 print *, ''
 print *, 'Setup for turbulence in a periodic box'
 print *, ''
!
!--boundaries
!
 call set_boundary(0.,1.,0.,1.,0.,1.)
!
!--units
!
 ! Molecular cloud conditions:
 !  L = 3 pc, rho = 1e-20 g/cm^3, c_s = 0.2 km/s
 call set_units(dist=3.0*pc, mass=1e-20*(3.0*pc)**3, time=3.0*pc/(0.2*km))
!
!--general parameters
!
 time = 0.
 hfact = 1.2
 gamma = 1.
!
!--setup particles
!
 if (id==master) then
    npartx = 64
    call prompt('Enter number of particles in x ',npartx,16)
 endif
 call bcast_mpi(npartx)
 deltax = dxbound/npartx

 if (id==master) then
    ilattice = 1
    call prompt('Select lattice type (1=cubic, 2=closepacked)',ilattice,1,2)
 endif

 call bcast_mpi(ilattice)
 if (id==master) then
    rhozero = 1.
    call prompt('Enter density (gives particle mass)',rhozero,0.)
 endif
 call bcast_mpi(rhozero)

 if (maxvxyzu < 4) then
    if (id==master) then
       polykset = 1.
       call prompt('Enter sound speed in code units (sets polyk)',polykset,0.)
    endif
    call bcast_mpi(polykset)
    polyk = polykset**2
    print*,' polyk = ',polyk
 else
    polyk = 0.
    polykset = 0.
 endif

 if (use_dust) then
    print *, ''
    print *, 'Setting up dusty turbulence:'
    print *, ''

    call set_dust_default_options()
    call set_dust_interactively()

    if (dust_method == 1) then
        use_dustfrac = .true.
        ndusttypes = ndustsmallinp
        ndustsmall = ndustsmallinp
        if (ndusttypes > 1) then
           call set_dustbinfrac(smincgs,smaxcgs,sindex,dustbinfrac(1:ndusttypes),grainsize(1:ndusttypes))
           grainsize(1:ndusttypes) = grainsize(1:ndusttypes)/udist
           graindens(1:ndusttypes) = graindens(1)/umass*udist**3
        else
           grainsize(1) = grainsizecgs/udist
           graindens(1) = graindenscgs/umass*udist**3
        endif
    elseif (dust_method == 2) then
        use_dustfrac = .false.
    endif
 endif

 if (mhd) then
    Bz_0 = 1.4142e-5
    print *, ''
    print *, 'Setting up MHD turbulence: (with uniform intial magnetic field in z-direction)'
    print *, ''
    if (id==master) call prompt('Enter initial magnetic field strength ',Bz_0)
 endif


 ! setup preferred values of .in file
 filename= trim(fileprefix)//'.in'
 inquire(file=filename,exist=iexist)
 if (.not. iexist) then
    tmax         = 1.00   ! run for 20 turbulent crossing times
    dtmax        = 0.0025
    nfulldump    = 5      ! output 4 full dumps per crossing time
    beta         = 4      ! legacy from Price & Federrath (2010), haven't checked recently if still required
 endif
 npart = 0
 npart_total = 0

 call set_unifdis(latticetype(ilattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                  deltax,hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)

 npartoftype(:) = 0
 npartoftype(igas) = npart
 print *, ' npart = ',npart,npart_total

 totmass = rhozero*dxbound*dybound*dzbound
 massoftype = totmass/npart_total
 print *, ' particle mass = ',massoftype(igas)

 if (use_dust .and. dust_method == 2) then
     call set_unifdis(latticetype(ilattice),id,master,xmin,xmax,ymin,ymax,zmin,zmax,&
                      deltax,hfact,npart,xyzh,periodic,nptot=npart_total,mask=i_belong)
     xyzh(1:3, npartoftype(igas):npart_total) = xyzh(1:3, npartoftype(igas):npart_total) + 0.5 * deltax
     npartoftype(idust) = npartoftype(igas)
     massoftype(idust)  = totmass * dust_to_gas / npartoftype(idust)
     massoftype(igas) = totmass * (1.0 - dust_to_gas) / npartoftype(igas)

     do i=npartoftype(igas)+1, npart_total
         call set_particle_type(i,idust)
     enddo

     print *, ' dust particle mass = ',massoftype(idust)

 endif

 do i=1,npartoftype(igas)
    call set_particle_type(i,igas)
    vxyzu(1:3,i) = 0.
    if (mhd) then
       Bxyz(:,i) = 0.
       Bxyz(3,i) = Bz_0
    endif

    if (use_dust .and. dust_method == 1) then
        if (ndusttypes > 1) then
            dustfrac(1:ndusttypes,i) = dust_to_gas*dustbinfrac(1:ndusttypes)
        else
            call set_dustfrac(dust_to_gas,dustfrac(:,i))
        endif
    endif
 enddo

 totmass = npartoftype(igas) * massoftype(igas) + npartoftype(idust) * massoftype(idust)
 print *, ' total mass = ', totmass

 if (mhd) ihavesetupB = .true.

end subroutine setpart

end module setup
