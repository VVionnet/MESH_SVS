module runsvs_mesh

    !> MESH modules.
    !*  mpi_module: Required for 'il1' and 'il2' indexing.
    !*  model_files_variables: Required for 'fls' object.
    !*  sa_mesh_common: Required for MESH variables and common routines.
    !*  model_dates: Required for 'ic' counter.
    use mpi_module
    use model_files_variables
    use sa_mesh_common
!todo: Replace 'cm' instance with 'vs' counterparts.
    use climate_forcing
    use model_dates

    !> SPS/SVS modules.
    !*  sfcbus_mod: For active variables 'vl', dictionary 'vd', and related indices.
    !*  sfc_options: Surface configuration options and variables.
    !*  svs_options: SVS configuration options and variables.
    use sfcbus_mod
    use sfc_options
    use svs_configs

    use str_mod, only: str_concat

    implicit none

    character(len=1024) :: msg_S

    !> SVS constants.
    type runsvs_mesh_constants
        integer :: NLANDCLASS = NCLASS
        real, dimension(NCLASS) :: Z0DAT = (/ &
            0.001, 0.001, 0.001, 1.75, 2.0, 1.0, 2.0, 3.0, 0.8, 0.1, &
            0.2, 0.2, 0.1, 0.1, 0.15, 0.15, 0.35, 0.25, 0.1, 0.25, &
            5.0, 0.1, 0.1, 0.1, 1.75, 0.5 /)
    end type

    !> SVS2 output
    integer :: iout_soil = 150
    integer :: iout_snow_bulk = 151
    integer :: iout_snow_profile = 152
    integer :: iout_snow_enbal = 153

    !> SVS1 output
    integer :: iout_svs1_soil = 160
    integer :: iout_svs1_snow = 161

    !> SVS variables names for I/O (direct variables).
    character(len = *), parameter, public :: VN_SVS_DEGLAT = 'DEGLAT'
    character(len = *), parameter, public :: VN_SVS_DEGLNG = 'DEGLNG'
    character(len = *), parameter, public :: VN_SVS_OBSERVED_FORCING = 'OBSERVED_FORCING'
    character(len = *), parameter, public :: VN_SVS_ZUSL = 'ZUSL'
    character(len = *), parameter, public :: VN_SVS_ZTSL = 'ZTSL'
    character(len = *), parameter, public :: VN_SVS_SIGMA_U = 'SIGMA_U'
    character(len = *), parameter, public :: VN_SVS_SIGMA_T = 'SIGMA_T'
    character(len = *), parameter, public :: VN_SVS_SLOP = 'SLOP'
    character(len = *), parameter, public :: VN_SVS_DRAINDENS = 'DRAINDENS'
    character(len = *), parameter, public :: VN_SVS_SOILTEXT = 'SOILTEXT'
    character(len = *), parameter, public :: VN_SVS_SCHMSOL = 'SCHMSOL'
    character(len = *), parameter, public :: VN_SVS_KHYD = 'KHYD'
    character(len = *), parameter, public :: VN_SVS_SAND = 'SAND'
    character(len = *), parameter, public :: VN_SVS_CLAY = 'CLAY'
    character(len = *), parameter, public :: VN_SVS_WSOIL = 'WSOIL'
    character(len = *), parameter, public :: VN_SVS_ISOIL = 'ISOIL'
    character(len = *), parameter, public :: VN_SVS_KTHERMAL = 'KTHERMAL'
    character(len = *), parameter, public :: VN_SVS_TGROUND = 'TGROUND'
    character(len = *), parameter, public :: VN_SVS_VF = 'VF'
    character(len = *), parameter, public :: VN_SVS_Z0V = 'Z0V'
    character(len = *), parameter, public :: VN_SVS_LNZ0 = 'LNZ0'
    character(len = *), parameter, public :: VN_SVS_TVEGE = 'TVEGE'
    character(len = *), parameter, public :: VN_SVS_WVEG = 'WVEG'
    character(len = *), parameter, public :: VN_SVS_TSNOW = 'TSNOW'
    character(len = *), parameter, public :: VN_SVS_SNODPL = 'SNODPL'
    character(len = *), parameter, public :: VN_SVS_SNODEN = 'SNODEN'
    character(len = *), parameter, public :: VN_SVS_SNOAL = 'SNOAL'
    character(len = *), parameter, public :: VN_SVS_WSNOW = 'WSNOW'
    character(len = *), parameter, public :: VN_SVS_TSNOWVEG = 'TSNOWVEG'
    character(len = *), parameter, public :: VN_SVS_SNVDP = 'SNVDP'
    character(len = *), parameter, public :: VN_SVS_SNVDEN = 'SNVDEN'
    character(len = *), parameter, public :: VN_SVS_SNVAL = 'SNVAL'
    character(len = *), parameter, public :: VN_SVS_WATPOND = 'WATPOND' !  For svs1 (with water ponding)
    character(len = *), parameter, public :: VN_SVS_WSNV = 'WSNV'
    character(len = *), parameter, public :: VN_SVS_TPSOIL = 'TPSOIL' ! For svs2  and svs1 (with soil freezing)
    character(len = *), parameter, public :: VN_SVS_TPSOILV = 'TPSOILV'! For svs2 only
    character(len = *), parameter, public :: VN_SVS_TPERM = 'TPERM' ! For svs2 only  and svs1 (with soil freezing)
    character(len = *), parameter, public :: VN_SVS_NSL = 'NSL' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWSCHEME = 'HSNOWSCHEME' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWDRIFT_CRO = 'HSNOWDRIFT_CRO' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWMETAMO = 'HSNOWMETAMO' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWRAD = 'HSNOWRAD' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWFALL = 'HSNOWFALL' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWCOND = 'HSNOWCOND' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWCOMP = 'HSNOWCOMP' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_HSNOWHOLD = 'HSNOWHOLD' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_LSNOWDRIFT_SUBLIM = 'LSNOWDRIFT_SUBLIM' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_SNOMA = 'SNOMA'
    character(len = *), parameter, public :: VN_SVS_SNVMA = 'SNVMA'
    character(len = *), parameter, public :: VN_SVS_SNOMA_SVS = 'SNOMA_ML'
    character(len = *), parameter, public :: VN_SVS_SNODEN_SVS = 'SNODEN_ML'
    character(len = *), parameter, public :: VN_SVS_SNOAGE_SVS = 'SNOAGE_ML'
    character(len = *), parameter, public :: VN_SVS_SNODIAMOPT_SVS = 'SNODOPT_ML'
    character(len = *), parameter, public :: VN_SVS_SNOSPHERI_SVS = 'SNOSPH_ML'
    character(len = *), parameter, public :: VN_SVS_SNOHIST_SVS = 'SNOHIST_ML'
    character(len = *), parameter, public :: VN_SVS_TSNOW_SVS = 'TSNOW_ML'
    character(len = *), parameter, public :: VN_SVS_WSNOW_SVS = 'WSNOW_ML'
    character(len = *), parameter, public :: VN_SVS_LOUT_SNOW_PROFILE = 'LOUT_SNOW_PROFILE' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_LOUT_SNOW_ENBAL = 'LOUT_SNOW_ENBAL' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_NPROFILE_DAY = 'NPROFILE_DAY' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_LSOIL_FREEZING_SVS1 = 'LSOIL_FREEZING_SVS1' ! For svs1 only 
    character(len = *), parameter, public :: VN_SVS_LWATER_PONDING_SVS1 = 'LWATER_PONDING_SVS1' ! For svs1 only 

    !> SVS variables names for I/O (modifiers/special conditions).
    character(len = *), parameter, public :: VN_SVS_SAND_N = 'SAND_N'
    character(len = *), parameter, public :: VN_SVS_CLAY_N = 'CLAY_N'
    character(len = *), parameter, public :: VN_SVS_WSOIL_N = 'WSOIL_N'
    character(len = *), parameter, public :: VN_SVS_ISOIL_N = 'ISOIL_N'
    character(len = *), parameter, public :: VN_SVS_TGROUND_N = 'TGROUND_N'
    character(len = *), parameter, public :: VN_SVS_VF_N = 'VF_N'
    character(len = *), parameter, public :: VN_SVS_Z0V_N = 'Z0V_N'
    character(len = *), parameter, public :: VN_SVS_TVEGE_N = 'TVEGE_N'
    character(len = *), parameter, public :: VN_SVS_TSNOW_N = 'TSNOW_N'
    character(len = *), parameter, public :: VN_SVS_TSNOWVEG_N = 'TSNOWVEG_N'
    character(len = *), parameter, public :: VN_SVS_TPSOIL_N = 'TPSOIL_N' ! For svs2  and svs1 (with soil freezing)
    character(len = *), parameter, public :: VN_SVS_TPSOILV_N = 'TPSOILV_N' ! For svs2 only


    !> SVS variables (for I/O).
    type runsvs_mesh_variables
        real, dimension(:), allocatable :: deglat
        real, dimension(:), allocatable :: deglng
        logical :: observed_forcing = .false.
        real, dimension(:), allocatable :: zusl
        real, dimension(:), allocatable :: ztsl
        real :: sigma_u = 0.995
        real :: sigma_t = 0.995
        real, dimension(:), allocatable :: slop
        real, dimension(:), allocatable :: draindens
        character(len = DEFAULT_FIELD_LENGTH) :: soiltext = 'NIL'
        character(len = DEFAULT_FIELD_LENGTH) :: schmsol = 'SVS'
        integer :: khyd = 6
        real, dimension(:, :), allocatable :: sand
        real, dimension(:, :), allocatable :: clay
        real, dimension(:, :), allocatable :: wsoil
        real, dimension(:, :), allocatable :: isoil
        real, dimension(:, :), allocatable :: tpsoil ! For svs2 and svs1 (with soil freezing)
        real, dimension(:, :), allocatable :: tpsoilv ! For svs2 only
        integer :: kthermal = 2
        real, dimension(:, :), allocatable :: tground
        real, dimension(:, :), allocatable :: vf
        real, dimension(:, :), allocatable :: z0v
        real, dimension(:), allocatable :: lnz0
        real, dimension(:, :), allocatable :: tvege
        real, dimension(:), allocatable :: wveg
        real, dimension(:, :), allocatable :: tsnow
        real, dimension(:), allocatable :: snodpl
        real, dimension(:), allocatable :: snoden
        real, dimension(:), allocatable :: snoal
        real, dimension(:), allocatable :: wsnow
        real, dimension(:, :), allocatable :: tsnowveg
        real, dimension(:), allocatable :: snvdp
        real, dimension(:), allocatable :: snvden
        real, dimension(:), allocatable :: snval
        real, dimension(:), allocatable :: watpond
        real, dimension(:), allocatable :: wsnv
        real, dimension(:), allocatable :: tperm ! For svs2 only
        integer :: nsl = 12 ! For svs2 only
        real, dimension(:,:), allocatable :: snoma_svs 
        real, dimension(:,:), allocatable :: snoden_svs 
        real, dimension(:,:), allocatable :: snodiamopt_svs 
        real, dimension(:,:), allocatable :: snoage_svs 
        real, dimension(:,:), allocatable :: snospheri_svs 
        real, dimension(:,:), allocatable :: snohist_svs 
        real, dimension(:,:), allocatable :: tsnow_svs
        real, dimension(:,:), allocatable :: wsnow_svs
        real, dimension(:,:), allocatable :: snomav_svs 
        real, dimension(:,:), allocatable :: snodenv_svs 
        real, dimension(:,:), allocatable :: snodiamoptv_svs 
        real, dimension(:,:), allocatable :: snoagev_svs 
        real, dimension(:,:), allocatable :: snospheriv_svs 
        real, dimension(:,:), allocatable :: snohistv_svs 
        real, dimension(:,:), allocatable :: tsnowv_svs
        real, dimension(:,:), allocatable :: wsnowv_svs
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowscheme = 'ES'
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowdrift_cro = 'ES'
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowmetamo = 'CI13'
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowrad = 'B92'
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowfall = 'VI12'
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowcond = 'Y81'
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowhold = 'B92'
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowcomp = 'B92'
        logical :: lsnowdrift_sublim = .true.
        logical :: lout_snow_profile = .false.
        logical :: lout_snow_enbal = .false.
        integer :: nprofile_day = 4 !
        logical :: lsoil_freezing_svs1 = .false.
        logical :: lwater_ponding_svs1 = .false.
    end type

    !* PROCESS_ACTIVE: Variable to enable SVS.
    type runsvs_mesh_container
        logical :: PROCESS_ACTIVE = .false.
        type(runsvs_mesh_constants) c
        type(runsvs_mesh_variables) vs
    end type

    type(runsvs_mesh_container), save, public :: svs_mesh

    private

    public &
        runsvs_mesh_init, runsvs_mesh_resume_states_seq, runsvs_mesh_within_tile, runsvs_mesh_save_states_seq, runsvs_mesh_finalize

    !> RPN/physics variables.
    !* phyinread_list_S: List of input variables.
    include "sfcinput.cdk"

    !> SPS/SVS internal variables.
    !* phy_bus: 2D variable bus to emulate a 'phybus' type bus (buffer variable).
    !* svs_bus: Variable bus for SVS (buffer variable).
    !* bus_length: Size of the first dimensions of the 'phy_bus' and 'svs_bus' variables.
    !* bus_ptr: Mapped index of SVS variables to surface variables from 'phy_bus'.
    !* time_dt: Model time-step in seconds (MESH: 'ic%dts').
    !* kount: Current time-step within 'time_dt', constant since SVS is called once per MESH time-step (constant: 1).
    !* trnch: 'Row number', used for GEM coupling (not used in SVS directly; constant: 1).
    !* ni (n): 'Running length', i.e., stride of variables, named 'n' in SVS (MESH: NML).
    !* ni (m): 'Horizontal dimension', named 'm' in SVS (passed same variable as for 'n' in 'sfc_main').
    !* nk: 'Vertical dimension', used to map temperature/humidity and momentum variables (not used in SVS directly; constant: 1).
    real, pointer, private :: phy_bus(:, :) => null()
    real, pointer, private :: svs_bus(:) => null()
    integer, private :: bus_length
    integer, allocatable, private :: bus_ptr(:)
    real, private :: time_dt = 0
    real, private :: lmo_winter = -1.
    integer :: kount_reset = 0
    integer, private :: kount = 0
    integer, parameter, private :: trnch = 1
    integer, private :: ni = 0
    integer, parameter, private :: nk = 1

    !> Constants.
    real, parameter, private :: deg2rad = acos(-1.0)/180.0, rad2deg = 180.0/acos(-1.0)

    contains

    subroutine runsvs_mesh_append_phyentvar(variable_name)

        !> Modules.
        use strings, only: lowercase

        !> Input/output variables.
        character(len = *), intent(in) :: variable_name

        !> Local variables.
        character(len = DEFAULT_FIELD_LENGTH) code

        !> Increment the variable count.
        phyinread_n = phyinread_n + 1

        !> Add the name.
        phyinread_list_S(phyinread_n) = trim(lowercase(variable_name))

        !> Check bounds.
        if (phyinread_n > PHYINREAD_MAX) then
            write(code, *) PHYINREAD_MAX
            call print_error( &
                "More variables than are supported (" // trim(adjustl(code)) // ") have been enabled in the surface input bus.")
            call program_abort()
        end if

    end subroutine

    subroutine phy_businit(ni, nk)

        !> For RPN/physics status.
        use phy_status, only: phy_error_L

        !> Modules.
   use cnv_options
        use phy_options
        use phybus

        !> Input/output variables.
        integer, intent(in) :: ni, nk

        !> Local variables.
        character(len = 6) :: nag, ntp, nmar, wwz, nuv, isss
   logical :: lcn_mpx, lcn_my2, lcn_p3i1, lcn_p3i2, lcn_p3i3, lcn_p3i4, lcn_none
   logical :: lbourg3d, lbourg
   logical :: lrslp
   logical :: lkfbe, lshal, lshbkf, lmid
   logical :: ladvtke, nadvtke, lmoistke
   logical :: lmoyhr, lmoyhrkf, lmoykfsh, lmoymid
   logical :: lgwdsm
   logical :: lccc2
   logical :: lghg, ltrigtau
   logical :: liuv
   logical :: lmoyhroz, lmoyhrgh, llinozout, llinghout, llinozage
   logical :: lcndsm
   logical :: lcons, lmoycons
   logical :: lhn_init, lsfcflx
   logical :: lsurfonly

        !> Includes.
        include "mcica.cdk"

        !> Variables from 'ens_perturb'.
        integer ptp_nc, spp_nc, ens_nc2d

        !> From 'phy_init'.
        ptp_nc = 0
        spp_nc = 0
        ens_nc2d = max(ptp_nc + spp_nc, 1)

   !# nagg is the dimension of aggregrated variables ('nagrege').
   write(nag,'(a,i2)') 'A*', nsurf + 1

   !# nipt is the number of tau/cloud top pressure bins in ISCCP histograms
   write(ntp,'(a,i2)') 'A*', ntau*nptop
   write(nuv,'(a,i2)') 'A*', RAD_NUVBRANDS

   !# nmar is the number of 2d Markov fields
   write(nmar,'(a,i2)') 'A*', ens_nc2d

   lcn_mpx  = (stcond(1:2) == 'MP')
   lcn_none = .not.lcn_mpx
   lcn_my2  = (stcond(1:6) == 'MP_MY2')
   lcn_p3i1 = (stcond == 'MP_P3' .and. p3_ncat >= 1)
   lcn_p3i2 = (stcond == 'MP_P3' .and. p3_ncat >= 2)
   lcn_p3i3 = (stcond == 'MP_P3' .and. p3_ncat >= 3)
   lcn_p3i4 = (stcond == 'MP_P3' .and. p3_ncat == 4)

   lbourg3d= (pcptype == 'BOURGE3D')
   lgwdsm  = (sgo_tdfilter > 0.)
   lrslp   = radslope
   ladvtke = advectke
   nadvtke = .not.ladvtke
   lmoyhr  = (moyhr > 0 .or. dynout)
   lkfbe   = any(convec == (/ &
        'BECHTOLD', &
        'KFC     ', &
        'KFC2    ', &
        'KFC3    '  &
        /))
   lshal   = (conv_shal /= 'NIL')
   lshbkf  = (conv_shal == 'BECHTOLD')
   lmid    = (conv_mid /= 'NIL')
   lmoyhrkf= (lmoyhr .and. lkfbe)
   lmoymid= (lmoyhr .and. lmid)
   lmoykfsh= (lmoyhr .and. lshbkf .and. bkf_lshalm)
   lbourg  = any(pcptype == (/&
        'BOURGE', &
        'NIL   '  &
        /))
   lmoistke= (fluvert == 'MOISTKE')
   lccc2   = (radia == 'CCCMARAD2')
   lghg    = (lccc2 .and. radghg_L)
   llinozage = (llinoz .and. age_linoz)              ! age of air tracer off 
   llinozout = (llinoz .and. out_linoz)
   llinghout = (llingh .and. out_linoz)
   lmoyhroz =(lmoyhr .and. llinoz .and. out_linoz)
   lmoyhrgh =(lmoyhr .and. llingh .and. out_linoz)
   lcndsm  = (cond_infilter > 0.)
   ltrigtau = (kfctrigtau > 0.)
   liuv    = (any(radia == (/&
        'CCCMARAD ', &
        'CCCMARAD2'  &
        /)) .and. kntraduv_S /= '')      

        wwz = '1'
        lsurfonly = (fluvert == 'SURFACE')
        if (lsurfonly) wwz = '0'
        isss = '0'
        if (tofd /= 'NIL') isss = '1'

   ! Activate energy budget diagnostics only if outputs are requested by the user
   lcons = .false.
   lmoycons = (lcons .and. lmoyhr)

   lhn_init = (lhn /= 'NIL')
   lsfcflx = (sfcflx_filter_order > 0)

#include "phyvar.hf"
        if (phy_error_L) then
            call print_error("An error occurred initializing the physics bus variable.")
            call program_abort()
        end if

    end subroutine

    subroutine runsvs_mesh_copy_vs_to_bus()

        !> For constants.
        use tdpack_const, only: omega

        !> Local variables.
        real, dimension(il1:il2) :: sumvfz0
        integer i

        !> Functions to derive indices by variable name in the bus.
#define a2(n,l) bus_ptr(vd%n%i)+l*ni
#define a1(n) a2(n,0)
#define z2(n,l) bus_ptr(vd%n%i)+l*ni+ni-1
#define z1(n) z2(n,0)

        !> Reset the bus variable.
        svs_bus = 0.0

        !> Assign momentum and thermodynamic levels if provided observations 'observed_forcing'.
        if (svs_mesh%vs%observed_forcing) then
            if (allocated(svs_mesh%vs%zusl)) svs_bus(a1(zusl):z1(zusl)) = svs_mesh%vs%zusl
            if (allocated(svs_mesh%vs%ztsl)) svs_bus(a1(ztsl):z1(ztsl)) = svs_mesh%vs%ztsl
        end if

        !> Assign variables (transforms from 'physimple_transforms3d').
        if (allocated(svs_mesh%vs%deglat)) svs_bus(a1(dlat):z1(dlat)) = svs_mesh%vs%deglat*deg2rad
        if (allocated(svs_mesh%vs%deglng)) then
            where (svs_mesh%vs%deglng < 0.0)
                svs_bus(a1(dlon):z1(dlon)) = (svs_mesh%vs%deglng + 360.0)*deg2rad
            elsewhere
                svs_bus(a1(dlon):z1(dlon)) = svs_mesh%vs%deglng*deg2rad
            end where
        end if

        ! Compute Coriolis factor
        svs_bus(a1(fcor):z1(fcor)) = 2. * omega *sin(svs_bus(a1(dlat):z1(dlat)))

        svs_bus(a1(z0):z1(z0)) = 0.0
        svs_bus(a1(mg):z1(mg)) = 1.0
        sumvfz0 = 0.0
        do i = 1199, 1174, -1
            svs_bus(a2(vegf, 1199 - i):z2(vegf, 1199 - i)) = svs_mesh%vs%vf(:, 1200 - i)
            if (allocated(svs_mesh%vs%z0v)) then
                svs_bus(a1(z0):z1(z0)) = svs_bus(a1(z0):z1(z0)) + svs_mesh%vs%vf(:, 1200 - i)*svs_mesh%vs%z0v(:, 1200 - i)
            else
                svs_bus(a1(z0):z1(z0)) = svs_bus(a1(z0):z1(z0)) + svs_mesh%vs%vf(:, 1200 - i)*svs_mesh%c%Z0DAT(1200 - i)
            end if
            sumvfz0 = sumvfz0 + svs_mesh%vs%vf(:, 1200 - i)
        end do
        where (sumvfz0 > 0.0)
            svs_bus(a1(z0):z1(z0)) = svs_bus(a1(z0):z1(z0))/sumvfz0
        end where
        if (allocated(svs_mesh%vs%lnz0)) then
            svs_bus(a1(z0en):z1(z0en)) = exp(svs_mesh%vs%lnz0)
            svs_bus(a1(z0mlanden):z1(z0mlanden)) = exp(svs_mesh%vs%lnz0)
        end if
!        svs_bus(a1(z0t):z1(z0t)) = svs_bus(a1(z0):z1(z0))
        if (allocated(svs_mesh%vs%slop)) svs_bus(a1(slop):z1(slop)) = svs_mesh%vs%slop
        if (allocated(svs_mesh%vs%draindens)) svs_bus(a1(draindens):z1(draindens)) = svs_mesh%vs%draindens
        if (soiltext == 'NIL') then
            do i = 1, nl_svs
                if (allocated(svs_mesh%vs%sand)) svs_bus(a2(sand, i - 1):z2(sand, i - 1)) = svs_mesh%vs%sand(:, i)
                if (allocated(svs_mesh%vs%clay)) svs_bus(a2(clay, i - 1):z2(clay, i - 1)) = svs_mesh%vs%clay(:, i)
            end do
            if (svs_mesh%vs%schmsol=='SVS') then
                call inisoili_svs(ni, trnch)
            else if (svs_mesh%vs%schmsol=='SVS2') then
                call inisoili_svs2(ni, trnch)
            endif
        else
            write(*,*) 'val stp',nl_stp,size( svs_mesh%vs%sand,2)

            do i = 1, nl_stp
                if (allocated(svs_mesh%vs%sand)) svs_bus(a2(sanden, i - 1):z2(sanden, i - 1)) = svs_mesh%vs%sand(:, i)
                if (allocated(svs_mesh%vs%clay)) svs_bus(a2(clayen, i - 1):z2(clayen, i - 1)) = svs_mesh%vs%clay(:, i)
            end do
        end if
        do i = 1, nl_svs
            if (allocated(svs_mesh%vs%wsoil)) svs_bus(a2(wsoil, i - 1):z2(wsoil, i - 1)) = svs_mesh%vs%wsoil(:, i)
            if (allocated(svs_mesh%vs%isoil)) svs_bus(a2(isoil, i - 1):z2(isoil, i - 1)) = svs_mesh%vs%isoil(:, i)
     
        end do
        do i = 1, svs_mesh%vs%kthermal
            if (allocated(svs_mesh%vs%tground)) svs_bus(a2(tground, i - 1):z2(tground, i - 1)) = svs_mesh%vs%tground(:, i)
        end do
        do i = 0, 1
            if (allocated(svs_mesh%vs%tvege)) svs_bus(a2(tvege, i):z2(tvege, i)) = svs_mesh%vs%tvege(:, i + 1)
        end do

        if(svs_mesh%vs%schmsol=='SVS2') then
           do i = 1, nl_svs
                if (allocated(svs_mesh%vs%tpsoil))  svs_bus(a2(tpsoil, i - 1):z2(tpsoil, i - 1)) = svs_mesh%vs%tpsoil(:, i)
                if (allocated(svs_mesh%vs%tpsoilv)) svs_bus(a2(tpsoilv, i - 1):z2(tpsoilv, i - 1)) = svs_mesh%vs%tpsoilv(:, i)
           end do
           if (allocated(svs_mesh%vs%tperm)) svs_bus(a1(tperm):z1(tperm)) = svs_mesh%vs%tperm
        endif

        if(svs_mesh%vs%schmsol=='SVS' .and. svs_mesh%vs%lsoil_freezing_svs1) then
           do i = 1, nl_svs
                if (allocated(svs_mesh%vs%tpsoil))  svs_bus(a2(tpsoil, i - 1):z2(tpsoil, i - 1)) = svs_mesh%vs%tpsoil(:, i)
           end do
           if (allocated(svs_mesh%vs%tperm)) svs_bus(a1(tperm):z1(tperm)) = svs_mesh%vs%tperm
        endif

        if(svs_mesh%vs%schmsol=='SVS' .and. svs_mesh%vs%lwater_ponding_svs1) then
           if (allocated(svs_mesh%vs%watpond)) svs_bus(a1(watpond):z1(watpond)) = svs_mesh%vs%watpond
        endif

        if (allocated(svs_mesh%vs%wveg)) svs_bus(a1(wveg):z1(wveg)) = svs_mesh%vs%wveg

        ! Snow initialisation
        if(svs_mesh%vs%schmsol=='SVS') then
          if (allocated(svs_mesh%vs%snodpl)) svs_bus(a1(snodpl):z1(snodpl)) = svs_mesh%vs%snodpl
          if (allocated(svs_mesh%vs%snoden)) svs_bus(a1(snoden):z1(snoden)) = svs_mesh%vs%snoden
          if (allocated(svs_mesh%vs%snoal)) svs_bus(a1(snoal):z1(snoal)) = svs_mesh%vs%snoal
          if (allocated(svs_mesh%vs%wsnow)) svs_bus(a1(wsnow):z1(wsnow)) = svs_mesh%vs%wsnow
          do i = 0, 1
              if (allocated(svs_mesh%vs%tsnow)) svs_bus(a2(tsnow, i):z2(tsnow, i)) = svs_mesh%vs%tsnow(:, i + 1)
          end do
          where (svs_bus(a1(snodpl):z1(snodpl)) == 0.0)
            svs_bus(a1(snoden):z1(snoden)) = 0.0
            svs_bus(a1(snoal):z1(snoal)) = 0.0
            svs_bus(a1(wsnow):z1(wsnow)) = 0.0
            svs_bus(a2(tsnow, 0):z2(tsnow, 0)) = 0.0
            svs_bus(a2(tsnow, 1):z2(tsnow, 1)) = 0.0
          end where
          if (allocated(svs_mesh%vs%snvdp)) svs_bus(a1(snvdp):z1(snvdp)) = svs_mesh%vs%snvdp
          if (allocated(svs_mesh%vs%snvden)) svs_bus(a1(snvden):z1(snvden)) = svs_mesh%vs%snvden
          if (allocated(svs_mesh%vs%snval)) svs_bus(a1(snval):z1(snval)) = svs_mesh%vs%snval
          if (allocated(svs_mesh%vs%wsnv)) svs_bus(a1(wsnv):z1(wsnv)) = svs_mesh%vs%wsnv
          do i = 0, 1
            if (allocated(svs_mesh%vs%tsnowveg)) svs_bus(a2(tsnowveg, i):z2(tsnowveg, i)) = svs_mesh%vs%tsnowveg(:, i + 1)
          end do
          where (svs_bus(a1(snvdp):z1(snvdp)) == 0.0)
            svs_bus(a1(snvden):z1(snvden)) = 0.0
            svs_bus(a1(snval):z1(snval)) = 0.0
            svs_bus(a1(wsnv):z1(wsnv)) = 0.0
            svs_bus(a2(tsnowveg, 0):z2(tsnowveg, 0)) = 0.0
            svs_bus(a2(tsnowveg, 1):z2(tsnowveg, 1)) = 0.0
          end where


        else if(svs_mesh%vs%schmsol=='SVS2') then
            do i = 1, svs_mesh%vs%nsl
                   if (allocated(svs_mesh%vs%snoage_svs)) svs_bus(a2(snoage_svs, i - 1):z2(snoage_svs, i - 1)) = svs_mesh%vs%snoage_svs(:, i)
                   if (allocated(svs_mesh%vs%snoden_svs)) svs_bus(a2(snoden_svs, i - 1):z2(snoden_svs, i - 1)) = svs_mesh%vs%snoden_svs(:, i)
                   if (allocated(svs_mesh%vs%snoma_svs)) svs_bus(a2(snoma_svs, i - 1):z2(snoma_svs, i - 1)) = svs_mesh%vs%snoma_svs(:, i)
                   if (allocated(svs_mesh%vs%snodiamopt_svs)) svs_bus(a2(snodiamopt_svs, i - 1):z2(snodiamopt_svs, i - 1)) = svs_mesh%vs%snodiamopt_svs(:, i)
                   if (allocated(svs_mesh%vs%snohist_svs)) svs_bus(a2(snohist_svs, i - 1):z2(snohist_svs, i - 1)) = svs_mesh%vs%snohist_svs(:, i)
                   if (allocated(svs_mesh%vs%snospheri_svs)) svs_bus(a2(snospheri_svs, i - 1):z2(snospheri_svs, i - 1)) = svs_mesh%vs%snospheri_svs(:, i)
                   if (allocated(svs_mesh%vs%tsnow_svs)) svs_bus(a2(tsnow_svs, i - 1):z2(tsnow_svs, i - 1)) = svs_mesh%vs%tsnow_svs(:, i)
                   if (allocated(svs_mesh%vs%wsnow_svs)) svs_bus(a2(wsnow_svs, i - 1):z2(wsnow_svs, i - 1)) = svs_mesh%vs%wsnow_svs(:, i)


                   if (allocated(svs_mesh%vs%snoagev_svs)) svs_bus(a2(snoagev_svs, i - 1):z2(snoagev_svs, i - 1)) = svs_mesh%vs%snoagev_svs(:, i)
                   if (allocated(svs_mesh%vs%snodenv_svs)) svs_bus(a2(snodenv_svs, i - 1):z2(snodenv_svs, i - 1)) = svs_mesh%vs%snodenv_svs(:, i)
                   if (allocated(svs_mesh%vs%snomav_svs)) svs_bus(a2(snomav_svs, i - 1):z2(snomav_svs, i - 1)) = svs_mesh%vs%snomav_svs(:, i)
                   if (allocated(svs_mesh%vs%snodiamoptv_svs)) svs_bus(a2(snodiamoptv_svs, i - 1):z2(snodiamoptv_svs, i - 1)) = svs_mesh%vs%snodiamoptv_svs(:, i)
                   if (allocated(svs_mesh%vs%snohistv_svs)) svs_bus(a2(snohistv_svs, i - 1):z2(snohistv_svs, i - 1)) = svs_mesh%vs%snohistv_svs(:, i)
                   if (allocated(svs_mesh%vs%snospheriv_svs)) svs_bus(a2(snospheriv_svs, i - 1):z2(snospheriv_svs, i - 1)) = svs_mesh%vs%snospheriv_svs(:, i)
                   if (allocated(svs_mesh%vs%tsnowv_svs)) svs_bus(a2(tsnowv_svs, i - 1):z2(tsnowv_svs, i - 1)) = svs_mesh%vs%tsnowv_svs(:, i)
                   if (allocated(svs_mesh%vs%wsnowv_svs)) svs_bus(a2(wsnowv_svs, i - 1):z2(wsnowv_svs, i - 1)) = svs_mesh%vs%wsnowv_svs(:, i)
            end do
            do i = 1, svs_mesh%vs%nsl
               where (svs_bus(a2(snoma_svs,i-1):z2(snoma_svs,i-1)) == 0.0)
                 svs_bus(a2(snoage_svs, i-1):z2(snoage_svs, i-1))  = 0. 
                 svs_bus(a2(snodiamopt_svs, i-1):z2(snodiamopt_svs, i-1))  = 0. 
                 svs_bus(a2(snospheri_svs, i-1):z2(snospheri_svs, i-1))  = 0. 
                 svs_bus(a2(snohist_svs, i-1):z2(snohist_svs, i-1))  = 0. 
                 svs_bus(a2(snoden_svs, i-1):z2(snoden_svs, i-1))  = 50.
                 svs_bus(a2(tsnow_svs, i-1):z2(tsnow_svs, i-1))  = 0. 
                 svs_bus(a2(wsnow_svs, i-1):z2(wsnow_svs, i-1))  = 0. 
              end where

              where (svs_bus(a2(snomav_svs,i-1):z2(snomav_svs,i-1)) == 0.0)
                 svs_bus(a2(snoagev_svs, i-1):z2(snoagev_svs, i-1))  = 0. 
                 svs_bus(a2(snodiamoptv_svs, i-1):z2(snodiamoptv_svs, i-1))  = 0. 
                 svs_bus(a2(snospheriv_svs, i-1):z2(snospheriv_svs, i-1))  = 0. 
                 svs_bus(a2(snohistv_svs, i-1):z2(snohistv_svs, i-1))  = 0. 
                 svs_bus(a2(snodenv_svs, i-1):z2(snodenv_svs, i-1))  = 50.
                 svs_bus(a2(tsnowv_svs, i-1):z2(tsnowv_svs, i-1))  = 0. 
                 svs_bus(a2(wsnowv_svs, i-1):z2(wsnowv_svs, i-1))  = 0. 
              end where

            end do



        end if


    end subroutine


    subroutine runsvs_mesh_init(shd, fls, cm)

        !> For RPN/physics status.
        use phy_status, only: phy_error_L

        !> For surface layer configuration.
        use sfclayer_mod, only: sl_put, SL_OK

        !> For rmnlib constant 'RMN_IS_OK'.
#include <rmnlib_basics.hf>

        type(ShedGridParams) :: shd
        type(fl_ids) :: fls
        type(clim_info) :: cm

        !> Local variables.
        character(len = DEFAULT_LINE_LENGTH) line
        character(len = DEFAULT_LINE_LENGTH) level
        character(len = DEFAULT_FIELD_LENGTH) code
!-        integer, allocatable :: rg_soil(:)
        integer :: j, i,  ierr, moyhr = 0

        !> Return if the process is not marked active.
        if (.not. svs_mesh%PROCESS_ACTIVE) then
            return
        else
            call print_new_section("RUNSVS is active: " // svs_mesh%vs%schmsol)
            call increase_tab()
        end if

        !> Check for required variables.
        ierr = 0
        if (.not. associated(vs%tile%fsin)) then
            call print_error("The driving variable '" // VN_FSIN // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. associated(vs%tile%flin)) then
            call print_error("The driving variable '" // VN_FLIN // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. associated(vs%tile%ta)) then
            call print_error("The driving variable '" // VN_TA // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. associated(vs%tile%qa)) then
            call print_error("The driving variable '" // VN_QA // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. associated(vs%tile%pres)) then
            call print_error("The driving variable '" // VN_PRES // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. associated(vs%tile%uv)) then
            call print_error("The driving variable '" // VN_UV // "' is not active or not associated with an input file.")
            ierr = 1
        end if
        if (.not. (associated(vs%tile%prern) .and. associated(vs%tile%presno)) .and. .not. associated(vs%tile%pre)) then
            call print_error( &
                "No driving variable for precipitation is active nor associated with an input file. The '" // VN_PRE // &
                "' variable or both the '" // VN_PRERN // "' and '" // VN_PRESNO // "' variables are required.")
            ierr = 1
        else if (associated(vs%tile%prern) .and. associated(vs%tile%presno) .and. associated(vs%tile%pre)) then
            call print_info( &
                "The '" // VN_PRERN // "' and '" // VN_PRESNO // "' variables are active. The '" // VN_PRE // &
                "' variable is also active but inputs on the field are not being used.")
        end if
        if (ierr /= 0) then
            call reset_tab()
            call print_error( &
                "The variables required to drive the model are not active or have not been associated with an input file.")
            call program_abort()
        end if

        !> Initialize surface options.
        ierr = sfc_options_init()
        if (.not. RMN_IS_OK(ierr)) then
            call print_error("An error occurred initializing the module.")
            call program_abort()
        end if

        !> Transfer the model time-step.
        time_dt = real(ic%dts)

        !> Set the number of active soil indices to the number of active land-tiles.
        ni = shd%lc%NML

        !> Set the surface scheme 
        schmsol = svs_mesh%vs%schmsol

        !> Note that model setup is equivalent to using external meteorological driving data ('offline' mode).
        atm_external = .true.

        !> Deactivate 'tplus' to use 'moins' variables.
        atm_tplus = .false.

        !> Deactivate 'radslope' to use 'flusolis'.
        radslope = .false.

        !> Deactivate 'svs_local_z0m' to use 'z0'.
!        svs_local_z0m = .false.

        !> Set the number of active surface layers to maximum of the surface IDs used in 'inisurf' checks.
        nsurf = indx_max !max(indx_soil, indx_water)

        !> User-set settings.
        use_photo = .true.
        use_eff_surf_tq = .true.
        sl_func_stab = 'BELJAARS91'
        sl_z0ref = .true.
        sl_lmin_soil = -1.0
        sl_lmin_glacier = 10.0
        sl_lmin_water = 10.0
        sl_lmin_seaice = 10.0
        read_emis = .false.
        limsnodp = .true.
        icemelt = .true.
        icelac = .false.
        diusst = 'FAIRALL'
        diusst_warmlayer = .true.
        diusst_coolskin = .true.
        diusst_warmlayer_lakes = .true.
        diusst_coolskin_lakes = .true.
        z0mtype = 'BELJAARS'
        z0ttype = 'DEACU12'
        salty_qsat = .true.
        urban_params_new = .true.
        kount_reset = 12

        !> Update the number of active surface layers for the physics bus.
!        nagrege = nsurf + 1

        !> Surface layer configuration.
        ierr = SL_OK
        if (ierr == SL_OK) ierr = sl_put('beta', beta)
        if (ierr == SL_OK) ierr = sl_put('rineutral', sl_rineutral)
        if (ierr == SL_OK) ierr = sl_put('tdiaglim', tdiaglim)
        if (ierr == SL_OK) ierr = sl_put('sl_stabfunc_stab', sl_func_stab)
        if (ierr == SL_OK) ierr = sl_put('sl_stabfunc_unstab', sl_func_unstab)
        if (ierr == SL_OK) ierr = sl_put('z0ref', sl_z0ref)
        if (ierr /= SL_OK) then
            call print_error("An error occurred configuring the surface layer module.")
            call program_abort()
        endif

        !> Transfer the soil dimension and soil mapping type.
        soiltext = svs_mesh%vs%soiltext
        if (soiltext == 'NIL') then
            nl_svs = shd%lc%IGND
            allocate(dl_svs(shd%lc%IGND), source = shd%lc%sl%zbot)
        else
            nl_svs = NL_SVS_DEFAULT
            allocate(dl_svs(NL_SVS_DEFAULT), source = DP_SVS_DEFAULT)
        end if

        !> Vegetation types (from 'sfc_nml').
        if (vf_type == "CCILCECO") then
            ntypel = 11
            ntypeh = 10
            allocate(vl_type(ntypel), vh_type(ntypeh))
            vl_type = (/10, 11, 12, 13, 14, 15, 16, 17, 20, 22, 23/)
            vh_type = (/4, 5, 6, 7, 8, 9, 18, 19, 25, 26/)
        else
            ntypel = 13
            ntypeh = 8
            allocate(vl_type(ntypel), vh_type(ntypeh))
            vl_type = (/10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 22, 23/)
            vh_type = (/4, 5, 6, 7, 8, 9, 25, 26/)
        end if

        ! Activate or not the soil freezing module in SVS1
        if(svs_mesh%vs%schmsol=='SVS') then
                lsoil_freezing_svs1 = svs_mesh%vs%lsoil_freezing_svs1
        endif

        ! Activate or not the water ponding module in SVS1
        if(svs_mesh%vs%schmsol=='SVS') then
                lwater_ponding_svs1 = svs_mesh%vs%lwater_ponding_svs1
        endif

        ierr =0
        ! Initialize number of snow layers (for multilayer snowpack schemes in SVS2)
        if(svs_mesh%vs%schmsol=='SVS2') then 
             nsl  = svs_mesh%vs%nsl
             hsnowscheme =  svs_mesh%vs%hsnowscheme
             hsnowdrift_cro = svs_mesh%vs%hsnowdrift_cro
             lsnowdrift_sublim = svs_mesh%vs%lsnowdrift_sublim
             hsnowcomp =  svs_mesh%vs%hsnowcomp
             hsnowcond =  svs_mesh%vs%hsnowcond
             hsnowrad =  svs_mesh%vs%hsnowrad
             hsnowfall =  svs_mesh%vs%hsnowfall
             hsnowhold =  svs_mesh%vs%hsnowhold
             hsnowmetamo =  svs_mesh%vs%hsnowmetamo

             if (.not.any(hsnowscheme == hsnowscheme_opt)) then
                 call str_concat(msg_S,hsnowscheme_opt ,', ')
                 call print_error('hsnowscheme = '//trim(hsnowscheme)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif
             if (.not.any(hsnowcomp == hsnowcomp_opt)) then
                 call str_concat(msg_S, hsnowcomp_opt,', ')
                 call print_error('hsnowcomp = '//trim(hsnowcomp)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif
             if (.not.any(hsnowmetamo == hsnowmetamo_opt)) then
                 call str_concat(msg_S, hsnowmetamo_opt,', ')
                 call print_error('hsnowmetamo = '//trim(hsnowmetamo)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif
             if (.not.any(hsnowdrift_cro == hsnowdrift_cro_opt)) then
                 call str_concat(msg_S, hsnowdrift_cro_opt,', ')
                 call print_error('hsnowdrift_cro = '//trim(hsnowdrift_cro)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif
             if (.not.any(hsnowcond == hsnowcond_opt)) then
                 call str_concat(msg_S, hsnowcond_opt,', ')
                 call print_error('hsnowcond = '//trim(hsnowcond)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif
             if (.not.any(hsnowrad == hsnowrad_opt)) then
                 call str_concat(msg_S, hsnowrad_opt,', ')
                 call print_error('hsnowrad = '//trim(hsnowrad)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif
             if (.not.any(hsnowfall == hsnowfall_opt)) then
                 call str_concat(msg_S, hsnowfall_opt,', ')
                 call print_error('hsnowfall = '//trim(hsnowfall)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif
             if (.not.any(hsnowhold == hsnowhold_opt)) then
                 call str_concat(msg_S, hsnowhold_opt,', ')
                 call print_error('hsnowhold = '//trim(hsnowhold)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif

            if (ierr /= 0) then
                call program_abort()
            end if

        endif

        write(*, nml = surface_cfgs)

        !> Initialize the physics bus.
        call phy_businit(ni, nk)

        !> Initialize the surface bus in the physics library.
        call sfc_businit(moyhr, ni, nk)
        if (phy_error_L) then
            call print_error("An error occurred initializing the surface bus variable.")
            call program_abort()
        end if

        !> Initialize the surface variable pointers.
        ierr = sfcbus_init()
        if (.not. RMN_IS_OK(ierr)) then
            call print_error("An error occurred initializing the surface variable pointers.")
            call program_abort()
        end if

        !> Inialize the empty variable list.
        phyinread_list_S = ''
        phyinread_n = 0

        !> Set the 'phyinread_list' list to emulate if the fields were read from file.
        call runsvs_mesh_append_phyentvar('vegf')
        if (allocated(svs_mesh%vs%lnz0)) then
            call runsvs_mesh_append_phyentvar('z0en')
            call runsvs_mesh_append_phyentvar('z0mlanden')
        end if
        call runsvs_mesh_append_phyentvar('slop')
        call runsvs_mesh_append_phyentvar('draindens')
        call runsvs_mesh_append_phyentvar('wsoil')
        call runsvs_mesh_append_phyentvar('isoil')
        call runsvs_mesh_append_phyentvar('tground')
        call runsvs_mesh_append_phyentvar('tvege')
        call runsvs_mesh_append_phyentvar('wveg')
        call runsvs_mesh_append_phyentvar('tsnow')
        call runsvs_mesh_append_phyentvar('snodpl')
        call runsvs_mesh_append_phyentvar('snoden')
        call runsvs_mesh_append_phyentvar('snoal')
        call runsvs_mesh_append_phyentvar('wsnow')
        call runsvs_mesh_append_phyentvar('tsnowveg')
        call runsvs_mesh_append_phyentvar('snvdp')
        call runsvs_mesh_append_phyentvar('snvden')
        call runsvs_mesh_append_phyentvar('snval')
        call runsvs_mesh_append_phyentvar('wsnv')
        if(svs_mesh%vs%schmsol=='SVS' .and. svs_mesh%vs%lsoil_freezing_svs1) then
            call runsvs_mesh_append_phyentvar('tpsoil')
            call runsvs_mesh_append_phyentvar('tperm')
        endif
        if(svs_mesh%vs%schmsol=='SVS' .and. svs_mesh%vs%lwater_ponding_svs1) then
            call runsvs_mesh_append_phyentvar('watpond')
        endif

        if(svs_mesh%vs%schmsol=='SVS2') then
             call runsvs_mesh_append_phyentvar('tpsoil')
             call runsvs_mesh_append_phyentvar('tpsoilv')
             call runsvs_mesh_append_phyentvar('tperm')

             call runsvs_mesh_append_phyentvar('snoden_svs')
             call runsvs_mesh_append_phyentvar('snoage_svs')
             call runsvs_mesh_append_phyentvar('snodiamopt_svs')
             call runsvs_mesh_append_phyentvar('snoma_svs')
             call runsvs_mesh_append_phyentvar('snospheri_svs')
             call runsvs_mesh_append_phyentvar('snohist_svs')
             call runsvs_mesh_append_phyentvar('tsnow_svs')
             call runsvs_mesh_append_phyentvar('wsnow_svs')

             call runsvs_mesh_append_phyentvar('snodenv_svs')
             call runsvs_mesh_append_phyentvar('snoagev_svs')
             call runsvs_mesh_append_phyentvar('snodiamoptv_svs')
             call runsvs_mesh_append_phyentvar('snomav_svs')
             call runsvs_mesh_append_phyentvar('snospheriv_svs')
             call runsvs_mesh_append_phyentvar('snohistv_svs')
             call runsvs_mesh_append_phyentvar('tsnowv_svs')
             call runsvs_mesh_append_phyentvar('wsnowv_svs')

        end if


        !> Update soil indices after the call to 'sfc_businit', which calls 'init_soil_text_levels'.
        if (soiltext == 'NIL') then

            !> Transfer soil variables.
            nl_stp = shd%lc%IGND
            allocate(weights(nl_svs, nl_stp))
            weights = 0.0
            do i = 1, nl_stp
                weights(i, i) = 1.0
            end do

            !> Overwrite the default levels set by the unknown 'soiltext' type.
            vd%sand%mul = nl_stp
            vl(vd%sand%i)%mul = nl_stp
            vd%clay%mul = nl_stp
            vl(vd%clay%i)%mul = nl_stp
        else

            !> Overwrite the default input level set by the unknown 'soiltext' type.
!-            nl_ste = shd%lc%IGND

            !> Add the 'sanden' and 'clayen' variables as inputs.
            call runsvs_mesh_append_phyentvar('sanden')
            call runsvs_mesh_append_phyentvar('clayen')

            !> Overwrite active soil layers for MESH.
            shd%lc%IGND = nl_svs
        end if
        khyd = svs_mesh%vs%khyd

        !> Required to activate the snow-related checks in 'coherence'.
        call runsvs_mesh_append_phyentvar('snodp')

        !> Create the range for soils.
!-        allocate(rg_soil(ni))
!-        rg_soil = 0
!-        do i = 1, ni
!-            rg_soil(i) = i
!-        end do

        !> Count the number of active 'soil' variables and build the surface pointer.
        allocate(bus_ptr(nvarsurf))
        bus_ptr = 1
        bus_length = 1
        do i = 1, nvarsurf

            !> Overrides to accommodate missing inialization of physics bus.
            vl(i)%niveaux = max(vl(i)%niveaux, 1)
            vl(i)%mul = max(vl(i)%mul, 1)
            vl(i)%mosaik = max(vl(i)%mosaik, 1)
print*,vl(i)%n,vl(i)%niveaux,vl(i)%mul,vl(i)%mosaik

            !> Increment the index count.
            bus_ptr(i) = bus_length
            bus_length = bus_length + vl(i)%niveaux*vl(i)%mul*vl(i)%mosaik*ni
        end do
        vd = transfer(vl, vd)
        bus_length = bus_length - 1

        !> Build the variable bus.
        if (bus_length > 1) then

            !> Allocate the bus.
            allocate(phy_bus(bus_length, trnch))
            phy_bus = 0.0

            !> Manually assign the surface pointer.
            do i = 1, nvarsurf
                if (associated(busptr(i)%ptr)) then
                    nullify(busptr(i)%ptr)
                end if
                busptr(i)%ptr => phy_bus(bus_ptr(i):(bus_ptr(i) + vl(i)%niveaux*vl(i)%mul*vl(i)%mosaik*ni - 1), :)
            end do

            !> Associate the 1D and 2D bus variables.
            svs_bus => phy_bus(:, trnch)
        else
            call print_error("There are not compatible surface tiles in the domain.")
            call program_abort()
        end if
        if (DIAGNOSEMODE) then
            write(code, *) nvarsurf
            call print_info(trim(adjustl(code)) // " internal surface variables are active.")
        end if

        !> Update bus variable.
        call runsvs_mesh_copy_vs_to_bus()

        !> Diagnostic summary of inputs at the first tile.
        if (DIAGNOSEMODE) then
            call reset_tab()
            call print_new_section('--------------------------------')
            call print_message('SVS DIAGNOSTICS')
            call print_message('--------------------------------')
            write(line, "('TILE:             ', i8)") 1
            call print_message(line)
            call print_message('--------------------------------')
            write(line, "('LATITUDE:         ', f10.1)") svs_bus(a1(dlat))*rad2deg
            call print_message(line)
            write(line, "('LONGITUDE:        ', f10.1)") svs_bus(a1(dlon))*rad2deg
            call print_message(line)
            call print_message('--------------------------------')
            write(line, "('ROUGHNESS LENGTH: ', f8.3)") svs_bus(a1(z0))
            call print_message(line)
            write(line, "('VEGETATION TEMP.: ', 2f8.3)") svs_bus(a1(tvege)), svs_bus(a1(tvege) + ni)
            call print_message(line)
            call print_message('VEGETATION COVER:')
            do i = 1199, 1174, -1
                write(line, "('% ', i5, '        ', f8.3)") i, svs_bus(a2(vegf, 1199 - i))*100.0
                call print_message(line)
            end do
            call print_message('--------------------------------')
            if (svs_mesh%vs%observed_forcing) then
                write(line, "('FORCING LEVEL:    ', (a))") 'height'
                call print_message(line)
                write(line, "(' THERMO. HEIGHT:   ', f8.3)") svs_bus(a1(ztsl))
                call print_message(line)
                write(line, "(' MOMENTUM HEIGHT:  ', f8.3)") svs_bus(a1(zusl))
                call print_message(line)
            else
                write(line, "('FORCING LEVEL:    ', (a))") 'sigma'
                call print_message(line)
                write(line, "(' THERMO. SIGMA:    ', f8.3)") svs_mesh%vs%sigma_t
                call print_message(line)
                write(line, "(' MOMENTUM SIGMA:   ', f8.3)") svs_mesh%vs%sigma_u
                call print_message(line)
            end if
            call print_message('--------------------------------')
            write(line, "('SLOPE:            ', f8.3)") svs_bus(a1(slop))
            call print_message(line)
            write(line, "('DRAIN.DENSITY     ', f8.3)") svs_bus(a1(draindens))
            call print_message(line)
            call print_message('--------------------------------')
            call print_message('SOIL MAPPING:')
            call print_message('DATABASE: ' // trim(soiltext))
            call print_message('WEIGHTS [METERS]:')
            do i = 1, nl_svs ! model layers
                write(line, "(' LAYER ', i3, ' DEPTH: ', f8.3)") i, dl_svs(i)
                call print_message(line)
                do j = 1, nl_stp ! database layers
                    if (soiltext == 'GSDE') then
                        write(line, "('  ', (a), ' DEPTH: ', f8.3, ' WEIGHT: ', f8.3)") 'DB', dl_gsde(j), weights(i, j)
                    else if (soiltext == 'SLC') then
                        write(line, "('  ', (a), ' DEPTH: ', f8.3, ' WEIGHT: ', f8.3)") 'DB', dl_slc(j), weights(i, j)
                    else if (soiltext == 'SOILGRIDS') then
                        write(line, "('  ', (a), ' DEPTH: ', f8.3, ' WEIGHT: ', f8.3)") 'DB', dl_soilgrids(j), weights(i, j)
                    else if (soiltext == 'NIL') then
                        write(line, "('  ', (a), ' DEPTH: ', f8.3, ' WEIGHT: ', f8.3)") 'DB', dl_svs(j), weights(i, j)
                    end if
                    call print_message(line)
                end do
            end do
            write(line, "('PERMEABLE LAYERS: ', i3)") khyd
            call print_message('SOIL TEXTURE:')
            call print_message('             % SAND    % CLAY')
            do i = 1, nl_svs ! model layers
                write(line, "(' LAYER ', i3, ': ', 999(f8.3, 2x))") i, svs_bus(a2(sand, i - 1)), svs_bus(a2(clay, i - 1))
                call print_message(line)
            end do
            call print_message('SOIL MOISTURE:')
            call print_message('             LIQUID    FROZEN')
            do i = 1, nl_svs ! permeable layers
                write(line, "(' LAYER ', i3, ': ', 999(f8.3, 2x))") i, svs_bus(a2(wsoil, i - 1)), svs_bus(a2(isoil, i - 1))
                call print_message(line)
            end do
            if(svs_mesh%vs%schmsol=='SVS') then
                write(line, "('SOIL TEMPERATURE: ', 2f8.3)") svs_bus(a1(tground)), svs_bus(a1(tground) + ni)
                call print_message(line)
                if(svs_mesh%vs%lsoil_freezing_svs1) then
                   call print_message('             Soil profile ')
                   do i = 1, nl_svs ! permeable layers
                      write(line, "(' LAYER ', i3, ': ', 999(f8.3, 1x))") i, svs_bus(a2(tpsoil, i - 1))
                      call print_message(line)
                   end do
                endif
            else if(svs_mesh%vs%schmsol=='SVS2') then
               call print_message('             Bare ground/low veg    High veg.')
               do i = 1, nl_svs ! permeable layers
                   write(line, "(' LAYER ', i3, ': ', 999(f8.3, 2x))") i, svs_bus(a2(tpsoil, i - 1)), svs_bus(a2(tpsoilv, i - 1))
                   call print_message(line)
               end do
            end if
            if(svs_mesh%vs%schmsol=='SVS2') then 
               call print_message('             Bare ground/low veg    High veg.')
               do i = 1, nsl ! snow
                   write(line, "(' LAYER ', i3, ': ', 999(f8.3, 2x))") i, svs_bus(a2(snoden_svs, i - 1)), svs_bus(a2(snodenv_svs, i - 1))
                   call print_message(line)
               end do
            end if

            call print_message('--------------------------------')
            call print_message('GROUND/LOW VEG. SNOW:')
            write(line, "(' SNOW TEMPERATURE:', 2f8.3)") svs_bus(a1(tsnow)), svs_bus(a1(tsnow) + ni)
            call print_message(line)
            write(line, "(' SNOW DEPTH:      ', 2f8.3)") svs_bus(a1(snodpl))
            call print_message(line)
            write(line, "(' SNOW DENSITY:    ', 2f8.3)") svs_bus(a1(snoden))
            call print_message(line)
            write(line, "(' SNOW ALBEDO:     ', 2f8.3)") svs_bus(a1(snoal))
            call print_message(line)
            write(line, "(' SNOW W/C:        ', 2f8.3)") svs_bus(a1(wsnow))
            call print_message(line)
            call print_message('HIGH VEG. SNOW:')
            write(line, "(' SNOW TEMPERATURE:', 2f8.3)") svs_bus(a1(tsnowveg)), svs_bus(a1(tsnowveg) + ni)
            call print_message(line)
            write(line, "(' SNOW DEPTH:      ', 2f8.3)") svs_bus(a1(snvdp))
            call print_message(line)
            write(line, "(' SNOW DENSITY:    ', 2f8.3)") svs_bus(a1(snvden))
            call print_message(line)
            write(line, "(' SNOW ALBEDO:     ', 2f8.3)") svs_bus(a1(snval))
            call print_message(line)
            write(line, "(' SNOW W/C:        ', 2f8.3)") svs_bus(a1(wsnv))
            call print_message(line)
            call print_message('--------------------------------')
        end if



ierr = 200
!!do i = 1, 26
!!write(code, *) i
!!open(ierr, file = 'output/2V_' // trim(adjustl(code)) // '.txt')
!!write(ierr, '(a)') '2V_' // trim(adjustl(code)) // ' '; ierr = ierr + 1
!!end do
!open(ierr, file = 'output/ACWF.txt'); write(ierr, '(a)') 'ACWF_0 '; ierr = ierr + 1
!open(ierr, file = 'output/AHFL.txt'); write(ierr, '(a)') 'AHFL_0 '; ierr = ierr + 1
!do i = 1, nl_svs
!write(code, *) i
!open(ierr, file = 'output/WSOL_' // trim(adjustl(code)) // '.txt')
!write(ierr, '(a)') 'WSOL_' // trim(adjustl(code)) // ' '; ierr = ierr + 1
!end do
!open(ierr, file = 'output/FL.txt'); write(ierr, '(a)') 'FL_0 '; ierr = ierr + 1
!open(ierr, file = 'output/FQ.txt'); write(ierr, '(a)') 'FQ_0 '; ierr = ierr + 1
!open(ierr, file = 'output/FV_1.txt'); write(ierr, '(a)') 'FV_1 '; ierr = ierr + 1
!open(ierr, file = 'output/EFLX.txt'); write(ierr, '(a)') 'EFLX_0 '; ierr = ierr + 1
!open(ierr, file = 'output/EG.txt'); write(ierr, '(a)') 'EG_0 '; ierr = ierr + 1
!open(ierr, file = 'output/ER.txt'); write(ierr, '(a)') 'ER_0 '; ierr = ierr + 1
!open(ierr, file = 'output/ETR.txt'); write(ierr, '(a)') 'ETR_0 '; ierr = ierr + 1
!open(ierr, file = 'output/FC_1.txt'); write(ierr, '(a)') 'FC_1 '; ierr = ierr + 1
!open(ierr, file = 'output/RSEF.txt'); write(ierr, '(a)') 'RSEF_0 '; ierr = ierr + 1
!open(ierr, file = 'output/Z0_1.txt'); write(ierr, '(a)') 'Z0_1 '; ierr = ierr + 1
!open(ierr, file = 'output/Z0HA.txt'); write(ierr, '(a)') 'Z0HA_0 '; ierr = ierr + 1
!open(ierr, file = 'output/Z0VH.txt'); write(ierr, '(a)') 'Z0VH_0 '; ierr = ierr + 1
!open(ierr, file = 'output/Z0VL.txt'); write(ierr, '(a)') 'Z0VL_0 '; ierr = ierr + 1
!open(ierr, file = 'output/ZH.txt'); write(ierr, '(a)') 'ZH_0 '; ierr = ierr + 1
!open(ierr, file = 'output/ZT_1.txt'); write(ierr, '(a)') 'ZT_1 '; ierr = ierr + 1
!open(ierr, file = 'output/ZU.txt'); write(ierr, '(a)') 'ZU_0 '; ierr = ierr + 1
!open(ierr, file = 'output/RSGR.txt'); write(ierr, '(a)') 'RSGR_0 '; ierr = ierr + 1
!open(ierr, file = 'output/RSVG.txt'); write(ierr, '(a)') 'RSVG_0 '; ierr = ierr + 1
!
!open(ierr, file = 'output/SNAL_0.txt'); write(ierr, '(a)') 'SNAL_0 '; ierr = ierr + 1
!open(ierr, file = 'output/SNDN_0.txt'); write(ierr, '(a)') 'SNDN_0 '; ierr = ierr + 1
!open(ierr, file = 'output/SNDP_0.txt'); write(ierr, '(a)') 'SNDP_0 '; ierr = ierr + 1
!open(ierr, file = 'output/SVAL_0.txt'); write(ierr, '(a)') 'SVAL_0 '; ierr = ierr + 1
!open(ierr, file = 'output/SVDN_0.txt'); write(ierr, '(a)') 'SVDN_0 '; ierr = ierr + 1
!open(ierr, file = 'output/SVDP_0.txt'); write(ierr, '(a)') 'SVDP_0 '; ierr = ierr + 1
!do i = 1, 2
!write(code, *) i
!open(ierr, file = 'output/TSN_' // trim(adjustl(code)) // '.txt')
!write(ierr, '(a)') 'TSN_' // trim(adjustl(code)) // ' '; ierr = ierr + 1
!end do
!do i = 1, 2
!write(code, *) i
!open(ierr, file = 'output/TSNV_' // trim(adjustl(code)) // '.txt')
!write(ierr, '(a)') 'TSNV_' // trim(adjustl(code)) // ' '; ierr = ierr + 1
!end do

    ! Prep SVS output files


    if(svs_mesh%vs%schmsol=='SVS2') then

       open(iout_soil, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_soil_hourly.csv', action = 'write')
       write(iout_soil, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       do j = 1, nl_svs
           write(level, FMT_GEN) j
           write(iout_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_ISOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_WSOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_TPSOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_TPSOILV) // '_' // trim(adjustl(level))
       end do
       write(iout_soil, *)


       open(iout_snow_bulk, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_snow_bulk_hourly.csv', action = 'write')
       write(iout_snow_bulk, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       write(iout_snow_bulk, FMT_CSV, advance = 'no') 'SNOMA', 'SNODP','SNODEN','SNOALB','WSNO','TSNO_SURF','RSNOW_AC','RAINRATE', 'SNOWRATE' 
       write(iout_snow_bulk, *)

       if(svs_mesh%vs%lout_snow_enbal) then 
           open(iout_snow_enbal, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_snow_enbal_hourly.csv', action = 'write')
           write(iout_snow_enbal, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
           write(iout_snow_enbal, FMT_CSV, advance = 'no') 'SNO_RNET', 'SNO_SWNET','SNO_LWNET','SNO_LE','SNO_H','SNO_BSUBL','SNO_G', 'SNO_HRAIN', 'SUBLIM_LV'
           write(iout_snow_enbal, *)
       endif

       if(svs_mesh%vs%lout_snow_profile) then 
          open(iout_snow_profile, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_snow_profile_hourly.csv', action = 'write')
          write(iout_snow_profile, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
          do j = 1, nsl
              write(level, FMT_GEN) j
              write(iout_snow_profile, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_SNOMA_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNODEN_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOAGE_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNODIAMOPT_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOSPHERI_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOHIST_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_TSNOW_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_WSNOW_SVS) // '_' // trim(adjustl(level)) 
          end do
          write(iout_snow_profile, *)
       endif


   else if(svs_mesh%vs%schmsol=='SVS') then

       open(iout_svs1_soil, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs1_soil_hourly.csv', action = 'write')
       write(iout_svs1_soil, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       do j = 1, nl_svs
           write(level, FMT_GEN) j
           write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_ISOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_WSOIL) // '_' // trim(adjustl(level))
       end do
       if(svs_mesh%vs%lsoil_freezing_svs1) then 
          do j = 1, nl_svs
             write(level, FMT_GEN) j
             write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_TPSOIL) // '_' // trim(adjustl(level))
         end do
       endif

       write(iout_svs1_soil, FMT_CSV, advance = 'no') 'TGROUND_1','TGROUND_2','TVEG_1','TVEG_2'
       write(iout_svs1_soil, *)

       open(iout_svs1_snow, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs1_snow_bulk_hourly.csv', action = 'write')
       write(iout_svs1_snow, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       write(iout_svs1_snow, FMT_CSV, advance = 'no') 'SNOMA', 'SNODP','SNODEN','SNOALB','WSNO','TSNO_1','TSNO_2'
       write(iout_svs1_snow, FMT_CSV, advance = 'no') 'SNVMA', 'SNVDP','SNVDEN','SNVALB','WSNV','TSNV_1','TSNV_2'
       write(iout_svs1_snow, *)


   endif





    end subroutine

    subroutine runsvs_mesh_resume_states_seq(fls, shd, resume_ts)

        !> MESH modules.
        !*  FLAGS: Required for 'RESUMEFLAG'.
!-        use FLAGS, only: RESUMEFLAG

        !> Input variables.
        type(fl_ids) fls
        type(ShedGridParams) shd
        !> Input variables (optional).
        logical, intent(in), optional :: resume_ts

        !> Return if the process is not marked active.
        if (.not. svs_mesh%PROCESS_ACTIVE) return

    end subroutine

    subroutine runsvs_mesh_copy_bus_to_vs()

        !> Local variables.
        integer i

        !> Functions to derive indices by variable name in the bus.
#define a2(n,l) bus_ptr(vd%n%i)+l*ni
#define a1(n) a2(n,0)
#define z2(n,l) bus_ptr(vd%n%i)+l*ni+ni-1
#define z1(n) z2(n,0)

        !> Allocate unallocated intermediary variables.
        if (.not. allocated(svs_mesh%vs%wsoil)) allocate(svs_mesh%vs%wsoil(ni, nl_svs))
        if (.not. allocated(svs_mesh%vs%isoil)) allocate(svs_mesh%vs%isoil(ni, nl_svs))
        if (.not. allocated(svs_mesh%vs%tground)) allocate(svs_mesh%vs%tground(ni, 2))
        if (.not. allocated(svs_mesh%vs%tvege)) allocate(svs_mesh%vs%tvege(ni, 2))
        if (.not. allocated(svs_mesh%vs%wveg)) allocate(svs_mesh%vs%wveg(ni))
        if (.not. allocated(svs_mesh%vs%snodpl)) allocate(svs_mesh%vs%snodpl(ni))
        if (.not. allocated(svs_mesh%vs%snoden)) allocate(svs_mesh%vs%snoden(ni))
        if (.not. allocated(svs_mesh%vs%snoal)) allocate(svs_mesh%vs%snoal(ni))
        if (.not. allocated(svs_mesh%vs%wsnow)) allocate(svs_mesh%vs%wsnow(ni))
        if (.not. allocated(svs_mesh%vs%tsnow)) allocate(svs_mesh%vs%tsnow(ni, 2))
        if (.not. allocated(svs_mesh%vs%snvdp)) allocate(svs_mesh%vs%snvdp(ni))
        if (.not. allocated(svs_mesh%vs%snvden)) allocate(svs_mesh%vs%snvden(ni))
        if (.not. allocated(svs_mesh%vs%snval)) allocate(svs_mesh%vs%snval(ni))
        if (.not. allocated(svs_mesh%vs%wsnv)) allocate(svs_mesh%vs%wsnv(ni))
        if (.not. allocated(svs_mesh%vs%tsnowveg)) allocate(svs_mesh%vs%tsnowveg(ni, 2))


        !> Assign variables.
        do i = 1, nl_svs
            svs_mesh%vs%wsoil(:, i) = svs_bus(a2(wsoil, i - 1):z2(wsoil, i - 1))
            svs_mesh%vs%isoil(:, i) = svs_bus(a2(isoil, i - 1):z2(isoil, i - 1))
        end do
        do i = 1, svs_mesh%vs%kthermal
            svs_mesh%vs%tground(:, i) = svs_bus(a2(tground, i - 1):z2(tground, i - 1))
        end do
        do i = 0, 1
            svs_mesh%vs%tvege(:, i + 1) = svs_bus(a2(tvege, i):z2(tvege, i))
        end do
        svs_mesh%vs%wveg = svs_bus(a1(wveg):z1(wveg))
        svs_mesh%vs%snodpl = svs_bus(a1(snodpl):z1(snodpl))
        svs_mesh%vs%snoden = svs_bus(a1(snoden):z1(snoden))
        svs_mesh%vs%snoal = svs_bus(a1(snoal):z1(snoal))
        svs_mesh%vs%wsnow = svs_bus(a1(wsnow):z1(wsnow))
        do i = 0, 1
            svs_mesh%vs%tsnow(:, i + 1) = svs_bus(a2(tsnow, i):z2(tsnow, i))
        end do
        svs_mesh%vs%snvdp = svs_bus(a1(snvdp):z1(snvdp))
        svs_mesh%vs%snvden = svs_bus(a1(snvden):z1(snvden))
        svs_mesh%vs%snval = svs_bus(a1(snval):z1(snval))
        svs_mesh%vs%wsnv = svs_bus(a1(wsnv):z1(wsnv))
        do i = 0, 1
            svs_mesh%vs%tsnowveg(:, i + 1) = svs_bus(a2(tsnowveg, i):z2(tsnowveg, i))
        end do

        if(svs_mesh%vs%schmsol=='SVS' .and. svs_mesh%vs%lsoil_freezing_svs1) then 
               if (.not. allocated(svs_mesh%vs%tpsoil)) allocate(svs_mesh%vs%tpsoil(ni, nl_svs))
               do i = 1, nl_svs
                  svs_mesh%vs%tpsoil(:, i) = svs_bus(a2(tpsoil, i - 1):z2(tpsoil, i - 1))
               end do
               if (.not. allocated(svs_mesh%vs%tperm)) allocate(svs_mesh%vs%tperm(ni))
               svs_mesh%vs%tperm = svs_bus(a1(tperm):z1(tperm ))
        end if

        if(svs_mesh%vs%schmsol=='SVS2') then 
            if (.not. allocated(svs_mesh%vs%tpsoil)) allocate(svs_mesh%vs%tpsoil(ni, nl_svs))
            if (.not. allocated(svs_mesh%vs%tpsoilv)) allocate(svs_mesh%vs%tpsoilv(ni, nl_svs))
            if (.not. allocated(svs_mesh%vs%tperm)) allocate(svs_mesh%vs%tperm(ni))

           ! if (.not. allocated(svs_mesh%vs%snoage_svs)) write(*,*) 'Alloc age'
           ! if (.not. allocated(svs_mesh%vs%snoma_svs)) write(*,*) 'Alloc mass'

            if (.not. allocated(svs_mesh%vs%snoage_svs)) allocate(svs_mesh%vs%snoage_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snodiamopt_svs)) allocate(svs_mesh%vs%snodiamopt_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snoma_svs)) allocate(svs_mesh%vs%snoma_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snoden_svs)) allocate(svs_mesh%vs%snoden_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snohist_svs)) allocate(svs_mesh%vs%snohist_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snospheri_svs)) allocate(svs_mesh%vs%snospheri_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%tsnow_svs)) allocate(svs_mesh%vs%tsnow_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%wsnow_svs)) allocate(svs_mesh%vs%wsnow_svs(ni,svs_mesh%vs%nsl))

            if (.not. allocated(svs_mesh%vs%snoagev_svs)) allocate(svs_mesh%vs%snoagev_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snodiamoptv_svs)) allocate(svs_mesh%vs%snodiamoptv_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snomav_svs)) allocate(svs_mesh%vs%snomav_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snodenv_svs)) allocate(svs_mesh%vs%snodenv_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snohistv_svs)) allocate(svs_mesh%vs%snohistv_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%snospheriv_svs)) allocate(svs_mesh%vs%snospheriv_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%tsnowv_svs)) allocate(svs_mesh%vs%tsnowv_svs(ni,svs_mesh%vs%nsl))
            if (.not. allocated(svs_mesh%vs%wsnowv_svs)) allocate(svs_mesh%vs%wsnowv_svs(ni,svs_mesh%vs%nsl))


            do i = 1, nl_svs
               svs_mesh%vs%tpsoil(:, i) = svs_bus(a2(tpsoil, i - 1):z2(tpsoil, i - 1))
               svs_mesh%vs%tpsoilv(:, i) = svs_bus(a2(tpsoilv, i - 1):z2(tpsoilv, i - 1))
            end do
            svs_mesh%vs%tperm = svs_bus(a1(tperm):z1(tperm ))
        
            do i = 1, svs_mesh%vs%nsl
               svs_mesh%vs%snoage_svs(:, i) = svs_bus(a2(snoage_svs, i - 1):z2(snoage_svs, i - 1))
               svs_mesh%vs%snoden_svs(:, i) = svs_bus(a2(snoden_svs, i - 1):z2(snoden_svs, i - 1))
               svs_mesh%vs%snospheri_svs(:, i) = svs_bus(a2(snospheri_svs, i - 1):z2(snospheri_svs, i - 1))
               svs_mesh%vs%snodiamopt_svs(:, i) = svs_bus(a2(snodiamopt_svs, i - 1):z2(snodiamopt_svs, i - 1))
               svs_mesh%vs%snoma_svs(:, i) = svs_bus(a2(snoma_svs, i - 1):z2(snoma_svs, i - 1))
               svs_mesh%vs%snohist_svs(:, i) = svs_bus(a2(snohist_svs, i - 1):z2(snohist_svs, i - 1))
               svs_mesh%vs%tsnow_svs(:, i) = svs_bus(a2(tsnow_svs, i - 1):z2(tsnow_svs, i - 1))
               svs_mesh%vs%wsnow_svs(:, i) = svs_bus(a2(wsnow_svs, i - 1):z2(wsnow_svs, i - 1))
               
               svs_mesh%vs%snoagev_svs(:, i) = svs_bus(a2(snoagev_svs, i - 1):z2(snoagev_svs, i - 1))
               svs_mesh%vs%snodenv_svs(:, i) = svs_bus(a2(snodenv_svs, i - 1):z2(snodenv_svs, i - 1))
               svs_mesh%vs%snospheriv_svs(:, i) = svs_bus(a2(snospheriv_svs, i - 1):z2(snospheriv_svs, i - 1))
               svs_mesh%vs%snodiamoptv_svs(:, i) = svs_bus(a2(snodiamoptv_svs, i - 1):z2(snodiamoptv_svs, i - 1))
               svs_mesh%vs%snomav_svs(:, i) = svs_bus(a2(snomav_svs, i - 1):z2(snomav_svs, i - 1))
               svs_mesh%vs%snohistv_svs(:, i) = svs_bus(a2(snohistv_svs, i - 1):z2(snohistv_svs, i - 1))
               svs_mesh%vs%tsnowv_svs(:, i) = svs_bus(a2(tsnowv_svs, i - 1):z2(tsnowv_svs, i - 1))
               svs_mesh%vs%wsnowv_svs(:, i) = svs_bus(a2(wsnowv_svs, i - 1):z2(wsnowv_svs, i - 1))
            end do


        endif

    end subroutine

    subroutine runsvs_mesh_within_tile(shd, fls, cm)

        !> For 'cmcdate_fromprint' function to convert from string date.
        use cmcdate_mod, only: cmcdate_fromprint

        !> For 'jdate_from_cmc' to convert to 'jdateo' for SVS.
        use mu_jdate_mod, only: jdate_from_cmc

        !> For RPN/physics status.
        use phy_status, only: phy_error_L

        !> For constants.
        use tdpack_const, only: rgasd, grav, cappa, tcdk

        use tdpack, only: fotvt

        type(ShedGridParams) :: shd
        type(fl_ids) :: fls
        type(clim_info) :: cm

        !> Local variable for string date (format: YYYYMMDD.hhmmss).
        character(len = 14) time_run_now

        !> Local variables.
        integer i, idateo, ierr, j, k, istat, nfreq

        real :: tt,hu
        real, dimension(il1:il2) :: tve
        !real, dimension(il1:il2) :: tt
        !real, dimension(il1:il2) :: hu

        !> Return if the process is not marked active.
        if (.not. svs_mesh%PROCESS_ACTIVE) return

        !> Update variables (equivalent to calls to 'phyput_input_param' and 'sfc_get_input_param').
        write(time_run_now, "(i4.4, 2i2.2, '.', 2i2.2)") ic%now%year, ic%now%month, ic%now%day, ic%now%hour, ic%now%mins
        idateo = cmcdate_fromprint(time_run_now)
        jdateo = jdate_from_cmc(idateo)

        !> Reset the initialization periodically (at least daily).
        if (ic%ts_count == 1 .or. (ic%now%hour == kount_reset .and. ic%now%mins == 0)) then
            call runsvs_mesh_copy_vs_to_bus()
            kount = 0
            call inichamp4(kount, trnch, ni, nk)
        end if

        !> Increment 'kount'.
        kount = kount + 1

    ! Write SVS hourly outputs
1010    format(9999(g15.7e2, ','))

        if(svs_mesh%vs%schmsol=='SVS2') then

           !if (ic%now%hour /= ic%next%hour) then !last time-step of hour
           if (ic%now%mins ==0) then! Full hour

              k=1 !>  Identity of the tile (offset relative to node-indexing).

              ! Write file containing soil outputs
              write(iout_soil, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              do i = 1, nl_svs
                 write(iout_soil, FMT_CSV, advance = 'no') &
                     busptr(vd%isoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch) , &
                     busptr(vd%wsoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch), &
                     busptr(vd%tpsoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch), &
                     busptr(vd%tpsoilv%i)%ptr(((i - 1)*ni + 1):i*ni, trnch)
              end do
              write(iout_soil, *)

              ! Write file containing bulk snow outputs
              write(iout_snow_bulk, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              write(iout_snow_bulk, FMT_CSV, advance = 'no') busptr(vd%snoma%i)%ptr(:, trnch),busptr(vd%snodpl%i)%ptr(:, trnch), &
                        busptr(vd%snoden%i)%ptr(:, trnch), busptr(vd%snoal%i)%ptr(:, trnch),busptr(vd%wsnow%i)%ptr(:, trnch), &                    
                        busptr(vd%tsnow_svs%i)%ptr(1:ni, trnch),busptr(vd%rsnows_acc%i)%ptr(:, trnch),  & 
                        busptr(vd%rainrate%i)%ptr(:, trnch),busptr(vd%snowrate%i)%ptr(:, trnch)
              write(iout_snow_bulk, *)

              if( svs_mesh%vs%lout_snow_enbal) then 
                 ! Write file containing snow energy balance outputs
                  write(iout_snow_enbal, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                  write(iout_snow_enbal, FMT_CSV, advance = 'no') busptr(vd%rnetsa%i)%ptr(:, trnch),busptr(vd%swnetsa%i)%ptr(:, trnch), &
                       busptr(vd%lwnetsa%i)%ptr(:, trnch), -1.0*busptr(vd%les%i)%ptr(:, trnch), -1.0*busptr(vd%hfluxsa%i)%ptr(:, trnch), &
                       busptr(vd%subldrifta%i)%ptr(:, trnch), -1.0*busptr(vd%gfluxsa%i)%ptr(:, trnch), busptr(vd%hpsa%i)%ptr(:, trnch), busptr(vd%esa%i)%ptr(:, trnch)
                  write(iout_snow_enbal, *)
              end if

              ! Write file containing  snow profile outputs
              ! The frequqncy depends on the value of nprofile_day that can be modifid by the user
              nfreq = 24/svs_mesh%vs%nprofile_day
              if ( modulo(ic%now%hour,nfreq) == 0 .and. svs_mesh%vs%lout_snow_profile ) then
                write(iout_snow_profile, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                do i = 1, nsl
                   write(iout_snow_profile, FMT_CSV, advance = 'no') &
                       busptr(vd%snoma_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch) , &
                       busptr(vd%snoden_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch), &
                       busptr(vd%snoage_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch), &
                       busptr(vd%snodiamopt_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch) , &
                       busptr(vd%snospheri_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch), &
                       busptr(vd%snohist_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch), &
                       busptr(vd%tsnow_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch), &
                       busptr(vd%wsnow_svs%i)%ptr(((i - 1)*ni + 1):i*ni, trnch)
                end do
                write(iout_snow_profile, *)

              end if

           end if

       else if(svs_mesh%vs%schmsol=='SVS') then

           if (ic%now%mins ==0) then! Full hour

              k=1 !>  Identity of the tile (offset relative to node-indexing).

              ! Write file containing soil outputs
              write(iout_svs1_soil, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              do i = 1, nl_svs
                 write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                     busptr(vd%isoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch) , &
                     busptr(vd%wsoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch)
              end do
              if(svs_mesh%vs%lsoil_freezing_svs1) then
                 do i = 1, nl_svs
                    write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                       busptr(vd%tpsoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch) 
                 end do
              end if 
              write(iout_svs1_soil, FMT_CSV, advance = 'no') busptr(vd%tground%i)%ptr(1:ni, trnch),busptr(vd%tground%i)%ptr((ni+1):2*ni, trnch),busptr(vd%tvege%i)%ptr(1:ni, trnch),busptr(vd%tvege%i)%ptr(ni+1:2*ni, trnch)
              write(iout_svs1_soil, *)

              ! Write file containing bulk snow outputs
              write(iout_svs1_snow, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              write(iout_svs1_snow, FMT_CSV, advance = 'no') busptr(vd%snoma%i)%ptr(:, trnch),busptr(vd%snodpl%i)%ptr(:, trnch), &
                        busptr(vd%snoden%i)%ptr(:, trnch), busptr(vd%snoal%i)%ptr(:, trnch),busptr(vd%wsnow%i)%ptr(:, trnch), &                    
                        busptr(vd%tsnow%i)%ptr(1:ni, trnch),busptr(vd%tsnow%i)%ptr((ni+1):2*ni, trnch)
              write(iout_svs1_snow, FMT_CSV, advance = 'no') busptr(vd%snvma%i)%ptr(:, trnch),busptr(vd%snvdp%i)%ptr(:, trnch), &
                        busptr(vd%snvden%i)%ptr(:, trnch), busptr(vd%snval%i)%ptr(:, trnch),busptr(vd%wsnv%i)%ptr(:, trnch), &                    
                        busptr(vd%tsnowveg%i)%ptr(1:ni, trnch),busptr(vd%tsnowveg%i)%ptr((ni+1):2*ni, trnch)
              write(iout_svs1_snow, *)


           end if




        end if


        !> Update 'lmin' if active (greater than zero).
        if (lmo_winter > 0.0) then
            if (ic%now%jday < 210) then

                !> Jun 15 -> 167.
                sl_lmin_soil = 1.0 + (lmo_winter - 1.0)*1.0/(1.0 + exp(0.3*(ic%now%jday - 167)))
            else

                !> Sep 15 -> 259.
                sl_lmin_soil = 1.0 + (lmo_winter - 1.0)*1.0/(1.0 + exp(-0.3*(ic%now%jday - 259)))
            end if
        end if

        !> Transfer driving variables.
        if (associated(vs%tile%prern) .and. associated(vs%tile%presno)) then
            busptr(vd%rainrate%i)%ptr(:, trnch) = vs%tile%prern/1000.0
            busptr(vd%snowrate%i)%ptr(:, trnch) = vs%tile%presno/1000.0
        else
            where (vs%tile%ta > tcdk)
                busptr(vd%rainrate%i)%ptr(:, trnch) = vs%tile%pre/1000.0
                busptr(vd%snowrate%i)%ptr(:, trnch) = 0.0
            elsewhere
                busptr(vd%rainrate%i)%ptr(:, trnch) = 0.0
                busptr(vd%snowrate%i)%ptr(:, trnch) = vs%tile%pre/1000.0
            end where
        end if
        busptr(vd%flusolis%i)%ptr(:, trnch) = vs%tile%fsin
        busptr(vd%fdsi%i)%ptr(:, trnch) = vs%tile%flin
        busptr(vd%tmoins%i)%ptr(:, trnch) = vs%tile%ta
        busptr(vd%humoins%i)%ptr(:, trnch) = vs%tile%qa
        busptr(vd%umoins%i)%ptr(:, trnch) = vs%tile%uv
        busptr(vd%vmoins%i)%ptr(:, trnch) = 0.0
        if (associated(vs%tile%uu) .and. associated(vs%tile%vv)) then
            busptr(vd%umoins%i)%ptr(:, trnch) = vs%tile%uu
            busptr(vd%vmoins%i)%ptr(:, trnch) = vs%tile%vv
        end if
        busptr(vd%pmoins%i)%ptr(:, trnch) = vs%tile%pres

        ! Compute virtual tempertaure needed when computing forcing height
        do i = il1,il2 
             tt = busptr(vd%tmoins%i)%ptr(i, trnch)
             hu = busptr(vd%humoins%i)%ptr(i, trnch)
             tve(i) = FOTVT(tt,hu)   
        end do  

        !> Update momentum and thermodynamic levels if using atmospheric forcing (not 'observed_forcing').
        if (.not. svs_mesh%vs%observed_forcing) then
            svs_bus(a1(zusl):z1(zusl)) = -rgasd/grav*log(svs_mesh%vs%sigma_u)*tve(:)  !*dat(ic%ts_count)
            svs_bus(a1(ztsl):z1(ztsl)) = -rgasd/grav*log(svs_mesh%vs%sigma_t)*tve(:)  !*dat(ic%ts_count)
        end if

        !> Required to replace the calculation in 'phystepinit'.
        busptr(vd%thetaa%i)%ptr(:, trnch) = svs_mesh%vs%sigma_t**(-cappa)*busptr(vd%tmoins%i)%ptr(:, trnch)



        !> Call SVS.
        if(svs_mesh%vs%schmsol=='SVS') then
             call svs(svs_bus, bus_length, bus_ptr, nvarsurf, time_dt, kount, trnch, ni, ni, nk)
        else if(svs_mesh%vs%schmsol=='SVS2') then
             call svs2(svs_bus, bus_length, bus_ptr, nvarsurf, time_dt, kount, trnch, ni, ni, nk)
        end if
        if (phy_error_L) then
            call print_error("An error occurred during the iteration of the SVS time-step.")
            call program_abort()
        end if

        !> Copy bus variable.
        call runsvs_mesh_copy_bus_to_vs()

!ierr = 200
!!do i = 1, 26
!!write(ierr, *) busptr(vd%vegf%i)%ptr(((i - 1)*ni + 1):i*ni, trnch); ierr = ierr + 1
!!end do
!write(ierr, *) busptr(vd%accevap%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%fvapliqaf%i)%ptr(:, trnch); ierr = ierr + 1
!do i = 1, nl_svs
!write(ierr, *) busptr(vd%wsoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch); ierr = ierr + 1
!end do
!write(ierr, *) busptr(vd%fl%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%fq%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%fv%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%eflux%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%eg%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%er%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%etr%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%fc%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%resaef%i)%ptr(:, trnch); ierr = ierr + 1
!
!
!write(ierr, *) busptr(vd%z0%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%z0ha%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%z0mvh%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%z0mvl%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%ztsl%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%z0t%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%zusl%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%resagr%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%resavg%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%snoal%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%snoden%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%snodpl%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%snval%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%snvden%i)%ptr(:, trnch); ierr = ierr + 1
!write(ierr, *) busptr(vd%snvdp%i)%ptr(:, trnch); ierr = ierr + 1
!do i = 1, 2
!write(ierr, *) busptr(vd%tsnow%i)%ptr(((i - 1)*ni + 1):i*ni, trnch); ierr = ierr + 1
!end do
!do i = 1, 2
!write(ierr, *) busptr(vd%tsnowveg%i)%ptr(((i - 1)*ni + 1):i*ni, trnch); ierr = ierr + 1
!end do


        !> Transfer variables.
        vs%tile%et = busptr(vd%wflux%i)%ptr(:, trnch)
        vs%tile%ovrflw = max(0.0, busptr(vd%runofftot%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch))/ic%dts
        do i = 1, khyd
            vs%tile%latflw(:, i) = max(0.0, busptr(vd%latflw%i)%ptr(((i - 1)*ni + 1):i*ni, trnch))/ic%dts
        end do
        vs%tile%drainsol = max(0.0, busptr(vd%watflow%i)%ptr((khyd*ni + 1):(khyd + 1)*ni, trnch))/ic%dts

        vs%tile%qacan = busptr(vd%qsurf%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch)
        vs%tile%lqwscan = busptr(vd%wveg%i)%ptr(:, trnch)
        vs%tile%tacan = busptr(vd%tsurf%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch)
        vs%tile%tcan = &
            busptr(vd%tvege%i)%ptr(1:ni, trnch)*0.25 + &
            busptr(vd%tvege%i)%ptr((ni + 1):, trnch)*0.25 + &
            busptr(vd%tsnowveg%i)%ptr(1:ni, trnch)*0.25 + &
            busptr(vd%tsnowveg%i)%ptr((ni + 1):, trnch)*0.25
        vs%tile%sno = busptr(vd%snoma%i)%ptr(:, trnch)
        vs%tile%albsno = &
            busptr(vd%snoal%i)%ptr(:, trnch)*0.5 + &
            busptr(vd%snval%i)%ptr(:, trnch)*0.5
        vs%tile%rhosno = ( &
            busptr(vd%snoro%i)%ptr(:, trnch)*0.5 + &
            busptr(vd%snvro%i)%ptr(:, trnch)*0.5)*900.0
        vs%tile%tsno = &
            busptr(vd%tsnow%i)%ptr(1:ni, trnch)*0.5 + &
            busptr(vd%tsnow%i)%ptr((ni + 1):, trnch)*0.5
        where (busptr(vd%snoma%i)%ptr(:, trnch) > 0.0)
            vs%tile%lqwssno = busptr(vd%wsnow%i)%ptr(:, trnch)
        elsewhere
            vs%tile%lqwssno = 0.0
        end where
        vs%tile%qevp = busptr(vd%fv%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch)
        vs%tile%qsens = busptr(vd%fc%i)%ptr(((indx_soil - 1)*ni + 1):indx_soil*ni, trnch)
        vs%tile%thicsol(:, 1) = busptr(vd%isoil%i)%ptr(1:ni, trnch)
        vs%tile%thlqsol(:, 1) = busptr(vd%wsoil%i)%ptr(1:ni, trnch)
        vs%tile%thlqsol(:, 2) = busptr(vd%wsoil%i)%ptr((ni + 1):2*ni, trnch)
        do i = 3, nl_svs
            vs%tile%thlqsol(:, i) = busptr(vd%wsoil%i)%ptr(((i - 1)*ni + 1):i*ni, trnch)
        end do
        vs%tile%tsol(:, 1) = busptr(vd%tground%i)%ptr(1:ni, trnch)
        do i = 2, 2
            vs%tile%tsol(:, i) = busptr(vd%tground%i)%ptr((ni + 1):, trnch)
        end do


    end subroutine

    subroutine runsvs_mesh_finalize(shd, fls)

        !> Input variables.
        type(ShedGridParams) shd
        type(fl_ids) fls

        !> Return if the process is not marked active.
        if (.not. svs_mesh%PROCESS_ACTIVE) return

    end subroutine

    subroutine runsvs_mesh_save_states_seq(fls, shd)

        !> Input variables.
        type(fl_ids) fls
        type(ShedGridParams) shd

        !> Return if the process is not marked active.
        if (.not. svs_mesh%PROCESS_ACTIVE) return

    end subroutine

end module
