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

!phybus-6.3:
    use phymem, only: phyvar, phymem_get_slabvars

    implicit none

    character(len=1024) :: msg_S

    !> SVS constants.
    type runsvs_mesh_constants
        integer :: NLANDCLASS = NCLASS
        real, dimension(NCLASS) :: Z0DAT = (/ &
            0.001, 0.001, 0.001, 1.75, 2.0, 1.0, 2.0, 3.0, 0.8, 0.1, &
            0.2, 0.2, 0.1, 0.1, 0.15, 0.15, 0.35, 0.25, 0.1, 0.25, &
            5.0, 0.1, 0.1, 0.1, 1.75, 0.5 /)
        real,dimension(NCLASS) :: Z0DAT_UP
    end type

    !> SVS2 output
    integer :: iout_soil = 150
    integer :: iout_snow_bulk = 151
    integer :: iout_snow_profile = 152
    integer :: iout_snow_enbal = 153
    integer :: iout_snow_bulk_vegh = 154
    integer :: iout_svs2_watbal = 155
    integer :: iout_svs2_restart = 159

    !> SVS1 output
    integer :: iout_svs1_soil = 160
    integer :: iout_svs1_snow = 161
    integer :: iout_svs1_watbal = 162
    integer :: iout_svs1_snow_enbal = 163

    real preacc_tot,wsoil_tot,isoil_tot,snow_tot,veg_tot     

    !> SVS variables names for I/O (direct variables).
    character(len = *), parameter, public :: VN_SVS_AGINGCOEF = 'AGINGCOEF' ! For svs2 only
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
    character(len = *), parameter, public :: VN_SVS_OC = 'OC'
    character(len = *), parameter, public :: VN_SVS_BULKSOIL = 'BULKSOIL'
    character(len = *), parameter, public :: VN_SVS_WSOIL = 'WSOIL'
    character(len = *), parameter, public :: VN_SVS_ISOIL = 'ISOIL'
    character(len = *), parameter, public :: VN_SVS_LATFL = 'LATFL'
    character(len = *), parameter, public :: VN_SVS_WATFL = 'WATFL'
    character(len = *), parameter, public :: VN_SVS_KTHERMAL = 'KTHERMAL'
    character(len = *), parameter, public :: VN_SVS_TGROUND = 'TGROUND'
    character(len = *), parameter, public :: VN_SVS_TGROUNDV = 'TGROUNDV' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_VF = 'VF'
    character(len = *), parameter, public :: VN_SVS_Z0V = 'Z0V'
    character(len = *), parameter, public :: VN_SVS_LNZ0 = 'LNZ0'
    character(len = *), parameter, public :: VN_SVS_TVEGE = 'TVEGE'
    character(len = *), parameter, public :: VN_SVS_TVEGEL = 'TVEGEL' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_TVEGEH = 'TVEGEH' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_TFL = 'TFL' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_WVEG = 'WVEG'
    character(len = *), parameter, public :: VN_SVS_WVEG_VL = 'WVEG_VL' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_WVEG_VH = 'WVEG_VH' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_WFL = 'WFL' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_WFL_ICE = 'WFL_ICE' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_SNCMA = 'SNCMA' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_TSNOW = 'TSNOW'
    character(len = *), parameter, public :: VN_SVS_SNODPL = 'SNODPL'
    character(len = *), parameter, public :: VN_SVS_SNODEN = 'SNODEN'
    character(len = *), parameter, public :: VN_SVS_SNOAL = 'SNOAL'
    character(len = *), parameter, public :: VN_SVS_WSNOW = 'WSNOW'
    character(len = *), parameter, public :: VN_SVS_TSNOWVEG = 'TSNOWVEG'
    character(len = *), parameter, public :: VN_SVS_SNVDP = 'SNVDP'
    character(len = *), parameter, public :: VN_SVS_SNVDEN = 'SNVDEN'
    character(len = *), parameter, public :: VN_SVS_SNVAL = 'SNVAL'
    character(len = *), parameter, public :: VN_SVS_WATPOND = 'WATPOND' !  For svs1 and svs2 (with water ponding)
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
    character(len = *), parameter, public :: VN_SVS_HSNOWRES = 'HSNOWRES' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LSNOWDRIFT_SUBLIM = 'LSNOWDRIFT_SUBLIM' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LSNOWAGING_VAR = 'LSNOWAGING_VAR' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_SNOMA = 'SNOMA'
    character(len = *), parameter, public :: VN_SVS_SNVMA = 'SNVMA'
    character(len = *), parameter, public :: VN_SVS_SNOMA_SVS = 'SNOMA_ML'
    character(len = *), parameter, public :: VN_SVS_SNODEN_SVS = 'SNODEN_ML'
    character(len = *), parameter, public :: VN_SVS_SNOAGE_SVS = 'SNOAGE_ML'
    character(len = *), parameter, public :: VN_SVS_SNODIAMOPT_SVS = 'SNODOPT_ML'
    character(len = *), parameter, public :: VN_SVS_SNOSPHERI_SVS = 'SNOSPH_ML'
    character(len = *), parameter, public :: VN_SVS_SNOHIST_SVS = 'SNOHIST_ML'
    character(len = *), parameter, public :: VN_SVS_SNOTYPE_SVS = 'SNOTYPE_ML'
    character(len = *), parameter, public :: VN_SVS_TSNOW_SVS = 'TSNOW_ML'
    character(len = *), parameter, public :: VN_SVS_WSNOW_SVS = 'WSNOW_ML'
    character(len = *), parameter, public :: VN_SVS_SNOMAV_SVS = 'SNOMAV_ML'
    character(len = *), parameter, public :: VN_SVS_SNODENV_SVS = 'SNODENV_ML'
    character(len = *), parameter, public :: VN_SVS_SNOAGEV_SVS = 'SNOAGEV_ML'
    character(len = *), parameter, public :: VN_SVS_SNODIAMOPTV_SVS = 'SNODOPTV_ML'
    character(len = *), parameter, public :: VN_SVS_SNOSPHERIV_SVS = 'SNOSPHV_ML'
    character(len = *), parameter, public :: VN_SVS_SNOHISTV_SVS = 'SNOHISTV_ML'
    character(len = *), parameter, public :: VN_SVS_SNOTYPEV_SVS = 'SNOTYPEV_ML'
    character(len = *), parameter, public :: VN_SVS_TSNOWV_SVS = 'TSNOWV_ML'
    character(len = *), parameter, public :: VN_SVS_WSNOWV_SVS = 'WSNOWV_ML'
    character(len = *), parameter, public :: VN_SVS_LOUT_SNOW_PROFILE = 'LOUT_SNOW_PROFILE' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LOUT_SNOW_ENBAL = 'LOUT_SNOW_ENBAL' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LOUT_SVS1_WATBAL = 'LOUT_SVS1_WATBAL ' ! For svs1 only 
    character(len = *), parameter, public :: VN_SVS_NPROFILE_DAY = 'NPROFILE_DAY' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LOUT_SVS2_WATBAL = 'LOUT_SVS2_WATBAL ' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LSOIL_FREEZING_SVS1 = 'LSOIL_FREEZING_SVS1' ! For svs1 only
    character(len = *), parameter, public :: VN_SVS_LMACROPORES_SVS = 'LMACROPORES_SVS' ! For svs1 only
    character(len = *), parameter, public :: VN_SVS_LPHASE_CHANGE_EFF_SVS1 = 'LPHASE_CHANGE_EFF_SVS1' ! For svs1 only
    character(len = *), parameter, public :: VN_SVS_LWATER_PONDING_SVS = 'LWATER_PONDING_SVS' ! For svs1 and svs2
    character(len = *), parameter, public :: VN_SVS_LUNIQUE_PROFILE_SVS2 = 'LUNIQUE_PROFILE_SVS2' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LBCHEAT_SVS2 = 'LBCHEAT_SVS2'! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LSNOW_INTERCEPTION_SVS2 = 'LSNOW_INTERCEPTION_SVS2' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_CANO_REF_FORCING = 'CANO_REF_FORCING' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LCANO_SVS2 = 'LCANO_SVS2' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LOUT_SNOW_VEGH = 'LOUT_SNOW_VEGH' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_VGH_DENS = 'VGH_DENS' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_HVEGLPOL = 'HVEGLPOL' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LWRITE_RESTART = 'LWRITE_RESTART' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_LFORLIT = 'LFORLIT' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_READ_OC = 'READ_OC' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_LREAD_RESTART = 'LREAD_RESTART' ! For svs2 only 
    character(len = *), parameter, public :: VN_SVS_LVAR_LMIN_STABLE = 'LVAR_LMIN_STABLE '
    character(len = *), parameter, public :: VN_SVS_LMO_WINTER = 'LMO_WINTER'  ! Used if LVAR_LMIN_STABLE == 'VAR'
    character(len = *), parameter, public :: VN_SVS_LMIN_STABLE = 'LMIN_STABLE' ! Used if LVAR_LMIN_STABLE == 'CST'

    !> SVS variables names for I/O (modifiers/special conditions).
    character(len = *), parameter, public :: VN_SVS_SAND_N = 'SAND_N'
    character(len = *), parameter, public :: VN_SVS_CLAY_N = 'CLAY_N'
    character(len = *), parameter, public :: VN_SVS_OC_N = 'OC_N'
    character(len = *), parameter, public :: VN_SVS_BULKSOIL_N = 'BULKSOIL_N'
    character(len = *), parameter, public :: VN_SVS_WSOIL_N = 'WSOIL_N'
    character(len = *), parameter, public :: VN_SVS_ISOIL_N = 'ISOIL_N'
    character(len = *), parameter, public :: VN_SVS_TGROUND_N = 'TGROUND_N'
    character(len = *), parameter, public :: VN_SVS_TGROUNDV_N = 'TGROUNDV_N'
    character(len = *), parameter, public :: VN_SVS_VF_N = 'VF_N'
    character(len = *), parameter, public :: VN_SVS_Z0V_N = 'Z0V_N'
    character(len = *), parameter, public :: VN_SVS_TVEGE_N = 'TVEGE_N'
    character(len = *), parameter, public :: VN_SVS_TVEGEL_N = 'TVEGEL_N'
    character(len = *), parameter, public :: VN_SVS_TVEGEH_N = 'TVEGEH_N'
    character(len = *), parameter, public :: VN_SVS_TFL_N = 'TFL_N'
    character(len = *), parameter, public :: VN_SVS_WFL_N = 'WFL_N'
    character(len = *), parameter, public :: VN_SVS_WFL_ICE_N = 'WFL_ICE_N'
    character(len = *), parameter, public :: VN_SVS_TSNOW_N = 'TSNOW_N'
    character(len = *), parameter, public :: VN_SVS_TSNOWVEG_N = 'TSNOWVEG_N'
    character(len = *), parameter, public :: VN_SVS_TPSOIL_N = 'TPSOIL_N' ! For svs2  and svs1 (with soil freezing)
    character(len = *), parameter, public :: VN_SVS_TPSOILV_N = 'TPSOILV_N' ! For svs2 only
    character(len = *), parameter, public :: VN_SVS_SNOMA_SVS_N = 'SNOMA_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNODEN_SVS_N = 'SNODEN_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNOAGE_SVS_N = 'SNOAGE_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNODIAMOPT_SVS_N = 'SNODOPT_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNOSPHERI_SVS_N = 'SNOSPH_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNOHIST_SVS_N = 'SNOHIST_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNOMAV_SVS_N = 'SNOMAV_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNODENV_SVS_N = 'SNODENV_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNOAGEV_SVS_N = 'SNOAGEV_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNODIAMOPTV_SVS_N = 'SNODOPTV_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNOSPHERIV_SVS_N = 'SNOSPHV_ML_N'
    character(len = *), parameter, public :: VN_SVS_SNOHISTV_SVS_N = 'SNOHISTV_ML_N'
    character(len = *), parameter, public :: VN_SVS_TSNOW_SVS_N = 'TSNOW_ML_N'
    character(len = *), parameter, public :: VN_SVS_WSNOW_SVS_N = 'WSNOW_ML_N'
    character(len = *), parameter, public :: VN_SVS_TSNOWV_SVS_N = 'TSNOWV_ML_N'
    character(len = *), parameter, public :: VN_SVS_WSNOWV_SVS_N = 'WSNOWV_ML_N'   

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
        real, dimension(:, :), allocatable :: bulksoil
        real, dimension(:, :), allocatable :: oc
        real, dimension(:, :), allocatable :: wsoil
        real, dimension(:, :), allocatable :: isoil
        real, dimension(:, :), allocatable :: tpsoil ! For svs2 and svs1 (with soil freezing)
        real, dimension(:, :), allocatable :: tpsoilv ! For svs2 only
        integer :: kthermal = 2
        real, dimension(:, :), allocatable :: tground
        real, dimension(:),    allocatable :: tgroundv
        real, dimension(:, :), allocatable :: vf
        real, dimension(:, :), allocatable :: z0v
        real, dimension(:), allocatable :: lnz0
        real, dimension(:, :), allocatable :: tvege
        real, dimension(:, :), allocatable :: tvegeh ! For svs2 only
        real, dimension(:), allocatable :: tvegel ! For svs2 only
        real, dimension(:), allocatable :: tfl ! For svs2 only
        real, dimension(:), allocatable :: wveg
        real, dimension(:), allocatable :: wveg_vl
        real, dimension(:), allocatable :: wveg_vh
        real, dimension(:), allocatable :: wfl  ! For svs2 only
        real, dimension(:), allocatable :: wfl_ice  ! For svs2 only
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
        real, dimension(:), allocatable :: sncma ! For svs2 only
        real, dimension(:), allocatable :: tperm ! For svs2 only
        real, dimension(:), allocatable :: agingcoef ! For svs2 only
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
        character(len = DEFAULT_FIELD_LENGTH) :: hsnowres = 'RIL'
        logical :: lsnowdrift_sublim = .true.
        logical :: lsnowaging_var = .false.
        logical :: lout_snow_profile = .false.
        logical :: lout_snow_enbal = .false.
        logical :: lout_snow_vegh = .false.
        logical :: lout_svs2_watbal = .false.
        logical :: lwrite_restart = .false.
        logical :: lforlit = .false.
        logical :: read_oc = .false.
        logical :: lread_restart = .false.
        logical :: lout_svs1_watbal = .false.
        integer :: nprofile_day = 4 !
        logical :: lsoil_freezing_svs1 = .false.
        logical :: lmacropores_svs = .false.
        logical :: lphase_change_eff_svs1 = .false.
        logical :: lwater_ponding_svs = .false.
        logical :: lunique_profile_svs2 = .true.
        character(len = DEFAULT_FIELD_LENGTH) :: lbcheat_svs2 = 'TPERM'
        logical :: lsnow_interception_svs2 = .false.
        character(len = DEFAULT_FIELD_LENGTH) :: cano_ref_forcing = 'FOR'
        logical :: lcano_svs2 = .false.
        real, dimension(:), allocatable :: vgh_dens
        real, dimension(:), allocatable :: hveglpol
        logical :: read_hveglpol = .true.

        character(len = DEFAULT_FIELD_LENGTH) :: LVAR_LMIN_STABLE = 'NON'
        real :: lmo_winter = -1.0 !
        real :: lmin_stable = 20. !
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
    integer :: kount_reset = 0
    integer, private :: kount = 0
    integer, parameter, private :: trnch = 1
    integer, private :: ni = 0
    integer, parameter, private :: nk = 1
!phybus-6.3:
    type(phyvar), pointer, contiguous :: pvars(:)

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
!        use phy_status, only: phy_error_L

        !> Modules.
   use bus_builder, only: bb_keylist, bb_n
   use wb_itf_mod
   use cnv_options
   use phy_options
   use phy_status, only: phy_error_L, PHY_OK, physeterror
   use phybusidx
   use ens_perturb, only: ptp_nc, spp_nc, ens_nc2d
   use microphy_utils, only: mp_phybusinit
   use phymem, only: phymem_init, phymem_add, phymem_find, phymem_alloc
   use phymem, only: nphyvars
   use phygridmap, only: phy_lcl_ni, phy_lcl_nj, phydim_ni, phydim_nj, phydim_nk
   use sfc_businit_mod, only: sfc_businit

        !> Input/output variables.
        integer, intent(in) :: ni, nk

        !> Local variables.
#include <rmnlib_basics.hf>
   include "surface.cdk"

   character(len=2), parameter :: E1 = 'e1'
   character(len=2), parameter :: D1 = 'd1'
   character(len=2), parameter :: P0 = 'p0'
   character(len=2), parameter :: P1 = 'p1'
   character(len=2), parameter :: V0 = 'v0'
   character(len=2), parameter :: U0 = 'u0'

   character(len=4), parameter :: LVLA = 'A'
   character(len=4), parameter :: LVLA4= 'A*4'
   character(len=4), parameter :: LVLE = 'E'
   character(len=4), parameter :: LVLM = 'M'
   character(len=4), parameter :: LVLM2= 'M*2'
   character(len=4), parameter :: LVLT = 'T'

   character(len=6)  :: nag, nmar, dwwz, nuv, psss, nccl
   !integer :: ier, iverb, nsurf, nextra, i
   integer :: ier, iverb, nextra, i
   logical :: lbourg3d, lbourg
   logical :: lkfbe, lshal, lshbkf, lmid
   logical :: lmoistke, lrpnint
   logical :: lmoyhr, lmoyhrkf, lmoykfsh, lmoymid
   logical :: lgwdsm, lgwd, ltofd
   logical :: lccc2
   logical :: lghg, ltrigtau, ltrigtauw
   logical :: liuv
   logical :: lmoyhroz, lmoyhrgh, llinozout, llinghout, llinozage
   logical :: lmoycons
   logical :: lhn_init, lsfcflx
   logical :: lsurfonly, lwindgust
   logical :: lpcp_frac, ladvzn, ls2, lmp, lcsun

   ier = phymem_init()
   if (.not.RMN_IS_OK(ier)) then
      call physeterror('phybusinit', 'Problem with phymem module init')
      return
   endif
   nagrege = nsurf + 1

   !# nagg is the dimension of aggregrated variables ('nagrege').
!   write(nag,'(a,i2)') 'A*', nsurf + 1
   write(nag,'(a,i2)') 'A*', nagrege

   !# nipt is the number of tau/cloud top pressure bins in ISCCP histograms
!   write(ntp,'(a,i2)') 'A*', ntau*nptop
   write(nuv,'(a,i2)') 'A*', RAD_NUVBRANDS
   write(nccl,'(a,i2)') 'A*', RAD_TCCL

   !# nmar is the number of 2d Markov fields
!phy_init:
      ptp_nc = 0
      ens_nc2d = max(ptp_nc + spp_nc, 1)
   write(nmar,'(a,i2)') 'A*', ens_nc2d

   ! Retrieve bus requirements for microphysics scheme
   if (mp_phybusinit() /= PHY_OK) then
      call physeterror('phybusinit', &
           'Cannot retrieve microphysics bus information')
      return
   endif

!   lcn_mpx  = (stcond(1:2) == 'MP')
!   lcn_none = .not.lcn_mpx
!   lcn_my2  = (stcond(1:6) == 'MP_MY2')
!   lcn_p3i1 = (stcond == 'MP_P3' .and. p3_ncat >= 1)
!   lcn_p3i2 = (stcond == 'MP_P3' .and. p3_ncat >= 2)
!   lcn_p3i3 = (stcond == 'MP_P3' .and. p3_ncat >= 3)
!   lcn_p3i4 = (stcond == 'MP_P3' .and. p3_ncat == 4)

   lbourg3d= (pcptype == 'BOURGE3D')
   lgwdsm  = (gwdrag /= 'NIL' .and. sgo_tdfilter > 0.)
   lgwd    = (gwdrag /= 'NIL')
   lmoyhr  = (moyhr > 0 .or. dynout)
   lkfbe   = any(convec == (/ &
        'BECHTOLD', &
        'KFC     ', &
        'KFC2    ' &
        /))
   lshal   = (conv_shal /= 'NIL')
   lshbkf  = (conv_shal == 'BECHTOLD')
   lmid    = (conv_mid /= 'NIL')
   lmoyhrkf= (lmoyhr .and. lkfbe)
   lmoymid= (lmoyhr .and. lmid)
   lmoykfsh= (lmoyhr .and. lshbkf .and. bkf_lshalm)
   lbourg  = any(pcptype == (/&
        'BOURGE ', &
        'NIL    ', &
        'SPS_W19', &
        'SPS_FRC', &
        'SPS_H13'  &
        /))
   lrpnint = (fluvert == 'RPNINT')
   lmoistke = (fluvert == 'MOISTKE')
   ladvzn  = (advectke .and. lrpnint)
   lccc2   = (radia == 'CCCMARAD2')
   lghg    = (lccc2 .and. radghg_L)
   ls2     = (stcond == 'S2')
   lcsun   = (stcond == 'CONSUN')
   lmp     = (stcond(1:3) == 'MP_')

   ltrigtau = (kfctrigtau > 0.)
   ltrigtauw = (deep_wavg .or. mid_wavg)

   dwwz = 'd1'
   lsurfonly = (fluvert == 'SURFACE')
   if (lsurfonly) dwwz = 'd0'
   psss = 'p0'
   if (tofd /= 'NIL') psss = 'p1'
   ltofd = (tofd /= 'NIL')
   lpcp_frac = lsurfonly .and. (pcptype == 'SPS_FRC')

   ! Activate energy budget diagnostics only if outputs are requested by the user
   ebdiag = .false.

   ! Activate ECMWF diagnostics only if outputs are requested by the user
   ecdiag = .false.

   ! Activate final-state screen-level diagnostics only if outputs are requested by the user
   fsdiag = .false.

   ! Activate wind gust estimate only if outputs are requested by the user
   lwindgust = .false.

   ! Activate lightning diagnostics only if outputs are requested by the user
   llight = .false.

   ! Activate refractivity diagnostics only if outputs are requested by the user
   lrefract = .false.

   ! Activate energy budget diagnostics only if outputs are requested by the user
   lcons = .false.
   lmoycons = (lcons .and. lmoyhr)

   lhn_init = (lhn /= 'NIL')
   lsfcflx = (sfcflx_filter_order > 0)

   etccdiag = .false.

   lhn_init = (lhn /= 'NIL')
   lsfcflx = (sfcflx_filter_order > 0)

   cmt_comp_diag = .false.

!phy_init
      phydim_ni = ni
      phydim_nj = 1

!?      call mapping2drivergrid()

      !# Establish physics v-grid information
      phydim_nk = nk

#include "phymkptr.hf"
#include "phyvar.hf"
        if (phy_error_L) then
            call print_error("An error occurred initializing the physics bus variable.")
            call program_abort()
        end if

   call sfc_businit(moyhr,ni,nk)
   if (phy_error_L) return

   !#NOTE: phymem_alloc must be done before any call to phymem_find
!phy_init
      nphyoutlist = -1
!      ier = wb_get_meta('itf_phy/PHYOUT', type1, sizeof1, nphyoutlist, options1)
!      if (.not.WB_IS_OK(ier)) nphyoutlist = -1
      nextra = 0
!      if (ptp_L) nextra = 3
      allocate(phyoutlist_S(max(1,nphyoutlist+nextra)))
      phyoutlist_S(:) = ' '
   if (debug_alldiag_L .or. nphyoutlist < 0) then
      ier = phymem_alloc(debug_mem_L, (/'*'/))
   else
      ier = phymem_alloc(debug_mem_L, phyoutlist_S)
   endif
   if (.not.RMN_IS_OK(ier)) &
        call physeterror('phybusinit', 'problem in phymem_alloc')
   if (phy_error_L .or. .not.RMN_IS_OK(ier)) return

#undef PHYMKPTR
#define PHYPTRGETIDX
#include "phymkptr.hf"
#include "phyvar.hf"
   if (phy_error_L) return   

   sigw = sigt

   if (lbourg3d) then
      fip = fip3d
      fneige = fneige3d
   endif

   if (qcmoinsmp > 0) qcmoins = qcmoinsmp
   if (qcphytdmp > 0) qcphytd = qcphytdmp
   if (qcplusmp > 0) qcplus = qcplusmp
   if (qrphytdmp > 0) qrphytd = qrphytdmp

    end subroutine

    subroutine runsvs_mesh_copy_vs_to_bus()

        !> For constants.
        use tdpack_const, only: omega

        !> Modules.
   use inisoili_svs_mod, only: inisoili_svs
   use inisoili_svs2_mod, only: inisoili_svs2

        ! For Crocus debug mode
        use mode_crodebug

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
                if(svs_mesh%vs%z0v(1, 1200 - i)>0.) then
                     svs_mesh%c%Z0DAT_UP(1200 - i) = svs_mesh%vs%z0v(1, 1200 - i)
                else
                     svs_mesh%c%Z0DAT_UP(1200 - i) = svs_mesh%c%Z0DAT(1200 - i)
                endif
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
                if (allocated(svs_mesh%vs%oc)) svs_bus(a2(oc, i - 1):z2(oc, i - 1)) = svs_mesh%vs%oc(:, i)
                if (allocated(svs_mesh%vs%bulksoil)) svs_bus(a2(bulksoil, i - 1):z2(bulksoil, i - 1)) = svs_mesh%vs%bulksoil(:, i)

            end do
            if (svs_mesh%vs%schmsol=='SVS') then
                call inisoili_svs(pvars, ni)
            else if (svs_mesh%vs%schmsol=='SVS2') then
                call inisoili_svs2(pvars, ni)
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
        if(svs_mesh%vs%schmsol=='SVS') then
           do i = 0, 1
             if (allocated(svs_mesh%vs%tvege)) svs_bus(a2(tvege, i):z2(tvege, i)) = svs_mesh%vs%tvege(:, i + 1)
           end do
           if (allocated(svs_mesh%vs%wveg)) svs_bus(a1(wveg):z1(wveg)) = svs_mesh%vs%wveg
        endif

        if(svs_mesh%vs%schmsol=='SVS2') then
           do i = 0, 1
            if (allocated(svs_mesh%vs%tvegeh)) svs_bus(a2(tvegeh, i):z2(tvegeh, i)) = svs_mesh%vs%tvegeh(:, i + 1)
           end do

           if (allocated(svs_mesh%vs%tgroundv)) svs_bus(a1(tgroundv):z1(tgroundv)) = svs_mesh%vs%tgroundv
           if (allocated(svs_mesh%vs%tvegel)) svs_bus(a1(tvegel):z1(tvegel)) = svs_mesh%vs%tvegel
           if (allocated(svs_mesh%vs%wveg_vh)) svs_bus(a1(wveg_vh):z1(wveg_vh)) = svs_mesh%vs%wveg_vh
           if (allocated(svs_mesh%vs%wveg_vl)) svs_bus(a1(wveg_vl):z1(wveg_vl)) = svs_mesh%vs%wveg_vl
           if (allocated(svs_mesh%vs%tfl)) svs_bus(a1(tfl):z1(tfl)) = svs_mesh%vs%tfl
           if (allocated(svs_mesh%vs%wfl)) svs_bus(a1(wfl):z1(wfl)) = svs_mesh%vs%wfl
           if (allocated(svs_mesh%vs%wfl_ice)) svs_bus(a1(wfl_ice):z1(wfl_ice)) = svs_mesh%vs%wfl_ice

           ! At initial step used tvege as the initial condition for tvegel and tvegeh
            !where (svs_bus(a2(tvegeh, i):z2(tvegeh, i)) ==0.)
            !         svs_bus(a2(tvegeh, i):z2(tvegeh, i))  = svs_bus(a2(tvege, i):z2(tvege, i))
            !end where
            !where (svs_bus(a2(tvegel, i):z2(tvegel, i)) ==0.)
            !         svs_bus(a2(tvegel, i):z2(tvegel, i))  = svs_bus(a2(tvege, i):z2(tvege, i))
            !end where

        endif

        ! Add option to use Crocus debug mode
        if(svs_mesh%vs%schmsol=='SVS2') then
                call init_crodebug(svs_mesh%vs%hsnowscheme)
        end if

        if(svs_mesh%vs%schmsol=='SVS2') then
           do i = 1, nl_svs
                if (allocated(svs_mesh%vs%tpsoil))  svs_bus(a2(tpsoil, i - 1):z2(tpsoil, i - 1)) = svs_mesh%vs%tpsoil(:, i)
                if (allocated(svs_mesh%vs%tpsoilv)) svs_bus(a2(tpsoilv, i - 1):z2(tpsoilv, i - 1)) = svs_mesh%vs%tpsoilv(:, i)
           end do
           if (allocated(svs_mesh%vs%tperm)) svs_bus(a1(tperm):z1(tperm)) = svs_mesh%vs%tperm
        endif

        if(svs_mesh%vs%schmsol=='SVS2') then
           if (allocated(svs_mesh%vs%vgh_dens)) svs_bus(a1(vgh_dens):z1(vgh_dens)) = svs_mesh%vs%vgh_dens
           if (allocated(svs_mesh%vs%hveglpol)) svs_bus(a1(hveglpol):z1(hveglpol)) = svs_mesh%vs%hveglpol
        endif

        if(svs_mesh%vs%schmsol=='SVS2' .and. svs_mesh%vs%lsnowaging_var) then
           if (allocated(svs_mesh%vs%agingcoef)) svs_bus(a1(agingcoef):z1(agingcoef)) = svs_mesh%vs%agingcoef              
        endif

        if(svs_mesh%vs%schmsol=='SVS' .and. svs_mesh%vs%lsoil_freezing_svs1) then
           do i = 1, nl_svs
                if (allocated(svs_mesh%vs%tpsoil))  svs_bus(a2(tpsoil, i - 1):z2(tpsoil, i - 1)) = svs_mesh%vs%tpsoil(:, i)
           end do
           if (allocated(svs_mesh%vs%tperm)) svs_bus(a1(tperm):z1(tperm)) = svs_mesh%vs%tperm
        endif

        if((svs_mesh%vs%schmsol=='SVS' .or. svs_mesh%vs%schmsol=='SVS2') .and. svs_mesh%vs%lwater_ponding_svs) then
           if (allocated(svs_mesh%vs%watpond)) svs_bus(a1(watpond):z1(watpond)) = svs_mesh%vs%watpond
        endif

        !if((svs_mesh%vs%schmsol=='SVS' .or. svs_mesh%vs%schmsol=='SVS2')) then
        !   lvar_lmin_stable = svs_mesh%vs%lvar_lmin_stable
        !endif

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

            if(svs_mesh%vs%lsnow_interception_svs2) then
                if (allocated(svs_mesh%vs%sncma))  svs_bus(a1(sncma):z1(sncma)) = svs_mesh%vs%sncma
            end if

            if(svs_mesh%vs%lcano_svs2) then
                if (allocated(svs_mesh%vs%sncma))  svs_bus(a1(sncma):z1(sncma)) = svs_mesh%vs%sncma
            end if
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
        use phy_status, only: phy_error_L, physeterror

        !> For surface layer configuration.
   use phymem
        use sfclayer, only: sl_put, SL_OK

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
!phybus-6.1        diusst_warmlayer = .true.
!phybus-6.1        diusst_coolskin = .true.
        diusst_warmlayer_lakes = .true.
!phybus-6.1        diusst_coolskin_lakes = .true.
        z0mtype = 'BELJAARS'
        z0ttype = 'DEACU12'
        salty_qsat = .true.
!phybus-6.1        urban_params_new = .true.
!phybus-6.3
        svs_urban_params = .true.
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

        ! Activate or not the macropores in SVS1
        if(svs_mesh%vs%schmsol=='SVS' .or. svs_mesh%vs%schmsol=='SVS2') then
                lmacropores_svs = svs_mesh%vs%lmacropores_svs
        endif

        ! Activate or not the phase change efficiency in SVS1
        if(svs_mesh%vs%schmsol=='SVS') then
                lphase_change_eff_svs1 = svs_mesh%vs%lphase_change_eff_svs1
        endif

        ! Activate or not the water ponding module in SVS1 and SVS2
        if(svs_mesh%vs%schmsol=='SVS' .or. svs_mesh%vs%schmsol=='SVS2') then
                lwater_ponding_svs = svs_mesh%vs%lwater_ponding_svs
        endif

        ! Activate or not the unique soil column in SVS2
        if(svs_mesh%vs%schmsol=='SVS2') then
                ! Activate or not the unique soil column in SVS2
                !lunique_profile_svs2 = svs_mesh%vs%lunique_profile_svs2
                ! Select the option used for the lower boundary condition of the heat diffusion in SVS2
                lbcheat_svs2 = svs_mesh%vs%lbcheat_svs2
        endif

        ! Activate or not snow interception in SVS2
        if(svs_mesh%vs%schmsol=='SVS2') then
                lsnow_interception_svs2 = svs_mesh%vs%lsnow_interception_svs2
        endif

        ! Choice of the canopy case study in SVS2, i.e., open to forest (under canopy), forest (under canopy), above canopy
        if(svs_mesh%vs%schmsol=='SVS2') then
                cano_ref_forcing = svs_mesh%vs%cano_ref_forcing
        endif

        if(svs_mesh%vs%schmsol=='SVS2') then
                lforlit = svs_mesh%vs%lforlit
                read_oc = svs_mesh%vs%read_oc
        endif
        ! Activate or not the canopy module for snow under the canopy
!?sps-a16        if(svs_mesh%vs%schmsol=='SVS2') then
!?sps-a16                lcano_svs2 = svs_mesh%vs%lcano_svs2
!?sps-a16        endif

        ! Activate or not the use of user-entered height of low veg.
        if(svs_mesh%vs%schmsol=='SVS2') then
                read_hveglpol = svs_mesh%vs%read_hveglpol
        endif

        ierr =0
        ! Initialize number of snow layers (for multilayer snowpack schemes in SVS2)
        if(svs_mesh%vs%schmsol=='SVS2') then
             nsl  = svs_mesh%vs%nsl
             hsnowscheme =  svs_mesh%vs%hsnowscheme
             hsnowdrift_cro = svs_mesh%vs%hsnowdrift_cro
             lsnowdrift_sublim = svs_mesh%vs%lsnowdrift_sublim
             lsnowaging_var = svs_mesh%vs%lsnowaging_var
             hsnowcomp =  svs_mesh%vs%hsnowcomp
             hsnowcond =  svs_mesh%vs%hsnowcond
             hsnowrad =  svs_mesh%vs%hsnowrad
             hsnowfall =  svs_mesh%vs%hsnowfall
             hsnowhold =  svs_mesh%vs%hsnowhold
             hsnowmetamo =  svs_mesh%vs%hsnowmetamo
             hsnowres =  svs_mesh%vs%hsnowres

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

             if (.not.any(hsnowres == hsnowres_opt)) then
                 call str_concat(msg_S, hsnowres_opt,', ')
                 call print_error('hsnowres = '//trim(hsnowres)//' is not a valid option. Choose among: '//trim(msg_S))
                 ierr = 1
             endif

            if (ierr /= 0) then
                call program_abort()
            end if

        endif

        write(*, nml = surface_cfgs)

        !> Initialize snowpack constants for Crocus and ES
        call ini_csts

        ! Update physical parameters for Crocus using values provided in MESH_parameter.txt
        !if(svs_mesh%vs%xvaging_noglacier>0.) then
        !       xvaging_noglacier=svs_mesh%vs%xvaging_noglacier
        !end if

        !> Initialize the physics bus.
        call phy_businit(ni, nk)

!physlb-6.3:
      nullify(pvars)
      ierr = phymem_get_slabvars(pvars, F_trnch = 1)
      if (.not.(RMN_IS_OK(ierr) .and. associated(pvars))) then
         call physeterror('physlb1', 'Problem getting slab vars pointers')
         return
      endif

        !> Initialize the surface bus in the physics library.
!moved:phy_businit
!        call sfc_businit(moyhr, ni, nk)
!        if (phy_error_L) then
!            call print_error("An error occurred initializing the surface bus variable.")
!            call program_abort()
!        end if

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
        if(svs_mesh%vs%schmsol=='SVS') then
           call runsvs_mesh_append_phyentvar('tvege')
           call runsvs_mesh_append_phyentvar('wveg')
        endif
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
        if( (svs_mesh%vs%schmsol=='SVS' .or. svs_mesh%vs%schmsol=='SVS2') .and. svs_mesh%vs%lwater_ponding_svs) then
            call runsvs_mesh_append_phyentvar('watpond')
        endif

        if(svs_mesh%vs%schmsol=='SVS2') then
             call runsvs_mesh_append_phyentvar('tvegeh')
             call runsvs_mesh_append_phyentvar('wveg_vh')
             call runsvs_mesh_append_phyentvar('tvegel')
             call runsvs_mesh_append_phyentvar('tfl')
             call runsvs_mesh_append_phyentvar('wfl')
             call runsvs_mesh_append_phyentvar('wfl_ice')

             call runsvs_mesh_append_phyentvar('wveg_vl')
             call runsvs_mesh_append_phyentvar('tgroundv')
             call runsvs_mesh_append_phyentvar('tpsoil')
             call runsvs_mesh_append_phyentvar('tpsoilv')
             call runsvs_mesh_append_phyentvar('tperm')
             call runsvs_mesh_append_phyentvar('vgh_dens')
             call runsvs_mesh_append_phyentvar('hveglpol')

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

             if(svs_mesh%vs%lsnow_interception_svs2) then
                 call runsvs_mesh_append_phyentvar('sncma')
             endif

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
            vd%oc%mul = nl_stp
            vl(vd%oc%i)%mul = nl_stp
            vd%bulksoil%mul = nl_stp
            vl(vd%bulksoil%i)%mul = nl_stp
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

!6.1:
! This macro only gets pointer address to pass as an argument to calling function --- to MANIPULATE VARIABLE, USE MACROS  MK... below
! def-6.1 (NAME2) busptr(vd%NAME2%i)%ptr(1,trnch)
! Assign local variables below
! def-6.1 (NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
! def-6.1 (NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)
!6.3:
! This macro only gets pointer address to pass as an argument to calling function --- to MANIPULATE VARIABLE, USE MACROS  MK... below
! def-6.3 (NAME2) pvars(vd%NAME2%idxv)%data(:)
! Assign local variables below
! def-6.3 (NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni) => pvars(vd%NAME2%idxv)%data(:)
! def-6.3 (NAME1,NAME2) nullify(NAME1); if (vd%NAME2%idxv > 0) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => pvars(vd%NAME2%idxv)%data(:)

            !> Manually assign the surface pointer.
!surface/copybus:
            do i = 1, nvarsurf
                if (vl(i)%niveaux <= 0 .or. .not.vl(i)%doagg_L) cycle
                if (vl(i)%idxv <= 0) cycle
                if (associated(pvars(vl(i)%idxv)%data)) then
                    nullify(pvars(vl(i)%idxv)%data)
                end if
                pvars(vl(i)%idxv)%data => phy_bus(bus_ptr(i):(bus_ptr(i) + vl(i)%niveaux*vl(i)%mul*vl(i)%mosaik*ni - 1), trnch)
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
            if(svs_mesh%vs%schmsol=='SVS') then
                write(line, "('VEGETATION TEMP.: ', 2f8.3)") svs_bus(a1(tvege)), svs_bus(a1(tvege) + ni)
                call print_message(line)
            else if(svs_mesh%vs%schmsol=='SVS2') then
                write(line, "('LOW VEGETATION TEMP.: ', 2f8.3)") svs_bus(a1(tvegel))
                call print_message(line)
                write(line, "('HIGH VEGETATION TEMP.: ', 2f8.3)") svs_bus(a1(tvegeh)), svs_bus(a1(tvegeh) + ni)
                call print_message(line)
            endif
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
            call print_message('             % SAND    % CLAY    % OC   %BULKSOIL')
            do i = 1, nl_svs ! model layers
                write(line, "(' LAYER ', i3, ': ', 999(f8.3, 3x))") i, svs_bus(a2(sand, i - 1)), svs_bus(a2(clay, i -1)),svs_bus(a2(oc, i - 1)), svs_bus(a2(bulksoil, i - 1))
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
               call print_message('             Soil profile ')
               do i = 1, nl_svs ! permeable layers
                   write(line, "(' LAYER ', i3, ': ', 999(f8.3, 1x))") i, svs_bus(a2(tpsoil, i - 1))
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

       ! Initialize variable for water mass balance
       preacc_tot = 0.

       open(iout_soil, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_soil_hourly.csv', action = 'write')
       write(iout_soil, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       do j = 1, nl_svs
           write(level, FMT_GEN) j
           write(iout_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_ISOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_WSOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_TPSOIL) // '_' // trim(adjustl(level))
       end do
       if(.not. svs_mesh%vs%lunique_profile_svs2) then
         do j = 1, nl_svs
           write(level, FMT_GEN) j
           write(iout_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_TPSOILV) // '_' // trim(adjustl(level))
         end do
       endif
       write(iout_soil, FMT_CSV, advance = 'no') 'TVEGL','TVEGH','TGROUND','TGROUNDV', 'TFL', 'WVEGL','WVEGH', 'WFL', 'WFL_ICE'
       write(iout_soil, *)


       open(iout_snow_bulk, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_snow_bulk_hourly.csv', action = 'write')
       write(iout_snow_bulk, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       write(iout_snow_bulk, FMT_CSV, advance = 'no') 'SNOMA', 'SNODP','SNODEN','SNOALB','WSNO','TSNO_SURF','RSNOW_AC','RAINRATE', 'SNOWRATE', 'SCF_LV'
       write(iout_snow_bulk, *)

       if(svs_mesh%vs%lout_snow_vegh) then
          open(iout_snow_bulk_vegh, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_snow_bulk_veg_hourly.csv', action = 'write')
          write(iout_snow_bulk_vegh, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
          write(iout_snow_bulk_vegh, FMT_CSV, advance = 'no') 'SNVMA',   &
          'SNVDP','SNVDEN','SNVALB','WSNV','TSNV_SURF','RSNV_AC','RAINRATE_VGH', 'SNOWRATE_VGH','WDRIFT_VGH'
          if(svs_mesh%vs%lsnow_interception_svs2) then
              write(iout_snow_bulk_vegh, FMT_CSV, advance = 'no') 'SNCMA'
          endif
          write(iout_snow_bulk_vegh, *)
      endif

       if(svs_mesh%vs%lout_snow_enbal) then
           open(iout_snow_enbal, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_snow_enbal_hourly.csv', action = 'write')
           write(iout_snow_enbal, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
           write(iout_snow_enbal, FMT_CSV, advance = 'no') 'SNO_RNET', 'SNO_SWNET','SNO_LWNET','SNO_LE','SNO_H','SNO_BSUBL','SNO_G', 'SNO_HRAIN', 'SUBLIM_LV'
           if(svs_mesh%vs%lout_snow_vegh) then
                   write(iout_snow_enbal, FMT_CSV, advance = 'no') 'LW_DCA','SW_DCA','TA_CA','HU_CA','WS_CA','SNV_LE','SNV_H' 
              if(svs_mesh%vs%lsnow_interception_svs2) then
                   write(iout_snow_enbal, FMT_CSV, advance = 'no') 'ESNC','ESNCAF','INTSNC','INTSNCAF','DRPSNC','DRPSNCAF','UNLSNC','UNLSNCAF', &
                   'MFSNC','MFSNCAF'
              endif
           endif
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
                            trim(VN_SVS_WSNOW_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOTYPE_SVS) // '_' // trim(adjustl(level))
          end do
          if(svs_mesh%vs%lout_snow_vegh) then
                do j = 1, nsl
                       write(level, FMT_GEN) j
                       write(iout_snow_profile, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_SNOMAV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNODENV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOAGEV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNODIAMOPTV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOSPHERIV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOHISTV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_TSNOWV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_WSNOWV_SVS) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_SNOTYPEV_SVS) // '_' // trim(adjustl(level))
               end do
          endif
          write(iout_snow_profile, *)
       endif

       if(svs_mesh%vs%lout_svs2_watbal) then
          open(iout_svs2_watbal, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs2_watbal_hourly.csv', action = 'write')
          write(iout_svs2_watbal, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
          write(iout_svs2_watbal, FMT_CSV, advance = 'no') 'PCP_AC','EVP_AC', 'LAT_AC', 'DRA_AC', 'ROF_AC','ROF_INS'
          write(iout_svs2_watbal, FMT_CSV, advance = 'no') 'WSOIL_TOT','ISOIL_TOT','SNOW_TOT','VEG_TOT'
          write(iout_svs2_watbal, *)
       endif

       write(*,*) 'restart',svs_mesh%vs%lwrite_restart
       if(svs_mesh%vs%lwrite_restart) then        
            open(iout_svs2_restart, file = './' // trim(fls%GENDIR_OUT) // '/' // 'restart_svs2.csv', action = 'write')
       endif

   else if(svs_mesh%vs%schmsol=='SVS') then

       ! Initialize variable for water balance
       preacc_tot = 0.     

       open(iout_svs1_soil, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs1_soil_hourly.csv', action = 'write')
       write(iout_svs1_soil, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       do j = 1, nl_svs
           write(level, FMT_GEN) j
           write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_ISOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_WSOIL) // '_' // trim(adjustl(level)), &
                            trim(VN_SVS_LATFL) // '_' // trim(adjustl(level))
       end do
       do j = 1, nl_svs+1
           write(level, FMT_GEN) j
           write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_WATFL) // '_' // trim(adjustl(level))
       end do
       if(svs_mesh%vs%lsoil_freezing_svs1) then
          do j = 1, nl_svs
             write(level, FMT_GEN) j
             write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                            trim(VN_SVS_TPSOIL) // '_' // trim(adjustl(level))
         end do
       endif

       write(iout_svs1_soil, FMT_CSV, advance = 'no') 'TGROUND_1','TGROUND_2','TVEG_1','TVEG_2','ALBSFC','FC','FV','FL','RNET_S','SATSFC'
       write(iout_svs1_soil, *)

       open(iout_svs1_snow, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs1_snow_bulk_hourly.csv', action = 'write')
       write(iout_svs1_snow, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
       write(iout_svs1_snow, FMT_CSV, advance = 'no') 'SNOMA', 'SNODP','SNODEN','SNOALB','WSNO','TSNO_1','TSNO_2'
       write(iout_svs1_snow, FMT_CSV, advance = 'no') 'SNVMA', 'SNVDP','SNVDEN','SNVALB','WSNV','TSNV_1','TSNV_2'
       write(iout_svs1_snow, *)

       if(svs_mesh%vs%lout_svs1_watbal) then 
          open(iout_svs1_watbal, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs1_watbal_hourly.csv', action = 'write')
          write(iout_svs1_watbal, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
          write(iout_svs1_watbal, FMT_CSV, advance = 'no') 'PCP_AC','EVP_AC', 'LAT_AC', 'DRA_AC', 'ROF_AC','ROF_INS'
          write(iout_svs1_watbal, FMT_CSV, advance = 'no') 'WSOIL_TOT','ISOIL_TOT','SNOW_TOT','VEG_TOT' 
          write(iout_svs1_watbal, *)
       endif

       if(svs_mesh%vs%lout_snow_enbal) then
           open(iout_svs1_snow_enbal, file = './' // trim(fls%GENDIR_OUT) // '/' // 'svs1_snow_enbal_hourly.csv', action = 'write')
           write(iout_svs1_snow_enbal, FMT_CSV, advance = 'no') 'YEAR', 'JDAY', 'HOUR', 'MINS'
           write(iout_svs1_snow_enbal, FMT_CSV, advance = 'no') 'SNO_RNET','SNO_LE','SNO_H','SNO_MELT_AC','SNO_MELTR_AC'
           if(svs_mesh%vs%lout_snow_vegh) then
                   write(iout_svs1_snow_enbal, FMT_CSV, advance = 'no') 'SNV_RNET','SNV_LE','SNV_H'
           endif
           write(iout_svs1_snow_enbal, *)
       endif       


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
        if(svs_mesh%vs%schmsol=='SVS') then
          if (.not. allocated(svs_mesh%vs%tvege)) allocate(svs_mesh%vs%tvege(ni, 2))
          if (.not. allocated(svs_mesh%vs%wveg)) allocate(svs_mesh%vs%wveg(ni))
        endif
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
        if(svs_mesh%vs%schmsol=='SVS') then
          do i = 0, 1
            svs_mesh%vs%tvege(:, i + 1) = svs_bus(a2(tvege, i):z2(tvege, i))
          end do
          svs_mesh%vs%wveg = svs_bus(a1(wveg):z1(wveg))
        endif
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

        if( (svs_mesh%vs%schmsol=='SVS' .or. svs_mesh%vs%schmsol=='SVS2') .and. svs_mesh%vs%lwater_ponding_svs) then 
               if (.not. allocated(svs_mesh%vs%watpond)) allocate(svs_mesh%vs%watpond(ni))
               svs_mesh%vs%watpond = svs_bus(a1(watpond):z1(watpond ))
        end if

        if(svs_mesh%vs%schmsol=='SVS2') then
            if (.not. allocated(svs_mesh%vs%tpsoil)) allocate(svs_mesh%vs%tpsoil(ni, nl_svs))
            if (.not. allocated(svs_mesh%vs%tpsoilv)) allocate(svs_mesh%vs%tpsoilv(ni, nl_svs))
            if (.not. allocated(svs_mesh%vs%tperm)) allocate(svs_mesh%vs%tperm(ni))

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
            if (.not. allocated(svs_mesh%vs%vgh_dens)) allocate(svs_mesh%vs%vgh_dens(ni))
            if (.not. allocated(svs_mesh%vs%hveglpol)) allocate(svs_mesh%vs%hveglpol(ni))

            do i = 1, nl_svs
               svs_mesh%vs%tpsoil(:, i) = svs_bus(a2(tpsoil, i - 1):z2(tpsoil, i - 1))
               svs_mesh%vs%tpsoilv(:, i) = svs_bus(a2(tpsoilv, i - 1):z2(tpsoilv, i - 1))
            end do

            svs_mesh%vs%tperm = svs_bus(a1(tperm):z1(tperm ))
            svs_mesh%vs%vgh_dens = svs_bus(a1(vgh_dens):z1(vgh_dens))
            svs_mesh%vs%hveglpol = svs_bus(a1(hveglpol):z1(hveglpol))

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

            if (.not. allocated(svs_mesh%vs%tvegel)) allocate(svs_mesh%vs%tvegel(ni))
            if (.not. allocated(svs_mesh%vs%tfl)) allocate(svs_mesh%vs%tfl(ni))
            if (.not. allocated(svs_mesh%vs%tvegeh)) allocate(svs_mesh%vs%tvegeh(ni,2))
            do i = 0, 1
                 svs_mesh%vs%tvegeh(:, i + 1) = svs_bus(a2(tvegeh, i):z2(tvegeh, i))
            end do
            svs_mesh%vs%tvegel = svs_bus(a1(tvegel):z1(tvegel ))
            svs_mesh%vs%tfl = svs_bus(a1(tfl):z1(tfl ))

            if (.not. allocated(svs_mesh%vs%wveg_vl)) allocate(svs_mesh%vs%wveg_vl(ni))
            if (.not. allocated(svs_mesh%vs%wveg_vh)) allocate(svs_mesh%vs%wveg_vh(ni))
            svs_mesh%vs%wveg_vl = svs_bus(a1(wveg_vl):z1(wveg_vl ))
            svs_mesh%vs%wveg_vh = svs_bus(a1(wveg_vh):z1(wveg_vh ))


            if (.not. allocated(svs_mesh%vs%wfl)) allocate(svs_mesh%vs%wfl(ni))
            if (.not. allocated(svs_mesh%vs%wfl_ice)) allocate(svs_mesh%vs%wfl_ice(ni))
            svs_mesh%vs%wfl = svs_bus(a1(wfl):z1(wfl ))
            svs_mesh%vs%wfl_ice = svs_bus(a1(wfl_ice):z1(wfl_ice ))

            if (.not. allocated(svs_mesh%vs%tgroundv)) allocate(svs_mesh%vs%tgroundv(ni))
            svs_mesh%vs%tgroundv = svs_bus(a1(tgroundv):z1(tgroundv ))

            if(svs_mesh%vs%lsnow_interception_svs2) then
               if (.not. allocated(svs_mesh%vs%sncma)) allocate(svs_mesh%vs%sncma(ni))
               svs_mesh%vs%sncma = svs_bus(a1(sncma):z1(sncma))
            end if

        endif

    end subroutine

    subroutine runsvs_mesh_within_tile(shd, fls, cm)

        !> For 'cmcdate_fromprint' function to convert from string date.
        use cmcdate_mod, only: cmcdate_fromprint

        !> For 'jdate_from_cmc' to convert to 'jdateo' for SVS.
        use mu_jdate_mod, only: jdate_from_cmc

        !> For RPN/physics status.
        use phy_status, only: phy_error_L, physeterror

        !> For constants.
        use tdpack_const, only: rgasd, grav, cappa, tcdk

        use tdpack, only: fotvt

        !> Modules.
        use inichamp, only: inichamp4
        use sfc_calcdiag, only: sfc_calcdiag3
        use svs2_mod, only: svs2
        use svs_mod, only: svs
        use aggveghigh_mod, only: aggveghigh
        use aggveglow_mod, only: aggveglow


        type(ShedGridParams) :: shd
        type(fl_ids) :: fls
        type(clim_info) :: cm

        type(dates_model) :: ts

        !> Local variable for string date (format: YYYYMMDD.hhmmss).
        character(len = 14) time_run_now

        !> Local variables.
        integer i, idateo, ierr, j, k, istat, nfreq,nstep_tot,nstep_now

        real :: tt,hu
        real, dimension(il1:il2) :: tve
        !real, dimension(il1:il2) :: tt
        !real, dimension(il1:il2) :: hu

        !> Return if the process is not marked active.
        if (.not. svs_mesh%PROCESS_ACTIVE) return

        !> Update variables (equivalent to calls to 'phyput_input_param' and 'sfc_get_input_param').


         ! > Get current and last time steps
         nstep_tot = jday_to_tsteps(ic%stop%year,ic%stop%jday,ic%stop%hour,ic%stop%mins,ic%dtmins )
         nstep_now = jday_to_tsteps(ic%now%year,ic%now%jday,ic%now%hour,ic%now%mins,ic%dtmins )
        


    ! Write SVS hourly outputs
1010    format(9999(g15.7e2, ','))

        if(svs_mesh%vs%schmsol=='SVS2') then

           ! Compute variable for water mass balance
           preacc_tot = preacc_tot +  1000.*sum(pvars(vd%rainrate%idxv)%data(:))*ic%dts +  1000.*sum(pvars(vd%snowrate%idxv)%data(:))*ic%dts
           call layer_thickness()
           wsoil_tot=0.
           isoil_tot=0.
           do j = 1, svs_mesh%vs%khyd
               wsoil_tot = wsoil_tot + 1000.0*pvars(vd%wsoil%idxv)%data(j)*delz(j) !mm
               isoil_tot = isoil_tot + 1000.0*pvars(vd%isoil%idxv)%data(j)*delz(j) !mm
           end do
           snow_tot = pvars(vd%snoma%idxv)%data(1)* (1.-pvars(vd%vegh%idxv)%data(1)) + & ! pvars(vd%svs_wtg%idxv)%data(5) +  & ! pvars(vd%vegl%idxv)%data(1) + &
                      pvars(vd%snvma%idxv)%data(1)*pvars(vd%vegh%idxv)%data(1)

           veg_tot =  pvars(vd%wveg_vl%idxv)%data(1)*pvars(vd%svs_wtg%idxv)%data(3) +  & ! Weight of vl
                      pvars(vd%wveg_vh%idxv)%data(1) *pvars(vd%vegh%idxv)%data(1) +  &
                      pvars(vd%sncma%idxv)%data(1) *pvars(vd%vegh%idxv)%data(1) 


           !if (ic%now%hour /= ic%next%hour) then !last time-step of hour
           if (ic%now%mins ==0) then! Full hour

              k=1 !>  Identity of the tile (offset relative to node-indexing).

              ! Write file containing soil outputs
              write(iout_soil, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              do i = 1, nl_svs
                 write(iout_soil, FMT_CSV, advance = 'no') &
                     pvars(vd%isoil%idxv)%data(((i - 1)*ni + 1):i*ni) , &
                     pvars(vd%wsoil%idxv)%data(((i - 1)*ni + 1):i*ni), &
                     pvars(vd%tpsoil%idxv)%data(((i - 1)*ni + 1):i*ni)
              end do
              if(.not. svs_mesh%vs%lunique_profile_svs2) then
                do i = 1, nl_svs
                   write(iout_soil, FMT_CSV, advance = 'no') &
                     pvars(vd%tpsoilv%idxv)%data(((i - 1)*ni + 1):i*ni)
                end do
              endif
              write(iout_soil, FMT_CSV, advance = 'no') pvars(vd%tvegel%idxv)%data(1:ni),pvars(vd%tvegeh%idxv)%data(1:ni), &
                      pvars(vd%tground%idxv)%data(1:ni) , pvars(vd%tgroundv%idxv)%data(1:ni), pvars(vd%tfl%idxv)%data(1:ni),&
                      pvars(vd%wveg_vl%idxv)%data(1), pvars(vd%wveg_vh%idxv)%data(1), pvars(vd%wfl%idxv)%data(1), pvars(vd%wfl_ice%idxv)%data(1)
              write(iout_soil, *)

              ! Write file containing bulk snow outputs
              write(iout_snow_bulk, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              write(iout_snow_bulk, FMT_CSV, advance = 'no') pvars(vd%snoma%idxv)%data(:),pvars(vd%snodpl%idxv)%data(:), &
                        pvars(vd%snoden%idxv)%data(:), pvars(vd%snoal%idxv)%data(:),pvars(vd%wsnow%idxv)%data(:), &
                        pvars(vd%tsnow_svs%idxv)%data(1:ni),pvars(vd%rsnows_acc%idxv)%data(:),  &
                        pvars(vd%rainrate%idxv)%data(:),pvars(vd%snowrate%idxv)%data(:),pvars(vd%PSNGRVL%idxv)%data(:)
              write(iout_snow_bulk, *)


              if( svs_mesh%vs%lout_snow_enbal) then
                 ! Write file containing snow energy balance outputs
                  write(iout_snow_enbal, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                  write(iout_snow_enbal, FMT_CSV, advance = 'no') pvars(vd%rnetsa%idxv)%data(:),pvars(vd%swnetsa%idxv)%data(:), &
                       pvars(vd%lwnetsa%idxv)%data(:), -1.0*pvars(vd%lfluxsa%idxv)%data(:), -1.0*pvars(vd%hfluxsa%idxv)%data(:), &
                       pvars(vd%subldrifta%idxv)%data(:), -1.0*pvars(vd%gfluxsa%idxv)%data(:), pvars(vd%hpsa%idxv)%data(:), pvars(vd%esa%idxv)%data(:)
                  if( svs_mesh%vs%lout_snow_vegh) then
                          write(iout_snow_enbal, FMT_CSV, advance = 'no') pvars(vd%lwca%idxv)%data(:),pvars(vd%swca%idxv)%data(:), &
                           pvars(vd%tca%idxv)%data(:), pvars(vd%qca%idxv)%data(:), pvars(vd%vca%idxv)%data(:),  &
                           -1.0*pvars(vd%lfluxsv%idxv)%data(:), -1.0*pvars(vd%HFLUXSV%idxv)%data(:)
                     if( svs_mesh%vs%lsnow_interception_svs2) then                    
                          write(iout_snow_enbal, FMT_CSV, advance = 'no')  -1.0*pvars(vd%esnc%idxv)%data(:), &
                           -1.0*pvars(vd%esncaf%idxv)%data(:),pvars(vd%intsnc%idxv)%data(:),pvars(vd%intsncaf%idxv)%data(:),          &
                           -1.0*pvars(vd%drpsnc%idxv)%data(:),-1.0*pvars(vd%drpsncaf%idxv)%data(:),-1.0*pvars(vd%unlsnc%idxv)%data(:),          &
                           -1.0*pvars(vd%unlsncaf%idxv)%data(:),-1.0*pvars(vd%mfsnc%idxv)%data(:),-1.0*pvars(vd%mfsncaf%idxv)%data(:) 
                     endif
                  endif
                  write(iout_snow_enbal, *)
              end if

              if( svs_mesh%vs%lout_snow_vegh) then
               ! Write file containing bulk snow outputs
                 write(iout_snow_bulk_vegh, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                 write(iout_snow_bulk_vegh, FMT_CSV, advance = 'no') pvars(vd%snvma%idxv)%data(:),pvars(vd%snvdp%idxv)%data(:), &
                        pvars(vd%snvden%idxv)%data(:), pvars(vd%snval%idxv)%data(:),pvars(vd%wsnv%idxv)%data(:), &
                        pvars(vd%tsnowv_svs%idxv)%data(1:ni),pvars(vd%rsnowsv_acc%idxv)%data(:),  &
                        pvars(vd%rainrate_vgh%idxv)%data(:),pvars(vd%snowrate_vgh%idxv)%data(:),pvars(vd%vca_drift%idxv)%data(:)
                 if( svs_mesh%vs%lsnow_interception_svs2) then
                      write(iout_snow_bulk_vegh, FMT_CSV, advance = 'no')  pvars(vd%sncma%idxv)%data(:)
                 endif
                 write(iout_snow_bulk_vegh, *)
              endif

              ! Write file containing  snow profile outputs
              ! The frequqncy depends on the value of nprofile_day that can be modifid by the user
              nfreq = 24/svs_mesh%vs%nprofile_day
              if ( modulo(ic%now%hour,nfreq) == 0 .and. svs_mesh%vs%lout_snow_profile ) then
                write(iout_snow_profile, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                do i = 1, nsl
                   write(iout_snow_profile, FMT_CSV, advance = 'no') &
                       pvars(vd%snoma_svs%idxv)%data(((i - 1)*ni + 1):i*ni) , &
                       pvars(vd%snoden_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snoage_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snodiamopt_svs%idxv)%data(((i - 1)*ni + 1):i*ni) , &
                       pvars(vd%snospheri_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snohist_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%tsnow_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%wsnow_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snotype_svs%idxv)%data(((i - 1)*ni + 1):i*ni)
                end do
                if( svs_mesh%vs%lout_snow_vegh) then
                  do i = 1, nsl
                     write(iout_snow_profile, FMT_CSV, advance = 'no') &
                       pvars(vd%snomav_svs%idxv)%data(((i - 1)*ni + 1):i*ni) , &
                       pvars(vd%snodenv_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snoagev_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snodiamoptv_svs%idxv)%data(((i - 1)*ni + 1):i*ni) , &
                       pvars(vd%snospheriv_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snohistv_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%tsnowv_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%wsnowv_svs%idxv)%data(((i - 1)*ni + 1):i*ni), &
                       pvars(vd%snotypev_svs%idxv)%data(((i - 1)*ni + 1):i*ni)
                  end do
                endif
                write(iout_snow_profile, *)

              end if

              ! Write file containing water balance information
              if(svs_mesh%vs%lout_svs2_watbal) then
                 write(iout_svs2_watbal, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                 write(iout_svs2_watbal, FMT_CSV, advance = 'no') preacc_tot,pvars(vd%accevap%idxv)%data(:)
                 write(iout_svs2_watbal, FMT_CSV, advance = 'no') pvars(vd%latflaf%idxv)%data(:),pvars(vd%drainaf%idxv)%data(:)
                 write(iout_svs2_watbal, FMT_CSV, advance = 'no') pvars(vd%runofftotaf%idxv)%data(1),pvars(vd%runofftot%idxv)%data(1)
                 write(iout_svs2_watbal, FMT_CSV, advance = 'no') wsoil_tot,isoil_tot,snow_tot,veg_tot
                 write(iout_svs2_watbal, *)
              end if
        end if

        if(nstep_now == nstep_tot-1 .and. svs_mesh%vs%lwrite_restart) then
          write(*,*) 'Write restart file'

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'wveg_vl',pvars(vd%wveg_vl%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'wveg_vh',pvars(vd%wveg_vh%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'SNCMA',pvars(vd%sncma%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tground', pvars(vd%tground%idxv)%data(:) 
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tgroundv', pvars(vd%tgroundv%idxv)%data(:) 
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tvegel', pvars(vd%tvegel%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tvegeh', pvars(vd%tvegeh%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tfl', pvars(vd%tfl%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'wfl', pvars(vd%wfl%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'wfl_ice', pvars(vd%wfl_ice%idxv)%data(:)
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tpsoil'
          do i = 1, nl_svs
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%tpsoil%idxv)%data(((i - 1)*ni + 1):i*ni)
          end do                
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'wsoil'
          do i = 1, nl_svs
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%wsoil%idxv)%data(((i - 1)*ni + 1):i*ni)
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'isoil'
          do i = 1, nl_svs
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%isoil%idxv)%data(((i - 1)*ni + 1):i*ni)
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snoma_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snoma_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snoden_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snoden_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snoage_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snoage_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snodopt_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snodiamopt_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snosph_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snospheri_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snohist_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snohist_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tsnow_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%tsnow_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'wsnow_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%wsnow_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snomav_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snomav_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snodenv_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snodenv_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snoagev_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snoagev_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snodoptv_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snodiamoptv_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snosphv_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snospheriv_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'snohistv_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%snohistv_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'tsnowv_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%tsnowv_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)

          write(iout_svs2_restart, FMT_GEN, advance = 'no') 'wsnowv_ml'
          do i = 1, nsl
               write(iout_svs2_restart, FMT_GEN, advance = 'no') pvars(vd%wsnowv_svs%idxv)%data(((i - 1)*ni + 1):i*ni) 
          end do  
          write(iout_svs2_restart, *)                  

        end if

       else if(svs_mesh%vs%schmsol=='SVS') then

           ! Compute variable for water mass balance
           preacc_tot = preacc_tot +  1000.*sum(pvars(vd%rainrate%idxv)%data(:))*ic%dts +  1000.*sum(pvars(vd%snowrate%idxv)%data(:))*ic%dts 
           call layer_thickness()
           wsoil_tot=0.
           isoil_tot=0.
           do j = 1, svs_mesh%vs%khyd
               wsoil_tot = wsoil_tot + 1000.0*pvars(vd%wsoil%idxv)%data(j)*delz(j) !mm
               isoil_tot = isoil_tot + 1000.0*pvars(vd%isoil%idxv)%data(j)*delz(j) !mm
           end do
           snow_tot = (pvars(vd%snoma%idxv)%data(1)+pvars(vd%wsnow%idxv)%data(1))*(1 - pvars(vd%vegh%idxv)%data(1)) +  &
                      (pvars(vd%snvma%idxv)%data(1)+pvars(vd%wsnv%idxv)%data(1))*pvars(vd%vegh%idxv)%data(1)
              
           veg_tot =  pvars(vd%wveg%idxv)%data(1) * (pvars(vd%vegh%idxv)%data(1) + pvars(vd%vegl%idxv)%data(1))


           if (ic%now%mins ==0) then! Full hour

              k=1 !>  Identity of the tile (offset relative to node-indexing).

              ! Write file containing soil outputs
              write(iout_svs1_soil, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              do i = 1, nl_svs
                 write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                     pvars(vd%isoil%idxv)%data(((i - 1)*ni + 1):i*ni) , &
                     pvars(vd%wsoil%idxv)%data(((i - 1)*ni + 1):i*ni), & 
                     pvars(vd%latflw%idxv)%data(((i - 1)*ni + 1):i*ni) 
              end do
              do i = 1, nl_svs+1
                 write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                     pvars(vd%watflow%idxv)%data(((i - 1)*ni + 1):i*ni) 
              end do
              if(svs_mesh%vs%lsoil_freezing_svs1) then
                 do i = 1, nl_svs
                    write(iout_svs1_soil, FMT_CSV, advance = 'no') &
                       pvars(vd%tpsoil%idxv)%data(((i - 1)*ni + 1):i*ni)
                 end do
              end if
              write(iout_svs1_soil, FMT_CSV, advance = 'no') pvars(vd%tground%idxv)%data(1:ni),pvars(vd%tground%idxv)%data((ni+1):2*ni), &
                      pvars(vd%tvege%idxv)%data(1:ni),pvars(vd%tvege%idxv)%data(ni+1:2*ni), &
                      pvars(vd%alvis%idxv)%data(1:ni),pvars(vd%fc%idxv)%data(1:ni),pvars(vd%fv%idxv)%data(1:ni), &
                      pvars(vd%fl%idxv)%data(1:ni),pvars(vd%rnet_s%idxv)%data(1:ni),0.0 !removed-6.3:pvars(vd%satsfc%idxv)%data(1:ni) 
              write(iout_svs1_soil, *)

              ! Write file containing bulk snow outputs
              write(iout_svs1_snow, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
              write(iout_svs1_snow, FMT_CSV, advance = 'no') pvars(vd%snoma%idxv)%data(:),pvars(vd%snodpl%idxv)%data(:), &
                        pvars(vd%snoden%idxv)%data(:), pvars(vd%snoal%idxv)%data(:),pvars(vd%wsnow%idxv)%data(:), &
                        pvars(vd%tsnow%idxv)%data(1:ni),pvars(vd%tsnow%idxv)%data((ni+1):2*ni)
              write(iout_svs1_snow, FMT_CSV, advance = 'no') pvars(vd%snvma%idxv)%data(:),pvars(vd%snvdp%idxv)%data(:), &
                        pvars(vd%snvden%idxv)%data(:), pvars(vd%snval%idxv)%data(:),pvars(vd%wsnv%idxv)%data(:), &
                        pvars(vd%tsnowveg%idxv)%data(1:ni),pvars(vd%tsnowveg%idxv)%data((ni+1):2*ni)
              write(iout_svs1_snow, *)

              ! Write file containing bulk snow outputs
              if(svs_mesh%vs%lout_svs1_watbal) then 
                 write(iout_svs1_watbal, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                 write(iout_svs1_watbal, FMT_CSV, advance = 'no') preacc_tot,pvars(vd%accevap%idxv)%data(:)
                 write(iout_svs1_watbal, FMT_CSV, advance = 'no') pvars(vd%latflaf%idxv)%data(:),pvars(vd%drainaf%idxv)%data(:)
                 write(iout_svs1_watbal, FMT_CSV, advance = 'no') pvars(vd%runofftotaf%idxv)%data(1),pvars(vd%runofftot%idxv)%data(1)
                 write(iout_svs1_watbal, FMT_CSV, advance = 'no') wsoil_tot,isoil_tot,snow_tot,veg_tot
                 write(iout_svs1_watbal, *)
              end if

              if( svs_mesh%vs%lout_snow_enbal) then
                 ! Write file containing snow energy balance outputs
                  write(iout_svs1_snow_enbal, FMT_CSV, advance = 'no') ic%now%year, ic%now%jday, ic%now%hour, ic%now%mins
                  write(iout_svs1_snow_enbal, FMT_CSV, advance = 'no') pvars(vd%rnetsa%idxv)%data(:),-1.0*pvars(vd%les%idxv)%data(:), &
                        -1.0*pvars(vd%hfluxsa%idxv)%data(:), pvars(vd%melts%idxv)%data(:),pvars(vd%meltsr%idxv)%data(:)
                  if( svs_mesh%vs%lout_snow_vegh) then
                          write(iout_svs1_snow_enbal, FMT_CSV, advance = 'no') pvars(vd%rnetsv%idxv)%data(:),-1.0*pvars(vd%lesv%idxv)%data(:), &
                        -1.0*pvars(vd%hfluxsv%idxv)%data(:)
                  endif
                  write(iout_svs1_snow_enbal, *)
              end if
           end if

        end if

        !> Reset the initialization periodically (at least daily).
        if (ic%ts_count == 1 .or. (ic%now%hour == kount_reset .and. ic%now%mins == 0)) then
            call runsvs_mesh_copy_vs_to_bus()
            kount = 0
            call inichamp4(pvars, kount, ni, nk)

            ! Update roughness length for high and low veg. to be consistent with information provided by the users in
            ! MESH_paramter.txt
            if (allocated(svs_mesh%vs%z0v)) then
                    call aggveghigh(pvars(vd%vegf%idxv)%data(:), log(svs_mesh%c%Z0DAT_UP), log(svs_mesh%c%Z0DAT_UP), pvars(vd%z0mvh%idxv)%data(:), &
                    pvars(vd%dlat%idxv)%data(:), ni, NCLASS)
                    pvars(vd%z0mvh%idxv)%data(:) = EXP(pvars(vd%z0mvh%idxv)%data(:))

                    call aggveglow(pvars(vd%vegf%idxv)%data(:), log(svs_mesh%c%Z0DAT_UP), log(svs_mesh%c%Z0DAT_UP), pvars(vd%z0mvl%idxv)%data(:), &
                    pvars(vd%dlat%idxv)%data(:), ni, NCLASS)
                    pvars(vd%z0mvl%idxv)%data(:) = EXP(pvars(vd%z0mvl%idxv)%data(:))
            endif

            ! Compute beginning date of simulation
            write(time_run_now, "(i4.4, 2i2.2, '.', 2i2.2)") ic%now%year, ic%now%month, ic%now%day, ic%now%hour, ic%now%mins
            idateo = cmcdate_fromprint(time_run_now)
            jdateo = jdate_from_cmc(idateo)

        end if


        !> Increment 'kount'.
        kount = kount + 1

        !# Options for minimum length of MO
        !# * 'NON'   :  lmo_winter is fixed to -1 and a minimum wind speed is imposed
        !# * 'CST'   :  a minimum constant length of MO is specified by the used (lmin_stable)
        !# * 'VAR'   :  a variable length of MO is used as a function of the season. It ranges between the value lmo_winter in wintertime (value provided by the user) and 1 in summer.
        if (svs_mesh%vs%lvar_lmin_stable .EQ. 'CST') then ! Use constant lmin = lmin_stable 
            sl_lmin_soil = svs_mesh%vs%lmin_stable
        else if (svs_mesh%vs%lvar_lmin_stable .EQ. 'NON') then
            sl_lmin_soil = -1.
        else
        !> Update 'lmin' if active (greater than zero).
            if (svs_mesh%vs%lmo_winter > 0.0) then
                if (ic%now%jday < 210) then

                    !> Jun 15 -> 167.
                    sl_lmin_soil = 1.0 + (svs_mesh%vs%lmo_winter - 1.0)*1.0/(1.0 + exp(0.3*(ic%now%jday - 167)))
                else

                    !> Sep 15 -> 259.
                    sl_lmin_soil = 1.0 + (svs_mesh%vs%lmo_winter - 1.0)*1.0/(1.0 + exp(-0.3*(ic%now%jday - 259)))
                end if
            end if
        end if

        !> Transfer driving variables.
        if (associated(vs%tile%prern) .and. associated(vs%tile%presno)) then
            pvars(vd%rainrate%idxv)%data(:) = vs%tile%prern/1000.0
            pvars(vd%snowrate%idxv)%data(:) = vs%tile%presno/1000.0
        else
            where (vs%tile%ta > tcdk)
                pvars(vd%rainrate%idxv)%data(:) = vs%tile%pre/1000.0
                pvars(vd%snowrate%idxv)%data(:) = 0.0
            elsewhere
                pvars(vd%rainrate%idxv)%data(:) = 0.0
                pvars(vd%snowrate%idxv)%data(:) = vs%tile%pre/1000.0
            end where
        end if
        pvars(vd%flusolis%idxv)%data(:) = vs%tile%fsin
        pvars(vd%fdsi%idxv)%data(:) = vs%tile%flin
        pvars(vd%tmoins%idxv)%data(:) = vs%tile%ta
        pvars(vd%humoins%idxv)%data(:) = vs%tile%qa
        pvars(vd%umoins%idxv)%data(:) = vs%tile%uv
        pvars(vd%vmoins%idxv)%data(:) = 0.0
        if (associated(vs%tile%uu) .and. associated(vs%tile%vv)) then
            pvars(vd%umoins%idxv)%data(:) = vs%tile%uu
            pvars(vd%vmoins%idxv)%data(:) = vs%tile%vv
        end if
        pvars(vd%pmoins%idxv)%data(:) = vs%tile%pres

        ! Compute virtual tempertaure needed when computing forcing height
        do i = il1,il2
             tt = pvars(vd%tmoins%idxv)%data(i)
             hu = pvars(vd%humoins%idxv)%data(i)
             tve(i) = FOTVT(tt,hu)
        end do

        !> Update momentum and thermodynamic levels if using atmospheric forcing (not 'observed_forcing').
        if (.not. svs_mesh%vs%observed_forcing) then
            svs_bus(a1(zusl):z1(zusl)) = -rgasd/grav*log(svs_mesh%vs%sigma_u)*tve(:)  !*dat(ic%ts_count)
            svs_bus(a1(ztsl):z1(ztsl)) = -rgasd/grav*log(svs_mesh%vs%sigma_t)*tve(:)  !*dat(ic%ts_count)
        end if


        !> Required to replace the calculation in 'phystepinit'.
        pvars(vd%thetaa%idxv)%data(:) = svs_mesh%vs%sigma_t**(-cappa)*pvars(vd%tmoins%idxv)%data(:)



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

!lacchr = (mod(step_driver-1, acchr) == 0)
!step_driver=moyhr=acchr=0 (reset every time-step, MESH handles aggregations)
   call sfc_calcdiag3(pvars, 0, 0, delt, kount, 0, ni)
   if (phy_error_L) return

        !> Copy bus variable.
        call runsvs_mesh_copy_bus_to_vs()

!ierr = 200
!!do i = 1, 26
!!write(ierr, *) pvars(vd%vegf%idxv)%data(((i - 1)*ni + 1):i*ni); ierr = ierr + 1
!!end do
!write(ierr, *) pvars(vd%accevap%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%fvapliqaf%idxv)%data(:); ierr = ierr + 1
!do i = 1, nl_svs
!write(ierr, *) pvars(vd%wsoil%idxv)%data(((i - 1)*ni + 1):i*ni); ierr = ierr + 1
!end do
!write(ierr, *) pvars(vd%fl%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%fq%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%fv%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni); ierr = ierr + 1
!write(ierr, *) pvars(vd%eflux%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%eg%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%er%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%etr%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%fc%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni); ierr = ierr + 1
!write(ierr, *) pvars(vd%resaef%idxv)%data(:); ierr = ierr + 1
!
!
!write(ierr, *) pvars(vd%z0%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni); ierr = ierr + 1
!write(ierr, *) pvars(vd%z0ha%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%z0mvh%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%z0mvl%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%ztsl%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%z0t%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni); ierr = ierr + 1
!write(ierr, *) pvars(vd%zusl%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%resagr%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%resavg%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%snoal%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%snoden%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%snodpl%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%snval%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%snvden%idxv)%data(:); ierr = ierr + 1
!write(ierr, *) pvars(vd%snvdp%idxv)%data(:); ierr = ierr + 1
!do i = 1, 2
!write(ierr, *) pvars(vd%tsnow%idxv)%data(((i - 1)*ni + 1):i*ni); ierr = ierr + 1
!end do
!do i = 1, 2
!write(ierr, *) pvars(vd%tsnowveg%idxv)%data(((i - 1)*ni + 1):i*ni); ierr = ierr + 1
!end do


        !> Transfer variables.
        vs%tile%et = pvars(vd%wflux%idxv)%data(:)
        vs%tile%ovrflw = max(0.0, pvars(vd%runofftot%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni))/ic%dts
        do i = 1, khyd
            vs%tile%latflw(:, i) = max(0.0, pvars(vd%latflw%idxv)%data(((i - 1)*ni + 1):i*ni))/ic%dts
        end do
        vs%tile%drainsol = max(0.0, pvars(vd%watflow%idxv)%data((khyd*ni + 1):(khyd + 1)*ni))/ic%dts

        vs%tile%qacan = pvars(vd%qsurf%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni)
!svs1? idxv=-1 with svs2
!print*,'wveg',vd%wveg%idxv
      if (vd%wveg%idxv > -1) then
        vs%tile%lqwscan = pvars(vd%wveg%idxv)%data(:)
      else
        vs%tile%lqwscan = 0.0
      end if
        vs%tile%tacan = pvars(vd%tsurf%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni)
!print*,'tvege',vd%tvege%idxv
!print*,'tsnowveg',vd%tsnowveg%idxv
      if (vd%tvege%idxv > -1) then
        vs%tile%tcan = &
            pvars(vd%tvege%idxv)%data(1:ni)*0.25 + &
            pvars(vd%tvege%idxv)%data((ni + 1):)*0.25 + &
            pvars(vd%tsnowveg%idxv)%data(1:ni)*0.25 + &
            pvars(vd%tsnowveg%idxv)%data((ni + 1):)*0.25
      else
        vs%tile%tcan = &
!            pvars(vd%tvege%idxv)%data(1:ni)*0.25 + &
!            pvars(vd%tvege%idxv)%data((ni + 1):)*0.25 + &
            pvars(vd%tsnowveg%idxv)%data(1:ni)*0.5 + &
            pvars(vd%tsnowveg%idxv)%data((ni + 1):)*0.5
      end if
        vs%tile%sno = pvars(vd%snoma%idxv)%data(:)
        vs%tile%albsno = &
            pvars(vd%snoal%idxv)%data(:)*0.5 + &
            pvars(vd%snval%idxv)%data(:)*0.5
        vs%tile%rhosno = ( &
            pvars(vd%snoro%idxv)%data(:)*0.5 + &
            pvars(vd%snvro%idxv)%data(:)*0.5)*900.0
        vs%tile%tsno = &
            pvars(vd%tsnow%idxv)%data(1:ni)*0.5 + &
            pvars(vd%tsnow%idxv)%data((ni + 1):)*0.5
        where (pvars(vd%snoma%idxv)%data(:) > 0.0)
            vs%tile%lqwssno = pvars(vd%wsnow%idxv)%data(:)
        elsewhere
            vs%tile%lqwssno = 0.0
        end where
        vs%tile%qevp = pvars(vd%fv%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni)
        vs%tile%qsens = pvars(vd%fc%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni)
        vs%tile%thicsol(:, 1) = pvars(vd%isoil%idxv)%data(1:ni)
        vs%tile%thlqsol(:, 1) = pvars(vd%wsoil%idxv)%data(1:ni)
        vs%tile%thlqsol(:, 2) = pvars(vd%wsoil%idxv)%data((ni + 1):2*ni)
        do i = 3, nl_svs
            vs%tile%thlqsol(:, i) = pvars(vd%wsoil%idxv)%data(((i - 1)*ni + 1):i*ni)
        end do
        vs%tile%tsol(:, 1) = pvars(vd%tground%idxv)%data(1:ni)
        do i = 2, 2
!6.3:tground doesn't have second dimension
!            vs%tile%tsol(:, i) = pvars(vd%tground%idxv)%data((ni + 1):)
            vs%tile%tsol(:, i) = pvars(vd%tground%idxv)%data(1:ni)
        end do

        ! Cumulate surface runoff for land surface tile
!         pvars(vd%runofftotaf%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni) =   &
!                           pvars(vd%runofftotaf%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni)   + &
!                           pvars(vd%runofftot%idxv)%data(((indx_soil - 1)*ni + 1):indx_soil*ni)

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
