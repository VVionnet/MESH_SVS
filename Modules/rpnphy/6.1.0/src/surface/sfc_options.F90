!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module sfc_options
   use, intrinsic :: iso_fortran_env, only: INT64
   implicit none
   public
   save
 
   !#
   integer, parameter :: CLASS_URB = 21 !# Class "urban"
   integer, parameter :: NCLASS    = 26 !# NUMBER OF CLASSES FOR NATURAL COVERS
   integer, parameter :: NCLASSURB = 12 !# NUMBER OF URBAN CLASSES
   integer, parameter :: NL = 3 !# nombre de niveaux dans la glace marine


   !# ----  FOR SVS LAND SCHEME ------
   integer, parameter :: MAX_NL_SVS = 50 !# maximum number of soil layers specified by user
   integer, parameter :: NL_SVS_DEFAULT = 7 ! default number of soil layers
   ! default depth of default soil layers in [METERS]
   real, dimension(NL_SVS_DEFAULT):: DP_SVS_DEFAULT =  (/ 0.05, 0.1, 0.2, 0.4, 1.0, 2.0, 3.0 /)
   !#----------------------------------

   real, parameter :: CRITEXTURE = 0.1
   real, parameter :: CRITLAC    = 0.01
   real, parameter :: CRITMASK   = 0.001
   real, parameter :: CRITSNOW   = 0.0001
   real, parameter :: CRITWATER  = 0.001
   real, parameter :: HIMIN      = 0.001
   real, parameter :: MINICEDP   = 0.05
   real, parameter :: N0RIB = 1.0E-5
   real, parameter :: SNOH0 = 0.1
   real, parameter :: VAMIN = 1.0e-4
   real, parameter :: Z0GLA = 0.0003

   !#
   logical           :: atm_external = .false.
   logical           :: atm_tplus  = .false.
   logical           :: climat     = .false.
   logical           :: cplocn     = .false.
   real              :: delt       = 0.
   integer(INT64)  :: jdateo     = 0
   logical           :: rad_off    = .false.
   logical           :: radslope   = .false.
   logical           :: update_alwater = .false.
   logical           :: z0veg_only = .false.
   logical           :: thermal_stress = .false.
   integer           :: kntveg     = -1
   logical           :: timings_L  = .false.
   integer           :: nphyoutlist = -1
   character(len=32), pointer :: phyoutlist_S(:) => NULL()

   !# Surface layer coefficients for exchange with PBL
   real :: bh91_a, bh91_b, bh91_c, bh91_d, d97_as, dg92_ci, l07_ah, l07_am

   !# Adjust surface temperature over snow after reading (coherency check)
   logical           :: adj_i0_snow = .true.
   namelist /surface_cfgs/ adj_i0_snow

   !# Prandtl number for neutral stability (initialized by SL module)
   real              :: beta        = 0.
   namelist /surface_cfgs/ beta

   !# Diurnal SST scheme
   !# * 'NIL    ' : No Diurnal SST scheme
   !# * 'FAIRALL' : #TODO: define
   character(len=16) :: diusst      = 'NIL'
   namelist /surface_cfgs/ diusst
   character(len=*), parameter :: DIUSST_OPT(2) = (/ &
        'NIL    ', &
        'FAIRALL'  &
        /)

   !# Diurnal SST scheme active coolskin if .true.
   logical           :: diusst_coolskin = .true.
   namelist /surface_cfgs/ diusst_coolskin

   !# Diurnal SST scheme active coolskin over freshwater lakes if .true.
   logical           :: diusst_coolskin_lakes = .true.
   namelist /surface_cfgs/ diusst_coolskin_lakes

   !# Diurnal SST scheme active warmlayer if .true.
   logical           :: diusst_warmlayer = .true.
   namelist /surface_cfgs/ diusst_warmlayer

   !# Diurnal SST scheme active warmlayer over freshwater lakes if .true.
   logical           :: diusst_warmlayer_lakes = .true.
   namelist /surface_cfgs/ diusst_warmlayer_lakes

   !# Depth of soil layers in [METERS] in SVS land surface scheme (schmsol=SVS)
   real :: dp_svs(MAX_NL_SVS) = -1.0
   namelist /surface_cfgs/ dp_svs

   !# Emissivity for ice (glacier and sea ice)
   !# * '_constant_' : A fixed floating point value used as a constant
   character(len=16) :: ice_emiss = '0.99'
   real              :: ice_emiss_const = -1.
   namelist /surface_cfgs/ ice_emiss

   !# Set water temperature of ice-covered lakes to 0C for points north of
   !# ice line if .true.
   !# needs an initialization file otherwise the model stops
   logical           :: icelac      = .false.
   namelist /surface_cfgs/ icelac

   !# Sea ice melting
   logical           :: icemelt     = .false.
   namelist /surface_cfgs/ icemelt

   !# Implicit surface fluxes if .true.; explicit fluxes if .false.
   logical           :: impflx      = .false.
   namelist /surface_cfgs/ impflx

   !# If .true. apply temporary fix to ISBA
   !# * timestep dependent KCOEF
   !# * No PSN factor for meting and freezing
   logical           :: isba_melting_fix = .false.
   namelist /surface_cfgs/ isba_melting_fix

   !# Use the vegetation-only roughness length to compute vegetation snow fraction
   logical           :: isba_snow_z0veg = .false.
   namelist /surface_cfgs/ isba_snow_z0veg

   !# Emissivity for bare soil (ISBA scheme only)
   !# * '_constant_' : A fixed floating point value used as a constant
   !# * 'CLIMATO'    : Value read from an input climatology file (EMIB)
   character(len=16) :: isba_soil_emiss = '0.95'
   real              :: isba_soil_emiss_const = -1.
   namelist /surface_cfgs/ isba_soil_emiss
   character(len=*), parameter :: ISBA_SOIL_EMISS_OPT(1) = (/ &
        'CLIMATO' &
        /)

   !# If .true., freeze precipitation reaching the ground in sub-zero conditions
   logical           :: isba_zr_freeze = .false.
   namelist /surface_cfgs/ isba_zr_freeze

   !# OBSOLETE, REPLACED by KHYD !!! WILL BE EVENTUALLY REMOVED
   !# Deepest active (permeable) soil layer in SVS land surface scheme (schmsol=SVS)
   integer           :: kdp    = -1
   namelist /surface_cfgs/ kdp

   !# Last/Deepest soil layer considered during the accumulation of
   !# lateral flow and drainage. Drainage is taken as the vertical flux
   !# leaving layer KHYD, and lateral flow as the sum of lateral flows from
   !# layers 1 to KHYD
   integer           :: khyd    = -1
   namelist /surface_cfgs/ khyd

   !# Number of snow layers in multi-layer snowpack scheme:
   !# In Crocus it refers to the maximal number of snow layers
   integer           :: nsl    = 12
   namelist /surface_cfgs/ nsl

   !# Vegetation field update frequency (units D,H,M,S,P)
   character(len=16) :: kntveg_S     = ''
   namelist /surface_cfgs/ kntveg_S

   !# Lead fraction for ice-covered lakes
   real              :: lake_leadfrac = 0.
   namelist /surface_cfgs/ lake_leadfrac

   !# Minimum fraction of leads in sea ice.&nbsp; Multiply ice fraction by (1.-leadfrac)
   real              :: leadfrac    = 0.03
   namelist /surface_cfgs/ leadfrac

   !# Limit snow depth to 10 cm for calculation of heat conductivity of snow
   !# over sea-ice and glacier if .true.
   logical           :: limsnodp    = .false.
   namelist /surface_cfgs/ limsnodp

   !# (coupling) fluxes over ocean are taken from ocean model if .true.
   logical           :: owflux      = .false.
   namelist /surface_cfgs/ owflux

   !# read-in land surface emissivity if .true.
   logical           :: read_emis     = .false.
   namelist /surface_cfgs/ read_emis

   !# read-in high vegetation roughness for SVS if .true.
   logical           :: read_z0vh     = .false.
   namelist /surface_cfgs/ read_z0vh

   !# Takes into account effect of ocean salinity on saturation specific
   !# humidity at ocean surface (boundary condition for LH flux calculation)
   logical           :: salty_qsat  = .false.
   namelist /surface_cfgs/ salty_qsat

   !# Land surface processes
   !# * 'NIL ' : No Land surface processes
   !# * 'ISBA' : Interaction Soil Biosphere Atmosphere (ISBA) land sfc scheme
   !# * 'SVS ' : Soil, Vegetation, and Snow (SVS) (Multibudget) land sfc scheme
   !# * 'SVS2' : Advanced version of the SVS land sfc scheme
   character(len=16) :: schmsol     = 'ISBA'
   namelist /surface_cfgs/ schmsol
   character(len=*), parameter :: SCHMSOL_OPT(4) = (/ &
        'NIL ', &
        'ISBA', &
        'SVS ', &
        'SVS2'  &
        /)

   !# Urban surface processes
   !# * 'NIL' : No Urban surface processes
   !# * 'TEB' : Town Energy Balance (TEB) urban scheme
   character(len=16) :: schmurb     = 'NIL'
   namelist /surface_cfgs/ schmurb
   character(len=*), parameter :: SCHMURB_OPT(2) = (/ &
        'NIL', &
        'TEB'  &
        /)

   !# Minimum Obukhov length (L) for glaciers
   real           :: sl_Lmin_glacier = -1.
   namelist /surface_cfgs/ sl_Lmin_glacier
   
   !# Minimum Obukhov length (L) for sea ice
   real           :: sl_Lmin_seaice = -1.
   namelist /surface_cfgs/ sl_Lmin_seaice
   
   !# Mimimum Obukhov length (L) for soil surfaces
   real           :: sl_Lmin_soil = -1.
   namelist /surface_cfgs/ sl_Lmin_soil

   !# Minimum Obukhov length (L) for water
   real           :: sl_Lmin_water = -1.
   namelist /surface_cfgs/ sl_Lmin_water
   
   !# Define bulk Ri values for near-neutral regime in the surface layer
   real           :: sl_rineutral = 0.
   namelist /surface_cfgs/ sl_rineutral

   !# Class of stability functions (stable case) to use in the surface layer
   !# * 'DELAGE97  ' : Use functions described by Delage (1997; BLM)
   !# * 'BELJAARS91' : Use functions described by Beljaars and Holtslag (1991; JAM)
   !# * 'LOCK07    ' : Use functions described by Lock (2007; Tech Report) employed at UKMO
   character(len=16) :: sl_func_stab = 'DELAGE97'
   namelist /surface_cfgs/ sl_func_stab
   character(len=*), parameter :: SL_FUNC_STAB_OPT(3) = (/ &
        'DELAGE97  ', &
        'BELJAARS91', &
        'LOCK07    '  &
        /)

   !# Class of stability functions (unstable case) to use in the surface layer
   !# * 'DELAGE92' : Use functions described by Delage and Girard (1992; BLM)
   !# * 'DYER74  ' : Use functions described by Dyer (1974; BLM)
   character(len=16) :: sl_func_unstab = 'DELAGE92'
   namelist /surface_cfgs/ sl_func_unstab
   character(len=*), parameter :: SL_FUNC_UNSTAB_OPT(2) = (/ &
        'DELAGE92', &
        'DYER74  '  &
        /)
 
   !# Use a reference roughness for surface layer calculations
   logical :: sl_z0ref = .false.
   namelist /surface_cfgs/ sl_z0ref

   !# Emissivity for snow
   !# * '_constant_' : A fixed floating point value used as a constant
   character(len=16) :: snow_emiss = '1.'
   real              :: snow_emiss_const = -1.
   namelist /surface_cfgs/ snow_emiss

   !#  Soil texture database/calculations for SVS land surface scheme
   !# * 'GSDE   '   : 8 layers of sand & clay info from Global Soil Dataset for ESMs (GSDE)
   !# * 'SLC    '   : 5 layers of sand & clay info from Soil Landscape of Canada (SLC)
   !# * 'SOILGRIDS' : 7 layers of sand & clay info from ISRIC ? World Soil Information
   character(len=16) :: soiltext    = 'GSDE'
   namelist /surface_cfgs/ soiltext
   character(len=*), parameter :: SOILTEXT_OPT(3) = (/ &
        'GSDE     ',  &
        'SLC      ',  &
        'SOILGRIDS' &
        /)

   !# If .true., SVS1 simulates soil freezing and thawing and its impact on hydrology
   logical           :: lsoil_freezing_svs1 = .false.
   namelist /surface_cfgs/ lsoil_freezing_svs1 

   !# If .true., SVS1 simulates water ponding at the surface
   logical           :: lwater_ponding_svs1 = .false.
   namelist /surface_cfgs/ lwater_ponding_svs1 


   !# Use snow albedo "I6" directly if .true.;
   !# Use snow age "XA" to calculate snow albedo if .false.
   logical           :: snoalb_anl  = .true.
   namelist /surface_cfgs/ snoalb_anl

   !# use dynamic calculation of z0h for bare ground + vegetation  for SVS if .true.
   logical           :: svs_dynamic_z0h     = .false.
   namelist /surface_cfgs/ svs_dynamic_z0h 

   
   !# use hrsurf based on soil texture for SVS if .true.
   logical           :: svs_hrsurf_sltext     = .false.
   namelist /surface_cfgs/ svs_hrsurf_sltext

   !# use local momentum (no snow) roughness for SVS if .true.
   logical           :: svs_local_z0m     = .false.
   namelist /surface_cfgs/ svs_local_z0m  
   



   
   !# Limit temperature inversions to 8K/40m in surface layer if .true.
   logical           :: tdiaglim    = .false.
   namelist /surface_cfgs/ tdiaglim

   !# OPTION FOR CALCULATION of AVERAGE LAND SURFACE TEMPERATURE AND HUMIDITY IN SVS
   !# .FALSE. :  Area-average only calculation for sfc T and Hum.
   !# .TRUE.  :  Option that uses effective surface temperature  and specific humidity instead
   !             of composite (area-averaged only) counterparts in surface flux calculations (D. Deacu)
   logical           :: use_eff_surf_tq    = .false.
   namelist /surface_cfgs/ use_eff_surf_tq

   !# OPTION TO USE PHOTOSYNTHESIS CODE FOR STOMATAL RESISTANCE in SVS
   logical           :: use_photo = .true.
   namelist /surface_cfgs/ use_photo
 
   !# Factor multiplying stomatal resistance in ISBA
   real              :: veg_rs_mult = 1.
   namelist /surface_cfgs/ veg_rs_mult

   !#  VF definitions and mapping in SVS 
   !# * 'CLASSIC' : Same VF definitions as ISBA
   !# * 'CCILCECO' : New VF definitions used in SVS only
   !# with geo. fields generated using CCILC 2015 + Ecobiomes
   character(len=16) :: vf_type    = 'CLASSIC'
   namelist /surface_cfgs/ vf_type
   character(len=*), parameter :: VFTYPE_OPT(2) = (/ &
        'CLASSIC  ',  &
        'CCILCECO '   &
        /)

   !# Emissivity for water
   !# * '_constant_' : A fixed floating point value used as a constant
   character(len=16) :: water_emiss = '1.'
   real              :: water_emiss_const = -1.
   namelist /surface_cfgs/ water_emiss

   !# Use directional roughness length if .true.
   logical           :: z0dir       = .false.
   namelist /surface_cfgs/ z0dir

   !# Constant value of thermal roughness length (m) applied over water within
   !# latitudinal band defined by z0tlat
   real              :: z0hcon      = 4.0e-5
   namelist /surface_cfgs/ z0hcon

   !# Minimum value of momentum roughness length (m)
   real              :: z0min       = 1.5e-5
   namelist /surface_cfgs/ z0min

   !# Momentum roughness length formulation over water
   !# * 'CHARNOCK' : #TODO: define
   !# * 'BELJAARS' : #TODO: define
   character(len=16) :: z0mtype     = 'CHARNOCK'
   namelist /surface_cfgs/ z0mtype
   character(len=*), parameter :: Z0MTYPE_OPT(2) = (/ &
        'CHARNOCK', &
        'BELJAARS'  &
        /)

   !# Roughness length for sea ice
   real              :: z0seaice    = 1.6e-4
   namelist /surface_cfgs/ z0seaice
   
   !# Thermal roughness length formulation over water
   !# * 'MOMENTUM' : Uses z0h = z0m (replaces key z0trdps300=.false.)
   !# * 'DEACU12'  : #TODO: define  (replaces key z0trdps300=.true.)
   !# * 'ECMWF'    : #TODO: define  (New formulation used by ECMWF)
   character(len=16) :: z0ttype     = 'MOMENTUM'
   namelist /surface_cfgs/ z0ttype
   character(len=*), parameter :: Z0TTYPE_OPT(3) = (/ &
        'MOMENTUM', &
        'DEACU12 ', &
        'ECMWF   '  &
        /)

   !# Thermal roughness length formulation over vegetation
   !# * 'FIXED' : Uses z0h = z0m 
   !# * 'ZILI95': evolves with u*
   character(len=16) :: z0tevol     = 'FIXED'
   namelist /surface_cfgs/ z0tevol
   character(len=*), parameter :: Z0TEVOL_OPT(2) = (/ &
        'FIXED ', &
        'ZILI95'  &
        /)

   !# Latitude (2 elements, in degrees) used to specify Z0T over water
   !# * If |lat| <= Z0TLAT(1) constant Z0T.
   !# * If |lat| >= Z0TLAT(2) Charnock's relation.
   !# * In between, linear interpolation is used.
   real              :: z0tlat(2)   = 0.
   namelist /surface_cfgs/ z0tlat

   !# Height at which to compute screen-level temperature (m)
   real              :: zt = 1.5
   namelist /surface_cfgs/ zt

   !# Height at which to compute anemomenter-level winds (m)
   real              :: zu = 10.
   namelist /surface_cfgs/ zu

   !# New urban surface parameters
   logical           :: urban_params_new = .false.
   namelist /surface_cfgs/ urban_params_new
 
   !# Adjust wind diagnostic in TEB with building height  if .true.
   logical           :: urb_diagwind = .false.
!!$   namelist /surface_cfgs/ urb_diagwind

   !# Adjust temperature diagnostic in TEB in the street  if .true.
   logical           :: urb_diagtemp = .false.
!!$   namelist /surface_cfgs/ urb_diagtemp


   ! Option for the snowpack representation in SVS-2
!                        'ES': Explicit snow scheme 
!                        'CRO': Crocus
   character(len=16) :: hsnowscheme = 'ES '
   namelist /surface_cfgs/ hsnowscheme
   character(len=*), parameter :: HSNOWSCHEME_OPT(2) = (/ &
        'ES ',  &
        'CRO'  &  
         /)

   ! Option for the snowdrift scheme for Crocus: Mechanical transformation of snow grain and compaction 
   !  and effect of wind  on falling snow properties
   !    'NONE': No snowdrift scheme
   !     For the 3 other options the mehanical transformation is activated and different options are
   !         available for the properties of falling snow:
   !       'DFLT': falling snow falls as purely dendritic
   !       'GA01': Gallee et al 2001
   !       'VI13': Vionnet et al 2013 (Default in SVS2)
   character(len=16) :: hsnowdrift_cro = 'VI13'
   namelist /surface_cfgs/ hsnowdrift_cro
   character(len=*), parameter :: HSNOWDRIFT_CRO_OPT(4) = (/ &
        'NONE',  &
        'DFLT',  &
        'GA01',  &  
        'VI13'  &  
         /)

   ! Option for the snowdrift scheme for ES: Mechanical transformation of snow grain and compaction 
   !    'DFLT': Snowdrift scheme activated  (Default in SVS2)
   !    'NONE': No snowdrift scheme
   character(len=16) :: hsnowdrift_es = 'DFLT'
   namelist /surface_cfgs/ hsnowdrift_es
   character(len=*), parameter :: HSNOWDRIFT_ES_OPT(2) = (/ &
        'NONE',  &
        'DFLT'  &
         /)

   !# Activate mass loss due to blowing snow sublimation in ES and Crocus (default TRUE in SVS2)
   logical           :: lsnowdrift_sublim = .true.
   namelist /surface_cfgs/ lsnowdrift_sublim

   ! Option for the metamorphism scheme for Crocus
   !    HSNOWMETAMO=B21 Correction of C13 to correctly handle the conversion from 
   !                     dendricity/sphericity/grain size to optical diameter/sphericity (Default in SVS2)
   !    HSNOWMETAMO=C13 Carmagnola et al 2014 
   !    HSNOWMETAMO=T07 Taillandier et al 2007
   !    HSNOWMETAMO=F06 Flanner et al 2006
   character(len=16) :: hsnowmetamo = 'B21'
   namelist /surface_cfgs/ hsnowmetamo
   character(len=*), parameter :: HSNOWMETAMO_OPT(4) = (/ &
        'C13',  &
        'T07',  &  
        'F06',  &  
        'B21'  &  
         /)

   ! Option for the radiative transfer scheme for Crocus
   !    HSNOWRAD=B92 Brun et al 1992  (Default in SVS-Cro)
   !    HSNOWRAD=T17 (Tuzet et al. 2017) (Libois et al. 2013) TARTES with impurities content scheme
   character(len=16) :: hsnowrad = 'B92'
   namelist /surface_cfgs/ hsnowrad
   character(len=*), parameter :: HSNOWRAD_OPT(2) = (/ &
        'B92',  &
        'T17'   &  
         /)

   ! Option for the falling snow density for Crocus
   !   HSNOWFALL=V12 Vionnet et al. 2012 from Brun et al. 1989 (Default in SVS2)
   !   HSNOWFALL=A76 Anderson et al. 1976
   !   HSNOWFALL=S02 Lehning el al. 2002
   !   HSNOWFALL=P75 Pahaut 1975
   !   HSNOWFALL=NZE Constant density 200 kg/m3 (defined snowcro.F90 )         
   character(len=16) :: hsnowfall = 'V12'
   namelist /surface_cfgs/ hsnowfall
   character(len=*), parameter :: HSNOWFALL_OPT(5) = (/ &
        'V12',  &
        'A76',  &  
        'S02',  &  
        'P75',  &  
        'NZE'  &  
         /)

   ! Option for the thermal conductivity scheme for Crocus
   !    HSNOWCOND=Y81 default Crocus from Yen et al. 1981 (Default in SVS2)
   !    HSNOWCOND=I02 ISBA_ES snow conductivity parametrization (Boone et al. 2002)
   !    HSNOWCOND=C11 Calonne et al. 2011 snow conductivity parametrization
   character(len=16) :: hsnowcond = 'Y81'
   namelist /surface_cfgs/ hsnowcond
   character(len=*), parameter :: HSNOWCOND_OPT(3) = (/ &
        'Y81',  &
        'I02',  &  
        'C11'   &  
         /)

   ! Option for the liquid water content scheme for Crocus
   !   HSNOWHOLD=B92 default Crocus from Brun et al. 1992 or Vionnet et al. 2012 (Default in SVS2)
   !   HSNOWHOLD=B02 ISBA_ES  parametrization (Boone et al. 2002)
   !   HSNOWHOLD=O04 CLM parametrization (Oleson et al 2004)
   !   HSNOWHOLD=SPK SNOWPACK parametrization (Lehning et al 2002)    
   character(len=16) :: hsnowhold = 'B92'
   namelist /surface_cfgs/ hsnowhold
   character(len=*), parameter :: HSNOWHOLD_OPT(4) = (/ &
        'B92',  &
        'B02',  &  
        'O04',  &  
        'SPK'   &   
         /)

   ! Option for the compaction scheme for Crocus
   !   HSNOWCOMP=B92 default Crocus from Brun et al. 1992 or Vionnet et al. 2012 (Default in SVS2)
   !   HSNOWCOMP=S14 use the settling param of Schleef et al. (2014) for fresh snow (less than 2 days) 
   !   HSNOWCOMP=T11 param of snow viscosity from Teufelsbauer (2011)   
   character(len=16) :: hsnowcomp = 'B92'
   namelist /surface_cfgs/ hsnowcomp
   character(len=*), parameter :: HSNOWCOMP_OPT(3) = (/ &
        'B92',  &
        'S14',  &  
        'T11'   &  
         /)

   ! Option for the turbulent fluxes in Crocus
!    HSNOWRES = 'DEF' = Default: Louis (ISBA: Noilhan and Mahfouf 1996)
!    HSNOWRES = 'RIL' = Limit Richarson number under very stable conditions (currently testing)
!    HSNOWRES = 'M98'  = Martin et Lejeune 1998 : older computation for turbulent fluxes coefficents in Crocus   
   character(len=16) :: hsnowres = 'RIL'
   namelist /surface_cfgs/ hsnowres
   character(len=*), parameter :: HSNOWRES_OPT(3) = (/ &
        'DEF',  &
        'RIL',  &  
        'M98'   &  
         /)


contains

   function sfc_options_init() result(F_istat)
      use sfclayer_mod, only: sl_get, SL_OK
      implicit none
      integer :: F_istat
#include <msg.h>
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.
      F_istat = sl_get('beta',beta)
      F_istat = min(sl_get('bh91_a',bh91_a),F_istat)
      F_istat = min(sl_get('bh91_b',bh91_b),F_istat)
      F_istat = min(sl_get('bh91_c',bh91_c),F_istat)
      F_istat = min(sl_get('bh91_d',bh91_d),F_istat)
      F_istat = min(sl_get('d97_as',d97_as),F_istat)
      F_istat = min(sl_get('dg92_ci',dg92_ci),F_istat)
      F_istat = min(sl_get('l07_ah',l07_ah),F_istat)
      F_istat = min(sl_get('l07_am',l07_am),F_istat)
      F_istat = min(sl_get('rineutral',sl_rineutral),F_istat)
      if (.not.RMN_IS_OK(F_istat)) &
           call msg(MSG_ERROR,'(sfc_options_init) cannot retrieve AS, BETA, CI or NEUTRAL')
      return
   end function sfc_options_init

end module sfc_options
