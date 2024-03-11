!SFX_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!SFX_LIC This is part of the SURFEX software governed by the CeCILL-C licence
!SFX_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!SFX_LIC for details. version 1.
SUBROUTINE INI_SURF_CSTS 
!     ##################
!
!!****  *INI_SURF_CSTS * - routine to initialize all surface parameter as
!!                         emissivity and albedo
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!      The physical constants are set to their default numerical values 
!!      or specified in namelist NAM_SURF_CSTS
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      B. Decharme       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    08/2009
!!      M Lafaysse 05/2014 : snow parameters
!!      B. Decharme    05/13 : Add NAM_SURF_REPROD_OPER for versions reproductibility
!!      P. Samuelsson 10/2014 MEB
!!      B. Decharme    01/16 : Update XCFFV
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_SURF_PAR,  ONLY : XUNDEF
!VV #ifndef CROCUS_EXT
!VV USE MODD_SURF_CONF, ONLY : CPROGNAME
!VV #endif
!
!VV USE MODD_WATER_PAR
!VV USE MODD_FLOOD_PAR
!VV USE MODD_MEB_PAR,   ONLY : XTAU_LW,                            &
!VV                            XRAGNC_FACTOR, XKDELTA_WR
USE MODD_SNOW_PAR,  ONLY : XEMISSN, XANSMIN, XANSMAX,          &
                           XAGLAMIN, XAGLAMAX, XHGLA,          &
                           XWSNV, XZ0SN, XZ0HSN,               &
                           X_RI_MAX,                           &
                           XTAU_SMELT,                         &
                           XALBICE1, XALBICE2, XALBICE3,       &
                           XRHOTHRESHOLD_ICE, XZ0ICEZ0SNOW,    &
                           XVAGING_NOGLACIER, XVAGING_GLACIER, &
                           XPERCENTAGEPORE,                    &
                           XPERCENTAGEPORE_FRZ,                &
                           XPERCENTAGEPORE_ICE,                &                           
                           LMEBREC,                            &
                           XANSFRACMEL, XTEMPANS, XANSMINMEB,  &
                           XIMPUR_WET, XIMPUR_DRY,          &
                           XPSR_SNOWMAK, XRHO_SNOWMAK,         &
                           XPTA_SEUIL, XPR_A, XPR_B, XPT,      &
                           XPP_D1, XPP_D2, XPP_D3, XPP_H1,     &
                           XPP_H2, XPP_H3, XPP_H4, XWT, XPTR , &
                           XTIMESNOWMAK,                       &
                           XPROD_SCHEME, XSM_END, XFREQ_GRO,   & !Grooming and Snowmaking option by P.Spandre 20160211
                           XSCAVEN_COEF
!
USE MODD_SNOW_METAMO, ONLY : XVVISC3
!
USE MODD_SURF_ATM, ONLY: XRIMAX,LALDTHRES,XCISMIN, XVMODMIN

!VV #ifndef CROCUS_EXT
!VV USE MODI_GET_LUOUT
!VV USE MODI_OPEN_NAMELIST
!VV USE MODI_CLOSE_NAMELIST
!VV USE MODE_POS_SURF
!VV !
 USE MODD_REPROD_OPER,  ONLY : XEVERG_RSMIN, XEVERG_VEG, &
                                    CDGAVG, CIMPLICIT_WIND,   &
                                    CQSAT, CCHARNOCK, CDGDIF
!VV USE MODI_TEST_NAM_VAR_SURF
!VV !
!VV #endif
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
USE PARKIND1  ,ONLY : JPRB
!
IMPLICIT NONE
!
INTEGER               :: ILUOUT    ! unit of output listing file
INTEGER               :: ILUNAM    ! namelist file  logical unit
LOGICAL               :: GFOUND    ! true if namelist is found
!
LOGICAL               :: LREPROD_OPER
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!VV #ifndef CROCUS_EXT
!VV NAMELIST/NAM_SURF_CSTS/ XEMISSN, XANSMIN, XANSMAX, XAGLAMIN, XAGLAMAX, &
!VV                         XALBWAT, XALBCOEF_TA96, XALBSCA_WAT, XEMISWAT, &
!VV                         XALBWATICE, XEMISWATICE, XHGLA, XWSNV, XCFFV,  &
!VV                         XZ0SN, XZ0HSN, XTAU_SMELT, XALBSEAICE,         &
!VV                         XZ0FLOOD, XALBWATSNOW,                         &
!VV                         LMEBREC,                                       &
!VV                         XANSFRACMEL, XTEMPANS, XANSMINMEB,             &
!VV                         XTAU_LW, XRAGNC_FACTOR
!VV !
!VV NAMELIST/NAM_SURF_SNOW_CSTS/ XZ0ICEZ0SNOW, XRHOTHRESHOLD_ICE,          &
!VV                              XALBICE1, XALBICE2, XALBICE3,             &
!VV                              XVAGING_NOGLACIER, XVAGING_GLACIER,       &
!VV                              XPERCENTAGEPORE, XPERCENTAGEPORE_FRZ,     &
!VV                              XVVISC3,X_RI_MAX,                         &
!VV                              XIMPUR_WET, XIMPUR_DRY, XPSR_SNOWMAK,     &
!VV                              XRHO_SNOWMAK, XPTA_SEUIL,                 &
!VV                              XPR_A, XPR_B, XPT, XTIMESNOWMAK,          &
!VV                              XPP_D1, XPP_D2, XPP_D3, XPP_H1,           &
!VV                              XPP_H2, XPP_H3, XPP_H4, XWT, XPTR ,       &
!VV                              XPROD_SCHEME, XSM_END, XFREQ_GRO,         &
!VV                              XSCAVEN_COEF
!VV !
!VV NAMELIST/NAM_REPROD_OPER/ LREPROD_OPER, XEVERG_RSMIN, XEVERG_VEG, &
!VV                           CDGAVG, CDGDIF, CIMPLICIT_WIND, CQSAT,  &
!VV                           CCHARNOCK
!VV #endif
!
!-------------------------------------------------------------------------------
!*       0. INIT
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('INI_SURF_CSTS',0,ZHOOK_HANDLE)
!VV !
!VV XALBWAT     = XUNDEF
!VV XALBSEAICE  = XUNDEF
!VV XALBWATICE  = XUNDEF
!VV XALBWATSNOW = XUNDEF
!VV XEMISWAT    = XUNDEF
!VV XEMISWATICE = XUNDEF
XEMISSN     = XUNDEF
!
!-------------------------------------------------------------------------------
!*       1. Default values
!-------------------------------------------------------------------------------
!
! Minimum and maximum values of the albedo of snow:
!
XANSMIN = 0.50 ! (-)
XANSMAX = 0.85 ! (-)
!
! Minimum and maximum values of the albedo of permanet snow/ice:
!
XAGLAMIN = 0.8 ! (-)
XAGLAMAX = 0.85 ! (-)
!
! Use recommended settings for snow albedo (FALSE = ISBA default)
! 
!VV LMEBREC=.FALSE.
!
! Fraction of maximum value of the albedo of snow that is reached for melting
! snow
!
XANSFRACMEL = 1.0 ! (-)
!
! Threeshold temperature above which the snow albedo starts to decrease 
!
XTEMPANS = 274.15 ! (K)
!
! Minimum value of the albedo of snow reached under canopy vegetation:
!
!VV XANSMINMEB = 0.30 ! (-)
!
! Height of aged snow in glacier case (allows Pn=1)
!
XHGLA    = 33.3 !(m)
! 
! Coefficient for calculation of snow fraction over vegetation
!
XWSNV = 5.0 !(-)
!
! Water direct albedo coefficient (option "TA96")
!
!VV XALBCOEF_TA96 =  0.037
!
! Water diffuse albedo
!
!VV XALBSCA_WAT =  0.06

! Coefficient for calculation of floodplain fraction over vegetation
!
!VV XCFFV = 4.0
!
! Roughness length of pure snow surface (m)
!
XZ0SN = 0.001
!
! Roughness length for heat of pure snow surface (m)
!
XZ0HSN = 0.0001
!
! Maximum Richardson number limit for very stable conditions over snow using the 'RIL' option
X_RI_MAX = 0.2
!X_RI_MAX = 0.026
!
! Snow Melt timescale with D95 (s): needed to prevent time step 
! dependence of melt when snow fraction < unity.
!
XTAU_SMELT = 300.
!
! Extinction coefficient for view factor for long-wave radiation 
!
!VV XTAU_LW = 0.5   ! -
!
! MEB resistance increase factor for canopy air sapce.
! If=1, then NO effect. It is generally >=1
! and is needed because the original parameterization
! does not account for extremely stable conditions,
! such as over a snowpack.
!
!VV XRAGNC_FACTOR= 200. ! -
!
! MEB maximum intercepted water fraction (on vegetation)
!
!VV XKDELTA_WR   = 0.25 ! -
!
! NAM_SURF_SNOW_CSTS
!
! Roughness length ratio between ice and snow
XZ0ICEZ0SNOW = 10.
!
!
! 3 bands spectral albedo for glacier ice (CROCUS)
! Default values from Lejeune et al 2009 (Zongo, Bolivia)
XALBICE1 = 0.38
XALBICE2 = 0.23
XALBICE3 = 0.08
!
! Gerbaux et al 2005 (Saint Sorlin)
! PALBICE1=0.23
! PALBICE2=0.16
! PALBICE3=0.05
!
! Options for MM snow production and grooming p.s 20160211
XPSR_SNOWMAK = 0.
XRHO_SNOWMAK = 0.
XPTA_SEUIL = 0.
XPR_A = 0.
XPR_B = 0.
XPT = 0.
XPP_D1 = 0.
XPP_D2 = 0.
XPP_D3 = 0.
XPP_H1 = 0.
XPP_H2 = 0.
XPP_H3 = 0.
XPP_H4 = 0.
XWT = 0.
XPTR = 0.
XTIMESNOWMAK = 0.
XPROD_SCHEME = 0.
XSM_END = (/0,0,0,0/)
XFREQ_GRO = 0
!
! Density threshold for ice detection kg.m-3
XRHOTHRESHOLD_ICE = 850.
!
! Parameters for ageing effect on albedo
XVAGING_NOGLACIER = 60.
XVAGING_GLACIER   = 900.

! percentage of the total pore volume to compute the max liquid water holding capacity   !Pahaut 1976
XPERCENTAGEPORE = 0.05
XPERCENTAGEPORE_FRZ = 1.
XPERCENTAGEPORE_ICE = 0.5
!
! Snow viscosity coefficient
XVVISC3= 0.023
!
! Roughness length for flood (m)
!
!VV XZ0FLOOD = 0.0002

!!! impurity value 
XIMPUR_DRY(1)=0. ! BC dry deposition at top of snowpack (g m-2 s-1)
XIMPUR_WET(1)=0.! BC Wet deposition of with precipitation (g m-2 s-1)

XIMPUR_DRY(2)=0. ! Dust dry deposition at top of snowpack (g m-2 s-1)
XIMPUR_WET(2)=0. ! Dust Wet deposition of with precipitation (g m-2 s-1)

XIMPUR_DRY(3:5)=0. ! Other types of impurities dry deposition at top of snowpack (g m-2 s-1)
XIMPUR_WET(3:5)=0. ! Other types of impurities deposition of with precipitation (g m-2 s-1)

! Scavenging coefficient of impurities

XSCAVEN_COEF=(/0.0,0.0,0.,0.,0./) 



! VV Add init needed for the externalized version of Crocus
XRIMAX=0.2
LALDTHRES = .FALSE. 
XCISMIN = 0. ! Not used when LALDTHRES = .FALSE.
XVMODMIN  = 0. ! Not used when LALDTHRES = .FALSE.


!VV #ifndef CROCUS_EXT
!VV !-------------------------------------------------------------------------------
!VV !
!VV ! * Reproductibility for SURFEX OPER
!VV !
!VV LREPROD_OPER = .FALSE. ! default
!VV !
!VV ! * Vegetation parameters for tropical forest
!VV !
!VV !XEVERG_RSMIN : old = 250. (Manzi 1993) but observations range 
!VV !               from 140 to 180. According to Delire et al. (1997) and 
!VV !               new tests over 6 local sites, 175. is recommended
!VV !               Should be the default after check with AROME/ALADIN
!VV !
!VV XEVERG_RSMIN = 175.  !Rsmin
!VV !
!VV !XEVERG_VEG : old = 0.99 (Manzi 1993) but according to Delire et al. (1997) and 
!VV !             new tests over 6 local sites, 1.0 is recommended because 0.99
!VV !             induces unrealistic bare soil evaporation for Tropical forest
!VV !             Should be the default after check with AROME/ALADIN
!VV !
!VV XEVERG_VEG   = 1.0  !Veg fraction
!VV !
!VV ! * Soil depth average
!VV !
!VV CDGAVG = 'INV'
!VV !
!VV ! * Soil depth with ISBA-DF
!VV !
!VV CDGDIF = 'ROOT'
!VV !
!VV ! * wind implicitation option
!VV !
!VV CIMPLICIT_WIND = 'NEW'

!
! * qsat computation
!
CQSAT = 'NEW'

!VV !
!VV ! * Charnock parameter
!VV !
!VV CCHARNOCK = 'NEW'
!VV !
!VV !-------------------------------------------------------------------------------
!VV !*       2. User values
!VV !-------------------------------------------------------------------------------
!VV !
!VV  CALL GET_LUOUT(CPROGNAME,ILUOUT)
!VV !    
!VV  CALL OPEN_NAMELIST(CPROGNAME,ILUNAM)
!VV !
!VV  CALL POSNAM(ILUNAM,'NAM_SURF_CSTS',GFOUND,ILUOUT)
!VV IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_SURF_CSTS)
!VV !
!VV IF(LMEBREC)THEN
!VV ! Fraction of maximum value of the albedo of snow that is reached for melting
!VV ! snow
!VV !
!VV   XANSFRACMEL = 0.85 ! (-)
!VV !
!VV ! Threeshold temperature above which the snow albedo starts to decrease 
!VV !
!VV   XTEMPANS = 268.15 ! (K)
!VV !
!VV ENDIF
!VV !
!VV  CALL POSNAM(ILUNAM,'NAM_SURF_SNOW_CSTS',GFOUND,ILUOUT)
!VV IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_SURF_SNOW_CSTS)
!VV !
!VV !-------------------------------------------------------------------------------
!VV !*       3. For Reproductibility
!VV !-------------------------------------------------------------------------------
!VV !
!VV  CALL POSNAM(ILUNAM,'NAM_REPROD_OPER',GFOUND,ILUOUT)
!VV IF (GFOUND) READ(UNIT=ILUNAM,NML=NAM_REPROD_OPER)
!VV !
!VV  CALL TEST_NAM_VAR_SURF(ILUOUT,'CDGAVG',CDGAVG,'ARI','INV')
!VV  CALL TEST_NAM_VAR_SURF(ILUOUT,'CDGDIF',CDGDIF,'SOIL','ROOT')
!VV  CALL TEST_NAM_VAR_SURF(ILUOUT,'CIMPLICIT_WIND',CIMPLICIT_WIND,'NEW','OLD')
!VV  CALL TEST_NAM_VAR_SURF(ILUOUT,'CQSAT',CIMPLICIT_WIND,'NEW','OLD')
!VV  CALL TEST_NAM_VAR_SURF(ILUOUT,'CCHARNOCK',CIMPLICIT_WIND,'NEW','OLD')
!VV !
!VV  CALL TEST_NAM_VAR_SURF(ILUOUT,'XEVERG_RSMIN',XEVERG_RSMIN,175.0,250.0)
!VV  CALL TEST_NAM_VAR_SURF(ILUOUT,'XEVERG_VEG',XEVERG_VEG,1.0,0.99) 
!VV !
!VV IF(LREPROD_OPER)THEN
!VV   XEVERG_RSMIN   = 250.
!VV   XEVERG_VEG     = 0.99
!VV   CDGAVG         = 'ARI'
!VV   CQSAT          = 'OLD'
!VV   CCHARNOCK      = 'OLD'
!VV ENDIF
!VV !
!VV ! Water global albedo (option "UNIF")
!VV !
!VV IF(XALBWAT==XUNDEF)THEN
!VV   IF(LREPROD_OPER)THEN
!VV     XALBWAT =  0.135
!VV   ELSE
!VV     XALBWAT =  0.065
!VV   ENDIF
!VV ENDIF
!VV !
!VV ! Sea ice albedo
!VV !
!VV IF(XALBSEAICE==XUNDEF)THEN
!VV   IF(LREPROD_OPER)THEN
!VV     XALBSEAICE =  0.85
!VV   ELSE
!VV     XALBSEAICE =  0.71
!VV   ENDIF
!VV ENDIF
!VV !
!VV ! water ice and snow albedo
!VV !
!VV IF(XALBWATICE==XUNDEF)THEN
!VV   IF(LREPROD_OPER)THEN
!VV     XALBWATICE =  0.85
!VV   ELSE
!VV     XALBWATICE =  0.40
!VV   ENDIF
!VV ENDIF
!VV !
!VV IF(XALBWATSNOW==XUNDEF)THEN
!VV   IF(LREPROD_OPER)THEN
!VV     XALBWATSNOW =  0.85
!VV   ELSE
!VV     XALBWATSNOW =  0.60
!VV   ENDIF
!VV ENDIF
!VV !                   
!VV ! Water emissivity
!VV !
!VV IF(XEMISWAT==XUNDEF)THEN
!VV   IF(LREPROD_OPER)THEN
!VV     XEMISWAT =  0.98
!VV   ELSE
!VV     XEMISWAT =  0.96
!VV   ENDIF
!VV ENDIF
!VV !
!VV ! Sea ice emissivity
!VV !
!VV IF(XEMISWATICE==XUNDEF)THEN
!VV   IF(LREPROD_OPER)THEN
!VV     XEMISWATICE =  1.0
!VV   ELSE
!VV     XEMISWATICE =  0.97
!VV   ENDIF
!VV ENDIF
!VV !
!VV !
!VV ! Snow emissivity:
!VV !
 IF(XEMISSN==XUNDEF)THEN
!VV   IF(LREPROD_OPER)THEN
!VV     XEMISSN =  1.0
!VV   ELSE
     XEMISSN =  0.99
!VV   ENDIF
 ENDIF
!VV !
!VV !-------------------------------------------------------------------------------
!VV !
!VV  CALL CLOSE_NAMELIST(CPROGNAME,ILUNAM)
!VV !
!VV #endif
IF (LHOOK) CALL DR_HOOK('INI_SURF_CSTS',1,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_SURF_CSTS 
