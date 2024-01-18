!auto_modi:spll_snowcro.D
MODULE MODI_SNOWCRO
INTERFACE
      SUBROUTINE SNOWCRO(HSNOWRES, TPTIME, OMEB, OGLACIER, HIMPLICIT_WIND,      &
                      PPEW_A_COEF, PPEW_B_COEF,                                 &
                      PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF, PPEQ_B_COEF,       &
                      PSNOWSWE,PSNOWRHO,PSNOWHEAT,PSNOWALB,                     &
                      PSNOWGRAN1,PSNOWGRAN2,PSNOWHIST,PSNOWAGE, PSNOWIMPUR,     &
                      PTSTEP,PPS,PSR,PRR,PPSN3L, PRSURF,                             &
                      PTA,PTG,PSW_RAD,PQA,PVMOD,PWIND_DRIFT, PLW_RAD, PRHOA,    &
                      PUREF,PEXNS,PEXNA,PDIRCOSZW,                              &
                      PZREF,PZ0,PZ0EFF,PZ0H,PALB,                               &
                      PSOILCOND,PD_G,                                           &
                      PSNOWLIQ,PSNOWTEMP,PSNOWDZ,                               &
                      PTHRUFAL,PGRNDFLUX,PEVAPCOR, PGFLXCOR,                    &
                      PSWNETSNOW,PSWNETSNOWS,PLWNETSNOW,                        &
                      PRNSNOW,PHSNOW,PGFLUXSNOW,                                &
                      PHPSNOW,PLES3L,PLEL3L,PEVAP,PSNDRIFT,PRI,                 &
                      PEMISNOW,PCDSNOW,PUSTAR,PCHSNOW,PSNOWHMASS,PQS,           &
                      PPERMSNOWFRAC,PZENITH,PANGL_ILLUM,PXLAT,PXLON,PBLOWSNW,   &
                      HSNOWDRIFT,OSNOWDRIFT_SUBLIM,OSNOW_ABS_ZENITH,            &
                      HSNOWMETAMO,HSNOWRAD,OATMORAD,P_DIR_SW, P_SCA_SW,         &
                      PSPEC_ALB, PDIFF_RATIO,PSPEC_TOT,PSNOWFLUX,PIMPWET,PIMPDRY,&
                      HSNOWFALL, HSNOWCOND,HSNOWHOLD,HSNOWCOMP,HSNOWZREF,       &
                      PSNOWMAK, OSNOWCOMPACT_BOOL, OSNOWMAK_BOOL, OSNOWTILLER,  &
                      OSELF_PROD, OSNOWMAK_PROP  )
USE MODD_TYPE_DATE_SURF, ONLY: DATE_TIME
IMPLICIT NONE
REAL, INTENT(IN)                       :: PTSTEP
TYPE(DATE_TIME), INTENT(IN)            :: TPTIME      ! current date and time
 CHARACTER(LEN=*), INTENT(IN)          :: HSNOWRES
LOGICAL, INTENT(IN)                    :: OMEB       ! True = coupled to MEB. This means surface fluxes ae IMPOSED
LOGICAL, INTENT(IN)                    :: OGLACIER   ! True = Over permanent snow and ice,
 CHARACTER(LEN=*), INTENT(IN)          :: HIMPLICIT_WIND   ! wind implicitation option
REAL, DIMENSION(:), INTENT(IN)         :: PPS, PTA, PSW_RAD, PQA, PVMOD, PWIND_DRIFT, PLW_RAD, PSR, PRR
REAL, DIMENSION(:,:), INTENT(IN)       :: P_DIR_SW, P_SCA_SW ! direct and diffuse spectral irradiance (W/m2/um)
REAL, DIMENSION(:,:), INTENT(IN)       :: PIMPWET, PIMPDRY  !Dry and wet deposit coefficient from Forcing File(g/m²/s)
REAL, DIMENSION(:), INTENT(IN)         :: PTG, PSOILCOND, PD_G, PPSN3L
REAL, DIMENSION(:), INTENT(IN)         :: PZREF, PUREF, PEXNS, PEXNA, PDIRCOSZW, PRHOA, PZ0, PZ0EFF, &
                                       PALB, PZ0H, PPERMSNOWFRAC
REAL, DIMENSION(:), INTENT(IN)         :: PPEW_A_COEF, PPEW_B_COEF,                   &
                                        PPET_A_COEF, PPEQ_A_COEF, PPET_B_COEF,      &
                                        PPEQ_B_COEF
REAL, DIMENSION(:), INTENT(OUT)      :: PSNOWALB
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWHEAT, PSNOWRHO, PSNOWSWE
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWGRAN1, PSNOWGRAN2, PSNOWHIST
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PSNOWAGE  ! Snow grain age
REAL, DIMENSION(:,:,:), INTENT(INOUT)  :: PSNOWIMPUR  ! Snow impurity content (g) (LOCATION,LAYER,NIMPUR)) Impur type :1/BC 2/Dust
REAL, DIMENSION(:,:), INTENT(OUT)    :: PSNOWTEMP
REAL, DIMENSION(:,:), INTENT(OUT)      :: PSNOWLIQ, PSNOWDZ
REAL, DIMENSION(:), INTENT(OUT)        :: PTHRUFAL,  PEVAPCOR,  PGFLXCOR
REAL, DIMENSION(:,:), INTENT(OUT)      :: PSPEC_ALB, PDIFF_RATIO,PSPEC_TOT  !! spectral albedo, diffuse to total irradiance ratio and total incoming spectral irradiance (npoints,nbands)
REAL, DIMENSION(:), INTENT(OUT)        ::   PSNOWFLUX
REAL, DIMENSION(:), INTENT(INOUT)      :: PGRNDFLUX
REAL, DIMENSION(:), INTENT(INOUT)      :: PRNSNOW, PHSNOW, PGFLUXSNOW, PLES3L, PLEL3L, &
                                          PHPSNOW, PCDSNOW, PUSTAR, PEVAP
REAL, DIMENSION(:), INTENT(OUT)        :: PSNDRIFT
REAL, DIMENSION(:), INTENT(INOUT)      :: PSWNETSNOW, PLWNETSNOW, PSWNETSNOWS
REAL, DIMENSION(:), INTENT(INOUT)      :: PCHSNOW,PRI
REAL, DIMENSION(:), INTENT(OUT)        :: PEMISNOW, PSNOWHMASS
REAL, DIMENSION(:), INTENT(OUT)        ::  PQS
REAL, DIMENSION(:), INTENT(IN)         :: PZENITH ! solar zenith angle
REAL, DIMENSION(:), INTENT(IN)         :: PANGL_ILLUM ! Effective illumination angle, Angle between the sun and the normal to the ground (=zenith if no slope) used in TARTES
REAL, DIMENSION(:), INTENT(IN)         :: PXLAT,PXLON ! LAT/LON after packing
REAL, DIMENSION(:,:), INTENT(INOUT)    :: PBLOWSNW !  Properties of deposited blowing snow (from Sytron or Meso-NH/Crocus)
 CHARACTER(4), INTENT(IN)              :: HSNOWDRIFT        ! Snowdrift scheme :
LOGICAL, INTENT(IN)                    :: OSNOWDRIFT_SUBLIM ! activate sublimation during drift
REAL, DIMENSION (:), INTENT(IN)        ::  PSNOWMAK        ! Snowmaking thickness (m)
LOGICAL, INTENT(IN)                    :: OSNOWCOMPACT_BOOL, OSNOWMAK_BOOL, OSNOWTILLER, &
                                          OSELF_PROD, OSNOWMAK_PROP
LOGICAL, INTENT(IN)                    :: OSNOW_ABS_ZENITH ! activate parametrization of solar absorption for polar regions
 CHARACTER(3), INTENT(IN)              :: HSNOWMETAMO, HSNOWRAD, HSNOWFALL, HSNOWCOND, HSNOWHOLD, HSNOWCOMP, HSNOWZREF
LOGICAL, INTENT(IN)                    :: OATMORAD ! activate atmotartes scheme
REAL, DIMENSION(:), INTENT(IN)         :: PRSURF
END SUBROUTINE SNOWCRO
END INTERFACE
END MODULE MODI_SNOWCRO
