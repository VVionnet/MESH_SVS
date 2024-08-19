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

      SUBROUTINE DRAG_SVS2 ( TGRS, TGRVS, TVGLS, TVGHS, TSV,  WD1, &
                              WR_VL, WR_VH,  THETAA, VMOD, VDIR, HU, RHOA, &
                              PS, RS, Z0, Z0LOC, Z0VG, WFC, WSAT, CLAY1,  &
                              SAND1, LAI_VL, LAI_VH, WRMAX_VL,WRMAX_VH,&
                              ZUSL, ZTSL, LAT, PSNVH, &
                              FCOR, Z0HA, VEGL, VEGH, &
                              VGH_DENS, Z0MVH,Z0MVL, Z0SNOW, Z0HSN,  VEGHEIGHT,  &
                              LAIVH, ZVCAN, FCANS,SNCMA, &
                              RESAGR,RESAGRV, RESA_VL, RESA_VH, RES_SNCA, RESA_SV, &
                              HUSURF,HUSURFGV, &
                              HRSURF,HRSURFGV, &
                              HV_VL, HV_VH, HVSN_VH, DEL_VL, DEL_VH,  &
                              Z0HBG, Z0HVL, Z0HVH,Z0HGV,  N)
      use tdpack
      use sfc_options
      use sfclayer_mod, only: sl_sfclayer,SL_OK
      use svs_configs
      use MODE_THERMOS
      use MODD_CSTS
      use canopy_csts, only: ZRALAI, RCHD, RADIUS_ICESPH, ZVENT, lres_snca,&
                             ZBETA
!
      implicit none
!!!#include <arch_specific.hf>

      INTEGER N
      REAL TGRS(N), TVGLS(N), TVGHS(N), WR_VL(N), WR_VH(N), THETAA(N), VMOD(N), VDIR(N), HU(N), CLAY1(N)
      REAL SAND1(N), PS(N), RS(N), Z0(N),Z0LOC(N), Z0VG(N), WFC(N,NL_SVS), WSAT(N,NL_SVS)
      REAL LAI_VL(N), LAI_VH(N),WRMAX_VH(N),WRMAX_VL(N), ZUSL(N), ZTSL(N), LAT(N)
      REAL FCOR(N), Z0HA(N), Z0HBG(N), Z0HVL(N), Z0HVH(N), Z0HGV(N), FCANS(N)
      REAL VGH_DENS(N), Z0MVH(N), Z0MVL(N), VEGHEIGHT(N), LAIVH(N), ZVCAN(N)
      REAL RESAGR(N),RESAGRV(N), RESA_VL(N), RESA_VH(N), RES_SNCA(N), RESA_SV(N)
      REAL HUSURF(N),HUSURFGV(N), HV_VL(N), HV_VH(N), HVSN_VH(N), DEL_VL(N), DEL_VH(N)
      REAL HRSURF(N),HRSURFGV(N), WD1(N), Z0SNOW(N), Z0HSN(N), PSNVH(N)
      REAL VEGH(N), VEGL(N), TGRVS(N), SNCMA(N), RHOA(N)
      REAL TSV(N,NSL)
!
!Author
!          S. Belair, M.Abrahamowicz, S.Z.Husain, N.Alavi, S.Zhang (June 2015)
!Revisions
! 001      Name (date) - Comment
!
!
!Object
!
!     Calculates the drag coefficients for heat and momentum transfers
!     over ground and vegetation(i.e., Ch and Cd).
!
!
!Method
!
!
!     1) computes hu, hv, and DEL
!
!     2) use this to find qsoil, the grid-averaged relative humidity
!        of the soil
!
!     3) find the transfer and resistance coefficients Ch, Cd, and Ra
!        Calculate the surface fluxes of heat, moisture,
!        and momentum over water surfaces.
!
!Arguments
!
!          - Input/Output -
! RESAGR    aerodynamical surface resistance for bare ground
! RESAGRV    aerodynamical surface resistance for bare ground
! RESA_VL    aerodynamical surface resistance for low vegetation
! RESA_VH    aerodynamical surface resistance for high vegetation
!
!          - Input -
! TGRS      skin (surface) temperature of bare ground
! TGRVS      skin (surface) temperature of ground below high veg.
! TVGLS      skin (surface) temperature of low vegetation
! TVGHS      skin (surface) temperature of high vegetation
! TSV       Temperature of the top snow layer in the high vegetation
! WD1       Soil volumetric water content (first level)
! WR_VL     water content retained by the low vegetation canopy
! WR_VH     water content retained by the high vegetation canopy
! THETAA    potential temperature at the lowest level
! VMOD      wind speed at the lowest level
! VDIR      wind direction at the lowest level
! HU        specific humidity of air at the lowest level
! PS        surface pressure
! RS        surface or stomatal resistance
! Z0        momentum roughness length (no snow)
! Z0LOC     local (no orography) land momentum roughness
! Z0VG      averaged vegetation-only momentum roughness
! WFC       volumetric water content at the field capacity
! WSAT      volumetric water content at saturation
! CLAY1     percentage of clay in first soil layer
! SAND1     percentage of sand in first soil layer
! LAI_VL    leaf Area Index of low vegetation
! LAI_VH    leaf Area Index of high vegetation
! WRMAX_VL  max volumetric water content retained on low vegetation
! WRMAX_VH  max volumetric water content retained on high vegetation
! ZTSL      reference height for temperature and humidity input
! ZUSL      reference height for wind input
! LAT       latitude
! FCOR      Coriolis factor
! Z0HA      AVERAGED Local roughness associated with exposed (no snow)
!           vegetation only (also no orography), for heat transfer
! VEGH     fraction of HIGH vegetation
! VEGL     fraction of LOW vegetation
! VGH_DENS   Density of trees in areas of high vegeation (m)
! VEGHEIGHT Height of trees in areas of high vegeation (m)
! ZVCAN   Wind speed within or above the canopy depending on CANO_REF_FORCING (m/s)
! FCANS        Canopy layer snowcover fractions
! SNCMA   Snow intercepted in high vegetation (kg m-2)
! RHOA   Air density (kg m-3)
! Z0SNOW  Roughness length for momentum of snow (m)
! Z0HSN Roughness length for heat of snow (m)
! PSNVH    fraction of HIGH vegetation covered by snow
!
!           - Output -
! HRSURF   relative humidity of the bare ground surface (1st layer)
! HUSURF    specific humidity of the bare ground surface
! HRSURFGV   relative humidity of the the snow-free ground below high veg (1st layer)
! HUSURFGV   specific humidity of the snow-free ground below high veg
! HV_VL        Halstead coefficient of low vegetation canopy
! HV_VH        Halstead coefficient of the high vegetation canopy
! HVSN_VH        Halstead coefficient of the high vegetation canopy accounting for intercepted snow
! DEL_VL    fraction of low veg. canopy covered by intercepted water
! DEL_VH    fraction of high veg canopy covered by intercepted water
! Z0HBG     Bare ground thermal roughness
! Z0HGV     Ground below high veg thermal roughness
! Z0HVL     LOW Vegetation thermal roughness
! Z0HVH     HIGH Vegetation thermal roughness
! ZRSURF     Aerodynamic surface resistance for soil under canopy (cf. Gouttevin et al. 2013)
! Z0MVH      local mom roughness length for high veg.
! Z0MVL      local mom roughness length for high veg.
! RES_SNCA  Resistance for sublimation of the intercepted snow in high canopy (

!
      INTEGER I, zopt

      real, dimension(n) :: temp, coef_vh, qsatgr, qsat_vl, qsat_vh, qsat_sv, &
           zqs_vl, zqs_vh, ctugr, ctugrv, ctuvh, ctuvl, wcrit_hrsurf, z0bg_n,ra,&
           z0gv_n, qsatgrv,wcrit_hrsurfgv, z0hg, zz0hgv, ZZ0HVH, ZZ0HVL, TSV_CORR
     real, dimension(n) :: ZUGV, ZTGV, ZZ0MGV, ZDH, QSATI_VH, VSUBL, Z0MVH_EFF, &
                           STABM_VH, ZUSTAR_VH, & ! Related to VGH
                           ZUSTAR_BELOW, ZRSURF, VMOD_BELOW, & ! Related to surface below VGH
                           LZZ0T_BELOW, STABT_BELOW, & ! Related to surface below VGH 
                           LZZ0M_BELOW, STABM_BELOW  ! Related to surface below VGH
     real :: ZFSURF
     REAL :: NU, MU, NR, DVAP
     REAL :: XI2,EXT2
!
!***********************************************************************
!
!
!
!
!
!------------------------------------------------------------------------
!
!
!
!         BARE GROUND LOCAL HEAT ROUGHNESS.  It is approximated by the
!         local momentum roughness of bare ground, times a scaling factor.

      DO I=1,N
         Z0BG_N(I) = Z0(I)
         Z0HG(I) = Z0M_TO_Z0H * Z0(I)
      END DO
!
!         GROUND BELOW HIGH VEG LOCAL HEAT ROUGHNESS.  It is approximated by the
!         local momentum roughness of bare ground, times a scaling factor.

      DO I=1,N
         Z0GV_N(I) = Z0(I) ! Value is the same as BG
         ZZ0HGV(I) = Z0M_TO_Z0H * Z0GV_N(I)
      END DO
!
!         HIGH AND LOW VEG LOCAL HEAT ROUGHNESS.  It is approximated by the
!         local momentum roughness of bare ground, times a scaling factor.

      DO I=1,N
!         Calculate an effective roughness length for the sparse canopy that accounts for VEG_DENS
!         as well as does a weighted average for the open canopy of the Z0 of snow and BGV (from Essery)    
         Z0MVH_EFF(I) =  VGH_DENS(I) * Z0MVH(I) + (1.-VGH_DENS(I))*(Z0SNOW(I)**PSNVH(I) * Z0GV_N(I)**(1.-PSNVH(I))) 
         ZZ0HVH(I) = Z0M_TO_Z0H * Z0MVH_EFF(I)
         ZZ0HVL(I) = Z0M_TO_Z0H * Z0MVL(I)
      END DO

!          Make sure that wind speed for the canopy is not equal to 0

      DO I=1,N
         ZVCAN(I) = MAX(ZVCAN(I), 0.1)
      END DO
!
!
!
!*       1.A     RELATIVE HYMIDITY OF THE EXPOSED BARE GROUND AND OF THE
!                GROUND BELOW HIGH VEG
!               -------------------------------------------------
!
!                        This relative humidity is related to
!                        the superficial soil moisture and the
!                        field capacity of the ground
!                        ** If the 1st soil layer is very shallow (under 5cm)
!                        might need to change the calc. to use a deeper layer
!
!                        Same value since only one soil profile is used
!
      if( svs_hrsurf_sltext ) then
         !use hrsurf formulation based on soil texture

         DO I=1,N
             ! Set soil water max for hrsurf calc. based on clay percentage)
            if ( clay1(i) .lt. 1.0 ) then
               wcrit_hrsurf(i) = wfc(i,1)
            else if ( clay1(i)  .lt. 40.0 ) then
               wcrit_hrsurf(i) =  (sand1(i)/((sand1(i)+clay1(i)))) * wfc(i,1) &
                    +             (clay1(i)/((sand1(i)+clay1(i)))) * wsat(i,1)
            else
               wcrit_hrsurf(i)= wsat(i,1)
            endif

            TEMP(I)   = PI*WD1(I)/WCRIT_HRSURF(I)
            HRSURF(I) = 0.5 * ( 1.-COS(TEMP(I)) )
            HRSURFGV(I) = HRSURF(I)
            WCRIT_HRSURFGV(I) = WCRIT_HRSURF(I)

         END DO
      else
         ! formulation based on field capacity
         DO I=1,N
            TEMP(I)   = PI*WD1(I)/WFC(I,1)
            HRSURF(I) = 0.5 * ( 1.-COS(TEMP(I)) )
            wcrit_hrsurf(i) = wfc(i,1)
            HRSURFGV(I) = HRSURF(I)
            WCRIT_HRSURFGV(I) = WCRIT_HRSURF(I)
         END DO
      endif
!
!
!*       1.B    SPECIFIC HUMIDITY OF THE EXPOSED BARE GROUND (HU)
!               -------------------------------------------------
!                         there is a specific treatment for dew
!                         (see Mahfouf and Noilhan, jam, 1991)
!
!                         first calculate the saturation vapor
!                         pressure and specific humidity
!
      DO I=1,N
        QSATGR(I) = FOQST( TGRS(I), PS(I) )
      END DO
!
!
      DO I=1,N
!
!
!                         when hu*qsat < qa, there are two
!                         possibilities
!
!                         low-level air is dry, i.e.,
!                         qa < qsat
!
!

        IF ( HRSURF(I)*QSATGR(I).LT.HU(I).AND.QSATGR(I).GT.HU(I) )&
                HRSURF(I) = HU(I) / QSATGR(I)

!
!
!                          b) low-level air is humid, i.e.,
!                          qa >= qsat
!
        IF ( HRSURF(I)*QSATGR(I).LT.HU(I).AND.QSATGR(I).LE.HU(I) )&
                  HRSURF(I) = 1.0

!
!                          for very humid soil (i.e., wg > wfc ),
!                          we take hu=1
!
        IF ( WD1(I).GT.WCRIT_HRSURF(I) ) HRSURF(I) = 1.0

!
      END DO
!
!                           Calculate specific humidity over ground
      DO I=1,N
        HUSURF(I) = HRSURF(I) * QSATGR(I)
      END DO
!
!
!
!*       1.C   SPECIFIC HUMIDITY OF THE GROUND BELOW HIGH VEGETATION (HU)
!               -------------------------------------------------
!
!                         there is a specific treatment for dew
!                         (see Mahfouf and Noilhan, jam, 1991)
!
          DO I=1,N
!
              IF(VEGH(I) > EPSILON_SVS) THEN ! High vegetation is present

!                         first calculate the saturation vapor
!                         pressure and specific humidity
                          QSATGRV(I) = FOQST( TGRVS(I), PS(I) )

!                         when hu*qsat < qa, there are two
!                         possibilities
!
!                         low-level air is dry, i.e.,
!                         qa < qsat
!
                  IF ( HRSURFGV(I)*QSATGRV(I).LT.HU(I).AND.QSATGRV(I).GT.HU(I) )&
                       HRSURFGV(I) = HU(I) / QSATGRV(I)
!
!                          b) low-level air is humid, i.e.,
!                          qa >= qsat
!
                  IF ( HRSURFGV(I)*QSATGRV(I).LT.HU(I).AND.QSATGRV(I).LE.HU(I) )&
                      HRSURFGV(I) = 1.0

!
!                          for very humid soil (i.e., wg > wfc ),
!                          we take hu=1
!
                  IF ( WD1(I).GT.WCRIT_HRSURFGV(I) ) HRSURFGV(I) = 1.0

!
!
!                           Calculate specific humidity over ground
                  HUSURFGV(I) = HRSURFGV(I) * QSATGRV(I)
!
            ELSE      ! No high vegetation
!
!               ! Use value from bare ground to avoid empty arrays
!
             HUSURFGV(I) = HUSURF(I)
!
            ENDIF
!
          END DO

!
!
!
!
!
!**     2.A     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR LOW VEGETATION
!*             ------------------------------------------------------------
!
!                         first calculate the saturation vapor
!                         pressure over low vegetation
!
!
      DO I=1,N
        QSAT_VL(I) = FOQST( TVGLS(I), PS(I) )
      END DO
!
!
!*                         then calculate the fraction of the foliage
!                          covered by intercepted water (DEL)
!
!
      DO I=1,N
!
!                          Calculate the maximum value of
!                          equivalent water content in the
!                          vegetation canopy
!
!                          calculate DEL
!
         IF(VEGL(I)>EPSILON_SVS) THEN   ! Low vegetation present in the grid cell
!
             DEL_VL(I) =   (WR_VL(I)/WRMAX_VL(I))**(2./3.)
!
             DEL_VL(I) = MIN(DEL_VL(I),0.25) ! from Boone et al. (2017), same as DEL_VH

         ELSE
             DEL_VL(I) = 0.
         ENDIF
!
      END DO
!
!
!                         calculate Hv based on previous time
!                         step resavg to use in flxsurf4
!
      DO I=1,N
!
!                         calculate Hv based on previous time
!                         step resavg to calculate specific
!                         humidity of low vegetation
!
         IF(VEGL(I)>EPSILON_SVS) THEN   ! Low vegetation present in the grid cell
             HV_VL(I) = 1. - MAX(0.,SIGN(1.,QSAT_VL(I)-HU(I)))&
                   *RS(I)*(1.-DEL_VL(I)) / (RESA_VL(I)+RS(I))

!      Atmospheric resistence for exchange between the the foliage
!      and the air within the canopy space (Dearrorff, 1978)
!        RA(I)  = 100.0 / (0.3*VMOD(I) + 0.3)  !
!       Equivalent of Hv defined with respect to the mean flow inside the canopy space
!        RPP(I) = 1. - MAX(0.,SIGN(1.,QSATVG(I)-QAF(I)))&
!                 *RS(I)*(1.-DEL(I)) / (RA(I)+RS(I))

!
!                         calculate specific humidity of vegetation
!
             ZQS_VL(I) = HV_VL(I) * QSAT_VL(I) + ( 1. - HV_VL(I) ) * HU(I)
!
         ELSE
             ZQS_VL(I) = 0.
             HV_VL(I) = 0.
         ENDIF
      END DO
!
!
!
!
      if ( svs_dynamic_z0h ) then
         zopt=9

         ! TO_DO NL: which values of Z0LOC, Z0HA to use here

         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGLS, ZQS_VL, Z0LOC, Z0HA, LAT, FCOR, z0mloc=z0loc, &
              optz0=zopt, L_min=sl_Lmin_soil, &
              coeft=CTUVL, z0t_optz0=Z0HVL )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
      else
        ! NL: Updated with roughn. lengths for VH instead of averaged roughn. length for veg
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGLS, ZQS_VL, Z0MVL, ZZ0HVL, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUVL )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif

         do i=1,n ! TO_DO NL: delete or update that with roughn. length for heat for VL?
            z0hvl(i)=zz0hvl(i)
         enddo
      endif

      DO I=1,N
         IF(VEGL(I)>EPSILON_SVS) THEN   ! Low vegetation present in the grid cell
             RESA_VL(I) = 1. / CTUVL(I)
         ELSE
             RESA_VL(I) = 1.
         ENDIF
      END DO
!
!
!
!
!**     2.B     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR HIGH VEGETATION
!*             ------------------------------------------------------------
!
!                         first calculate the saturation vapor
!                         pressure over high vegetation
!
!
      DO I=1,N
        QSAT_VH(I) = FOQST( TVGHS(I), PS(I) )
      END DO
      QSATI_VH(:) = QSATI( TVGHS(:), PS(:) )
!
!
!*                         then calculate the fraction of the foliage
!                          covered by intercepted water (DEL_VH)    
!
!
      DO I=1,N
!
!                          Calculate the maximum value of
!                          equivalent water content in the
!                          vegetation canopy
!
!
!                          calculate DEL_VH
!
         IF(VEGH(I) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell
             COEF_VH(I) = 1. + 2.*LAI_VH(I)
!
             IF (HU(I) .GT. QSAT_VH(I)) THEN ! Condensation
                DEL_VH(I) = 1.
             ELSE
             	!DEL_VH(I) =   MIN(WR_VH(I),WRMAX_VH(I)) / &
                !   ( (1.-COEF_VH(I))*MIN(WR_VH(I),WRMAX_VH(I)) + COEF_VH(I)*WRMAX_VH(I) )

!               From Deardoff (1978), also used in Gouttevin et al. 2015 and same equation as for FCANS
!               Avoids using COEFF_VH as we do not know where it comes from
                DEL_VH(I) = (WR_VH(I) / WRMAX_VH(I))**(2./3.)

             	DEL_VH(I) = MIN(DEL_VH(I),0.25) ! from Boone et al. (2017)
             ENDIF

         ELSE
             COEF_VH(I) = 1.
             DEL_VH(I) = 0.
         ENDIF
!
      END DO
!
!
!
!                         calculate HV_VH based on previous time
!                         step resa_vh to use in flxsurf4
!
      DO I=1,N
      !
!                         calculate HV_VH based on previous time
!                         step RESA_VH to calculate specific 
      !                         humidity of vegetation
      !
         IF(VEGH(I) .GE. EPSILON_SVS) THEN   ! High vegetation present in the grid cell
               HV_VH(I) = 1. - MAX(0.,SIGN(1.,QSAT_VH(I)-HU(I)))&
                  *RS(I)*(1.-DEL_VH(I)) / (RESA_VH(I)+RS(I))

      !      Atmospheric resistance for exchange between the the foliage
      !      and the air within the canopy space (Dearrorff, 1978)
      !        RA(I)  = 100.0 / (0.3*VMOD(I) + 0.3)  !
      !       Equivalent of Hv defined with respect to the mean flow inside the canopy space
      !        RPP(I) = 1. - MAX(0.,SIGN(1.,QSATVG(I)-QAF(I)))&
      !                 *RS(I)*(1.-DEL(I)) / (RA(I)+RS(I))

      !
      !                         calculate specific humidity of high vegetation
      !

             ZQS_VH(I) = (1.-FCANS(I)) * (HV_VH(I) * QSAT_VH(I) + ( 1. - HV_VH(I) ) * HU(I)) &
                         + FCANS(I) * QSATI_VH(I)

      !
         ELSE
               ZQS_VH(I) = 0.
               HV_VH(I) = 0.
         ENDIF
      END DO
   !
   !
!          Calculate aerodynamic resistance of high vegetation
   !
   !
      if ( svs_dynamic_z0h ) then
         zopt=9
         ! TO_DO NL: Which values of Z0LOC and Z0HA to use here
         ! Which height to use? Should be above the canopy
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGHS, ZQS_VH, Z0LOC, Z0HA, LAT, FCOR, z0mloc=z0loc, &
              optz0=zopt, L_min=sl_Lmin_soil, &
              coeft=CTUVH, z0t_optz0=Z0HVH )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
      else

         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGHS, ZQS_VH, Z0MVH_EFF, ZZ0HVH, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUVH, stabm=STABM_VH )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif

         do i=1,n 
            z0hvh(i)=zz0hvh(i)
           ! ustar above the canopy used in the aero resistances for turbulent fluxes
           ! stability momentum is used for ustar, stabm_gv = lzzo + am (L. 655 in sfclayer.F90)
           ZUSTAR_VH(I) = VMOD(I) * KARMAN / STABM_VH(I)
         enddo
      endif

      DO I=1,N
         IF(VEGH(I) .GT. EPSILON_SVS) THEN   ! High vegetation present in the grid cell
           RESA_VH(I) = 1. / CTUVH(I)

!          TODO: should VEG_DENS impact the aero resistance for a sparse canopy?
!          IF (VGH_DENS(I) .GT. EPSILON_SVS) THEN
!             ! From Mazzotti et al. (2020). Decreases the resistance (higher flux) as VGH_DENS decreases (more open)
!             RESA_VH(I) = RESA_VH(I) / VGH_DENS(I)**0.5
!          ENDIF
         ELSE
             RESA_VH(I) = 1.
         ENDIF
      END DO
!
!
!
!
!**     2.C     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR BARE GROUND
!*             ---------------------------------------------------------------
!
!                      *************************************
!                      DANS LE CALL DE FLXSURF, LE Z0H DEVRAIT ETRE
!                      REMPLACER PAR UN Z0H_LOCAL JUSTE POUR LE SOL NU
!                      *************************************
!
      if ( svs_dynamic_z0h ) then
         zopt=9
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TGRS, HUSURF, Z0LOC, Z0HG, LAT, FCOR, optz0=zopt ,&
              z0mloc=z0loc, L_min=sl_Lmin_soil, &
              coeft=CTUGR, z0t_optz0=Z0HBG )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error returned by sl_sfclayer()')
            return
         endif

      else
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TGRS, HUSURF, Z0BG_N, Z0HG, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUGR )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error returned by sl_sfclayer()')
            return
         endif
         
         do i=1,N
            z0hbg(i)=z0hg(i)
         enddo

      endif


      DO I=1,N
        RESAGR(I) = 1. / CTUGR(I)
      END DO
!
!
!
!
!
!
!
!
!**     2.D     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR GROUND
!                 BELOW HIGH VEG.
!
!*             ---------------------------------------------------------------
!
!                      *************************************
!                      DANS LE CALL DE FLXSURF, LE Z0H DEVRAIT ETRE
!                      REMPLACER PAR UN Z0H_LOCAL JUSTE POUR LE SOL NU
!                      *************************************
!

      if ( svs_dynamic_z0h ) then
         zopt=9
         ! TO_DO NL: check the roughness lengths used
         i = sl_sfclayer( THETAA, HU, ZVCAN, VDIR, ZUSL, ZTSL, &
              TGRVS, HUSURFGV, Z0LOC, ZZ0HGV, LAT, FCOR, optz0=zopt ,&
              z0mloc=z0loc, L_min=sl_Lmin_soil, &
              coeft=CTUGRV, z0t_optz0=Z0HGV )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error returned by sl_sfclayer()')
            return
         endif

      else

         IF (CANO_REF_FORCING == 'ABV') THEN ! Reference height above the canopy. In this case, z0 should be the canopy roughness lengths and the heights above canopy

             DO I=1,N
                 ! Determine the surface resistance using the momentum stability calculated above
                 ZFSURF = 1. + ZRALAI * (1. - EXP(- LAIVH(I) * VGH_DENS(I)))

                 ! The heat roughness length should be the one at the surface below canopy
                 ZRSURF(I) = LOG(Z0MVH_EFF(I) / ZZ0HGV(I)) / (ZUSTAR_VH(I) * KARMAN) * ZFSURF 

                 ! Estimate ustar below canopy/above surface from ZRSURF and roughness lengths of surface (see SNOWPACK code)
                 ZUSTAR_BELOW(I) = 1./(ZRSURF(I)*KARMAN) * LOG(VEGHEIGHT(I)/ZZ0HGV(I))
                 
                 VMOD_BELOW(I) = ZUSTAR_BELOW(I)/KARMAN * LOG(VEGHEIGHT(I)/Z0GV_N(I))


             ENDDO


             ! Calculate the atm stability below the canopy
             ! Temperature and height of vegetation are used instead of atm temperature and zref
             i = sl_sfclayer( TVGHS, ZQS_VH, VMOD_BELOW, VDIR, VEGHEIGHT, VEGHEIGHT, &
                  TGRVS, HUSURFGV, Z0GV_N, ZZ0HGV, LAT, FCOR, &
                  L_min=sl_Lmin_soil, &
                  stabm=STABM_BELOW, lzz0m=LZZ0M_BELOW, &
                  stabt=STABT_BELOW, lzz0t=LZZ0T_BELOW )

             do i=1,N
                 ! Apply atm stability to RSURF 
                 ZUSTAR_BELOW(I) = VMOD_BELOW(I) * KARMAN  / (LOG(VEGHEIGHT(I)/Z0GV_N(I)) + STABM_BELOW(I)-LZZ0M_BELOW(I))
                 ZRSURF(I) = 1./(ZUSTAR_BELOW(I)*KARMAN) * (LOG(VEGHEIGHT(I)/ZZ0HGV(I)) + STABT_BELOW(I)-LZZ0T_BELOW(I)) 
                 CTUGRV(I) = 1. / (RESA_VH(I) + ZRSURF(I))

                 Z0HGV(i) = ZZ0HGV(i)
             enddo

         ELSE ! O2F or FOREST

             i = sl_sfclayer( THETAA, HU, ZVCAN, VDIR, ZUSL, ZTSL, &
                  TGRVS, HUSURFGV, Z0GV_N, ZZ0HGV, LAT, FCOR, &
                  L_min=sl_Lmin_soil, &
                  coeft=CTUGRV )

             if (i /= SL_OK) then
                call physeterror('drag_svs', 'error returned by sl_sfclayer()')
                return
             endif


             do i=1,N
                Z0HGV(i) = ZZ0HGV(i)
             enddo
         ENDIF

      endif


      DO I=1,N
        IF(VEGH(I) .GE. EPSILON_SVS) THEN
            RESAGRV(I) = 1. / CTUGRV(I)
         ELSE
            RESAGRV(I) = 1.
         ENDIF
      END DO
!
!
!
!
!
!
!
!
!
!
!**     2.E     SURFACE AERO RESISTANCE FOR TURBULENT FLUXES FOR SNOW
!                 BELOW HIGH VEG.
!
!*             ---------------------------------------------------------------
!
!                      *************************************
!                      DANS LE CALL DE FLXSURF, LE Z0H DEVRAIT ETRE
!                      REMPLACER PAR UN Z0H_LOCAL JUSTE POUR LE SOL NU
!                      *************************************
!


     DO I=1,N

         IF (TSV(I,1) .LT. EPSILON_SVS) THEN ! If there is no snow, TSV == 0 and it creates issues
            TSV_CORR(I) = 273.15
         ELSE
            TSV_CORR(I) = TSV(I,1)
         ENDIF
     
     ENDDO

     QSAT_SV(:) = QSATI( TSV_CORR(:), PS(:) )

     IF (CANO_REF_FORCING == 'ABV') THEN ! Reference height above the canopy. 

         DO I=1,N

             ! Determine the surface resistance using the momentum stability calculated above
             ZFSURF = 1. + ZRALAI * (1. - EXP(- LAIVH(I) * VGH_DENS(I)))

             ! The heat roughness length should be the one at the surface below canopy
             ZRSURF(I) = LOG(Z0MVH_EFF(I) / Z0HSN(I)) / (ZUSTAR_VH(I) * KARMAN) * ZFSURF

             ! Estimate ustar below canopy/above surface from ZRSURF and roughness lengths of surface (see SNOWPACK code)
             ZUSTAR_BELOW(I) = 1./(ZRSURF(I)*KARMAN) * LOG(VEGHEIGHT(I)/Z0HSN(I))
             
             VMOD_BELOW(I) = ZUSTAR_BELOW(I)/KARMAN * LOG(VEGHEIGHT(I)/Z0SNOW(I))

         ENDDO

         ! Calculate the atm stability below the canopy
         ! Temperature and height of vegetation are used instead of atm temperature and zref
         i = sl_sfclayer( TVGHS, ZQS_VH, VMOD_BELOW, VDIR, VEGHEIGHT, VEGHEIGHT, &
              TSV_CORR, QSAT_SV, Z0SNOW, Z0HSN, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              stabm=STABM_BELOW, lzz0m=LZZ0M_BELOW, & 
              stabt=STABT_BELOW, lzz0t=LZZ0T_BELOW )

         do i=1,N
            IF(VEGH(I) .GE. EPSILON_SVS) THEN

                 ! Apply atm stability to RSURF
                 ZUSTAR_BELOW(I) = VMOD_BELOW(I) * KARMAN  / (LOG(VEGHEIGHT(I)/Z0SNOW(I)) + STABM_BELOW(I)-LZZ0M_BELOW(I))
                 ZRSURF(I) = 1./(ZUSTAR_BELOW(I)*KARMAN) * (LOG(VEGHEIGHT(I)/Z0HSN(I)) + STABT_BELOW(I)-LZZ0T_BELOW(I)) 

                 ! TODO: VV and NL, should the VEG_DEN impact the aero resistance for a sparse canopy
                 RESA_SV(I) = RESA_VH(I) + ZRSURF(I)

            ELSE
                 ! Set to 0 in that case. If it's 0, aero resistance for snow is calculated in surface_aero_cond
                 RESA_SV(I) = 0. 
            ENDIF
         ENDDO


     ELSE ! O2F or FOREST 

         DO i=1,N
            ! Set to 0 in that case. If it's 0, aero resistance for snow is calculated in surface_aero_cond
            RESA_SV(I) = 0.
         ENDDO

     ENDIF



!*       3.     Resistance of the snow intercepted in the high vegetation (RES_SNCA)
!               ---------------------------------------------------------------
!

      DO I=1,N

         IF (FCANS(I) .GT. 0. .AND. lres_snca) THEN

            ! Fraction of the entire forest height [-]
            XI2 = 1.-ZVENT

            ! Canopy wind speed extinction coefficient [-]
            ! Ellis et al (2010) (EL10) refers to Eagleson (2002) to justify the formulation of this coefficient
            EXT2 = ZBETA * LAI_VH(I)

            ! Computation of ventilation wind speed of intercepted snow derived from above-caopny wind speed  [Eq 8 in EL10]
            ! Estimated within canopy wind speed at fraction XI2 of the entire tree height [Eq 8 in EL10]  [m s-1]
            VSUBL(I) = ZVCAN(I) * EXP(-1. * EXT2 * XI2)

            ! Sutherland's equation for kinematic viscosity
            MU = 1.8325e-5*416.16/( THETAA(I)+120.)*(THETAA(I)/296.16)*SQRT(THETAA(I)/296.16)/RHOA(I)

            ! Compute Reynolds Number     [-]
            NR = 2.0 * RADIUS_ICESPH * ZVCAN(I) / MU

            ! Compute the Nusselt Number  [-]
            NU = 1.79 + 0.606 * SQRT(NR)

            ! Compute diffusivity of water vapour in air [m2 s-1]
            DVAP = 2.063e-5 * (TVGHS(I)/273.15)**(1.75)

            ! Resistance for snow within the canopy
            RES_SNCA(I) = 2.*917.*RADIUS_ICESPH**2 /(3.*0.02*FCANS(I)**(-0.4/0.67)*SNCMA(I)*NU*DVAP)

         ELSE
            RES_SNCA(I) = 0.
         ENDIF


      END DO
!


!*       4.     HALSTEAD COEFFICIENT (RELATIVE HUMIDITY OF THE VEGETATION) (HV)
!               ---------------------------------------------------------------
!
!
!                          Here we calculate "true" new HV based on updated
!                          resavg
!
      DO I=1,N


        HV_VL(I) = 1. - MAX(0.,SIGN(1.,QSAT_VL(I)-HU(I)))&
                 *RS(I)*(1.-DEL_VL(I)) / (RESA_VL(I)+RS(I))

        HV_VH(I) = 1. - MAX(0.,SIGN(1.,QSAT_VH(I)-HU(I)))&
                 *RS(I)*(1.-DEL_VH(I)) / (RESA_VH(I)+RS(I))

        ! Account for intercepted snow in the high vegetation
        HVSN_VH(I) = FCANS(I) * (CHLF + CHLC) *  RESA_VH(I) / (RESA_VH(I) + RES_SNCA(I)) &
                    + (1.-FCANS(I)) * HV_VH(I) * CHLC


      END DO
!


      RETURN
      END
