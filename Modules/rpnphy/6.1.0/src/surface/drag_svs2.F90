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

      SUBROUTINE DRAG_SVS2 ( TGRS, TGRVS, TVGLS, TVGHS, WD1, &
                              WR_VL, WR_VH,  THETAA, VMOD, VDIR, HU,  &
                              PS, RS, Z0, Z0LOC, Z0VG, WFC, WSAT, CLAY1,  &
                              SAND1, LAI_VL, LAI_VH, WRMAX_VL,WRMAX_VH,&
                              ZUSL, ZTSL, LAT, &
                              FCOR, Z0HA, VEGL, VEGH, CLUMPING,&
                              VGH_DENS, Z0MVH,Z0MVL, VEGHEIGHT,  &
                              LAIVH, ZVCAN, &
                              RESAGR,RESAGRV, RESA_VL, RESA_VH, &
                              HUSURF,HUSURFGV, &
                              HRSURF,HRSURFGV, &
                              HV_VL, HV_VH, DEL_VL, DEL_VH,  &
                              Z0HBG, Z0HVL, Z0HVH,Z0HGV,  N)
      use tdpack
      use sfc_options
      use sfclayer_mod, only: sl_sfclayer,SL_OK
      use svs_configs
!
      implicit none
!!!#include <arch_specific.hf>

      INTEGER N
      REAL TGRS(N), TVGLS(N), TVGHS(N), WR_VL(N), WR_VH(N), THETAA(N), VMOD(N), VDIR(N), HU(N), CLAY1(N)
      REAL SAND1(N), PS(N), RS(N), Z0(N),Z0LOC(N), Z0VG(N), WFC(N,NL_SVS), WSAT(N,NL_SVS)
      REAL LAI_VL(N), LAI_VH(N),WRMAX_VH(N),WRMAX_VL(N), ZUSL(N), ZTSL(N), LAT(N)
      REAL FCOR(N), Z0HA(N), Z0HBG(N), Z0HVL(N), Z0HVH(N), Z0HGV(N)
      REAL VGH_DENS(N), Z0MVH(N), Z0MVL(N), VEGHEIGHT(N), LAIVH(N), ZVCAN(N)
      REAL RESAGR(N),RESAGRV(N), RESA_VL(N), RESA_VH(N)
      REAL HUSURF(N),HUSURFGV(N), HV_VL(N), HV_VH(N), DEL_VL(N), DEL_VH(N)
      REAL HRSURF(N),HRSURFGV(N), WD1(N)
      REAL VEGH(N), VEGL(N), TGRVS(N)
      REAL CLUMPING
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
! CLUMPING   coefficient to convert LAI to effective LAI
! VGH_DENS   Density of trees in areas of high vegeation (m)
! VEGHEIGHT Height of trees in areas of high vegeation (m)
! ZVCAN   Wind speed within or above the canopy depending on CANO_REF_FORCING (m/s)
!
!           - Output -
! HRSURF   relative humidity of the bare ground surface (1st layer)
! HUSURF    specific humidity of the bare ground surface
! HRSURFGV   relative humidity of the the snow-free ground below high veg (1st layer)
! HUSURFGV   specific humidity of the snow-free ground below high veg
! HV_VL        Halstead coefficient (i.e., relative humidity of the
!           low vegetation canopy)
! HV_VH        Halstead coefficient (i.e., relative humidity of the
!           high vegetation canopy)
! DEL_VL    fraction of low veg. canopy covered by intercepted water
! DEL_VH    fraction of high veg canopy covered by intercepted water
! Z0HBG     Bare ground thermal roughness
! Z0HGV     Ground below high veg thermal roughness
! Z0HVL     LOW Vegetation thermal roughness
! Z0HVH     HIGH Vegetation thermal roughness
! ZRSURF     Aerodynamic surface resistance for soil under canopy (cf. Gouttevin et al. 2013)
! Z0MVH      local mom roughness length for high veg.
! Z0MVL      local mom roughness length for high veg.

!
      INTEGER I, zopt

      real, dimension(n) :: temp, coef_vl, coef_vh, qsatgr, qsat_vl, qsat_vh, &
           zqs_vl, zqs_vh, ctugr, ctugrv, ctuvg, wcrit_hrsurf, z0bg_n,ra,&
           z0gv_n, qsatgrv,wcrit_hrsurfgv, z0hg, zz0hgv, ZZ0HVH, ZZ0HVL
     real, dimension(n) :: ZUGV, ZTGV, ZZ0MGV, ZDH

     real :: ZRSURF, LZZ0, LZZ0T, RESAGRV_NEUTRAL, ZUSTAR, ZFSURF

     REAL, PARAMETER :: ZRALAI = 3.! Parameter for excess resistance introduced by canopy between surface and ref level (cf Table 1, Gouttevin et al. 2015)
     REAL, PARAMETER :: ZRCHD = 0.67    ! Ratio of displacement height to canopy height
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
         Z0GV_N(I) = 1.0 ! Value to be modified
         ZZ0HGV(I) = Z0M_TO_Z0H * Z0GV_N(I)
      END DO
!
!         HIGH AND LOW VEG LOCAL HEAT ROUGHNESS.  It is approximated by the
!         local momentum roughness of bare ground, times a scaling factor.

      DO I=1,N
         ZZ0HVH(I) = Z0M_TO_Z0H * Z0MVH(I)
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
!**     2.A     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR BARE GROUND
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
!**     2.B     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR GROUND
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
              TGRVS, HUSURF, Z0LOC, ZZ0HGV, LAT, FCOR, optz0=zopt ,&
              z0mloc=z0loc, L_min=sl_Lmin_soil, &
              coeft=CTUGRV, z0t_optz0=Z0HGV )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error returned by sl_sfclayer()')
            return
         endif

      else

         IF (CANO_REF_FORCING == 'ABV') THEN ! Reference height above the canopy. In this case, z0 should be the canopy roughness lengths and the heights above canopy
            DO I=1,N
                 ! WARNING NL: Might need to be updated following conversation with Stephane B. and Maria A.
                 ZUGV(I) = ZUSL(I) + VEGHEIGHT(I)
                 ZTGV(I)  = ZTSL(I) + VEGHEIGHT(I)
                 ZDH(I) = VEGHEIGHT(I)*ZRCHD
            ENDDO


             ! This CTUGRV should be for high vegetation
             i = sl_sfclayer( THETAA, HU, ZVCAN, VDIR, ZUGV, ZTGV, &
                  TGRVS, HUSURF, Z0MVH, ZZ0HVH, LAT, FCOR, &
                  L_min=sl_Lmin_soil, &
                  coeft=CTUGRV )

             if (i /= SL_OK) then
                call physeterror('drag_svs', 'error returned by sl_sfclayer()')
                return
             endif

             DO I=1,N

                 ! Compute aerodymanical resistance with stability atm correction
                 RESAGRV(I) = 1. / CTUGRV(I)

                 ! ustar above the canopy used in the aero resistances for turbulent fluxes
                 ZUSTAR = ZVCAN(I) * KARMAN / LOG((ZUGV(I)-ZDH(I))/Z0MVH(I))
                 ZFSURF = 1. + ZRALAI * (1. - EXP(-CLUMPING * LAIVH(I) * VGH_DENS(I)))
                 ZRSURF = LOG(Z0MVH(I) / ZZ0HGV(I)) / (ZUSTAR * KARMAN) * ZFSURF ! The heat roughness length should be the one at the surface below canopy

                 ! Compute aerodymanical resistance for neutral stability
                 ! cf sfclayer_mod (L.588-605 for computation of lzz0 and lzz0t (z0ref== T)
                 LZZ0 = LOG((ZUGV(I) + Z0MVH(I)) / Z0MVH(I))
                 LZZ0T = LOG((ZTGV(I) + ZZ0HVH(I)) / ZZ0HVH(I))
                 RESAGRV_NEUTRAL = 1. / (VMOD(I) * KARMAN * KARMAN / (LZZ0 * LZZ0T))

                ! Apply stability correction  to ZRSURF
                ZRSURF = ZRSURF * RESAGRV(I) / RESAGRV_NEUTRAL

                CTUGRV(I) = 1. / (RESAGRV(I) + ZRSURF)
             ENDDO

             do i=1,N
                Z0HGV(i) = ZZ0HGV(i)
             enddo
         ELSE

             i = sl_sfclayer( THETAA, HU, ZVCAN, VDIR, ZUSL, ZTSL, &
                  TGRVS, HUSURF, Z0GV_N, ZZ0HGV, LAT, FCOR, &
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
        IF(VEGH(I)>EPSILON_SVS) THEN
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

!**     2.C     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR LOW VEGETATION
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
             COEF_VL(I) = 1. + 2.*LAI_VL(I)
!

             DEL_VL(I) =   MIN(WR_VL(I),WRMAX_VL(I)) / &
              ((1.-COEF_VL(I))*MIN(WR_VL(I),WRMAX_VL(I)) +COEF_VL(I)*WRMAX_VL(I) )
!
             DEL_VL(I) = MIN(DEL_VL(I),0.1)

         ELSE
             COEF_VL(I) = 1.
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
              coeft=CTUVG, z0t_optz0=Z0HVL )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
      else
        ! NL: Updated with roughn. lengths for VH instead of averaged roughn. length for veg
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGLS, ZQS_VL, Z0MVL, ZZ0HVL, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUVG )

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
             RESA_VL(I) = 1. / CTUVG(I)
         ELSE
             RESA_VL(I) = 1.
         ENDIF
      END DO
!
!
!**     2.D     SURFACE TRANSFER COEFFICIENTS FOR HEAT (CH) FOR HIGH VEGETATION
!*             ------------------------------------------------------------
!
!                         first calculate the saturation vapor
!                         pressure over vegetation
!
!
      DO I=1,N
        QSAT_VH(I) = FOQST( TVGHS(I), PS(I) )
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
!
!                          calculate DEL
!
         IF(VEGH(I)>EPSILON_SVS) THEN   ! High vegetation present in the grid cell
             COEF_VH(I) = 1. + 2.*LAI_VH(I)
!
             DEL_VH(I) =   MIN(WR_VH(I),WRMAX_VH(I)) / &
                   ( (1.-COEF_VH(I))*MIN(WR_VH(I),WRMAX_VH(I)) + COEF_VH(I)*WRMAX_VH(I) )
!
             DEL_VH(I) = MIN(DEL_VH(I),0.1)

         ELSE
             COEF_VH(I) = 1.
             DEL_VH(I) = 0.
         ENDIF
!
      END DO
!
!
!
!

!
!
!                         calculate Hv based on previous time
!                         step resavg to use in flxsurf4
!
      DO I=1,N
      !
      !                         calculate Hv based on previous time
      !                         step resavg to calculate specific
      !                         humidity of vegetation
      !
         IF(VEGH(I)>EPSILON_SVS) THEN   ! High vegetation present in the grid cell
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
               ZQS_VH(I) = HV_VH(I) * QSAT_VH(I) + ( 1. - HV_VH(I) ) * HU(I)
      !
         ELSE
               ZQS_VH(I) = 0.
               HV_VH(I) = 0.
         ENDIF
      END DO
   !
   !
   !
   !
      if ( svs_dynamic_z0h ) then
         zopt=9
         ! TO_DO NL: Which values of Z0LOC and Z0HA to use here
         ! Which height to use? Should be above the canopy
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGHS, ZQS_VH, Z0LOC, Z0HA, LAT, FCOR, z0mloc=z0loc, &
              optz0=zopt, L_min=sl_Lmin_soil, &
              coeft=CTUVG, z0t_optz0=Z0HVH )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif
      else
         ! NL: Updated with roughn. lengths for VH instead of averaged roughn. length for veg
         i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
              TVGHS, ZQS_VH, Z0MVH, ZZ0HVH, LAT, FCOR, &
              L_min=sl_Lmin_soil, &
              coeft=CTUVG )

         if (i /= SL_OK) then
            call physeterror('drag_svs', 'error 2 returned by sl_sfclayer()')
            return
         endif

         do i=1,n ! TO_DO NL: delete or update that with roughn. length for heat for VL?
            z0hvh(i)=zz0hvh(i)
         enddo
      endif

      DO I=1,N
         IF(VEGH(I)>EPSILON_SVS) THEN   ! High vegetation present in the grid cell
             RESA_VH(I) = 1. / CTUVG(I)
         ELSE
             RESA_VH(I) = 1.
         ENDIF
      END DO


!*       3.     HALSTEAD COEFFICIENT (RELATIVE HUMIDITY OF THE VEGETATION) (HV)
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

      END DO
!
      RETURN
      END
