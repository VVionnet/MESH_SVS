!-------------------------------------- LICENCE BEGIN ------------------------
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

      SUBROUTINE SNOW_INTERCEPTION_SVS2 (DT, TVEG, WIND_TOP, ISWR,      &
                     RR, SR, WVEG, WVEG_MAX, SNCMA, ESUBSNWC,SUBSNWC_CUM, LAIVH, VEGH, HM_CAN,   &
                     VGH_DENS, SCAP, RR_VEG, SR_VEG,  FCANS, N)


      use tdpack
      use sfclayer_mod, only: sl_sfclayer,SL_OK
      use sfc_options
      use svs_configs
      USE MODE_THERMOS
      USE MODD_CSTS
      implicit none

      INTEGER N


      REAL SNCMA(N), RR(N), SR(N),  RR_VEG(N), SR_VEG(N), ISWR(N)
      REAL LAIVH(N), TVEG(N), VEGH(N),WIND_TOP(N),VGH_DENS(N)
      REAL ESUBSNWC(N), SUBSNWC_CUM(N), SCAP(N), FCANS(N), HM_CAN(N)
      REAL WVEG(N), WVEG_MAX(N), TAIR(N)
      REAL DT

!
!
!Author
!             V. Vionnet (April 2023)
!
!Object
!             Module to simulate interception of snow by high vegetation and
!             sublimation and unloading of intercepted snow using the
!             approach proposed in the FSM2 model
!
!Reference    Hedstrom and Pomeroy (1998) for the snow interception module (referred below as HP98)
!
!Arguments

!             - Input/Output (Prognostic variables)
!
! SNCMA       mass of intercepted snow in the canopy [kg/m2]

!             - Input (Forcing) -

! RR          liquid precipitation rate at the top of the high-vegetation [kg/m2/s]
! SR          solid  precipitation rate at the top of the high-vegetation [kg/m2/s]
! TVEG           Canopy temperature within the high vegetation [K]
! WIND_TOP    Wind speed at canopy top [m/s]
! ISWR        Solar radiation incident at the top of the high-vegetation [W m-2]
! ESUBSNWC      sublimation rate from incercepted snow [kg m-2 s-1]
!
!             - Input (other parameters)
! DT          time step [s]
! LAIVH       Vegetation leaf area index for HIGH vegetation only [m2/m2]
! VEGH        fraction of HIGH vegetation [0-1]
! VGH_DENS    tree cover density in  HIGH vegetation [0-1]
! SCAP       Vegetation layer snow capacities (kg m-2)
! HM_CAN    Heat mass for the high vegetation layer (J K-1 m-2)
!
!             - Output
!
! SUBSNWC_CUM   Cumulated mass loss due to sublimation of incercepted snow [kg m-2]
! RR_VEG       liquid precipitation rate below high-vegetation [kg/m2/s]
! SR_VEG       solid  precipitation rate below high-vegetation [kg/m2/s]
! FCANS        Canopy layer snowcover fractions from FSM2

      INTEGER I

      REAL TUNL  ! Unloading time constant

      REAL, DIMENSION(N) :: DIRECT_SNOW ! Snow mass directly transferred through the canopy (kg/m^2)
      REAL, DIMENSION(N) :: DRIP_CPY ! Mass of liquid water dripping below the canopy (kg/m^2)
      REAL, DIMENSION(N) :: INTCPT !  Canopy interception (kg/m^2)
      REAL, DIMENSION(N) :: NET_SNOW ! Snow mass sent to the snowpack below the canopy  (kg/m^2)



      REAL  SNW_UNLOAD !  Unloaded mass in solid form (kg/m^2)
      REAL, DIMENSION(N) :: UNLOAD !  Unloaded mass (kg/m^2)
      REAL, DIMENSION(N) :: HM_CAN_INI !  Unloaded mass (kg/m^2)
      REAL MELT


     ! Parameters used in the snow unloading code
     !
      REAL, PARAMETER  ::  TCNC = 240*3600. ! Canopy unloading time scale for cold snow (s)
      REAL, PARAMETER ::  TCNM = 48*3600. ! Canopy unloading time scale for melting snow (s)

     !
     ! 0. Initialise variables
     !
     ! Heat mass of the canopy without the snow
     HM_CAN_INI(:) = HM_CAN(:)

     !
     ! 1. Evolution of snow mass intercepted by the canopy
     !

      ! Initialize snowfall mass directly transfered through the canopy
      DIRECT_SNOW(:) = 0.
      NET_SNOW(:) = 0.

      DRIP_CPY(:) = 0.


      DO I=1,N
        IF (VEGH(I).GE.EPSILON_SVS) THEN
            IF(SNCMA(I)> 0. .OR. SR(I)>0.) THEN  ! Snow is present on the canopy or occurrence of snowfall


                    !!!!!!!!
                    ! Interception Code (Simple version of HP98)
                    !!!!!!!!

                    ! Compute Canopy interception, check Pomeroy VGH_DENS position
                    INTCPT(I) = (SCAP(I) - SNCMA(I))*(1. -exp(-VGH_DENS(I) * SR(I)*DT/SCAP(I)))

                    ! Update intercepted snow mass
                    SNCMA(I) = SNCMA(I) + INTCPT(I)

                    ! Update the amount of snow that falls below high vegetation
                    DIRECT_SNOW(I)  = SR(I)*DT - INTCPT(I)

                    ! Remove mass of intercepted lost by sublimation
                    ! Sublimation calculated in the previous time step so should be
                    ! applied to snow already on canopy

                    IF(ESUBSNWC(I)>0) THEN
                        IF (SNCMA(I)>0) THEN
                            SNCMA(I) = SNCMA(I) - ESUBSNWC(I)*DT
                        ENDIF
                    ELSE
                        IF (TVEG(I) < 273.15) THEN
                            SNCMA(I) = SNCMA(I) - ESUBSNWC(I)*DT
                        ENDIF
                    ENDIF

                    SNCMA(I) = MAX(0.,SNCMA(I))

                    SUBSNWC_CUM(I) = SUBSNWC_CUM(I) +  ESUBSNWC(I)*DT

                    !!!!!!!!
                    ! Unloading Code
                    !!!!!!!!

                    ! Unloading time constant
                    IF(TVEG(I) >= 273.15) THEN ! Melting conditions
                        TUNL = TCNM
                    ELSE
                        TUNL = TCNC
                    ENDIF

                    ! Compute unloaded snow mass
                    UNLOAD(I) = SNCMA(I)  * DT / TUNL

                    ! Update intercepted snow mass
                    SNCMA(I)  = SNCMA(I) - UNLOAD(I)

                    NET_SNOW(I) = DIRECT_SNOW(I) + UNLOAD(I)


                    !!!!!!!!
                    ! Melting of snow canopy
                    !!!!!!!!

                    HM_CAN(I) = HM_CAN_INI(I) + CPI * SNCMA(I)

                    IF(TVEG(I)>273.15 .AND. SNCMA(I) > 0) THEN

                        !MELT = CPI * SNCMA(I) * (TVEG(I) - 273.15) / CHLF
                        !MELT = HM_CAN(I) * (TVEG(I) - 273.15) / CHLF
                        MELT = 5.56E-6 * DT * (SNCMA(I)/SCAP(I)) * (TVEG(I) - 273.15)

                        IF (MELT  > SNCMA(I)) THEN
                            MELT = SNCMA(I)
                        ENDIF

                        DRIP_CPY(I) = MELT

                        !TVEG(I) = TVEG(I) - CHLF * MELT / HM_CAN(I)
                        !TVEG(I) = TVEG(I) - CHLF * MELT / (CPI * SNCMA(I))
                        TVEG(I) = TVEG(I) -   CHLF * MELT / HM_CAN(I)

                        SNCMA(I) = SNCMA(I) - MELT

                    ELSE
                        DRIP_CPY(I) = 0.    ! Unload considered as solid
                    ENDIF


                    ! Canopy layer snowcover fractions Boone
                    FCANS(I) = 0.89*(SNCMA(I)/SCAP(I))**0.3/(1.+EXP(-4.7*(SNCMA(I)/SCAP(I)-0.45)))

                    ! Update canopy heat mass with new intercepted snow
                    HM_CAN(I) = HM_CAN_INI(I) + CPI * SNCMA(I)

            ENDIF
         ENDIF
       END DO



      DO I=1,N

         IF (VEGH(I).GE.EPSILON_SVS) THEN

            ! Rain is not intercepted by vegetation when snow is present under high veg.
            ! Therefore rain under high veg accounts for total rain above high veh plus dripping
            RR_VEG(I) = RR(I)+DRIP_CPY(I)/DT

            ! Snowfall rate below high vegetation
            SR_VEG(I) = NET_SNOW(I)/DT

         ELSE
            ! no high veg, no modification of rainfall and snowfall

            RR_VEG(I) = RR(I)
            SR_VEG(I) = SR(I)

         ENDIF

      ENDDO

      END SUBROUTINE SNOW_INTERCEPTION_SVS2
