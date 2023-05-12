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

      SUBROUTINE SNOW_INTERCEPTION_SVS2 (T, HU,PS,VMOD,     & 
                     RR, SR, SNCMA, ESUBSNWC,SUBSNWC_CUM, LAIVH, VEGH,    &
                     RR_VEG, SR_VEG, DT,  N)


      use tdpack
      use sfclayer_mod, only: sl_sfclayer,SL_OK
      use sfc_options
      use svs_configs
      USE MODE_THERMOS 
      USE MODD_CSTS
      implicit none

      INTEGER N


      REAL SNCMA(N), RR(N), SR(N),  RR_VEG(N), SR_VEG(N)
      REAL PS(N), RHOA(N), HU(N)
      REAL LAIVH(N), T(N), VEGH(N),VMOD(N)
      REAL ESUBSNWC(N), SUBSNWC_CUM(N)
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
! T           air temperature at the forcing level (above high-vegetation) [K]
! HU          Specific humidity of air at the model lowest level [kg kg-1]
! VMOD        wind speed at the forcing level (above high-vegetation) [m/s]
! PS          Surface pressure [Pa]

!             - Input (other parameters)
! DT          time step [s]
! LAIVH       Vegetation leaf area index for HIGH vegetation only [m2/m2]
! VEGH        fraction of HIGH vegetation [0-1]


!             - Output 
!
! ESUBSNWC      latent heat of sublimation from incercepted snow [W m-2]
! SUBSNWC_CUM   Cumulated mass loss due to sublimation of incercepted snow [kg m-2]
! RR_VEG       liquid precipitation rate below high-vegetation [kg/m2/s]
! SR_VEG       solid  precipitation rate below high-vegetation [kg/m2/s]      

      INTEGER I

      REAL TUNL  ! Unloading time constant
      REAL, DIMENSION(N) :: SCAP 
      REAL, DIMENSION(N) :: DIRECT_SNOW ! Snow mass directly transferred through the canopy (kg/m^2)
      REAL, DIMENSION(N) :: DRIP_CPY ! Mass of liquid water dripping below the canopy (kg/m^2)
      REAL, DIMENSION(N) :: FVEG   !  Canopy cover fraction
      REAL, DIMENSION(N) :: INTCPT !  Canopy interception (kg/m^2)
      REAL, DIMENSION(N) :: NET_SNOW ! Snow mass sent to the snowpack below the canopy  (kg/m^2)
      REAL, DIMENSION(N) :: PQSAT ! Specific humidity at saturation [kg kg-1]
      REAL   SUB_CPY ! Intercepted snow mass lost by sublimation (kg/m^2)
      REAL  SNW_UNLOAD !  Unloaded mass in solid form (kg/m^2)      
      REAL, DIMENSION(N) :: UNLOAD !  Unloaded mass (kg/m^2)
      REAL, DIMENSION(N) :: VSUBL ! Ventilation wind speed used in the calcultation of sublimation

     
     ! Parameters used in the snow interception code
     !
      REAL, PARAMETER ::  CVAI = 4.4 ! Intercepted snow capacity per unit lai (kg/m^2) from Essery et a. (2003)
                  ! Note that HP98 are using a temperature-dependant param for falling snow density that could be considered later
                  ! The maximum canopy snow interception load in HP98 also depends on the type of trees. 

     ! Parameters used in the snow sublimation code
     !
      REAL, PARAMETER ::   CSUB =  0.002/3600.  ! Coefficient of transfert in Lundquist et al (2021) (kg s m**-3 Pa**-1)
      
      REAL   ::   EPS  ! Ratio of molecular weights of water and dry air

     ! Parameters used in the snow unloading code
     !
      REAL ::  TCNC ! Canopy unloading time scale for cold snow (s)
      REAL ::  TCNM  ! Canopy unloading time scale for melting snow (s)

     !
     ! 0. Initialize parameters
     !

      EPS =  XRD/XRV
      TCNC = 240*3600.
      TCNM = 48*3600.
     !
     ! 1. Configure the options used in the code
     !

     
     ! 2. Evolution of snow mass intercepted by the canopy
     !

      ! Initialize snowfall mass directly transfered through the canopy
      DIRECT_SNOW(:) = 0.

      ! Compute specific humidity at saturation
      PQSAT(:) = QSATI( T(:), PS(:))


      DO I=1,N
         IF(SNCMA(I)> 0. .OR. SR(I)>0.) THEN  ! Snow is present on the canopy or occurrence of snowfall

                !write(*,*) 'In inter', SNCMA(1), SR(1)*DT
                !!!!!!!!
                ! Interception Code (Simple version of HP98)
                !!!!!!!!

                ! Compute canopy snow interception capacity (kg m-2) (Eq 12 in HP98)
                ! CVAI may depends on type of trees and density of falling snow (see comment above) 
                SCAP(I) =  CVAI * LAIVH(I)

                ! Compute horizontal Canopy coverage 
                ! TODO: To be replaced by input from external database
                 FVEG(I) = 0.29 * LOG(LAIVH(I))+ 0.55
                 FVEG(I) = MIN(1., MAX(0., FVEG(I)))   

                !write(*,*) 'Veg prop',  LAIVH(I), FVEG(I)          

                ! Compute Canopy interception
                INTCPT(I) = (SCAP(I) - SNCMA(I))*(1 -exp(-FVEG(I)*SR(I)*DT/SCAP(I)))

                ! Update intercepted snow mass 
                SNCMA(I) = SNCMA(I) + INTCPT(I)

                ! Update the amount of snow that falls below high vegetation
                DIRECT_SNOW(I)  = SR(I)*DT - INTCPT(I)

                !write(*,*) 'After inter', SNCMA(1), INTCPT(1)

                !!!!!!!!
                ! Sublimation Code  (Lundquist et al, 2021)
                !!!!!!!!     

                VSUBL(I) = 0.7*VMOD(I) ! Initial and gross assumption. 
                                   ! Need to be revised! 
          

                IF( SNCMA(I)>EPSILON_SVS .AND. T(I) < 273.15) THEN
                   ! Mass of intercepted snow loss due to sublimation
                   ! Eq. 5 in Lundquist et al. (2021)
                   ! Only account for sublimation so far
                   SUB_CPY  = MAX(0.,CSUB * VSUBL(I) * PS(I) / EPS *(PQSAT(I)-HU(I)) * DT)
                  ! write(*,*) 'Subl',PQSAT(I), HU(I),SUB_CPY
                ENDIF

                ! Remove mass of intercepted lost by sublimation
                IF(SUB_CPY>SNCMA(I)) THEN
                   SUB_CPY = SNCMA(I)
                   SNCMA(I) = 0.
                ELSE
                   SNCMA(I) = SNCMA(I) - SUB_CPY
                ENDIF

                SUBSNWC_CUM(I) = SUBSNWC_CUM(I) +  SUB_CPY
                ESUBSNWC(I) = SUB_CPY/DT 


                !write(*,*) 'After Subl', SNCMA(1), SUB_CPY

                !!!!!!!!
                ! Unloading Code 
                !!!!!!!!

                ! Unloading time constant
                IF(T(I) > 273.15) THEN ! Melting conditions
                   TUNL = TCNM
                ELSE
                   TUNL = TCNC
                ENDIF

                ! Compute unloaded snow mass 
                UNLOAD(I) = SNCMA(I)  * DT / TUNL

                ! Update intercepted snow mass
                SNCMA(I)  = SNCMA(I) - UNLOAD(I)

               !!!!!!!!
               ! Update snowfall and rainfall rate sent to snow below high-veh 
               !!!!!!!! 
               IF(T(I)>=273.15) THEN 
                    DRIP_CPY(I) =  UNLOAD(I) ! Unload considered as liquid
                    SNW_UNLOAD = 0. 
                ELSE
                    DRIP_CPY(I) = 0.    ! Unload considered as solid
                    SNW_UNLOAD = UNLOAD(I)
               ENDIF   

                ! Update the total snow mass sent to the snowpack below
                ! high veg
                NET_SNOW(I) = DIRECT_SNOW(I) + UNLOAD(I)

                !write(*,*) 'After Unload', SNCMA(1), UNLOAD(1),T(1)
                !write(*,*) 'Drip', DRIP_CPY(1),'Solid', SNW_UNLOAD

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
