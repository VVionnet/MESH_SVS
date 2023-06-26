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

      SUBROUTINE SNOW_INTERCEPTION_SVS2 (T, HU,PS,WIND_TOP, ISWR,RHOA,      & 
                     RR, SR, SNCMA, ESUBSNWC,SUBSNWC_CUM, LAIVH, VEGH,    &
                     VGH_DENS, RR_VEG, SR_VEG, DT,  N)


      use tdpack
      use sfclayer_mod, only: sl_sfclayer,SL_OK
      use sfc_options
      use svs_configs
      USE MODE_THERMOS 
      USE MODD_CSTS
      implicit none

      INTEGER N


      REAL SNCMA(N), RR(N), SR(N),  RR_VEG(N), SR_VEG(N), ISWR(N)
      REAL PS(N), RHOA(N), HU(N)
      REAL LAIVH(N), T(N), VEGH(N),WIND_TOP(N),VGH_DENS(N)
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
! T           air temperature within the high vegetation [K]
! HU          Specific humidity of air within the high vegetation [kg kg-1]
! WIND_TOP    Wind speed at canopy top [m/s]
! PS          Surface pressure [Pa]
! ISWR        Solar radiation incident at the top of the high-vegetation [W m-2]

!             - Input (other parameters)
! DT          time step [s]
! LAIVH       Vegetation leaf area index for HIGH vegetation only [m2/m2]
! VEGH        fraction of HIGH vegetation [0-1]
! VGH_DENS    tree cover density in  HIGH vegetation [0-1]


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
      REAL, DIMENSION(N) :: INTCPT !  Canopy interception (kg/m^2)
      REAL, DIMENSION(N) :: NET_SNOW ! Snow mass sent to the snowpack below the canopy  (kg/m^2)
      REAL, DIMENSION(N) :: PQSAT ! Specific humidity at saturation [kg kg-1]
      REAL, DIMENSION(N) :: PPSAT ! Vapor presure at saturation [Pa]

      REAL, DIMENSION(N) :: SNCMAT ! Mass of intercepted snow at the endof the time step [kg/m2]


      REAL   SUB_CPY ! Intercepted snow mass lost by sublimation (kg/m^2)
      REAL  SNW_UNLOAD !  Unloaded mass in solid form (kg/m^2)      
      REAL, DIMENSION(N) :: UNLOAD !  Unloaded mass (kg/m^2)
      REAL, DIMENSION(N) :: VSUBL ! Ventilation wind speed used in the calcultation of sublimation


      CHARACTER(LEN=4) :: HDIFU  ! Option to compute the molecular diffusivity of water vapour in air
      CHARACTER(LEN=4) HSUBL_CANO ! Option used for the snow canopy sulimation

      REAL XI2,EXT2, WINDEXT2, UVENT, EVAP, LAMBDA, SSTAR, A1, B1,  LS, J
      REAL C1, SIGMA, SVDENS, SUB_RATE,SUB_POT, MPM, CE
      REAL NU, NR, MU,DVAP
     
     ! Parameters used in the snow interception code
     !
      REAL, PARAMETER ::  CVAI = 4.4 ! Intercepted snow capacity per unit lai (kg/m^2) from Essery et a. (2003)
                  ! Note that HP98 are using a temperature-dependant param for falling snow density that could be considered later
                  ! The maximum canopy snow interception load in HP98 also depends on the type of trees. 

     ! Parameters used in the snow sublimation code
     !
     ! LU21             
      REAL, PARAMETER ::   CSUB =  0.002/3600.  ! Coefficient of transfert in Lundquist et al (2021) (kg s m**-3 Pa**-1)

     ! HP98
      REAL, PARAMETER ::     ALPHA = 5   ! Shape parameter used when computing the averahe sublimation rate
      REAL, PARAMETER ::     GAMA = 1.15 ! Parameter used in the computation of the exponential wind profile in the canopy
      REAL, PARAMETER ::    ZVENT = 0.75 ! Ratio between ventilation wind speed height and tree height [-]
      REAL, PARAMETER :: RADIUS_ICESPH = 5e-4 ! Radius of single 'ideal' ice shpere [m]
      REAL, PARAMETER :: ALBEDO_ICESPH = 0.8  ! Albedo of single 'ideal' ice shpere [-] 
      REAL, PARAMETER :: KS = 0.0114          ! Snow shape coefficient for jack pine. Taken from Pomeroy et al. (1998)
      REAL, PARAMETER :: FRACT = 0.37 ! Fractal dimension of intercepted snow [-]


      REAL   EPS ! Ratio of molecular weights of water and dry air

     ! Parameters used in the snow unloading code
     !
      REAL, PARAMETER  ::  TCNC = 240*3600. ! Canopy unloading time scale for cold snow (s)
      !REAL, PARAMETER ::  TCNM = 48*3600. ! Canopy unloading time scale for melting snow (s)
      REAL, PARAMETER ::  TCNM = 6.*3600. ! Canopy unloading time scale for melting snow (s)


     ! Physical constant used in the code 
     ! TODO: Use the constants from GEM

      REAL, SAVE :: MWAT = 0.0180153 ! Molecular weight of water [kg mol-1]
      REAL, SAVE :: RGAZ = 8.314  !  Universal gas constant [J mol-1 K-1]
      REAL, SAVE :: CICE = 2.102e3; !  Heat capacity of ice [J kg-1 K-1]

     !
     ! 0. Initialize parameters
     !

      EPS =  XRD/XRV
     !
     ! 1. Configure the options used in the code
     !

      HDIFU = 'HP98' ! Selection the formulation used to compute the molecular  diffusivity  of  water  vapour  in  air
                      !    'HP98': formulation used in HP98
                      !    'HP13': formulation used in Harder and Pomeroy (2013)
                      !    'TM66': original formulation used in Thorpe and Mason (1966) including the contribution of air pressure
                      !    'PU98': formulation used in Pruppacher et al.(1998)

      HSUBL_CANO = 'NONE'  ! Select the approach used to compute the mass loss due sublimation of intercepted snow 
                      ! 'LU21': simple approach from Lundquist et al. ! (2021)
                      ! 'HP98': formulation proposed by Hedstrom and Pomeroy (1998) 
                      ! 'NONE': no mass loss due to sublimation of intercepted snow

     
     ! 2. Evolution of snow mass intercepted by the canopy
     !

      ! Initialize snowfall mass directly transfered through the canopy
      DIRECT_SNOW(:) = 0.

      ! Compute specific humidity and vapor pressure at saturation
      PQSAT(:) = QSATI( T(:), PS(:))
      PPSAT(:) = PSAT(T(:))


      DO I=1,N
         IF(SNCMA(I)> 0. .OR. SR(I)>0.) THEN  ! Snow is present on the canopy or occurrence of snowfall

                write(*,*) 'In inter', SNCMA(1), SR(1)*DT
                !!!!!!!!
                ! Interception Code (Simple version of HP98)
                !!!!!!!!

                ! Compute canopy snow interception capacity (kg m-2) (Eq 12 in HP98)
                ! CVAI may depends on type of trees and density of falling snow (see comment above) 
                SCAP(I) =  CVAI * LAIVH(I)

                ! Compute Canopy interception
                INTCPT(I) = (SCAP(I) - SNCMA(I))*(1 -exp(-VGH_DENS(I)*SR(I)*DT/SCAP(I)))

                ! Update intercepted snow mass 
                SNCMA(I) = SNCMA(I) + INTCPT(I)

                ! Update the amount of snow that falls below high vegetation
                DIRECT_SNOW(I)  = SR(I)*DT - INTCPT(I)

                write(*,*) 'After inter', SNCMA(1), INTCPT(1)

                !!!!!!!!
                ! Sublimation Code  (Lundquist et al, 2021)
                !!!!!!!!     

                IF(HSUBL_CANO == 'LU21') THEN
                    VSUBL(I) = 0.7*WIND_TOP(I) ! Initial and gross assumption. 
                                   ! Need to be revised! 
          

                   IF( SNCMA(I)>EPSILON_SVS .AND. T(I) < 273.15) THEN
                      ! Mass of intercepted snow loss due to sublimation
                      ! Eq. 5 in Lundquist et al. (2021)
                      ! Only account for sublimation so far
                      SUB_CPY  = MAX(0.,CSUB * VSUBL(I) * PS(I) / EPS *(PQSAT(I)-HU(I)) * DT)
                     ! write(*,*) 'Subl',PQSAT(I), HU(I),SUB_CPY
                   ENDIF

                ELSE IF(HSUBL_CANO == 'HP98') THEN

                   ! Fraction of the entire forest height [-]
                   XI2 = 1-ZVENT

                   ! Canopy wind speed extinction coefficient [-]
                   ! Ellis et al (2010) (EL10) refers to Eagleson (2002) to justify the formulation of this coefficient
                   EXT2 = GAMA * LAIVH(I) 

                   ! Computation of ventilation wind speed of intercepted snow derived from above-caopny wind speed  [Eq 8 in EL10]
                   ! Estimated within canopy wind speed at fraction XI2 of the entire tree height [Eq 8 in EL10]  [m s-1]
                   VSUBL(I) = WIND_TOP(I) * EXP(-1. * EXT2 * XI2)

                    ! Vaport density at saturation 
                    SVDENS = PPSAT(I) * MWAT / (RGAZ * T(I))

                   ! Derive partial vapor pressure from speficic humidity and pressure [Pa]
                   EVAP = HU(I) * PS(I) / (0.622+0.378 * HU(I) ) 

                   ! Sutherland's equation for kinematic viscosity
                   MU=1.8325e-5*416.16/( T(I)+120)*(T(I)/296.16)*SQRT(T(I)/296.16)/RHOA(I)

                   ! Compute thermal conductivity of air [J m-1 s-1 K-1]
                   LAMBDA = 0.000063 * T(I) + 0.00673

                   ! Compute latent heat of sublimation [J kg-1]
                   ! We use here the formulation of Roges and Yau (1989) as in Harder and Pomeroy (2013)
                   ! Note that the canopy module in CRHM used a constant value
                   ! TODO: Crocus approach should be used for consistency
                   IF(T(I)<273.15) THEN
                      LS =  1000.0 * (2834.1 - 0.29 *(T(I)-273.15) - 0.004*(T(I)-273.15)**2.)
                   ELSE 
                      LS = 1000.0 * (2501.0 - (2.361 * (T(I)-273.15)))
                   ENDIF    

                   ! Compute diffusivity of water vapour in air [m2 s-1]
                   IF(HDIFU=='HP98') THEN
                      DVAP = 2.063e-5 * (T(I)/273.15)**(-1.75)
                   ELSE IF(HDIFU=='HP13') THEN
                      DVAP = 2.063e-5 * (T(I)/273.15)**1.75
                   ELSE IF(HDIFU=='TM66') THEN
                     DVAP = 2.063e-5 *   (T(I)/273.15)**1.75 * (PS(I)/101325)
                    ELSE IF(HDIFU=='PU98') THEN
                      DVAP =  0.211E-4 * (T(I)/273.15)**1.94 * (101325/PS(I))
                   ENDIF

                   ! Compute Reynolds Number     [-]
                   NR  = 2.0 * RADIUS_ICESPH * UVENT / MU

                   ! Compute the Nusselt Number  [-]
                   NU = 1.79 + 0.606 * SQRT(NR)

                   ! Incoming shortwave radiation to the ideal snow particle [W m-2]
                   SSTAR = PI* RADIUS_ICESPH**2. * (1.0 - ALBEDO_ICESPH) * ISWR(I)

                   ! Compute the term used in equation for the sublimation rate of an ice sphere
                   ! from Thorpe and Masson (1966)
                   A1  = LAMBDA * T(I) *NU
                   B1 = ( LS * MWAT /(RGAZ * T(I)) ) -1.0 
                   J = B1/A1
                   C1  = 1.0 / (DVAP * NU * SVDENS )

                   ! Compute water vapour deficit with respect to  ice  [-]
                   SIGMA = EVAP/PPSAT(I) - 1.0  

                   !write(*,*) 'LS, SVDENS',LS, SVDENS

                   ! Compute the mean mass of a snow particle assuming that the size distribution of 
                   ! intercepted snow particle follow a gamma distribution with alpha  = 5
                   ! and a mean radius equals to RADIUS_ICESPH
                   ! This approach is used in PBSM3D (see EQ 23 in Pomeroy et al (1993))
                   ! However, it is never described in the papers about the canopy module in CRHM (J. Pomeroy, personal
                   ! communication, 2022)
                   MPM = 4.0 / 3.0 * PI * 917. * RADIUS_ICESPH**3. * (1.0 + 3.0 / ALPHA + 2.0 / ALPHA**2.)

                   ! Sublimation rate coefficient of single 'ideal' ice sphere [s-1]
                  SUB_RATE = (2.0 *PI * RADIUS_ICESPH *SIGMA - SSTAR *J)/ ( LS *J + C1) / MPM  


                   ! Compute the intercepted snow exposure coefficient
                   ! 
                   IF ((SNCMA(I)/SCAP(I)) <=  0.0) THEN
                      CE = 0.07
                   ELSE
                     CE = KS * (SNCMA(I)/SCAP(I))**(-1.0*FRACT) 
                   ENDIF

                   ! Calculate 'potential' canopy sublimation [s-1]
                   SUB_POT = SUB_RATE * CE

                   ! Limit sublimation to canopy snow available and take sublimated snow away from canopy snow at timestep start
                   SUB_CPY = MAX(0.,-SNCMA(I)*SUB_POT*DT)  ! Ensure that only sublimation is computed (neglect solid condensation)

                ELSE IF(HSUBL_CANO == 'NONE') THEN
                  SUB_CPY = 0. 
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


                ! Update the snow mass intercepted in the canopy
                !SNCMA(I) = SNCMAT(I)

                write(*,*) 'After Unload', SNCMA(1), UNLOAD(1),T(1)
                write(*,*) 'Drip', DRIP_CPY(1),'Solid', SNW_UNLOAD

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
