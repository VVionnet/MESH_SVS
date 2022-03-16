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
      SUBROUTINE SNOW_CANOPY_SVS (T, T2M, HU,  RHOA,PS,  U10M, V10M, ISWR,    & 
                     RR, SR, SNCMA, ESUBSNWC,SUBSNWC_CUM, LAIVH,HTVH,VEGH,    &
                     RR_VEG, SR_VEG, DT,  N)


      use tdpack
      use sfclayer_mod, only: sl_sfclayer,SL_OK
      use sfc_options
      use svs_configs
      implicit none

      INTEGER N


      REAL SNCMA(N), RR(N), SR(N),  RR_VEG(N), SR_VEG(N)
      REAL HTVH(N), PS(N), RHOA(N), HU(N),ISWR(N)
      REAL LAIVH(N), T(N), T2M(N), VEGH(N),U10M(N), V10M(N)
      REAL ESUBSNWC(N), SUBSNWC_CUM(N)
      REAL DT
!
!
!Author
!             V. Vionnet (March 2022)
!
!Object
!             Module to simulate interception of snow by high vegetation and
!             sublimation and unloading of intercepted snow
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
! T2M         2-m air temperature [K]
! HU          Specific humidity of air at the model lowest level [kg kg-1]
! RHOA        air density at the forcing level [kg m-3]
! U10M        U-component of wind at 10 m [m/s]
! V10M        V-component of wind at 10 m [m/s]
! ISWR        Solar radiation incident at the top of the high-vegetation [W m-2]
! PS          Surface pressure [Pa]

!             - Input (other parameters)
! DT          time step [s]
! LAIVH       Vegetation leaf area index for HIGH vegetation only [m2/m2]
! HTVH        Vegetation height for HIGH vegetation only [m]
! VEGH        fraction of HIGH vegetation [0-1]


!             - Output 
!
! ESUBSNWC      latent heat of sublimation from incercepted snow [W m-2]
! SUBSNWC_CUM   Cumulated mass loss due to sublimation of incercepted snow [kg m-2]
! RR_VEG       liquid precipitation rate below high-vegetation [kg/m2/s]
! SR_VEG       solid  precipitation rate below high-vegetation [kg/m2/s]




      INTEGER I

      REAL, DIMENSION(N) ::    RHOS, LSTAR,DIRECT_SNOW, UFHT, I1, CC
      REAL, DIMENSION(N) ::    TI, DRIP_CPY, SNW_UNLOAD, NET_SNOW
      
      REAL XI2,EXT2, WINDEXT2, UVENT, ESI, EVAP, LAMBDA, SSTAR, A1, B1,  LS, J
      REAL C1, SIGMA, SVDENS, SUB_RATE,SUB_POT,SUB_CPY, ALPHA, MPM, CE
      REAL NU, NR, MU,DVAP
      
     ! Parameters used in the snow interception code
     !
      REAL   SBAR  

     ! Parameters used in the sublimation and unlodading code
     !
      REAL  GAMA , ZVENT, RADIUS_ICESPH,ALBEDO_ICESPH, KS, FRACT
      REAL UNLOAD_T, UNLOAD_T_WATER, CUNLOAD
    

     ! Physical constant used in the code 

      REAL, SAVE :: MWAT = 0.0180153 ! Molecular weight of water [kg mol-1]
      REAL, SAVE :: RGAZ = 8.314  !  Universal gas constant [J mol-1 K-1]
      REAL, SAVE :: CICE = 2.102e3; !  Heat capacity of ice [J kg-1 K-1]

     ! Options in the code
      CHARACTER(LEN=4) :: HDIFU


     !
     ! 0. Configure the options used in the code
     !
      HDIFU = 'HP98' ! Selection the formulation used to compute the molecular  diffusivity  of  water  vapour  in  air
                      !    'HP98': formulation used in HP98
                      !    'HP13': formulation used in Harder and Pomeroy (2013)
                      !    'TM66': original formulation used in Thorpe and Mason (1966) including the contribution of air pressure
                      !    'PU98': formulation used in Pruppacher et al.(1998)

     !
     ! 1. Initialize Parameters used in the code
     !

     ! Parameters used in the snow interception code

      SBAR = 6.6  ! Maximum canopy snow interception load, [kg m-2]
                  ! 6.6 for pine; 5.9 for Spruce. (Hedstrom and Pomeroy, 1998)

     ! Parameters used in the sublimation and unlodading code

      GAMA = 1.15
      ZVENT = 0.75  ! Ratio between ventilation wind speed height relative and tree height [-]
      RADIUS_ICESPH = 5e-4 ! Radius of single 'ideal' ice shpere [m]
      ALBEDO_ICESPH = 0.8  ! Albedo of single 'ideal' ice shpere [-] 
      KS = 0.0114          ! Snow shape coefficient for jack pine. Taken from Pomeroy et al. (1998)
      FRACT = 0.37 ! Fractal dimension of intercepted snow [-]
      UNLOAD_T = 274.15 ! if ice-bulb temp >= unload_t : canopy snow is unloaded as snow
      UNLOAD_T_WATER = 277.15 ! if ice-bulb temp >= unload_t_water : canopy snow is unloaded as liquid water
      CUNLOAD = -1./(3.5*24*3600)*LOG(0.678)   ! Unloading mass rate (used in exponential decay of intercepted snow, Eq 14 in HP98)
                                               ! The value for CUNLOAD is taken from Mahat and Tarboton (WRR, 2013) 

     !
     ! 2. Evolution of snow mass intercepted by the canopy
     !

      IF(SNCMA(1)>0 .or. SR(1)>0.) THEN
          write(*,*) 'In Cano',SNCMA(1),SR(1)*DT
      ENDIF

      ! Initialize snowfall mass directly transfered through the canopy
      DIRECT_SNOW(:) = 0.

      ! Initialize canopy drip
      DRIP_CPY(:) = 0.


      DO I=1,N
          IF(SNCMA(I)> 0. .OR. SR(I)>0.) THEN  ! Snow is present on the canopy or occurrence of snowfall

                !!!!!!!!
                ! Interception Code 
                !!!!!!!!

                ! Compute density of falling snow (kg m-3) (Eq 13 in HP98)
                 !RHOS(I)  =  67.92 + 51.25 * exp((T2M(I)-273.15) / 2.59) 
                 RHOS(I)  =  67.92 + 51.25 * exp((T(I)-273.15) / 2.59)      ! CHOICE OF TEMPERATURE NEED TO BE REVISED

        
                ! Compute the maximal intercepted snow load (kg m-2) (Eq 12 in HP98)
                ! SBAR may depend on the type of trees (see comment above)
                 LSTAR(I) =  SBAR * (0.27 *46.0 / RHOS(I)) * LAIVH(I) 

                 write(*,*) 'RHO,LSTAR',RHOS(I),LSTAR(I)

                 ! Update intercepted snow mass if case of changes in the maximal intercepted snow mass 
                 IF(SNCMA(I) > LSTAR(I)) THEN
                     DIRECT_SNOW(I) = SNCMA(I) - LSTAR(I)  ! To added to the snowfall below high vegetation. 
                     SNCMA(I) = LSTAR(I)     
                 ENDIF


                 ! Wind at the top of the canopy (m/s) 
                 ! Taken so far equals to the 10 m wind speed. 
                 ! NEED TO BE REVISED TO ACOUNT FOR ACTUAL CANOPY HEIGHT 
                  UFHT(I)  = (U10M(I)**2+V10M(I)**2)**0.5

                 ! Compute horizontal Canopy coverage
                  CC(I) = 0.29 * LOG(LAIVH(I))+ 0.55
                  CC(I) = MIN(1., MAX(0., CC(I)))

                  WRITE(*,*) 'T,UFHT, CC, LAI',T(I),  UFHT(I), CC(I),LAIVH(I)

                  ! Compute amount of snow interecpted in the canopy
                  ! See the discussion pp 1618-1619 in HP98
                  ! For wind speed below 1 m/s the canopy-leaf contact areas per unit ground area, CL is supposed to be:
                  !   CL  = CC
                  ! For wind speed larger than 1 m/s, CL  = 1 
                  ! The origin of the condition SR(I)/LSTAR(I) <50. is not clear
                  I1(I) = 0. 
                  IF( SR(I)>0. .AND. SR(I)/LSTAR(I) <50.) THEN 
                        IF(UFHT(I)<1.) THEN
                             I1(I)  = (LSTAR(I) -SNCMA(I))* (1.0 - EXP(-CC(I) * SR(I)*DT/LSTAR(I))) 
                        ELSE
                             I1(I)  = (LSTAR(I) -SNCMA(I))* (1.0 - EXP(-SR(I)*DT/LSTAR(I)))
                        ENDIF
                        IF(I1(I)<0) I1(I) = 0.
                        
                        SNCMA(I) = SNCMA(I)+I1(I)

                        ! Update the amount of snow that falls below the high vegetation
                        DIRECT_SNOW(I)  = DIRECT_SNOW(I)+ SR(I)*DT - I1(I)
                  ENDIF

                  write(*,*) 'Af inter',SNCMA(I),I1(I),DIRECT_SNOW(I) 


                !!!!!!!!
                ! Sublimation Code 
                !!!!!!!!

                ! Fraction of the entire forest height [-]
                XI2 = 1-ZVENT

                ! Canopy wind speed extinction coefficient [-]
                ! EL10 refers to Eagleson (2002) to justify the formulation of this coefficient
                EXT2 = GAMA * LAIVH(I) 

                ! Computation of ventilation wind speed of intercepted snow derived from above-caopny wind speed  [Eq 8 in EL10]
                ! Estimated within canopy wind speed at fraction XI2 of the entire tree height [Eq 8 in EL10]  [m s-1]
                UVENT = UFHT(I) * EXP(-1 * EXT2 * XI2)

                ! Saturation vapor pressure with respect to ice (Buick, 1981) [Pa]
                ESI = 611.15 * EXP(22.452 * (T(I)-273.15) / (T(I)-273.15 + 272.55))

                ! Vaport density at saturation 
                SVDENS = ESI * MWAT / (RGAZ * T(I))

                ! Derive partial vapor pressure from speficic humidity and pressure [Pa]
                EVAP = HU(I) * PS(I) / (0.622+0.378 * HU(I) ) 

                ! Sutherland's equation for kinematic viscosity
                MU=1.8325e-5*416.16/( T(I)+120)*(T(I)/296.16)*SQRT(T(I)/296.16)/RHOA(I)

                ! Compute thermal conductivity of air [J m-1 s-1 K-1]
                LAMBDA = 0.000063 * T(I) + 0.00673

                ! Compute latent heat of sublimation [J kg-1]
                ! We use here the formulation of Roges and Yau (1989) as in Harder and Pomeroy (2013)
                ! Note that the canopy module in CRHM used a constant value
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
                SIGMA = EVAP/ESI - 1.0  

                write(*,*) 'LS, SVDENS',LS, SVDENS

                ! Compute the mean mass of a snow particle assuming that the size distribution of 
                ! intercepted snow particle follow a gamma distribution with alpha  = 5
                ! and a mean radius equals to RADIUS_ICESPH
                ! This approach is used in PBSM3D (see EQ 23 in Pomeroy et al (1993))
                ! However, it is never described in the papers about the canopy module in CRHM
                ALPHA = 5.0
                MPM = 4.0 / 3.0 * PI * 917. * RADIUS_ICESPH**3. * (1.0 + 3.0 / ALPHA + 2.0 / ALPHA**2.)

                ! Sublimation rate coefficient of single 'ideal' ice sphere [s-1]
                SUB_RATE = (2.0 *PI * RADIUS_ICESPH *SIGMA - SSTAR *J)/ ( LS *J + C1) / MPM  
                ! 

                ! Compute the intercepted snow exposure coefficient
                ! 
                IF ((SNCMA(I)/LSTAR(I)) <=  0.0) THEN
                    CE = 0.07
                ELSE
                    CE = KS * (SNCMA(I)/LSTAR(I))**(-1.0*FRACT) 
                ENDIF

                ! Calculate 'potential' canopy sublimation [s-1]
                SUB_POT = SUB_RATE * CE

                ! Limit sublimation to canopy snow available and take sublimated snow away from canopy snow at timestep start
                SUB_CPY = MAX(0.,-SNCMA(I)*SUB_POT*DT)  ! Ensure that only sublimation is computed (neglect solid condensation)
 
                write(*,*) 'UVENT, NR, NU, MU',UVENT, NR, NU, MU
                write(*,*) 'SIGMA',SIGMA
                write(*,*) 'SUB_CPY',SUB_CPY

                ! Remove mass of intercepted lost by sublimation
                IF(SUB_CPY>SNCMA(I)) THEN
                   SUB_CPY = SNCMA(I)
                   SNCMA(I) = 0.
                ELSE
                   SNCMA(I) = SNCMA(I) - SUB_CPY
                ENDIF

                SUBSNWC_CUM(I) = SUBSNWC_CUM(I) +  SUB_CPY
                ESUBSNWC(I) = SUB_CPY/DT!*LS

                write(*,*) 'Af sub',SNCMA(I),ESUBSNWC(I),SUBSNWC_CUM(I)

                !!!!!!!!
                ! Unloading Code 
                !!!!!!!!

                ! Compute ice bulb temperature Ti [K]
                TI(I) =   T(I) - SUB_POT * LS / CICE 

                IF(TI(I) > UNLOAD_T_WATER) THEN ! All intercepted snow unload as liquid water over 6 hours 
                    DRIP_CPY(I) = SNCMA(I) * DT/ (6.*3600) ! Unload in liquid Form (drip)
                    SNCMA(I)  = SNCMA(I) - DRIP_CPY(I)  

                    write(*,*) 'Drip unload',DRIP_CPY(I)

                ELSE IF(TI(I) >=  UNLOAD_T) THEN  ! All intercepted snow unload as snow over 6 hours 
                    SNW_UNLOAD(I) = SNCMA(I) * DT/ (6.*3600) ! Unload in solid Form
                    SNCMA(I)  = SNCMA(I) - SNW_UNLOAD(I)

                   write(*,*) 'Warm unload',SNW_UNLOAD(I)


                ELSE     ! Cold conditions: unloading assumed to follow an exponential decay (HP98) 
                    
                   SNW_UNLOAD(I) = SNCMA(I) * CUNLOAD * DT

                   write(*,*) 'Cold unload',SNW_UNLOAD(I)
                   IF(SNW_UNLOAD(I)> SNCMA(I)) THEN
                        SNW_UNLOAD(I) = SNCMA(I)
                        SNCMA(I)  = 0.
                   ELSE
                        SNCMA(I) = SNCMA(I)-SNW_UNLOAD(I)
                   ENDIF
                ENDIF

                write(*,*) 'Af unload',SNCMA(I) 
              

                NET_SNOW(I) = DIRECT_SNOW(I) + SNW_UNLOAD(I)


             ENDIF
      ENDDO   

      DO I=1,N

         IF (VEGH(I).GE.EPSILON_SVS) THEN

            ! Rain is not intercepted by vegetation when snow is present under high veg. 
            ! Therefore rain under high veg accounts for total rain above high veh plus dripping 
            RR_VEG(I) = RR(I)+DRIP_CPY(I)/DT  

            ! Snowfall rate below high vegetation
            SR_VEG(I) = DIRECT_SNOW(I)/DT

             write(*,*) 'Rain veg',RR_VEG(I),RR(I)
             write(*,*) 'Snow veg',SR_VEG(I),SR(I)

         ELSE
            ! no high veg, no modification of rainfall and snowfall

            RR_VEG(I) = RR(I)
            SR_VEG(I) = SR(I)

         ENDIF
      ENDDO

      END

