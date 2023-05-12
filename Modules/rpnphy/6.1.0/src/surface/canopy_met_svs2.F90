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

      SUBROUTINE CANOPY_MET_SVS2 (T, HU,PS,VMOD, ISW, ILW, TVEG, ZU,ZT, &
                                 ZENITH, ISW_CAN,ILW_CAN,               &
                                 LAIVH,VEGTRANS,SKYVIEW,EMISVH,N )

     use tdpack
      implicit none

      INTEGER N

      REAL T(N), HU(N), PS(N), VMOD(N)
      REAL ISW(N), ILW(N),ZENITH(N) 
      REAL TVEG(N,2)
      REAL ZU(N), ZT(N)
      REAL VEGTRANS(N),SKYVIEW(N), EMISVH(N), LAIVH(N)
      REAL ISW_CAN(N), ILW_CAN(N)


!
!
!Author
!             V. Vionnet (Mai 2023)
!
!Object
!             Module to simulate the impact of forest of near surface
!             meteorological variables in SVS2:
!                 - radiative fluxes
!                 - wind speed
!                 - air temperature and relative humidity
!

      
      INTEGER I


      REAL KT ! Sky clearness index
      REAL DFRAC ! Fractrion of diffuse radiation 
      REAL FVEG
      REAL DIR_SW,SCA_SW ! Incoming direct and diffuse radiation
      REAL TDIF, TDIR ! Canopy transmittivity for direct and diffuse radiation

      LOGICAL LRAD_SVS ! Use the formulation of SVS for the impact of forest on radiative fluxes 
                       ! Else use the FSM approach

     
     ! Parameters used in the radaitive code
      REAL, PARAMETER :: KEXT = 0.5    ! Vegetation light extinction coefficient              

     !
     ! 0. Initialize parameters and options
     !

      LRAD_SVS = .TRUE.


     !
     ! xxxxx. Compute impact of forest on incoming radiative fluxes 
     !
       !write(*,*) 'Above Can',ISW(1),ILW(1)
       !write(*,*) 'Zenith',ISW(1),ZENITH(1)*180/3.14159

      IF(LRAD_SVS) THEN
           DO I=1,N
               ! Vegtrans computed for high vegetation in vegi_svs2 and derived from Sicart et al (2004)
               ISW_CAN(I) = ISW(I) * VEGTRANS(I)
               ! Skyview computed for high vegetation in vegi_svs2 and derived from Verseghy et al. (1993),
               ! Use deep temperature from the force restore
               ILW_CAN(I) = SKYVIEW(I) * ILW(I) + (1 - SKYVIEW(I))*EMISVH(I) * STEFAN * TVEG(I,2)**4.

               !write(*,*) 'Can SVS',ISW_CAN(I),ILW_CAN(I)
           ENDDO
      ELSE
           DO I=1,N
              ! Compute sky clearness
              IF(ZENITH(I)>0) THEN
                  KT = ISW(I)/ (1367.*COS(ZENITH(I)))
              ELSE
                  KT = 0.
              ENDIF

              ! Compute fraction of diffuse radiation 
              IF(KT>0.8) THEN
                   DFRAC  = 0.165
              ELSE IF(KT>0.22) THEN
                   DFRAC = 0.95 - 0.16*KT + 4.39*KT**2. - 16.64*KT**3. +12.34*KT**4. 
              ELSE
                    DFRAC = 1 - 0.09*KT
              ENDIF

              SCA_SW = DFRAC*ISW(I)  ! Diffuse incoming SW
              DIR_SW = (1.-DFRAC)*ISW(I) ! Direct incoming SW

              FVEG =  0.29 * LOG(LAIVH(I))+ 0.55

              TDIF = EXP(-1.6*KEXT*FVEG)
              IF(ZENITH(I)>0) THEN
                 TDIR = EXP(-KEXT*LAIVH(I)/COS(ZENITH(I)))
              ELSE
                 TDIR = TDIF
              ENDIF

              ! Subcanopy shortwave radiation
              ISW_CAN(I) = TDIF * SCA_SW + TDIR *DIR_SW

               ! Subcanopy longwave radiation 
               ! Use forcing air temperature as a proxy
              ILW_CAN(I) = TDIF * ILW(I) + (1-TDIF)*STEFAN *T(I)**4

           ENDDO

      ENDIF




      END SUBROUTINE CANOPY_MET_SVS2
