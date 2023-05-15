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

      SUBROUTINE CANOPY_MET_SVS2 (T, HU,PS,VMOD, VDIR,  ISW, ILW, TVEG, ZU,ZT,          &
                                 FCOR, LAT, SNVDP,TSNOWV_SVS, TPSOILV,  &
                                 SUNCOS,ISW_CAN,ILW_CAN,VMOD_CAN,T_CAN,HU_CAN,   &
                                 HREF_SUBCANO, Z0MVH, & 
                                 VEGH, LAIVH,VEGTRANS,SKYVIEW,EMISVH,VEGHEIGHT, N )

     use tdpack
     use sfclayer_mod,   only : sl_prelim,sl_sfclayer,SL_OK
     USE MODE_THERMOS 
      implicit none

      INTEGER N

      REAL T(N), HU(N), PS(N), VMOD(N), VDIR(N)
      REAL ISW(N), ILW(N),SUNCOS(N) 
      REAL TVEG(N,2),TSNOWV_SVS(N),SNVDP(N),TPSOILV(N)
      REAL ZU(N), ZT(N),Z0MVH(N),FCOR(N),LAT(N)
      REAL VEGTRANS(N),SKYVIEW(N), EMISVH(N), LAIVH(N),VEGHEIGHT(N),VEGH(N)
      REAL ISW_CAN(N), ILW_CAN(N)
      REAL VMOD_CAN(N), T_CAN(N),HU_CAN(N),HREF_SUBCANO(N)
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

      
      INTEGER I,J


      REAL KT ! Sky clearness index
      REAL DFRAC ! Fractrion of diffuse radiation 
      REAL FVEG
      REAL DIR_SW,SCA_SW ! Incoming direct and diffuse radiation
      REAL TDIF, TDIR ! Canopy transmittivity for direct and diffuse radiation

      LOGICAL LRAD_SVS ! Use the formulation of SVS for the impact of forest on radiative fluxes 
                       ! Else use the FSM approach

      REAL ISW_CAN2, ILW_CAN2

      REAL VMOD_TOP ! Wind speed at canopy top
      REAL DH ! Displacement height
      REAL ZUREF ! Height of teh wind speed forcing abve the ground
      LOGICAL LFORCING_ABOVE ! True if forcing height is given above the canopy  (forcing from NWP system or reanalysis)
                          ! False is forcing height is the absolute height above the ground (mdoel driven by observed forcing) 

      REAL, DIMENSION(N) :: TSREF !  
      REAL, DIMENSION(N) :: QSREF !  
     ! Parameters used in the radiative code
      REAL, PARAMETER :: KEXT = 0.5    ! Vegetation light extinction coefficient              

      ! Parameters used in the wind code
      REAL, PARAMETER :: RCHD = 0.67    ! Ratio of displacement height to canopy height
      REAL, PARAMETER :: HSUBCANO = 2. ! Sub canopy reference height for wind, tair and hu  
      REAL, PARAMETER :: WCAN = 2.5        ! Canopy wind decay coefficient      

     !
     ! 0. Initialize parameters and options
     !

      LRAD_SVS = .TRUE.

      LFORCING_ABOVE  = .FALSE.

     !
     !  xxxx. Compute impact of forest on wind speed
     !  
     ! 
      DO I=1,N
         

         IF(VEGH(I)>0) THEN   ! High veegtation present in the grid cell

           IF(ZU(I) > HSUBCANO)    THEN
                 HREF_SUBCANO(I) =  HSUBCANO
           ELSE 
                 HREF_SUBCANO(I) =  ZU(I)
           ENDIF

           IF(LFORCING_ABOVE) THEN
             ZUREF =  VEGHEIGHT(I) + ZU(I)
           ELSE    
             ZUREF = ZU(I)
           ENDIF

           write(*,*) 'HForc', ZU(I),ZUREF,VEGHEIGHT(I) 

           IF(ZUREF> VEGHEIGHT(I)) THEN  ! Wind speed frocing above the canopy
             ! Displacement height 
             DH = VEGHEIGHT(I)*RCHD 

             ! Wind speed at canopy top, assuming logarithmic profile above canopy, dense canopy 
             VMOD_TOP = VMOD(I) * LOG((VEGHEIGHT(I)-DH)/Z0MVH(I))/LOG((ZUREF-DH)/Z0MVH(I))

             ! Wind speed at canopy base height, dense canopy, assuming exp profile between canopy top and canopy base height
             VMOD_CAN(I) = VMOD_TOP*EXP(WCAN*(HSUBCANO/VEGHEIGHT(I)-1))

             write(*,*) 'Wind',VMOD(I) ,VMOD_TOP, VMOD_CAN(I)
           ELSE IF(ZUREF> HSUBCANO) THEN ! Wind speed frocing in the canopy
             ! Wind speed at canopy base height, dense canopy, assuming exp profile between canopy top and canopy base height
             VMOD_CAN(I) = VMOD(I)*EXP(WCAN*(HSUBCANO/ZUREF-1.))

           ELSE   ! Wind speed frocing below the canopy bottom
             VMOD_CAN(I) = VMOD(I)
           ENDIF
         ELSE  
             ! No high vegetation present 
             HREF_SUBCANO(I) =  ZU(I)    
             VMOD_CAN(I) = VMOD(I)
         ENDIF    
      ENDDO


     !i = sl_sfclayer(xt_canyon,xq_canyon,xu_canyon,zvdir,xbld_height/2.0,xbld_height/2.0,xt_road(:,1),xq_canyon, &
     !   xz0_road,xz0_road/10.0,zdlat,zfcor,optz0=8,hghtm_diag=zu,hghtt_diag=zt,      &
     !   t_diag=ztdiag)
     !
     !  xxxx. Compute impact of forest on temperature and hunidity
     !  
     ! 
      DO I=1,N         
         IF(SNVDP(I)>0.) THEN 
               TSREF(I) =   TSNOWV_SVS(I)
         ELSE
               TSREF(I) =   TPSOILV(I)
         ENDIF                 
      ENDDO

      write(*,*) 'Snow',SNVDP(1),TSNOWV_SVS(1)

      QSREF(:) = QSAT(TSREF(:),PS(:))

       J = SL_SFCLAYER(T,HU,VMOD,VDIR,ZU,ZT,TSREF,QSREF, &  
              Z0MVH,Z0MVH/300,LAT,FCOR,hghtt_diag=HSUBCANO,t_diag=T_CAN)   

     write(*,*) 'TCAN', T(1), TSREF(1), T_CAN(1)


     !
     !  xxxx. Compute impact of forest on incoming radiative fluxes 
     !
       !write(*,*) 'Above Can',ISW(1),ILW(1)
       !write(*,*) 'Zenith',SUNCOS(1)

      !IF(LRAD_SVS) THEN
           DO I=1,N
               ! Vegtrans computed for high vegetation in vegi_svs2 and derived from Sicart et al (2004)
               ISW_CAN(I) = ISW(I) * VEGTRANS(I)
               ! Skyview computed for high vegetation in vegi_svs2 and derived from Verseghy et al. (1993),
               ! Use deep temperature from the force restore
               ILW_CAN(I) = SKYVIEW(I) * ILW(I) + (1 - SKYVIEW(I))*EMISVH(I) * STEFAN * TVEG(I,2)**4.

           ENDDO
      !ELSE
           DO I=1,N
              ! Compute sky clearness
              ! SUNCOS: cos of zenithal angle (max when sun is at its
              ! highest above the horizon).
              IF(SUNCOS(I)>0) THEN
                  KT = ISW(I)/ (1367.*SUNCOS(I))
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

              TDIF = EXP(-1.6*KEXT*LAIVH(I))
              IF(SUNCOS(I)>0) THEN
                 TDIR = EXP(-KEXT*LAIVH(I)/SUNCOS(I))
              ELSE
                 TDIR = TDIF
              ENDIF

              ! Subcanopy shortwave radiation
              ISW_CAN(I) = TDIF * SCA_SW + TDIR *DIR_SW
              !ISW_CAN2 = TDIF * SCA_SW + TDIR *DIR_SW

               ! Subcanopy longwave radiation 
               ! Use forcing air temperature as a proxy
               ILW_CAN(I) = TDIF * ILW(I) + (1-TDIF)*STEFAN *T(I)**4
              !ILW_CAN2 = TDIF * ILW(I) + (1-TDIF)*STEFAN *T(I)**4

               !write(*,*) 'Dif',ISW(1),DIR_SW,SCA_SW
               !write(*,*) 'KT',KT, DFRAC
               !write(*,*) 'Trans',VEGTRANS(1),TDIR,TDIF
               !write(*,*) 'Can SW',ISW(1),ISW_CAN(1),ISW_CAN2
               !write(*,*) 'Can LW',ILW(1),ILW_CAN(1),ILW_CAN2 
           ENDDO

      !ENDIF




      END SUBROUTINE CANOPY_MET_SVS2
