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

      SUBROUTINE EBUDGET_SVS2_ONEPROFILE_SKIN(TSA, WD, WF , &
                   TGRS,TGRD,TVGLS,TVGHS,TVGHD, TP, TPERM,  &                       
                   PGRNDFLUX, PGRNDFLUXV, DT, VMOD, VDIR, LAT, & 
                   RG, ALVG, LAI, GAMVEG, ALVL, ALVH,  & 
                   ALGR,EMGR, & 
                   RAT, THETAA, FCOR, ZUSL, ZTSL, HU, PS, &  
                   RHOA, WTA, Z0, Z0LOC, Z0H, & 
                   HRSURF, HV_VL, HV_VH, DEL_VL, DEL_VH, RS, & 
                   CG,CVP, EMISVL, EMISVH, PSNG, &  
                   RESAGR, RESA_VL, RESA_VH, RESASA, RESASV, RESAGRV, &
                   RNETSN, HFLUXSN,LESLNOFRAC, LESNOFRAC, ESNOFRAC, & 
                   ALPHAS, &  
                   TSNS, & 
                   RNETSV, HFLUXSV,LESLVNOFRAC, LESVNOFRAC, ESVNOFRAC,  &
                   ALPHASV, & 
                   TSVS, & 
                   VEGH, VEGL, VGHEIGHT, PSNVH,PSNVHA, PSURFVHA, & 
                   SKYVIEW,SKYVIEWA, &
                   SOILHCAP, SOILCOND, &  
                   RR,WR_VL,WR_VH,SNM,SVM, &
                   VTRA, ALBT, & 
                   RNET, HFLUX, LE, LEG, LEVL, LEVH, LES,LESV, LEGV,  & 
                   LER_VL, LETR_VL,LER_VH, LETR_VH, &
                   EG, ER_VL, ETR_VL ,ER_VH, ETR_VH, GFLUX, EFLUX, & 
                   BM, FQ, BT, RESAEF, & 
                   LEFF, & 
                   FTEMP, FVAP, ZQS, FRV, & 
                   ALFAT, ALFAQ, ILMO, HST, TRAD, N,    &
                   QVEG, QGV, QGR, &
                   TAF, QAF, VAF, & 
!                   RGVG,FIVG,IRGV,IRVG, &
!                   HGV,LGV, &
!                   FGRV, RNGV,grflux, &
                    RPP,Z0HA)

      use tdpack
      use sfclayer_mod, only: sl_sfclayer,SL_OK
      use sfc_options
      use svs_configs
      use svs2_tile_configs

      implicit none
!!!#include <arch_specific.hf>


      INTEGER N!, NSL
      
      REAL TSA(N),  DT, VMOD(N)
      REAL VDIR(N), LAT(N), PGRNDFLUX(N), PGRNDFLUXV(N)
      REAL TGRS(N), TGRD(N), TVGLS(N), TVGHS(N), TVGHD(N)
      REAL TSNS(N,NSL), TSVS(N,NSL)
      REAL TP(N,NL_SVS), TPERM(N)
      REAL WD(N,NL_SVS), WF(N,NL_SVS)   
      REAL RG(N), ALVG(N), ALVL(N), ALVH(N), RAT(N), THETAA(N), FCOR(N)
      REAL ZUSL(N), ZTSL(N)
      REAL HU(N), PS(N), RHOA(N), WTA(N,svs2_tilesp1), Z0(N)
      REAL Z0LOC(N), Z0H(N)
      REAL HV_VL(N), HV_VH(N), DEL_VL(N), DEL_VH(N),  RS(N)
      REAL CG(N), CVP(N),  PSNG(N), EMISVL(N), EMISVH(N)
      REAL LAI(N), GAMVEG(N), ALGR(N), EMGR(N)
      REAL RNET(N), HFLUX(N), LE(N), ALPHAS(N)
      REAL ALBT(N)
      REAL LEG(N), LEVL(N), LEVH(N), LER_VL(N), LETR_VL(N),LEGV(N), GFLUX(N)
      REAL LER_VH(N), LETR_VH(N)
      REAL EFLUX(N), BM(N), FQ(N), BT(N), LES(N)
      REAL FTEMP(N), FVAP(N), ER_VL(N), ETR_VL(N),ER_VH(N), ETR_VH(N) 
      REAL LEFF(N), ZQS(N), FRV(N)
      REAL EG(N), HRSURF(N)
      REAL RESAGR(N), RESA_VL(N), RESA_VH(N), RESASA(N), RESASV(N), RESAEF(N)
      REAL RNETSN(N), HFLUXSN(N),LESLNOFRAC(N),LESNOFRAC(N), ESNOFRAC(N)
      REAL RNETSV(N), HFLUXSV(N),LESLVNOFRAC(N),LESVNOFRAC(N), ESVNOFRAC(N)
      REAL ALPHASV(N)
      REAL ALFAT(N), ALFAQ(N), LESV(N)
      REAL VEGH(N), VEGL(N), PSNVH(N), PSNVHA(N),PSURFVHA(N)
      REAL SKYVIEW(N),SKYVIEWA(N), ILMO(N), HST(N), TRAD(N), VTRA(N)
      REAL WR_VL(N), WR_VH(N), RR(N),SNM(N),SVM(N)
      REAL VGHEIGHT(N)
      REAL SOILHCAP(N,NL_SVS),SOILCOND(N,NL_SVS)

!  ajout temporaire pour tests
      REAL VAF(N), TAF(N), QAF(N), QVEG(N), QGV(N), QGR(N)
!     REAL RGVG(N), FIVG(N), IRGV(N), IRVG(N), HGV(N), LGV(N), FGRV(N), grflux(N)
      REAL RNGV(N)
      REAL RPP(N)
      REAL Z0HA(N), RESAGRV(N)



!
!Author
!          S. Belair et al. (January 2016)
!Revisions
! 001      Bug fixes: M. Abrahamowicz, S. Z. Husain, N. Gauthier, E. Gaborit, 
!          V. Vionnet, D. Deacu      
!
!Object
!
!     Calculates the evolution of the surface and deep-soil temperature
!     (i.e., Ts and T2), as well as all the surface fluxes.
!
!
!!!!  METHOD
!!    ------
!
!     1- find the grid-averaged albedo, emissivity, and roughness length
!     2- compute the za, zb, and zc terms involved in the numerical
!        resolution of the equations for Ts and T2.
!     3- find Ts(t) and T2(t).
!     4- derive the surface fluxes.!
!
!Arguments
!
!
!          - Input/Output -
! TGRS      (bare) ground temperature -- S for "skin"
! TGRD      mean ground temperature -- D for "deep"
! TVGLS     low vegetation temperature -- S for "skin"
! TVGHS     high vegetation temperature -- S for "skin"
! TVGHD     mean high vegetation temperature -- D for "deep"      
! TS        surface  temperature (new) as seen from ground
! Z0H       agg. thermal roughness length for land surface
!           Output only when svs_dynamic_z0h=.true.
! TP(:,:)   soil temperature profile 

!          - Input -
! VMOD      module of the low-level wind
! VDIR      wind direction at the lowest level
! LAT       latitude
! WD        soil water content 
! WF        frozen soil water
! DT        timestep
! RG        global radiation (downward solar)
! ALVG      AVERAGED surface albedo associated with vegetation type
! LAI       AVERAGED vegetation leaf area index
! GAMVEG    AVERAGED parameter related to the vegetation height
! ALGR      albedo of bare ground (soil)
! EMGR      emissivity of bare ground (soil) 
! EMISVL    emissivity of low vegetation
! EMISVH    emissivity of high vegetation
! ALVL      albedo of low vegetation
! ALVH      albedo of high vegetation
! RAT       atmospheric radiation incident on the ground (NIR)
! THETAA    air potential temperature at the lowest level
! FCOR      Coriolis factor
! ZTSL      reference height for temperature and humidity input
! ZUSL      reference height for wind input
! HU        specific humidity of air at the lowest level
! PS        surface pressure
! RHOA      air density near the surface
! WTA       Weights for SVS surface types as seen from SPACE
! Z0        momentum roughness length (no snow)
! Z0LOC     local land momentum roughness length (no orography) 
! HRSURF    relative humidity of the bare ground surface (1st soil layer)
! HV_VL     Halstead coefficient (relative humidity of low veg. canopy)
! HV_VH     Halstead coefficient (relative humidity of high veg. canopy)
! DEL_VL    portion of the low veg. leaves covered by water
! DEL_VH    portion of the high veg. leaves covered by water
! RS        stomatal resistance
! CG        soil thermal coefficient
! CVP       AVERAGED vegetation thermal coefficient (with LAI effect)
! PSNG      fraction of bare ground or low veg. covered by snow
! RESAGR    aerodynamical surface resistance for bare ground
! RESA_VL    aerodynamical surface resistance for low vegetation
! RESA_VH    aerodynamical surface resistance for high vegetation
! RESASA    aerodynamical surface resistance for snow on bare ground/low veg
! RESASV    aerodynamical surface resistance for snow under high veg.
! RNETSN    net radiation over snow 
! HFLUXSN   sensible heat flux over snow 
! ESNOFRAC  water vapor flux from the snow surface (kg/m2/s)
! LESLNOFRAC latent heat flux of evaporation from the snow surface (W/m2)
! LESSNOFRAC latent heat flux of sublimation from the snow surface (W/m2)
! ALPHAS    albedo of snow
! TSNS      snow temperature at time t+dt (update in snow_svs.F90)
!
! RNETSV    net radiation over snow-under-vegetation
! HFLUXSV   sensible heat flux over snow-under-vegetation
! ESVNOFRAC water vapor flux from the snow-under-vegetation  (kg/m2/s)
! LESLVNOFRAC latent heat flux of evaporation from the snow-under-vegetation (W/m2)
! LESSVNOFRAC latent heat flux of sublimation from the snow-under-vegetation (W/m2)
! ALPHASV   albedo of snow-under-veg
! TSVS      snow-under-veg temperature at time t+dt (update in snow_svs.F90)
!
! VEGH     fraction of HIGH vegetation
! VEGL     fraction of LOW vegetation
! PSNVH    fraction of HIGH vegetation covered by snow
! PSNVHA   fraction of HIGH vegetation covered by snow as seen from
!          the ATMOSPHERE 
! PSURFVHA  fraction of surface seen through sparse high vegetation  
! RR       Liquid precipitation rate at the surface in [mm/s]
! WR_VL    Water retained by low vegetation
! WR_VH    Water retained by high vegetation
! SVM      snow water equivalent (SWE) for snow under high veg [kg/m2]
!
!           - Output -
! ALBT      total surface albedo (snow + vegetation + bare ground)
! RNET      net radiation
! HFLUX     sensible heat flux

! LE        latent heat flux
! LEG       latent heat flux over bare ground
! LEVL      latent heat flux over low vegetation
! LEVH      latent heat flux over high vegetation
! LES       latent heat flux over snow
! LESV      latent heat flux over snow-under-veg
! LER_VL       direct latent heat flux from low vegetation leaves
! LETR_VL      evapotranspiration latent heat flux from low veg
! LER_VH       direct latent heat flux from high vegetation leaves
! LETR_VH      evapotranspiration latent heat flux from high veg

! EG         evaporation rate (no fraction) over bare ground (1st soil layer)
! ER_VL        direct evaporation rate (no fraction) from low veg. leaves
! ETR_VL       evapotranspiration rate from low veg.  (no faction) 
! ER_VH       direct evaporation rate (no fraction) from high veg. leaves
! ETR_VH       evapotranspiration rate from high veg. (no faction)  
!
! GFLUX     ground flux
! EFLUX     water vapor flux
! DWATERDT  net tendency of melting-freezing of soil water
! TS        surface  temperature (new) as seen from ground
! TD        mean soil temperature
! TSA       surface  temperature (new) as seen from space
! LEFF      effective latent heat
! FTEMP      land sfc avg turbulent surface flux of temperature
! FVAP       land sfc avg turbulent surface flux of vapor
! ILMO       land sfc avg (1/length of Monin-Obukov) 
! HST        land sfc avg height of the boundary layer 
! FRV        land sfc average friction velocity 
! BM         homogeneous boundary condition term in the
!            diffusion equation for U and V
! BT         homogeneous boundary condition term in the
!            diffusion equation for T and Q
! FQ         land sfc average momentum flux
! RESAEF     effective aerodynamic resistance for land sfc
! ALFAT      inhomogeneous boundary term in the diffusion equation for Theta
! ALFAQ      inhomogeneous boundary term in the diffusion equation for Q
! ZQS       area-averaged specific  humidity of a model tile
! TRAD      averaged radiative temperature

!
!
      INTEGER I,zopt, K, I_FLUX
!
!
      REAL EMISSN, EMSOIL, KCOEF, RHOW
      REAL BFREEZ, RAIN1, RAIN2
      REAL ABARK

      REAL SOILCOND1, SOILCOND2    

      REAL SURF_ENBAL_VEG 

      REAL BETAA 
      DATA BETAA/1.0/                           ! fully-implicit time scheme
 
!     a mettre dans automatic plus tard
      REAL, DIMENSION(NL_SVS)   :: DZ       ! Local variable 
      REAL, DIMENSION(NL_SVS)   :: DELZZ    ! (m) thickness of the soil between 2 temperature levels
      REAL, DIMENSION(N,NL_SVS) :: SOILCD   ! W/(m K) two layers averaged soil heat conductivity 
      REAL, DIMENSION(N,NL_SVS) :: A2, B2, C2, D2, A4, B4, C4, D4   

!
!     MULTIBUDGET VARIABLES 
!     GR:ground, SN:snow, VG:vegetation, AG: aggregated
!     VGL: low vegetation, VGH: high vegetation        
       real, dimension(n) :: a3h, b3h, c3h, zhvgl,zhvgh, freezfrac, emvg, &
            alvglai, zqsatgr, zdqsatgr, zqsatvgl, zdqsatvgl, zqsatgrt, &
            zqsatvglt,zqsatvght,zqsatvgh, zdqsatvgh, &
            rnetgr, rnetvgl,rnetvgh, hfluxgr, hfluxvgl, hfluxvgh, &
            roragr, roravgl,roravgh,    &
            zqsatsno, tgrst, tgrdt, tvgst, tvgdt, esf, esvf, evlf, evhf, &
            tvglst,tvgldt,tvghst,tvghdt,   &            
            egf, ev_vl,ev_vh, zqsatsnv, legnofrac,legvnofrac,  &
            cmu, cm, ctu, vmod_lmin

       real, dimension(n) ::  ZQSATGRV, ZDQSATGRV, LOWVEG, HIGHVEG,  &
                               ZQSGRV,ZQSVGH, ZA,Z0TEMP,Z0HG,TGRVS,              &
                               TA4FLX,QA4FLX,ZU4FLX,ZT4FLX,VIT,Z0M4FLX,Z0H4FLX, &
                               CTUGRV,RORAGRV,TPVPRIM,TAFPRIM,DIFTEMP,ZQSATGRVT

       real, dimension(n) :: ABG, BBG, SKINCOND_BG,AVL, BVL,SKINCOND_VL, &
                      SURF_ENBAL_HVEGB,SURF_ENBAL_HVEGD, FRACBG, FRACVL,& 
                      FRACVH, FRACSNG, FRACSNVH       

       real,dimension(n, svs2_tilesp1) :: wtg                     

!************************************************************************
!
!
!
!                                THE FOLLOWING SHOULD BE PUT IN 
!                                A COMMON COMDECK
!
      EMISSN = 0.97
      EMSOIL = 0.94
      RHOW   = 1000.  
      KCOEF  = 1.E-6
      BFREEZ = 4.
!                                Albedo of Bark (S. Wang, Ecological Modelling, 2005)
      ABARK  = 0.15
!
!
!
!!       1.     GRID-AVERAGED ALBEDO, EMISSIVITY, AND ROUGHNESS LENGTH
!       ------------------------------------------------------
!                          (considering snow surfaces)
!
      DO I=1,N
!
!
!                               Calculate grid-averaged albedo 
!
!                               Using PSNVHA
!        
         ALBT(I)  = AG_SVS2( WTA(I,indx_svs2_bg), WTA(I,indx_svs2_vl), &
                      WTA(I,indx_svs2_vh),WTA(I,indx_svs2_sn),    &
                      WTA(I,indx_svs2_sv), WTA(I,indx_svs2_gv), &
                      ALGR(I),ALVL(I),ALVH(I),            &
                      ALPHAS(I),ALPHASV(I), ALGR(I)) ! VV TO BE MODIFIED 
!
! 
!                               Recalculate vegetation-only albedo to take LAI
!                               effect into account, can only consider it 
!                               if have high Vegetation (i.e. Gamveg greater
!                               or equal to 0.01). In the high vegetation case,
!                               weight albedo of leaves vs. bark according to LAI
!
            IF(GAMVEG(I).GE.0.01)THEN 
               ALVGLAI(I) = MIN( LAI(I)      , LAI0 )   * ALVG(I) / LAI0 &  
                          + MAX( LAI0-LAI(I) , 0.0    ) * ABARK  / LAI0
            ELSE
               ALVGLAI(I) = ALVG(I)
            ENDIF
!
      END DO
!

!
!*      1.bis      FRACTIONS OF SVS2 TYPES AS SEEN FROM GROUND
!               --------------------------

      CALL WEIGHTS_SVS2(VEGH,VEGL,PSNG,PSNVH,PSURFVHA,'GROUND',N,WTG)
!           
!
!      
!        2.     LATENT HEAT COEFFICIENTS - CORRECTION DUE TO FREEZING
!               AND MELTING OF SOIL WATER
!               -----------------------------------------------------
!
!                               Using the fraction of frozen water
!                               in the soil, calculate the "effective"
!                               latent heat of evaporation/sublimation
!
        DO I=1,N
           FREEZFRAC(I) = WF(I,1) / (WD(I,1)+WF(I,1)+EPSILON_SVS)
           LEFF(I)      = FREEZFRAC(I)      * (CHLC+CHLF)  &
                         + (1.-FREEZFRAC(I)) *  CHLC
        END DO
!
!
!
!
!!       3.A.  COEFFICIENTS FOR THE TIME INTEGRATION OF 
!!               BARE GROUND SKIN SURFACE TEMPERATURE      
!               --------------------------------------------
!
!
       DO I=1,N
!       
!         Thermodynamic functions used in the linearisation of the
!         latent heat flux
          ZQSATGR(I)  = FOQST( TGRS(I),PS(I) )
          ZDQSATGR(I) = FODQS( ZQSATGR(I),TGRS(I) )
!
!         Function zrsra used in the computation of the sensible heat
!         flux
!
          RORAGR(I) = RHOA(I) / RESAGR(I)
!
!         Skin conductivity for bare ground
!           
          SKINCOND_BG(I) =  2 * SOILCOND(I,1)/DELZ(1)
!       
          ABG(I) =  SKINCOND_BG(I) + 4.*EMGR(I)*STEFAN*(TGRS(I)**3)  &
               + RORAGR(I)*CPD + RORAGR(I)*LEFF(I)*HRSURF(I)*ZDQSATGR(I)
!
          BBG(I) = SKINCOND_BG(I)*TP(I,1) + (1.-ALGR(I))*RG(I) +EMGR(I)*RAT(I) &
               + 3.*EMGR(I)*STEFAN*(TGRS(I)**4) +RORAGR(I)*CPD*THETAA(I) &
               + RORAGR(I)*LEFF(I)*HRSURF(I)*ZDQSATGR(I)*TGRS(I) &
               - RORAGR(I)*LEFF(I)*(HRSURF(I)*ZQSATGR(I)-HU(I))         
!          
!         Update bare soil skin surface temperature
!       
          TGRST(I) = BBG(I)/ABG(I)

       END DO        
!!            TGRD AT TIME 'T+DT': VV To be changed 
!               -----------------
!
      DO I=1,N
        TGRDT(I) = (TGRD(I) + DT*TGRST(I)/86400.) /   &  
                      (1.+DT/86400.)
      END DO

!
!
!!       3B.     COEFFICIENTS FOR THE TIME INTEGRATION OF 
!!               LOW and HIGH VEGETATION TEMPERATURE      
!               --------------------------------------------
!               Skin temperature for low vegetation
!                  
!		FORCE RESTORE SCHEME FOR HIGH VEGETATION ONLY
!		Note that the ground thermal coefficient is still included when 
!		computing the vegetation thermal coefficient. This will have to be revised! 
!               --------------------------------------------
!         CALCULATE ONLY IF VEGETATION NON-ZERO PRESENT, OTHERWISE USE BARE GROUND TO IMPOSE DEFAULT PHYSICAL VALUE
!
!               --------------------------------------------
!        3B1. Calculation for low vegetation       
!               --------------------------------------------
       DO I=1,N      
          IF ( (VEGL(I)*(1-PSNG(I))).GE.EPSILON_SVS ) THEN
             ! EXPOSED LOW VEGETATION PRESENT
!
!         Thermodynamic functions used in the linearisation of the
!         latent heat flux             
             
             ZQSATVGL(I)  = FOQST( TVGLS(I),PS(I) )
             ZDQSATVGL(I) = FODQS( ZQSATVGL(I),TVGLS(I) )
!
!         Function zrsra used in the computation of the sensible heat
!         flux             
!
             RORAVGL(I) = RHOA(I) / RESA_VL(I)
!
!         Skin conductivity for low vegetation
!           
          SKINCOND_VL(I) =  10.
!             
!       
          AVL(I) =  SKINCOND_VL(I) + 4.*EMISVL(I)*STEFAN*(TVGLS(I)**3)  &
               + RORAVGL(I)*CPD + RORAVGL(I)*CHLC * HV_VL(I)*ZDQSATVGL(I)
!
          BVL(I) = SKINCOND_VL(I)*TP(I,1) + (1.-ALVL(I))*RG(I) +EMISVL(I)*RAT(I) &
               + 3.*EMISVL(I)*STEFAN*(TVGLS(I)**4) +RORAVGL(I)*CPD*THETAA(I) &
               + RORAVGL(I)*CHLC * HV_VL(I)*ZDQSATVGL(I)*TGRS(I) &
               - RORAVGL(I)*CHLC * HV_VL(I)*(ZQSATVGL(I)-HU(I))              
!
!          
!         Update low vegetation skin surface temperature
!       
          TVGLST(I) = BVL(I)/AVL(I)  
!          

          ELSE
             ! NO LOW VEGETATION -- USE BARE GROUND VALUES or ZERO to fill arrays to avoid numerical errors
             ZQSATVGL(I)  =  ZQSATGR(I)
             ZDQSATVGL(I) = ZDQSATGR(I)
             RORAVGL(I) = RORAGR(I)
             TVGLST(I) = TGRST(I)
          ENDIF


       ENDDO
!       
!               --------------------------------------------
!        3B2. Calculation for high vegetation       
!               --------------------------------------------      
!
       DO I=1,N      
          IF ( VEGH(I).GE.EPSILON_SVS ) THEN
             ! HIGH VEGETATION PRESENT
!
!                            Thermodynamic functions
!
             
             ZQSATVGH(I)  = FOQST( TVGHS(I),PS(I) )
             ZDQSATVGH(I) = FODQS( ZQSATVGH(I),TVGHS(I) )
!
!                              function zrsra      
!
             RORAVGH(I) = RHOA(I) / RESA_VH(I)
!
!                              Fraction of high vegetation to total vegetation
!    
!             FRACH(I) =  MIN(   (VEGH(I)*PSNVH(I))/(VEGH(I)+VEGL(I)*(1-PSNG(I))) , 1.0) 
!
!                                        terms za, zb, and zc for the
!                                              calculation of tvgs(t)
             A3H(I) = 1. / DT + CVP(I) *  & 
                    (4. * EMISVH(I) * STEFAN * (TVGHS(I)**3)  &  
                    +  RORAVGH(I) * ZDQSATVGH(I) * CHLC * HV_VH(I) &  
                    +  RORAVGH(I) * CPD )  & 
                    + 2. * PI / 86400.
!
             B3H(I) = 1. / DT + CVP(I) *   &
                    (3. * EMISVH(I) * STEFAN * (TVGHS(I)** 3)  &   
                    + RORAVGH(I) * ZDQSATVGH(I) * CHLC* HV_VH(I) )
           
!
             C3H(I) = 2. * PI * TVGHD(I) / 86400. &   
                     + CVP(I) *  & 
                     ( RORAVGH(I) * CPD * THETAA(I)  &  
                     + RG(I) * (1. - ALVGLAI(I)) + EMISVH(I)*RAT(I)  & 
                     - RORAVGH(I)  & 
                     * CHLC * HV_VH(I) * (ZQSATVGH(I)-HU(I)) )


             TVGHST(I) =  ( TVGHS(I)*B3H(I) + C3H(I) ) / A3H(I)

          ELSE
             ! NO VEGETATION -- USE BARE GROUND VALUES or ZERO to fill arrays to avoid numerical errors
             ZQSATVGH(I)  =  ZQSATGR(I)
             ZDQSATVGH(I) = ZDQSATGR(I)
             RORAVGH(I) = RORAGR(I)
!             FRACH(I) = 0.0
             TVGHST(I) = TGRST(I)
          ENDIF


       ENDDO


!!           TVGD AT TIME 'T+DT': Mean temperature from FR scheme for
!                      high veg
!               -----------------
!
      DO I=1,N
!                   Note that as an added precaution,
!                   we set the vegetation temperature to
!                   that of the ground, when no vegetation is present
!
         IF(VEGH(I).ge.EPSILON_SVS)THEN
            TVGHDT(I) = (TVGHD(I) + DT*TVGHST(I)/86400.) / (1.+DT/86400.)
         ELSE
            TVGHDT(I) = TGRDT(I)
         ENDIF
            
      END DO

!
!
!
!       4.     SURFACE ENERGYE BALANCE FOR GROUND BELOW HIGH VEG
!               --------------------------------------------
!               --------------------------------------------
       DO I=1,N
          IF(VEGH(I).ge.EPSILON_SVS) THEN   

!             first calculate the saturation vapor
!             pressure over ground under canopy
              ZQSATGRV(I)  = FOQST( TP(I,1), PS(I) )
              ZDQSATGRV(I) = FODQS( ZQSATGRV(I),TP(I,1) )
              ZQSATVGHT(I)  = FOQST( TVGHST(I), PS(I) )

              ! Derive specific humidity at the surface of the ground
              ! below high vegetation and over high vegetation 
              ZQSGRV(I) = HRSURF(I) * ZQSATGRV(I)
              ZQSVGH(I) = RPP(I) * ZQSATVGHT(I) + ( 1. - RPP(I) ) *QAF(I) ! VV TO BE MODIFIED

              ! Compute meteorology in the canopy air space 
              VAF(I) = 0.3 * VMOD(I)   !high veg
              TAF(I) = (0.3 * THETAA(I) + 0.6 *TVGHST(I) + 0.1 * TP(I,1))   ! Temp. inside canopy
              QAF(I) = (0.3 * HU(I) + 0.6 * ZQSVGH(I) + 0.1 * ZQSGRV(I))    ! Hum.  inside canopy

              ! Compute height and roughness used on the computation
              ! of the fluxes 
              ZA(I) = 2. * VGHEIGHT(I) / 3.   ! reference height in canopy for flux
              ZA(I) = MAX(2.0,ZA(I))
              Z0TEMP(I) = 1.0                 ! bare ground local momentum roughness
              Z0HG(I)=0.2                     ! bare ground local heat roughness

              ! Derive surface temperature of the soil below high
              ! vegetation. 
              TGRVS(I) = TP(I,1) 

              ! Select variable used to compute the fluxes 
              TA4FLX(I)  = TAF(I)
              QA4FLX(I)  = QAF(I)
              ZU4FLX(I)  = ZA(I)
              ZT4FLX(I)  = ZA(I)
              VIT(I)     = VAF(I)
              Z0M4FLX(I) = Z0TEMP(I)
              Z0H4FLX(I) = Z0HG(I)

              i_flux = sl_sfclayer( TA4FLX, QA4FLX, VIT, VDIR, ZU4FLX, ZT4FLX, &
                    TGRVS, ZQSGRV, Z0M4FLX, Z0H4FLX, LAT, FCOR, &
                    L_min=sl_Lmin_soil,coeft=CTUGRV )

              if (i_flux /= SL_OK) then
                print*, 'Abort. ebud_svs2 bec of err in sl_sfclayer()'
              stop
              endif 

              ! Compute aerodymanical resistance 
              RESAGRV(I) = 1. / CTUGRV(I)

              !   function zrsra ground under veg
              RORAGRV(I) = RHOA(I) / RESAGRV(I)

              !drag coef from CLASS for high veg only

              TPVPRIM(I) = TP(I,1)*(1.0 + 0.61 * ZQSGRV(I) )
              TAFPRIM(I) = TAF(I)*(1.0 + 0.61 * QAF(I) )
              DIFTEMP(I) = TPVPRIM(I) - TAFPRIM(I)
        
              IF (DIFTEMP(I).GT.1.0 )THEN
                 RORAGRV(I) = RHOA(I) * (0.0019 * (TPVPRIM(I) - TAFPRIM(I))**(1./3.) )
              ELSEIF (DIFTEMP(I).GT.0.001 .and. DIFTEMP(I).LE.1.0) THEN
                 RORAGRV(I) = RHOA(I) * (0.0019 * (TPVPRIM(I) - TAFPRIM(I)) )
              ELSE
                 RORAGRV(I) = 0.0
              ENDIF

             ! Contribution of the surface energy budget for high
             ! vegetation in the coefficient used in the tridiagnonal
             ! matrix (see next section)
             ! Contribution to coef D
             SURF_ENBAL_HVEGD(I) = ( VTRA(I)*(1.-ALGR(I))*RG(I) + SKYVIEWA(I)*EMGR(I)*RAT(I) &
                            + 3.*EMGR(I)*STEFAN*(TP(I,1)**4) &
                            +(1.-SKYVIEWA(I))*EMGR(I)*EMISVH(I)*STEFAN*(TVGHDT(I)**4) &
                            + RORAGRV(I)*CPD*TAF(I) &
                            + RORAGRV(I)*LEFF(I)*HRSURF(I)*ZDQSATGRV(I)*TP(I,1) &
                            - RORAGRV(I)*LEFF(I)*(HRSURF(I)*ZQSATGRV(I)-QAF(I)) )

             ! Contribution to coef B
             SURF_ENBAL_HVEGB(I) = 4.*EMGR(I)*STEFAN*(TP(I,1)**3)  +  RORAGRV(I)*CPD &
                             +  RORAGRV(I)*LEFF(I)*HRSURF(I)*ZDQSATGRV(I) 


            ELSE     
              ! 
              ! High vegetation not present in the grid cell       
              SURF_ENBAL_HVEGD(I)= 0.       
              SURF_ENBAL_HVEGB(I)= 0.       
          ENDIF
       END DO                

!
!
!       5.     COEFFICIENTS FOR THE TIME INTEGRATION OF TP(I,K)
!               SOIL TEMPERATURE
!               --------------------------------------------
!               --------------------------------------------

!       Interfacial Soil thermal conductivity
!       Inverse-weighted arithmetic mean of the soil thermal conductivity 
!       at the interface between two consecutive layers
!
       DO I=1,N
!         Thermal conductivity - Mean over two layers
          DO K=1,NL_SVS-1
             SOILCOND1 = SOILCOND(I,K)
             SOILCOND2 = SOILCOND(I,K+1)

            ! Inverse-weighted arithmetic mean of the soil thermal conductivity 
            ! at the interface between two consecutive layers
            SOILCD(I,K) = (DELZ(K)+ DELZ(K+1))/( DELZ(K)/SOILCOND1 + DELZ(K+1)/SOILCOND2)
          END DO
!         special case for deepest soil layers
          SOILCD(I,NL_SVS) = SOILCOND(I,NL_SVS)
       END DO

       ! Compute thickness of the soil layer between 2 temperature levels
       DO K=1,NL_SVS-1
          DELZZ(K) = (DELZ(K) + DELZ(K+1)) / 2.
       END DO
       DELZZ(NL_SVS) = DELZ(NL_SVS)
!  
!
!      Compute fractions of the respective type of cover
!   
!       DO I=1,N
!          ! Fraction of snow-free bare ground   
!          FRACBG(I) = (1-PSNG(I))*(1.-VEGH(I)-VEGL(I))  
!          ! Fraction of snow-free low vegetation  
!          FRACVL(I) = (1-PSNG(I))*VEGL(I)
!          ! Fraction of snow-free ground below high vegetation 
!          FRACVH(I) = (1-PSNVH(I))*VEGH(I)
!          ! Fraction of snow covered bg/lv
!          FRACSNG(I) = PSNG(I)*(1.-VEGH(I))
!          ! Fraction of snow covered ground below high vegetation 
!          FRACSNVH(I) = PSNVH(I)*VEGH(I)
!       END DO 

!
!
!       6.     SOLVE HEAT EQUATION IN THE SOIL
!               TRIDIAGONAL SYSTEM
!               --------------------------------------------
!               --------------------------------------------
!
!                              coefficients A, B, C and term D for the
!                              matrix inversion for the calculation of TP(t) 
!
!
       DO K=2,NL_SVS-1

          DO I=1,N

             A2(I,K) = (-BETAA) * (DT * SOILCD(I,K-1)) / &
                    (SOILHCAP(I,K) * DELZ(K) * DELZZ(K-1))

             C2(I,K) = (-BETAA) * (DT * SOILCD(I,K)) / &
                    (SOILHCAP(I,K) * DELZ(K) * DELZZ(K))

             B2(I,K) = 1.0 - A2(I,K) - C2(I,K)
             
             D2(I,K) = ( (1.-BETAA)*DT / (SOILHCAP(I,K)*DELZ(K)) ) * &  
                    (  SOILCD(I,K) / DELZZ(K) * (TP(I,K+1) - TP(I,K)) &
                    +  SOILCD(I,K-1) / DELZZ(K-1) * (TP(I,K-1) - TP(I,K)) ) &
                    +  TP(I,K)
          END DO
          
       END DO                      
!
!
!                 Add the upper boundary condition
!                 Accounting for the 6 types of land surface in SVS2: 
!                     - snow-free exposed bare ground WTG_2
!                     - snow-free low vegetation WTG_3
!                     - snow-free high vegetation WTG_4
!                     - snow-covered bare ground and low vegetation  WTG_5
!                     - snow below high-vegetation WTG_6  
!                     - snow free ground below high vegetation WTG_7     
!
       DO I=1,N
          B2(I,1) = DELZ(1)*SOILHCAP(I,1)/DT +           &
                    + WTG(I,indx_svs2_bg) * SKINCOND_BG(I)         &   
                    + WTG(I,indx_svs2_vl)  * SKINCOND_VL(I)         &                             
                    + WTG(I,indx_svs2_gv) * SURF_ENBAL_HVEGB(I)    &
                    +  BETAA*SOILCD(I,1)/DELZZ(1)
          
          C2(I,1) = (-BETAA)*SOILCD(I,1)/DELZZ(1)

          A2(I,1) = 0.0

           
          D2(I,1) = WTG(I,indx_svs2_bg)  * SKINCOND_BG(I) * TGRST(I)  &
                  + WTG(I,indx_svs2_vl)  * SKINCOND_VL(I) * TVGLST(I)  &
                  + WTG(I,indx_svs2_gv)  * SURF_ENBAL_HVEGD(I) &
                  + WTG(I,indx_svs2_sn)  * PGRNDFLUX(I)   &
                  + WTG(I,indx_svs2_sv)  * PGRNDFLUXV(I)  & 
                  + (1.-BETAA)*(SOILCD(I,1)/DELZZ(1))*(TP(I,2)-TP(I,1)) &
                  + ( DELZ(1)*SOILHCAP(I,1)/DT )*TP(I,1)

       END DO
!
!
!
!                               add the lower boundary condition

       DO I=1,N


          IF(LBCHEAT_SVS2=='TPERM') THEN
            !  Prescribed T at bottom 

            A2(I,NL_SVS) = -BETAA * DT * SOILCD(I,NL_SVS-1) / &
                    (SOILHCAP(I,NL_SVS) * DELZ(NL_SVS)*DELZZ(NL_SVS-1))


            B2(I,NL_SVS) = 1. + ( BETAA*DT*SOILCD(I,NL_SVS) ) / ( SOILHCAP(I,NL_SVS)*DELZ(NL_SVS)*DELZZ(NL_SVS) ) + &
                     ( BETAA*DT*SOILCD(I,NL_SVS-1) ) / ( SOILHCAP(I,NL_SVS)*DELZ(NL_SVS)*DELZZ(NL_SVS-1) )

            C2(I,NL_SVS) = 0.0

            D2(I,NL_SVS) = TP(I,NL_SVS) +  ( DT/(SOILHCAP(I,NL_SVS)*DELZ(NL_SVS)) ) * ( ( SOILCD(I,NL_SVS)/DELZZ(NL_SVS) ) * TPERM(I) &
               - ( (1. - BETAA) * SOILCD(I,NL_SVS-1)/DELZZ(NL_SVS-1) ) * (TP(I,NL_SVS)-TP(I,NL_SVS-1)) &
               - ( (1. - BETAA) * SOILCD(I,NL_SVS)/DELZZ(NL_SVS) ) * TP(I,NL_SVS) )

          ELSE IF(LBCHEAT_SVS2=="0FLUX") THEN
            ! flux-zero at bottom of the soil column

            A2(I,NL_SVS) = -(BETAA) * DT * SOILCD(I,NL_SVS-1) /  &
                    (SOILHCAP(I,NL_SVS) * DELZ(NL_SVS) * DELZZ(NL_SVS-1))

            B2(I,NL_SVS) = 1. - A2(I,NL_SVS)

            C2(I,NL_SVS) = 0.0

            D2(I,NL_SVS) = TP(I,NL_SVS) - ( (1. - BETAA) * DT * SOILCD(I,NL_SVS-1) / &
                    (SOILHCAP(I,NL_SVS) * DELZ(NL_SVS)*DELZZ(NL_SVS-1)) ) * TP(I,NL_SVS) &
                   +( (1. - BETAA) * DT * SOILCD(I,NL_SVS-1) / &
                    (SOILHCAP(I,NL_SVS) * DELZ(NL_SVS)*DELZZ(NL_SVS-1)) ) * TP(I,NL_SVS-1)

         ENDIF

       END DO
!
!                             matrix inversion to solve evoluation of TP
!
!
      CALL DIFUVD2(TP, A2, B2, C2, D2, D2, N, N, NL_SVS)
!
!
!!       7.     FLUX CALCULATIONS
!               -----------------
!
!

      DO I=1,N
!                                            recalculate the qsat functions
!
        ZQSATGRT(I)  = FOQST(  TGRST(I)  ,  PS(I)   )
        ZQSATVGLT(I) = FOQST(  TVGLST(I)  ,  PS(I)  )
        ZQSATVGHT(I) = FOQST(  TVGHST(I)  ,  PS(I)  )
        ZQSATGRVT(I) = FOQST(  TP(I,1)   ,  PS(I)   )
!
      ENDDO
!                     
!
!
      DO I=1,N
!
!                                            ---------------
!                                            NET RADIATION 
!                                            ---------------
!
!                                            Net radiation over exposed bare ground
!
        RNETGR(I) = (1. - ALGR(I)) * RG(I) + EMGR(I) *&  
                 (RAT(I) - STEFAN * (TGRST(I)** 4))
!
!                                            Net radiation over low vegetation
!
        RNETVGL(I) = (1. - ALVL(I)) * RG(I) + EMISVL(I) *&  
                 (RAT(I) - STEFAN * (TVGLST(I)** 4))! VV TO BE MODIFIED

!                                            Net radiation over high vegetation
!
        RNETVGH(I) = (1. - ALVGLAI(I)) * RG(I) + EMISVH(I) *&  
                 (RAT(I) - STEFAN * (TVGHST(I)** 4))! VV TO BE MODIFIED
!
!                                            Net radiation over snow-free ground 
!                                            below high veg TO BE ADDED
        
!
!                                            AGGREGATED net radiation (including snow) 

!                                    
        RNET(I) = AG_SVS2( WTA(I,indx_svs2_bg), WTA(I,indx_svs2_vl), &
                      WTA(I,indx_svs2_vh),WTA(I,indx_svs2_sn), &
                      WTA(I,indx_svs2_sv), WTA(I,indx_svs2_gv), &
                      RNETGR(I),RNETVGL(I),RNETVGH(I), &
                      RNETSN(I),RNETSV(I),RNETGR(I) ) ! VV TO BE MODIFIED
!        
!
!                                            ---------------
!                                            SENSIBLE HEAT FLUX 
!                                            ---------------
!
!
!                                            Sensible heat flux from the ground
!
        HFLUXGR(I) = RHOA(I) * CPD * (TGRST(I) - THETAA(I)) / RESAGR(I)
!
!                                            Sensible heat flux from the low vegetation
!
        HFLUXVGL(I) = RHOA(I) * CPD * (TVGLST(I) - THETAA(I)) / RESA_VL(I)
!
!                                            Sensible heat flux from the high vegetation
!
        HFLUXVGH(I) = RHOA(I) * CPD * (TVGHST(I) - THETAA(I)) / RESA_VH(I)


!
!                                             AGGREGATED sensible heat flux (including snow)
        HFLUX(I) = AG_SVS2( WTA(I,indx_svs2_bg), WTA(I,indx_svs2_vl), &
                      WTA(I,indx_svs2_vh),WTA(I,indx_svs2_sn), &
                      WTA(I,indx_svs2_sv), WTA(I,indx_svs2_gv), &
                      HFLUXGR(I),HFLUXVGL(I),HFLUXVGH(I), &
                      HFLUXSN(I),HFLUXSV(I),HFLUXGR(I))  ! VV TO BE MODIFIED
                              
!
!
!                                            ---------------
!                                            LATENT HEAT FLUXES 
!                                            ---------------
!        ------------------
!        EXPOSED BARE GROUND
!                                            Latent heat of evaporation from
!                                            the exposed bare ground
!
        LEGNOFRAC(I) = RHOA(I) * LEFF(I) * (HRSURF(I)* ZQSATGRT(I) - HU(I)) / RESAGR(I)
        LEG(I) = WTA(I,indx_svs2_bg) * LEGNOFRAC(I)          
!
!
!                                            Water vapor flux from ground
        EGF(I) = WTA(I,indx_svs2_bg) * (HRSURF(I)* ZQSATGRT(I) - HU(I)) / RESAGR(I)

!
!                                            Latent heat of evaporation from
!                                            the ground under vegetation
!                                            TO BE MODIFIED                                          
!
!VV : Need to carefully check the closue of the water mass balance
!
        LEGVNOFRAC(I) = RORAGRV(I) * LEFF(I) * (HRSURF(I)* ZQSATGRVT(I) - HU(I))
        LEGV(I) = (VEGH(I)*(1.-PSNVH(I))+VEGL(I)*(1.-PSNG(I))) * LEGVNOFRAC(I)
!
!
!                                            Evaporation rate from ground (for hydro_svs.ftn)
!
        EG(I) = RHOA(I)*(HRSURF(I)* ZQSATGRT(I) - HU(I)) / RESAGR(I)

!
!        ------------------
!        SNOW FREE LOW VEGETATION --- SET FLUXES TO ZERO if NO SNOW FREE LOW VEGETATION
!        
          IF(VEGL(I)*(1.-PSNG(I)).ge.EPSILON_SVS) THEN
!
!            Check if if qsat> HU                  
             ZHVGL(I) = MAX(0.0 , SIGN(1.,ZQSATVGLT(I) - HU(I)))
!             
!            Transpiration rate (for hydro_svs.ftn)
             ETR_VL(I) = RHOA(I)*ZHVGL(I)*(1. - DEL_VL(I))*(ZQSATVGLT(I) -HU(I))/(RESA_VL(I) + RS(I))
!
!            Latent heat of transpiration (with fraction) 
             LETR_VL(I) = WTA(I,indx_svs2_vl)* CHLC * ETR_VL(I)
!
!            Evapotranspiration rate from low veg. (for hydro_svs.ftn) (no fraction)
             EV_VL(I) =  RHOA(I)*HV_VL(I) * (ZQSATVGLT(I) - HU(I)) / RESA_VL(I)
!
             !  EV is limited to WR/DT+RR+ETR to avoid negative WR in hydro_svs when direct evaporation exceeds rainrate
             !  When snow is present, rain falls through vegetation to snow bank... so is not considered in evaporation... This is to conserve water budget.
             IF( SNM(I).GE.CRITSNOWMASS) THEN
                ! snow over low veg is present, rain falls directly to snow
                EV_VL(I) = MIN (EV_VL(I),(WR_VL(I)/DT+ETR_VL(I)))
             ELSE 
                ! no snow present, all rain is considered evaporation
                EV_VL(I) = MIN (EV_VL(I),(WR_VL(I)/DT+ETR_VL(I)+RR(I)))
             ENDIF
        
        ELSE
           ! NO SNOW FREE LOW VEGETATION --- SET FLUXES TO ZERO
           ZHVGL(I) = 0.0
           LETR_VL(I) = 0.0
           ETR_VL(I) = 0.0
           EV_VL(I) = 0.0
        ENDIF

!       Water vapor flux from low vegetation (including fraction)
        EVLF(I) =  WTA(I,indx_svs2_vh)  * EV_VL(I)/ RHOA(I)

!       Latent heat of evaporation from low vegetation (including fraction)
        LEVL(I) = RHOA(I) * CHLC * EVLF(I)
!
!       Direct evapo. rate from low veg. (no. fraction) (for hydro_svs.ftn)
        ER_VL(I) = EV_VL(I) - ETR_VL(I)
!
!       Latent heat of direct evaporation  (including fraction)
        LER_VL(I)  = LEVL(I) - LETR_VL(I)
!
!        ------------------
!        HIGH VEGETATION --- SET FLUXES TO ZERO if NO HIGH VEGETATION
!        
          IF(VEGH(I) .ge.EPSILON_SVS) THEN
!
!            Check if if qsat> HU
             ZHVGH(I) = MAX(0.0 , SIGN(1.,ZQSATVGHT(I) - HU(I)))
!
!            Transpiration rate (for hydro_svs.ftn) (no fraction)
             ETR_VH(I) = RHOA(I)*ZHVGH(I)*(1. - DEL_VH(I))*(ZQSATVGHT(I)-HU(I))/(RESA_VH(I) + RS(I))
!
!            Latent heat of transpiration (with fraction 
             LETR_VH(I) = WTA(I,indx_svs2_vh)* CHLC * ETR_VH(I)
!             
!            Evapotranspiration rate from low veg. (for hydro_svs.ftn)
             EV_VH(I) =  RHOA(I)*HV_VH(I) * (ZQSATVGHT(I) - HU(I)) / RESA_VH(I)

             !  EV is limited to WR/DT+RR+ETR to avoid negative WR in hydro_svs when direct evaporation exceeds rainrate
             !  When snow is present, rain falls through vegetation to snow bank... so is not considered in evaporation... This is to conserve water budget.
             IF( SVM(I).GE.CRITSNOWMASS) THEN
                ! snow is present below high veg, rain falls directly to snow
                EV_VH(I) = MIN (EV_VH(I),(WR_VH(I)/DT+ETR_VH(I)))
             ELSE 
                ! no snow present, all rain is considered evaporation
                EV_VH(I) = MIN (EV_VH(I),(WR_VH(I)/DT+ETR_VH(I)+RR(I)))
             ENDIF
        
        ELSE
           ! NO HIGH VEGETATION --- SET FLUXES TO ZERO
           ZHVGL(I) = 0.0
           LETR_VH(I) = 0.0
           ETR_VH(I) = 0.0
           EV_VH(I) = 0.0
        ENDIF

!       Water vapor flux from low vegetation (including fraction)
        EVHF(I) =  WTA(I,indx_svs2_vh)  * EV_VH(I)/ RHOA(I)

!       Latent heat of evaporation from low vegetation (including fraction)
        LEVH(I) = RHOA(I) * CHLC * EVHF(I)
!
!       Direct evapo. rate from low veg. (no. fraction) (for hydro_svs.ftn)
        ER_VH(I) = EV_VH(I) - ETR_VH(I)
!
!       Latent heat of direct evaporation  (including fraction)
        LER_VH(I)  = LEVH(I) - LETR_VH(I)
!        
!        ------------------
!        SNOW ABOVE BARE GROUND AND LOW VEGETATION
!
!                           Calculate latent heat snow weighted
!                           by grid-cell snow-coverage fraction
!       
        LES(I)  =  WTA(I,indx_svs2_sn) *  (LESLNOFRAC(I) + LESNOFRAC(I))
        ESF(I)  =  WTA(I,indx_svs2_sn) *  ESNOFRAC(I)
!
!        ------------------
!        SNOW BELOW HIGH VEGETATION      
!                           Same for snow-under-vegetation
!
        LESV(I) =  WTA(I,indx_svs2_sv)   *  (LESLVNOFRAC(I) + LESVNOFRAC(I))
        ESVF(I) =  WTA(I,indx_svs2_sv)   *  ESVNOFRAC(I)
!
!                           Total latent heat of evaporation
!                           (Including snow contribution) VV NEED to INCLUDE GROUND BELOW HIGH VEG
!
        LE(I) = LEG(I) + LEVL(I) + LEVH(I) + LES(I) + LESV(I) 
!
!                                            Total water vapor flux
!                                            (Including snow contribution)
!       VV NEED to INCLUDE GROUND BELOW HIGH VEG
        
        EFLUX(I) = EGF(I) + EVLF(I) +  EVHF(I) + ESF(I) + ESVF(I)
!
!                                            Heat flux into the ground
!
        GFLUX(I) = RNET(I) - HFLUX(I) - LE(I)
!
      ENDDO
!
!
!
!
!*       8.     NEW "land-tile"-AVERAGED QUANTITIES  (FRV, ZQS, BM, FQ)
!               -----------------------------------------------  
!
!
!
      DO I=1,N
!            Re-calculate snow saturation humidity
!
        ZQSATSNO(I)=FOQST( TSNS(I,1), PS(I) )
        ZQSATSNV(I)=FOQST( TSVS(I,1), PS(I) )
        !
      ENDDO

      IF ( .NOT. use_eff_surf_tq  ) THEN
        !  Area-average weighted mean calculation for land sfc temperature and humidity

        !
        DO i=1,n
           ! Calculate land-tile-averaged specific humidity
!
           ZQS(I) =     WTA(I,indx_svs2_bg)      *    HRSURF(I)        * ZQSATGRT(I) & 
                      + WTA(I,indx_svs2_sn)                            * ZQSATSNO(I) &  
                      + WTA(I,indx_svs2_sv)                            * ZQSATSNV(I) & 
                      + WTA(I,indx_svs2_gv)      *    HRSURF(I)        * ZQSATGRT(I) &  ! VV TO BE MODIFIED
                      + WTA(I,indx_svs2_vl)      *     HV_VL(I)        *ZQSATVGLT(I) &
                      + WTA(I,indx_svs2_vl)      *(1.-HV_VL(I))        *HU(I)        &
                      + WTA(I,indx_svs2_vh)      *     HV_VH(I)        *ZQSATVGHT(I) &
                      + WTA(I,indx_svs2_vh)      *(1.-HV_VH(I))       * HU(I)
!              
!
!              Calculate land-tile-averaged surface temperature
!              i.e., aggregate skin temperatures of diff. surfaces 
!
         TSA(I) = AG_SVS2( WTA(I,indx_svs2_bg), WTA(I,indx_svs2_vl), &
                      WTA(I,indx_svs2_vh),WTA(I,indx_svs2_sn), &
                      WTA(I,indx_svs2_sv), WTA(I,indx_svs2_gv), &
                      TGRST(I),TVGLST(I),TVGHST(I), &
                      TSNS(I,1),TSVS(I,1),TGRST(I) ) ! VV TO BE MODIFIED

        ENDDO
              
     else
        ! Consider aerodynamic resistance in calc. of effective ("land-tile-averaged") land sfc 
        ! temperature and humidity (see svs_configs.ftn90)

!
!        Calculate effective areodynamic resistance  ! VV TO BE MODIFIED
!         
        DO i=1,n
           RESAEF(I) = 1. / ( WTA(I,indx_svs2_bg)/RESAGR(I) +WTA(I,indx_svs2_vl)/RESA_VL(I) + &
                              WTA(I,indx_svs2_vh)/RESA_VH(I)+WTA(I,indx_svs2_gv)/RESAGR(I) + &
                              WTA(I,indx_svs2_sn)/RESASA(I) + WTA(I,indx_svs2_sv)/RESASV(I) )
!
!          Calculate effective land sfc specific humdity
!         
           ZQS(I) =     WTA(I,indx_svs2_bg)*    HRSURF(I)  * ZQSATGRT(I)/RESAGR(I) & 
                      + WTA(I,indx_svs2_sn)                * ZQSATSNO(I)/RESASA(I) &  
                      + WTA(I,indx_svs2_sv)                * ZQSATSNV(I)/RESASV(I) & 
                      + WTA(I,indx_svs2_gv)*    HRSURF(I)  * ZQSATGRT(I)/RESAGR(I) &  ! VV TO BE MODIFIED
                      +(WTA(I,indx_svs2_vl)*     HV_VL(I) * ZQSATVGLT(I)          &
                      + WTA(I,indx_svs2_vl)*(1.-HV_VL(I))  * HU(I))/RESA_VL(I)     &
                      +(WTA(I,indx_svs2_vh)*     HV_VH(I)  * ZQSATVGHT(I)         &
                      + WTA(I,indx_svs2_vh)*(1.-HV_VH(I))*HU(I))/RESA_VH(I)

!
!          Calculate effective land sfc temperature
! 
         TSA(I) = RESAEF(I) * AG_SVS2( WTA(I,indx_svs2_bg), WTA(I,indx_svs2_vl), &
                      WTA(I,indx_svs2_vh),WTA(I,indx_svs2_sn),  &
                      WTA(I,indx_svs2_sv), WTA(I,indx_svs2_gv), &
                      TGRST(I)/RESAGR(I),TVGLST(I)/RESA_VL(I),TVGHST(I)/RESA_VH(I), &
                      TSNS(I,1)/RESASA(I),TSVS(I,1)/RESASV(I),TGRST(I)/RESAGR(I) ) ! VV TO BE MODIFIED

        ENDDO
     endif


!    Calculated averaged radiative Temperature

     DO i=1,n

         TRAD(I) = AG_SVS2( WTA(I,indx_svs2_bg), WTA(I,indx_svs2_vl), &
                      WTA(I,indx_svs2_vh),WTA(I,indx_svs2_sn),  &
                      WTA(I,indx_svs2_sv), WTA(I,indx_svs2_gv), &
                      TGRST(I)**4,TVGLST(I)**4,TVGHST(I)**4, &
                      TSNS(I,1)**4,TSVS(I,1)**4,TGRST(I)**4 ) ! VV TO BE MODIFIED
        
        TRAD(I)=TRAD(I)**(1./4.)
     ENDDO
!
!
!                                             Calculate surface layer transfer coefficients
!                                             and fluxes using ZQS and TSA. 
!                                             Here want FRV (average friction velocity for 
!                                             grid-cell, used in diasurf2.ftn), CMU and 
!                                             CTU(Turbulent transfer coeff. for thermodynamics) 
!
!

! DO we still need FRV (UE), ILMO, HST ???? !!!! as output
!
!
     if( svs_dynamic_z0h ) then
        zopt=9
        i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
                         TSA, ZQS, Z0, Z0H, LAT, FCOR, &
                         optz0=zopt, z0mloc=Z0LOC, &
                         L_min=sl_Lmin_soil,spdlim=vmod_lmin, &
                         coefm=CMU, coeft=CTU, &
                         flux_t=FTEMP, flux_q=FVAP, &
                         ilmo=ILMO, ue=FRV, h=HST, &
                         z0t_optz0=Z0H )


        if (i /= SL_OK) then
           call physeterror('ebudget_svs', 'error returned by sl_sfclayer()')
           return
        endif
     else
      
        i = sl_sfclayer( THETAA, HU, VMOD, VDIR, ZUSL, ZTSL, &
                         TSA, ZQS, Z0, Z0H, LAT, FCOR, &
                         L_min=sl_Lmin_soil,spdlim=vmod_lmin, &
                         coefm=CMU, coeft=CTU, &
                         flux_t=FTEMP, flux_q=FVAP, &
                         ilmo=ILMO, ue=FRV, h=HST )


        if (i /= SL_OK) then
           call physeterror('ebudget_svs', 'error returned by sl_sfclayer()')
           return
        endif
        
   endif


      DO I=1,N
!
!
!        TERMS SAME FOR IMPLICIT AND EXPLICIT FORMULATION
        
         CM(I)  = CMU(I) / FRV(I)

         if ( sl_Lmin_soil > 0.) then
            ! use wind module consistent with imposed Monin-Obukhov length
            bm(i) = vmod_lmin(i) * (cm(i)**2)  

            fq(i) = rhoa(i) * bm(i) * vmod_lmin(i)

         else
            bm(i) = vmod(i) * (cm(i)**2)  

            fq(i) = rhoa(i) * bm(i) * vmod(i)
         endif



      ENDDO

      
      ! -------  EXPLICIT FORMULATION
      IF (.NOT.IMPFLX) THEN
         DO I=1,N   
 
         !   inhomogeneous boundary term in the diffusion equation for Theta
         !   alfat = ctu * (ta -t_surf) = -ftemp
            ALFAT(I)   =  -FTEMP(I)
!            inhomogeneous boundary term in the diffusion equation for Q
         !   alfaq = ctu * (qa -q_surf) = -fvap
            ALFAQ(I)   =  -FVAP(I)

         ! homogeneous boundary condition term in the
         ! diffusion equation for T and Q          
            BT(I) = 0.0

         ENDDO
      ENDIF


      ! -------  IMPLICIT FORMULATION
      IF (IMPFLX) THEN
         DO I=1,N   
      
         !  inhomogeneous boundary term in the diffusion equation for Theta
            ALFAT(I)   =  -CTU(I)  * TSA(I)

         ! inhomogeneous boundary term in the diffusion equation for Q
            ALFAQ(I)   =  -CTU(I)  * ZQS(I)

         ! homogeneous boundary condition term in the
         ! diffusion equation for T and Q          
            BT(I) = CTU(I)

         ENDDO
      ENDIF

!
!*       9.     UPDATE TEMPERATURE VARIABLES
!              -----------------------------
!
      DO I=1,N
        TGRS(I)   = TGRST(I)
        TGRD(I)   = TGRDT(I)
        TVGLS(I)   = TVGLST(I)
        TVGHS(I)   = TVGHST(I)
        TVGHD(I)   = TVGHDT(I)      

      ENDDO
!
!
      RETURN
      END
