!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------

subroutine inisoili_svs2(ni, trnch)
   use sfcbus_mod
   use svs_configs
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>

   integer ni, trnch

   !@Author  Maria Abrahamowicz, Stephane Belair , Vincent Fortin (20xx)
   !@Object  Compute soil properties for given soil texture. Compute these properties on 
   !         native levels of database providing soil texture, and then map properties 
   !         unto SVS levels.
   !@Arguments
   !             - Input -
   ! NI          longueur d'une tranche horizontale

   integer :: i, k, kk
   REAL b, usb, fb, crit1_wfcint, crit2_wfcint, soillayer_depth, beta_soc
   REAL :: depth_sapric = 1., & ! Depth (m) where the soil properties reach sapric values (Decharm et al. 2016)
           depth_fibric = 0.01  ! Depth (m) where the soil properties start to depart from fibric values (Decharm et al. 2016)
   
   ! "geo" variables are on the levels of the geophysical soil texture datbase
   REAL, dimension(ni,nl_stp) :: wsat_geo, wwilt_geo, wfc_geo, b_geo, psisat_geo, &
           ksat_geo, wfcint_geo, fb_geo ,quartz_geo,rhosoil_geo,conddry_geo,condsld_geo, soilhcapz_dry
   ! 100% soil organic content variables 
   REAL, dimension(ni,nl_stp) :: wsat_soc, wwilt_soc, wfc_soc, b_soc, psisat_soc, &
           ksat_soc, wfcint_soc, fb_soc ,conddry_soc,condsld_soc, hcap_soc
   real, pointer, dimension(:) :: zcgsat, zgrkef, zdraindens, zslop

   ! variables on the levels of SVS
   real, pointer, dimension(:,:) :: zbcoef, zclay, zfbcof, zksat, zpsisat, zsand, zwfc, zwfcint, zwsat, zwwilt 
   real, pointer, dimension(:,:) :: zconddry, zcondsld, zquartz, zrhosoil, zsoilhcapz_dry, zsoc


  
#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zcgsat, cgsat)
   MKPTR1D(zdraindens, draindens)
   MKPTR1D(zgrkef, grkef)
   MKPTR1D(zslop, slop)


   MKPTR2D(zbcoef, bcoef)
   MKPTR2D(zclay, clay)
   MKPTR2D(zconddry, conddry)
   MKPTR2D(zcondsld, condsld)  
   MKPTR2D(zfbcof, fbcof)
   MKPTR2D(zksat, ksat)
   MKPTR2D(zpsisat , psisat)
   MKPTR2D(zquartz, quartz)  
   MKPTR2D(zrhosoil, rhosoil)  
   MKPTR2D(zsand, sand)
   MKPTR2D(zwfc, wfc)
   MKPTR2D(zwfcint, wfcint)
   MKPTR2D(zwsat, wsat)
   MKPTR2D(zwwilt , wwilt)
   MKPTR2D(zsoilhcapz_dry , soilhcapz_dry)
   MKPTR2D(zsoc , soc)
   !call subroutine to compute layer thicknesses
   call layer_thickness()

   ! calculate soil parameters on native GEO layers, and then map them unto model layers. 
   ! calculate weights to be used in phybusinit.... because here... we are
   ! re-doing the calculation for each row of domain...
   ! but the weights are the same !

   ! Compute averaged values for sand and clay for each layer when SOILGRID is used
   if ( soiltext == "SOILGRIDS" ) then   
      do i=1,ni
         do k=1,(nl_stp-1)
            zsand(i,k) = (zsand(i,k)+zsand(i,k+1))/2.
            zclay(i,k) = (zclay(i,k)+zclay(i,k+1))/2.
         enddo
      enddo
   endif

   !     Computer soil properties for GEO layers
   do i=1,ni
      do k=1,nl_stp
         wsat_geo  (i,k)  =  -0.00126   * zsand(i,k) + 0.489
         wwilt_geo (i,k)  =  37.1342e-3 * sqrt(max(1.,zclay(i,k)))
         wfc_geo   (i,k)  =  89.0467e-3 * max(1.,zclay(i,k))**0.3496
         psisat_geo(i,k)  =  0.01 * ( 10.0**(-0.0131 * zsand(i,k) + 1.88) )
         ksat_geo  (i,k)  =  ( 10.0**(0.0153 * zsand(i,k) - 0.884) ) * 7.0556E-6

         b                 =  0.137 * zclay(i,k)  + 3.501
         b_geo     (i,k)  =  b
         usb               =  1./b
         fb                =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
         fb_geo(i,k)      =  fb
         ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
         ! Ensure that wc at fc stays between wilting point and saturation

         crit1_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*(wsat_geo(i,k)/wwilt_geo(i,k)*fb)**b
         crit2_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*fb**b

         if (abs(zslop(i)).gt.crit1_wfcint) then
            wfcint_geo(i,k) = wwilt_geo(i,k)        
         elseif (abs(zslop(i)).lt.crit2_wfcint) then
            wfcint_geo(i,k) = wsat_geo(i,k) 
         elseif (zslop(i).ne.0.0) then
            wfcint_geo(i,k) = wsat_geo(i,k) * fb * &
                 ( psisat_geo(i,k)/ABS(zslop(i)) *2. * zdraindens(i) )**usb
         else
            wfcint_geo(i,k) = wfc_geo(i,k)
         endif
         
         ! Compute soil thermal properties for heat diffusion

!       Quartz content (ref : NL95 & PL98)): 
        quartz_geo(i,k)  = 0.038 + 0.0095*zsand(i,k)

!       Soil dry density (PL98):
        rhosoil_geo(i,k) = (1.0-wsat_geo(i,k))*2700.

!       Soil solid conductivity:
        if (quartz_geo(i,k).gt.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (2.0**(1.0-quartz_geo(i,k))) 
        endif
        if (quartz_geo(i,k).le.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (3.0**(1.0-quartz_geo(i,k)))
        endif

!       Soil dry conductivity:
        conddry_geo(i,k) = (0.135*rhosoil_geo(i,k) + 64.7) / &
                        (2700. - 0.947*rhosoil_geo(i,k))   

!       Soil heat capacity
        if ((zsand(i,k)+zclay(i,k)) .gt. epsilon_svs) then
            ! From Lawrence and Slater 2008
            soilhcapz_dry(i,k) = (2.128*zsand(i,k) + 2.385*zclay(i,k))/(zsand(i,k)+zclay(i,k))*1.E6
        else    
            ! Use previous formulation in SVS - constant value
            soilhcapz_dry(i,k) = 2700. * 733.
        endif


      enddo
   enddo
   ! "Map" GEO soil properties unto model soil layers
   Do i = 1 , ni
      Do k = 1, nl_svs
         do kk = 1 , nl_stp
            
            zwsat  (i,k)  = zwsat  (i,k) + wsat_geo  (i,kk)  * weights( k , kk)
            zwwilt (i,k)  = zwwilt (i,k) + wwilt_geo (i,kk)  * weights( k , kk)
            
            zwfc   (i,k)  = zwfc   (i,k) + wfc_geo   (i,kk)  * weights( k , kk)
            zbcoef (i,k)  = zbcoef (i,k) + b_geo     (i,kk)  * weights( k , kk)
            zfbcof (i,k)  = zfbcof (i,k) + fb_geo    (i,kk)  * weights( k , kk)
            zpsisat(i,k)  = zpsisat(i,k) + psisat_geo(i,kk)  * weights( k , kk)
            zksat  (i,k)  = zksat  (i,k) + ksat_geo  (i,kk)  * weights( k , kk)
            zwfcint(i,k)  = zwfcint(i,k) + wfcint_geo(i,kk)  * weights( k , kk)
            zconddry  (i,k)  = zconddry  (i,k) + conddry_geo  (i,kk)  * weights( k , kk)
            zcondsld  (i,k)  = zcondsld  (i,k) + condsld_geo  (i,kk)  * weights( k , kk)
            zquartz   (i,k)  = zquartz   (i,k) + quartz_geo   (i,kk)  * weights( k , kk)
            zrhosoil  (i,k)  = zrhosoil  (i,k) + rhosoil_geo  (i,kk)  * weights( k , kk)
            zsoilhcapz_dry  (i,k)  = zsoilhcapz_dry  (i,k) + soilhcapz_dry  (i,kk)  * weights( k , kk)

         enddo
      enddo
      ! compute thermal coeff. 
      ! for 1st model layer only --- here simply use 1st GEO soil texture !!! Do not map !
      zcgsat (i)  = ( -1.557e-2 * zsand(i,1) &
           -  1.441e-2 * zclay(i,1) + 4.7021 ) * 1.E-6 
      
      ! Compute effective parameter for watdrain
      zgrkef(i)   = 2.* zdraindens(i) * zslop(i)

   enddo
   !     Computer 100% soil organic content properties for typical peat soil profile (Decharme et al. 2016)
   do i=1,ni
      do k=1,nl_stp

         if (k .EQ. 1) then
            soillayer_depth = delz(k) * 0.5
         else
            soillayer_depth = soillayer_depth + (delz(k) + delz(k-1))*0.5
         endif

         beta_soc = log(0.845/0.93)/log(depth_sapric/depth_fibric) ! See Table 1 in Decharme et al. 2016
         wsat_soc  (i,k)  =  0.93 * (soillayer_depth/depth_fibric)**(beta_soc)
         wsat_soc     (i,k)  =  min(max( wsat_soc(i,k), 0.845), 0.93) ! Bound the value to min and max

         beta_soc = log(0.222/0.073)/log(depth_sapric/depth_fibric)
         wwilt_soc (i,k)  =  0.073 * (soillayer_depth/depth_fibric)**(beta_soc)
         wwilt_soc     (i,k)  =  min(max( wwilt_soc(i,k), 0.073), 0.222)

         beta_soc = log(0.719/0.369)/log(depth_sapric/depth_fibric)
         wfc_soc   (i,k)  =  0.369 * (soillayer_depth/depth_fibric)**(beta_soc)
         wfc_soc     (i,k)  =  min(max( wfc_soc(i,k), 0.369), 0.719)

         beta_soc = log(0.0101/0.0103)/log(depth_sapric/depth_fibric)
         psisat_soc(i,k)  =  0.0103 * (soillayer_depth/depth_fibric)**(beta_soc)
         psisat_soc     (i,k)  =  min(max( psisat_soc(i,k), 0.0101), 0.0103)

         beta_soc = log((1.E-7)/(2.8E-4))/log(depth_sapric/depth_fibric)
         ksat_soc  (i,k)  = 2.8E-4 * (soillayer_depth/depth_fibric)**(beta_soc)
         ksat_soc     (i,k)  =  min(max( ksat_soc(i,k), 1.E-7), 2.8E-4)

         beta_soc = log(12./2.7)/log(depth_sapric/depth_fibric)
         b_soc     (i,k)  =  2.7 * (soillayer_depth/depth_fibric)**(beta_soc)
         b_soc     (i,k)  =  min(max( b_soc(i,k), 2.7), 12.)

         hcap_soc(i,k) = 2.5E6

         condsld_soc(i,k) =  0.25 

         conddry_soc(i,k) = 0.05   

      enddo
   enddo

   ! model soil layers incluting soil organic content
   Do i = 1 , ni
      Do k = 1, nl_svs


         ! Arithmetic mean

         zwsat  (i,k)  = (1.-zsoc(i,k)/100.) * zwsat  (i,k) + zsoc(i,k)/100. * wsat_soc  (i,k)
         zwwilt (i,k)  = (1.-zsoc(i,k)/100.) * zwwilt  (i,k) + zsoc(i,k)/100. * wwilt_soc  (i,k)
         zwfc   (i,k)  = (1.-zsoc(i,k)/100.) * zwfc  (i,k) + zsoc(i,k)/100. * wfc_soc  (i,k)
         zpsisat(i,k)  = (1.-zsoc(i,k)/100.) * zpsisat  (i,k) + zsoc(i,k)/100. * psisat_soc  (i,k)
         zsoilhcapz_dry (i,k) = (1.-zsoc(i,k)/100.) * zsoilhcapz_dry  (i,k) + zsoc(i,k)/100. * hcap_soc  (i,k)

         ! Geometric mean
         zksat  (i,k)  = zksat(i,k)**(1.-zsoc(i,k)/100.) * ksat_soc(i,k)**(zsoc(i,k)/100.)
         zconddry  (i,k)  = zconddry(i,k)**(1.-zsoc(i,k)/100.) * conddry_soc  (i,k)**(zsoc(i,k)/100.)
         zcondsld  (i,k)  = zcondsld(i,k)**(1.-zsoc(i,k)/100.) * condsld_soc(i,k)**(zsoc(i,k)/100.)

         ! bcoeff and fbcoeff, Arithmetic mean
         b                 =  (1.-zsoc(i,k)/100.) * zbcoef  (i,k) + zsoc(i,k)/100. * b_soc  (i,k)
         zbcoef     (i,k)  =  b
         usb               =  1./b
         fb                =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
         zfbcof(i,k)      =  fb

         ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
         ! Ensure that wc at fc stays between wilting point and saturation

         crit1_wfcint   = 2.*zdraindens(i)*zpsisat(i,k)*(zwsat(i,k)/zwwilt(i,k)*fb)**b
         crit2_wfcint   = 2.*zdraindens(i)*zpsisat(i,k)*fb**b

         if (abs(zslop(i)).gt.crit1_wfcint) then
            zwfcint(i,k) = zwwilt(i,k)        
         elseif (abs(zslop(i)).lt.crit2_wfcint) then
            zwfcint(i,k) = zwsat(i,k) 
         elseif (zslop(i).ne.0.0) then
            zwfcint(i,k) = zwsat(i,k) * fb * &
                 ( zpsisat(i,k)/ABS(zslop(i)) *2. * zdraindens(i) )**usb
         else
            zwfcint(i,k) = zwfc(i,k)
         endif

      enddo

   enddo

   return
 end subroutine inisoili_svs2
