import pandas as pd
import math

#Load the soil_profile.txt file
soil_profile = pd.read_csv('soil_profile.txt', sep=r'\s+',index_col = 0)

#Constants
CHLF = 334000    #J K-1     #Latent heat of fusion for water
TRPL = 273.16   #K          #Triple point of water
GRAV = 9.80616  #m s-2      #Gravitational constant

wsol_profile = []
isol_profile = []
for lay in range(len(soil_profile)):
    
    SAND = soil_profile['Sand_%'].iloc[lay]
    CLAY = soil_profile['Clay_%'].iloc[lay]
    TSOIL = soil_profile['TSOIL_K'].iloc[lay]
    TWC = soil_profile['Tot_wat_cont'].iloc[lay]

    #Soil texture-based parameters
    WSAT = -0.00126*SAND+0.489
    PSISAT = -0.01*(10**(-0.0131*SAND+1.88))
    b_coef = 0.137*CLAY+3.501

    #Make sure that the total water content does not exceed water content at saturation
    if (TWC > WSAT):
        TWC = WSAT

    #Matric potential at TSOIL
    PSIMAX = min(PSISAT, CHLF*(TSOIL-TRPL)/(GRAV*TSOIL))

    #Max liquid water content at TSOIL
    WORK = PSIMAX/PSISAT
    WORKLOG = math.log(WORK)/b_coef
    WSOLMAX = WSAT*math.exp(-WORKLOG)

    #Liquid and Ice content
    WSOL = float(min(TWC,WSOLMAX))
    ISOL = float(TWC-WSOL)

    wsol_profile.append(WSOL)
    isol_profile.append(ISOL)

wsol_isol_profile = pd.concat([pd.DataFrame(wsol_profile),pd.DataFrame(isol_profile)],axis= 1)
wsol_isol_profile.index = soil_profile.index
     
soil_profile = pd.concat([soil_profile, wsol_isol_profile], axis = 1)
soil_profile.columns.values[-2:] = ['wsoil', 'isoil']
soil_profile = soil_profile.T

soil_profile.to_csv('soil_profile_balanced.txt', sep=' ', index=True)



