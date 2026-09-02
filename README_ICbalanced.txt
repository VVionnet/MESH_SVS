!Steps to estimate balanced liquid-ice partitionning to initialize a SVS2 under subfreezing conditions

files:
 - soil_profile.txt (soil profile inpu file)
 - balance_IC_subfreezing.py (script)
 *Note that the soil_profile.txt should be located in the same repo balance_IC_subfreezing.py 

1. Change the soil_profile.txt file with your soil_profile.
    - The header of this file is the following:
        - Layer:          Layer index   
        - Thickness_m:    Layer thickness (m)
        - Depth_m:        Depth of the layer's lower boundary (m)
        - Sand_%:         Layer's sand percentage (%)
        - Clay_%:         Layer's clay percentage (%)
        - TSOIL_K         Layer's temperature (K)
        - Tot_wat_cont    Layer's Total water content (-)
    
    - This file expects numeric values separated by spacings

2. Run the balance_IC_subfreezing.py script.
    - This script computes the balanced liquid-ice partitionning of a soil layer at subfreezing conditions based on the Clausius-Clapeyron freezing-point depression and the Clapp-Hornberger/Brooks-Corey retention curve, as integrated in SVS2.

    - The script outputs the file soil_profile_balanced.txt, a transposed soil_profile.txt matrix with two new rows at the end (wsoil and isoil).

3. Copy the rows wsoil and isoil into the MESH_parameters.txt file of your MESH-SVS experiment.
