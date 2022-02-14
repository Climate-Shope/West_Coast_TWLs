West Coast Total Water Levels Scripts used to develop the data for the publication 
"Characterizing Storm-Induced Coastal Change Hazards Along the United States West Coast"
by James B. Shope, Li H. Erikson, Patrick L. Barnard, Curt D. Storlazzi, Katherine Serafin,
Kara Doran, Hilary Stockdon, Borja Reguero, Fernando Mendez, Sonia Castanedo, Alba Cid, 
Laura Cagigal, and Peter Ruggiero. 2022

This readme provides guidance on the included subfolders, their intended use, and the
order in which to run the scripts to regionally project total water levels (TWLs). 

The intent of these scripts is to  provide an overview of the methodology to project
extreme total water levels along a coastline. Limited data are provided, and it is assumed
the user will be able to supplement environmental sensitivity index, topography, tidal, 
and nearshore wave data that is not hosted here. 

---------------------------------------------------------------------------------------------

Included folders:

1. _Template: This is the primary folder containing scripts to calculate the TWLs. The scripts
are presented in a modular fashion to allow offloading of processing on multiple machines. The provided
example was used to calculate hourly TWLs and extreme TWLs in Douglas County Oregon. 


2. _STEP_Functions: These are a suite of utility functions used to calculate specific wave parameters,
limited geographic processing, and other miscellaneous functions to aid in data preprocessing. It is not
recommended that the user modify the functions within this folder.

3. Extreme_Value_Scripts: A suite of functions to aid in the calculation of extreme value fits to determine 
return level events. It is not recommended that the user modify these scripts. 

4. TWL_Slope_Temp: Contains functions to calculate beach slopes regionally, correct extracted parameters
of cross shore elevation profiles, and specific functions called in _Template scripts to determine the 
appropriate runup methodology, calculate the runup, and extract relevant cross shore elevation profile 
parameters (such as a dune toe or crest). 


------------------------------------------------------------------------------------------------

Calculating TWLs:

To calculate the TWLs for the provided example, the user will need to run the codes in the _Template and
TWL_Slope_Temp/Beach_Slope_Scripts folders in the following order. Note, if you have all the requisite data
to call all of the scripts from start to finish in a batch file, it is not recommended to do so. Some of these
processes can take days to weeks to compute on a singular machine. Finally, some optional codes are included that
are located in TWL_Slope_Temp\Correct_Profiles to provide some morphology specific fixes to common problems with
automatically calculated runup methods, toe locaitons, and TWL calculations associated with rock bench profiles. 
Depending on the region in question, these optional scripts may not be needed by the user. In the provided example,
most locations did not require these fixes. 



Step 1 cast transects  
-This assumes a shoreline shapefile that can be read by the matlab mapping toolbox. A shoreline for Douglas County in  
X and Y coordinates (UTM10), a shapefile of the shore, and coordinates of transect locations is provided. 

Order:  
M0_Smooth_Shoreline  
M1a_Cast_Transects_From_DEM  
M1b_Order_line_segements  
M1c_compare_remove_points  
M1d_Generate_Updated_Transects  

Step 2 interpolate elevation profiles  
-This step takes rasters of wave data and topography and extracts the relevant elevation and wave parameters along
the computed transect locations. Additonally, this prepares location point values to be used with NOAA's vDatum
tool to compute Mean High Water (MHW) and Mean Sea Level (MSL) approximations for each location.  

Order:  
M2a_Calculate_Offshore_Points  
M2b_interp_SWAN  
M2c_Interp_Profiles_a  
M2c_Interp_Profiles_b  
M2c_Interp_Profiles_c  
M2d_Interp_Minor_Profiles_a  
M2d_Interp_Minor_Profiles_b  
M2e_Prep_for_vdatum  
M2e_Calculate_MHW  
M2e_Calculate_MSL  
M2f_Profile_Island_Filter  
M2g_Minor_Porfile_Island_Filter  

Step 3 generate downscaled wave time series  
-Computes a timeseries of wave conditions from a representative sample following Camus and others 2011. These are
split into A and B processes to allow for offloading of the time consuming process to other machines. In past
experience, this process often takes several days on one machine. 

Order:  
M3a_reconstruct_timeseries_nearshore_A  
M3a_reconstruct_timeseries_nearshore_B  

Step 4 Archive shoreline environment  
-Assigns a NOAA Environmental Sensitivity Index shoreline type to each transect based on the intersection
of NOAA ESI shapefiles and transect XY coordinates. 

Order:  
M4a_Form_ESI_Intersect  
M4b_Get_Unique_ESI  

Step 5 Extract relevant morphology elevations and calculate TWLs  
-The profile analysis and update will extract the necessary components of cross shore geometry
to empirically calculate TWLs at each time step and provide some corrections. Please note that
this process is very regionally specific, and the user may need to provide different thresholds
for profile simplification and key elevations. The troubleshoot functions correct for profiles 
with abnormally large TWL values and profiles where NaN values have occurred. 

Order:  
M5a_profile_morphology_analysis  
M5b2_Update_2A_and_TAW_Morphology  
**   
Optional fixes to profile morphology analysis if determined necessary for the region  
TWL_Slope_Temp\Correct_Profiles\Fix_Toe_Onshore  
TWL_Slope_Temp\Correct_Profiles\Fix_Runup_Method  
**  
M5c_TWL_Calculation  
**  
Optional fixes to rocky bench profile TWL calculations if determined necessary for the region  
TWL_Slope_Temp\Correct_Profiles\Fix_2A_Profiles  
**  
M5d_TWL_Troubleshoot  
M5e_TWL_Troubleshoot_NaNs  


Step 6 Extreme value analysis  
-Calculates extreme, return interval events based on an annual GEV or a peaks over threshold
GPD fit method by determining which process better fits the data. 

Order:  
M6a_extreme_value_analysis_year_max_annual_Daily  


Step 7 Calculate first overtopping regime  
-Compared calculated TWLS to onshore geometry to determine the probability of the profile experiencing
collision, overtopping, or inundation for a given storm event. 

Order: (These are iterative corrections to Overtopping Regimes and Runup values)  
M8a_Overtopping_Regime_Iteration1 <--------- First attempt at determining Overtopping Regime  


Step 8 Reassessment of TWL calculations  
-Provides a logical check on TWL and overtopping errors and recalculates profiles if an error is identified. 

Order:  
M7a_Tweak_Runup_Iteration1 <--------- First pass to correct Runup errors  
M8b_Overtopping_Regime_Iteration2 <--------- Second pass to correct Overtopping Regime errors  


Step 9 Final reassessment of TWL calculations  
-Provides a second logical check on TWL and overtopping errors and recalculates profiles if an error is identified.  

Order:  
M7b_Tweak_Runup_Iteration2 <--------- Second pass to correct Runup errors  
M8c_Overtopping_Regime_Iteration3 <--------- Last pass to correct Overtopping Regime errors from first 2  


Step 10 Generate Output Text files  

Order:  
M9_Generate_RP_txt_Files  


Step 11 Correct Small Errors in Return Period Text Files  

Order:  
M10a_Update_Probabilities_for_overwash <--------- Corrects where Overwash probabilities are erroneous  
M10b_Outlier_Analysis <--------- Performs a simple outlier analysis to remove transects with excessively high TWL estimates  
____________________________________________________________  

Other files in Template Folder include those generated by
the above scripts for temporary use and can be ignored when
running the above. The one exception is the "Douglas_TS_Input"
file that is for use with M3a_reconstruct_timeseries_nearshore_B.
This file can be offloaded to another computer for wave 
time series computation.   

----------------------------------------------------------------------------------------------------------------------------  

Beach Slopes  

The above script sequence assumes that the user already has a series of beach slopes for use. If not, the following
will calculate the beach slopes from the elevation profiles. In sequence, these would be computed between steps 4 and 5.  

Step 1 Navigate to ...TWL_Slope_Temp/Beach_Slope_Scripts  
-These scripts are hosted and calculated externally from the _Template scripts  



Step 2  
-Calculate initial beach slopes and generate a regional moving average slope for each profile. Note that the created folders
use the name Stockdon for the slopes to indicate that they are intended for use with the Stockdon and others (2006) runup 
methodology. Gets a baseline. 

Order:  
Get_BeachSlope  
Get_Major_Slope_Averages  



Step 3  
-Redo with an updated iteration and provide a methodology to fill in beach slopes where they cannot be empirically
calculated from the elevation profiles, such as along sea cliffs.   

Order:  
Get_BeachSlope_Gen_Regional  
Get_Major_Slope_Averages_Limit_Slopes  
Fill_BeachSlope_Limit_Slope  

