
% OVEARCHING DATA STRUCTURE ----------------------------------------------

directories.TRANSECTS   = [ROOT filesep '01_DATA' filesep 'Transects']; 
directories.DEM         = [ROOT filesep '01_DATA' filesep 'DEM']; 

% OVEARCHING RESULTS STRUCTTURE -------------------------------------------

if exist('NAME_GRD_transects','var') 
    directories.profiles            = [ROOT filesep '03_Results\',NAME_GRD_transects,'\interp_profiles_wave_model' ]; 
    directories.dir_swan_interp     = [ROOT filesep '03_Results\',NAME_GRD_transects,'\interp_swan'                   ]; 
    directories.timeseries          = [ROOT filesep '03_Results\',NAME_GRD_transects,'\time_series_reconstructed',      ]; 
    
    directories.extreme_analysis= [ROOT filesep '03_Results\',NAME_GRD_transects,'\extreme_value_analysis',           ]; 
    directories.runup           = [ROOT filesep '03_Results\',NAME_GRD_transects,'\runup_series',                   ];     
    directories.slope           = [ROOT filesep '03_Results\',NAME_GRD_transects,'\slope_values',                   ]; 
    directories.minor           = [ROOT filesep '03_Results\',NAME_GRD_transects,'\minor_transects',                   ]; 

end

