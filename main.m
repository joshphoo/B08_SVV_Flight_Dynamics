% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
% |                         Kut-Simulatie main.m                          | 
% |                                 2019                                  | 
% |                             Brakke Willy                              | 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

% Start
% Clear workspace
clearvars

%----------------------------PROGRAM---------------------------------------

%----------------------Read Excell Data------------------------------------
run('excell_data_reader_ref.m')
% run('excell_data_reader_real.m')
%----------------------Calculate Mach number-------------------------------
run('machnumber.m')
%----------------------Obtain the Mass Balance-----------------------------
run('MassBalance.m')
%----------------------Obtain Engine Data----------------------------------
run('Engine_part.m')
% %----------------------Obtain Aerodynamic coefficients-------------------
run('Aerodynamic_Coefficients.m')
% %----------------------Run Flightdata Graphmaker-------------------------
run('flightdata_graphmaker.m')
% %----------------------Run Eigenmotions----------------------------------
run('Leanne.m')
