% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
% |                       PH-LAB Simulatie main.m                         | 
% |                                 2019                                  | 
% |                             Brakke Willy                              | 
% +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 

% Start
% Clear workspace
clearvars

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%! ZORG DAT JE DE JUISTE EXCELL_DATA_READER.M EN FLIGHTDATA_GRAPHMAKER.M  ! 
%! GEBRUIKT: REAL VOOR WERKELIJKE TESTVLUCHT, REF VOOR REFERENTIEDATA.    !
%! COMMENT DE ANDERE(NIET DELETEN)                                        !
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%----------------------------PROGRAM---------------------------------------

%----------------------Read Excell Data------------------------------------
% run('excell_data_reader_ref.m')
run('excell_data_reader_real.m')
%----------------------Calculate Mach number-------------------------------
run('machnumber.m')
%----------------------Obtain the Mass Balance-----------------------------
run('MassBalance.m')
%----------------------Obtain Engine Data----------------------------------
run('Engine_part.m')
% %----------------------Obtain Aerodynamic coefficients-------------------
run('Aerodynamic_Coefficients.m')
% %----------------------Run Flightdata Graphmaker-------------------------
% run('flightdata_graphmaker.m')
run('real_flightdata_graphmaker.m')     %Enige verschil is dat deze een min heeft voor de asym. angles
% %----------------------Run Eigenmotions----------------------------------
run('Leanne.m')
%run('Leanne_real.m')
