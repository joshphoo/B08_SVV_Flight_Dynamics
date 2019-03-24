% --------- Main File
% Some extra info
% Start

%----------------------------PROGRAM---------------------------------------

%----------------------Read Excell Data------------------------------------
%run('excell_data_reader_ref.m')
run('excell_data_reader_real.m')
%----------------------Obtain parameters from Marloes----------------------
run('parametersmarloes.m')
%----------------------Obtain the Mass Balance-----------------------------
run('MassBalance.m')
%----------------------Obtain Engine Data-------------------------------------
run('Engine_part.m')
%----------------------Obtain Aerodynamic coefficients---------------------
run('Aerodynamic_Coefficients.m')
%----------------------Run Flightdata Graphmaker---------------------------
run('flightdata_graphmaker.m')
