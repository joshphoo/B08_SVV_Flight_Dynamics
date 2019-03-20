% --------- Main File
% Some extra info
% Start

%----------------------------PROGRAM---------------------------------------

%----------------------Read Excell Data------------------------------------
run('excell_data_reader_ref.m')
run('excell_data_reader_real.m')
%----------------------Obtain parameters from Marloes----------------------
run('parametersmarloes.m')
%----------------------Run Engine Part-------------------------------------
run('Engine_part.m')
%----------------------Run Flightdata Graphmaker---------------------------
run('flightdata_graphmaker.m')
